# Load package and initialization ---------------------------------------------
suppressMessages(library(tidyverse))
suppressMessages(library(phyloseq))
suppressMessages(library(ape))
suppressMessages(library(readxl))
suppressMessages(library(writexl))
suppressMessages(library(plotly))
suppressMessages(library(htmlwidgets))
suppressMessages(library(ggrepel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(scatterplot3d))

distance <- 'bray'
set.seed(100)


# Load sample frame and quantitative data -------------------------------------
sampleFrame <- read.table("sample.txt", header = T, sep = '\t')
sampleFrame <- sampleFrame[, c('Sample', 'Group')]
quantFrame <- read.table("infile.txt", header = T, sep = '\t', quote = '', check.names = F, fileEncoding = "UTF-8", comment.char = '')

quantFrame <- extractSampleData(quantFrame, sampleFrame$Sample)

# Set index as row names, and delete feature with all same data
quantFrame <- quantFrame %>% column_to_rownames(var = colnames(quantFrame)[1])
groupNum <- length(unique(sampleFrame$Group))

# Calculate number of group contains more than 3 sample
groupCount <- sampleFrame %>%
              group_by(Group) %>%
              summarize(count = n())


# Ordination analysis --------------------------------------------------------
## Build phyloseq object
sampleIndexFrame <- sampleFrame %>% column_to_rownames(var = 'Sample')
quantObj <- otu_table(quantFrame, taxa_are_rows = TRUE)
sampleObj <- sample_data(sampleIndexFrame)

phyloseqObj <- phyloseq(quantObj, sampleObj)

## Calc distance
distanceObj <- distance(phyloseqObj, distance)

## Ordination analysis
ordinationList <- ordinate(phyloseqObj,
                           method = "PCoA",
                           distance = distanceObj)
## Write score results
allComponentFrame <- as.data.frame(ordinationList$vectors) %>%
                         rownames_to_column(var = 'Sample')
colnames(allComponentFrame) <- gsub('Axis.', 'PCoA', colnames(allComponentFrame))
## Add group information
allComponentFrame <- allComponentFrame %>%
                     inner_join(sampleFrame, by = "Sample")


# Plotting static 2D plot ----------------------------------------------------
p <- ggplot(allComponentFrame, aes(PCoA1, PCoA2))
p <- p + geom_point(aes(color = Group), size = 3)

suppressMessages(library(vegan))
adonisResult <- vegan::adonis2(distanceObj ~ Group, sampleFrame, 
                            permutations = 999, by = NULL)
R2 <- round(as.numeric(adonisResult$R2[1]), 3)
pvalue <- formatPvalue(adonisResult$`Pr(>F)`[1])
statFlag <- paste0('Adonis R2: ', R2, ', P-value: ', pvalue)
eigVector <- round(ordinationList$values$Relative_eig * 100, 2)
p <- p + labs(subtitle = statFlag,
                  x = paste0("PCoA1 (", eigVector[1], "%)"),
                  y = paste0("PCoA2 (", eigVector[2], "%)"))
p <- p + labs(title = paste0('2D ', 'PCoA', ' Plot'))
p <- p + theme_bw()
p <- p + theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 23),
               plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, size = 21),
               legend.position = "right",
               legend.box = "horizontal",
               text = element_text(size = 16),
               axis.text = element_text(size = 18), axis.title = element_text(size = 20),
               legend.text = element_text(size = 18), legend.title = element_text(size = 20),
               # 删除次要网格线
               panel.grid.minor = element_blank())

# Add ellipse if group contain more than 3 samples
threeRepeatGroups <- filter(groupCount, count == 3)
### Plot ellipse for group with 3 samples
suppressMessages(library(car))
validFrame <- allComponentFrame %>%
                  filter(Group %in% threeRepeatGroups$Group) %>%
                  mutate(Group = factor(Group, levels = unique(Group)))
# Use car::dataEllipse to calculate coordinate of ellipse
set.seed(100)
ellipseList <- dataEllipse(validFrame[, 2], validFrame[, 3], col = groupColors,
                               validFrame$Group, levels = 0.95, robust = TRUE)
ellipseList <- lapply(names(ellipseList),
                          function(x){data.frame(ellipseList[[x]], Group = x)})
ellipseFrame <- bind_rows(ellipseList)
ellipseFrame <- data.frame(Sample = ellipseFrame$Group, ellipseFrame)
# Draw ellipse
p <- p + geom_polygon(data = ellipseFrame, aes(x, y, fill = Group),
                          alpha = 0.2, show.legend = FALSE)
### Plot ellipse for group with more than 3 samples
p <- p + stat_ellipse(aes(fill = Group), show.legend = FALSE,
                      geom = "polygon", alpha = 0.2)
# Save image to file
static2DPDF.ellipse <- "PCoA_ellipse.pdf"
static2DPNG.ellipse <- "PCoA_ellipse.png"
# Save pdf format
width <- 9
height <- 6
cairo_pdf(file = static2DPDF.ellipse, width = width, height = height)
cp <- dev.cur()
png(static2DPNG.ellipse, type = 'cairo', width = width,
    height = height, units = "in", res = 300)
dev.control("enable")
print(p)
dev.copy(which = cp)
graphics.off()

