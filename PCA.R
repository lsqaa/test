# Load package and initialization ---------------------------------------------
suppressMessages(library(tidyverse))
suppressMessages(library(readxl))
suppressMessages(library(writexl))
suppressMessages(library(plotly))
suppressMessages(library(htmlwidgets))
suppressMessages(library(ggrepel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(scatterplot3d))

# Load sample frame and quantitative data -------------------------------------
sampleFrame <- read.table("sample.txt", header = T, sep = '\t')
quantFrame <- read.table("infile.txt", header = T, sep = '\t', quote = '', check.names = F, fileEncoding = "UTF-8", comment.char = '')
quantFrame <- extractSampleData(quantFrame, sampleFrame$Sample)

# Set index as row names, and delete feature with all same data
quantFrame <- quantFrame %>% column_to_rownames(var = colnames(quantFrame)[1])

# Set colors according to group number
groupNum <- length(unique(sampleFrame$Group))


# PCA analysis ----------------------------------------------------------------
quantFrame <- t(quantFrame) %>% scale %>% t
pca <- prcomp(t(as.matrix(quantFrame)), center = FALSE)
pca <- summary(pca)
# Convert all components and variance proportion of PCA from matrix to dataFrame
allComponentFrame <- pca$x %>%
                     as.data.frame %>%
                     rownames_to_column(var = 'Sample') %>%
                     as_tibble
varianceRatioFrame <- pca$importance %>%
                      as.data.frame %>%
                      rownames_to_column(var = 'Description')

# Data preparation for plotting -----------------------------------------------
# Add group information
allComponentFrame <- allComponentFrame %>%
                     inner_join(sampleFrame, by = "Sample")
# Extract proportion of variance of each component
pc1VarianceRatio <- pca$importance[2, 1] * 100
if (ncol(pca$x) >= 2){
    pc2VarianceRatio <- pca$importance[2, 2] * 100
}
if (ncol(pca$x) >= 3){
    pc3VarianceRatio <- pca$importance[2, 3] * 100
}
# Calculate number of group contains more than 3 sample
groupCount <- sampleFrame %>%
              group_by(Group) %>%
              summarize(count = n())


# Plotting static 2D PCA ------------------------------------------------------
pdf(NULL)
p <- ggplot(data = allComponentFrame, aes(x = PC1, y = PC2, label = Sample))
p <- p + geom_point(aes(color = Group), size = 1) +
         geom_point(aes(color = Group), size = 4, alpha = 0.5)
p <- p + labs(title = '2D PCA Plot',
              x = paste0("PC1 (", round(pc1VarianceRatio, 2), "%)"),
              y = paste0("PC2 (", round(pc2VarianceRatio, 2), "%)"))
p <- p + theme_bw()
p <- p + theme(plot.title = element_text(hjust = 0.5))
p <- p + theme(legend.position = "right",
               legend.box = "horizontal", text = element_text(size = 13))

# Add ellipse if group contain more than 3 samples
threeRepeatGroups <- filter(groupCount, count == 3)
### Plot ellipse for group with 3 samples
pdf(NULL)
suppressMessages(library(car))
validFrame <- allComponentFrame %>%
                  filter(Group %in% threeRepeatGroups$Group) %>%
                  mutate(Group = factor(Group, levels = unique(Group)))
# Use car::dataEllipse to calculate coordinate of ellipse
set.seed(100)
ellipseList <- dataEllipse(validFrame$PC1, validFrame$PC2, col = pcaColors,
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
static2DPDF.ellipse <- "PCA_ellipse.pdf"
static2DPNG.ellipse <- "PCA_ellipse.png"
# Save pdf format
graphics.off()
cairo_pdf(file = static2DPDF.ellipse, width =  9, height = 6)
cp <- dev.cur()
png(static2DPNG.ellipse, width = 9, height = 6, units = "in", res = 300)
dev.control("enable")
print(p)
# Copy png from pdf device
dev.copy(which = cp)
graphics.off()

