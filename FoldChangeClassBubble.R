# Load packages -----------------------------------------------------------------------------------
suppressMessages(library(tidyverse))
suppressMessages(library(Cairo))
suppressMessages(library(RColorBrewer))
suppressMessages(library(readxl))
suppressMessages(library(plotly))
suppressMessages(library(htmlwidgets))


# Parse command line arguments --------------------------------------------------------------------
log2FC.cutoff <- log2(1) %>% abs

# Load data and transform
data4plot.tbl <- read_tsv("infile.txt",
                             col_types = cols(),
                             locale = locale(encoding = 'UTF-8'),
                             na = c('', 'N/A'))
data4plot.tbl <- mutate_at(data4plot.tbl, c('Log2FC', 'VIP'), as.double) %>%
                     arrange(desc(Log2FC)) %>%
                     transmute(Index = Index,
                               Compounds = Compounds,
                               Class = .data[[class.name]],
                               Log2FC = Log2FC,
                               VIP = VIP)

# Set bubble fill color palette -------------------------------------------------------------------
color.pal <- strsplit("#ffd92f,#4daf4a,#377eb8,#984ea3", split = ',')[[1]]

class.num <- data4plot.tbl$Class %>% unique %>% length
class.color <- colorRampPalette(color.pal)(class.num)


# Plot static bubble ------------------------------------------------------------------------------
# Set x-axis limits
log2FC.vec <- data4plot.tbl$Log2FC %>% abs %>% as.vector
log2FC.vec <- log2FC.vec[!is.infinite(log2FC.vec)]
if (length(log2FC.vec) == 0) {
    xlim <- 16
} else {
    xlim <- max(log2FC.vec, na.rm = TRUE) %>% ceiling
}

# Plot PDF
p <- ggplot(data = data4plot.tbl,
            aes(x = Log2FC, y = Class, fill = Class, size = VIP,
                index= Index, compound = Compounds))
p <- p + geom_point(shape = 21, alpha = .7)
p <- p + scale_fill_manual(values = class.color)
p <- p + theme_bw()
p <- p + theme(panel.grid = element_blank(),
               panel.border = element_blank(),
               legend.text = element_text(size = 22),
               axis.line = element_line(color = 'black', size = .5, lineend = 'square'))
p <- p + scale_y_discrete(limits = sort(unique(data4plot.tbl$Class), decreasing = TRUE))
p <- p + scale_x_continuous(n.breaks = 7, limits = c(-xlim, xlim))
# p <- p + geom_vline(xintercept = -1, color = '#525252', linetype = 'longdash', size = .5, alpha = .5)
# p <- p + geom_vline(xintercept = 1, color = '#525252', linetype = 'longdash', size = .5, alpha = .5)
p <- p + geom_vline(xintercept = 0, color = '#525252', linetype = 'longdash', size = .5, alpha = .5)
static <- p + labs(x = expression(paste(Log[2], '(Fold Change)')),
                   size = 'Variable Importance in Projection (VIP)')
static <- static + guides(size = guide_legend(order = 1), fill = guide_legend(order = 2,override.aes = list(size = 3), ncol = 1)) + theme(text = element_text(size = 30))
output.pdf <- 'Log2FC_Class_Scatter.pdf'
CairoPDF(file = output.pdf, width = 21, height = 13)
print(static)
graphics.off()

output.png <- 'Log2FC_Class_Scatter.png'
png(file = output.png, width = 21, height = 13, res = 300, units = 'in')
print(static)
graphics.off()
