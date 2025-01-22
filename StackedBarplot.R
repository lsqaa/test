# Function defination -------------------------------------

stackedBarplot <- function(dataFrame){
    sampleNames <- colnames(dataFrame)
    # Convert to longer table for ggplot2
    dataFrame <- dataFrame %>% rownames_to_column(var = firstColumn)
    longFrame <- pivot_longer(dataFrame, - all_of(firstColumn),
                              names_to = 'Sample', values_to = 'value') %>%
                 as.data.frame()
    ### Put Others on top
    longFrame <- rbind(filter(longFrame, get(firstColumn) %in% 'Others'),
                       filter(longFrame, ! get(firstColumn) %in% 'Others'))
    longFrame[, 1] <- factor(longFrame[, 1], levels = unique(longFrame[, 1]))
    longFrame <- transform(longFrame,
                           Sample = factor(Sample, levels = sampleNames))
    # Plot
    p <- ggplot(longFrame, mapping = aes(Sample, value, fill = get(firstColumn))) + 
         geom_bar(stat='identity', position='stack') + 
         labs(x = NULL, y = 'Relative abundance') + 
         theme_classic() + 
         theme(text = element_text(size = 16),
               axis.text = element_text(color = 'black', size = 18),
               axis.title = element_text(size = 20),
               legend.text = element_text(size = 18),
               legend.title = element_text(size = 20)) + 
         theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
         coord_cartesian(expand = FALSE)
    # Adjust width
        width <- ncol(dataFrame) * 0.28 + 
             max(str_length(dataFrame[, 1])) * 0.1 + 3
        width <- max(width, 6)
    # Adjust height
        height <- nrow(dataFrame) * 0.19 + 
                max(str_length(colnames(dataFrame)[-1])) * 0.06 + 2
        height <- max(height, 7)
    # Out file
    pdfFigure <- "barplot.pdf"
    pngFigure <- "barplot.png"
    # Save pdf format
    cairo_pdf(file = pdfFigure, family = args$fontfamily, width = width, height = height)
    cp <- dev.cur()
    png(pngFigure, family = args$fontfamily, type = 'cairo', width = width, 
        height = height, units = "in", res = 300)
    dev.control("enable")
    print(p)
    dev.copy(which = cp)
    graphics.off()
}

# Load package and initialization -------------------------------------------
suppressMessages(library(ggplot2))
suppressMessages(library(readxl))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))

# Main code, data import and visialization -----------------------------------
# Load relative abundance file
sampleFrame <- read.table("sample.txt", header = T, sep = '\t')
abundanceFrame <- read.table("infile.txt", header = T, sep = '\t', quote = '', check.names = F, fileEncoding = "UTF-8", comment.char = '')

firstColumn <- colnames(abundanceFrame)[1]
abundanceFrame <- column_to_rownames(abundanceFrame, var = firstColumn)
abundanceFrame <- abundanceFrame[, sampleFrame$Sample, drop = FALSE]

 # Rename low content features to Others
 featureNames <- rownames(abundanceFrame)
 sortedNames <- names(sort(rowSums(abundanceFrame), decreasing = TRUE))
 sortedNames <- sortedNames[sortedNames != 'Others']
 top10Name <- sortedNames[1:10]
 featureNames[!(featureNames %in% top10Name)] <- 'Others'
 # Merge all features of Others 
 abundanceFrame <- rownames_to_column(abundanceFrame, var = firstColumn)
 abundanceFrame[, 1] <- featureNames
 abundanceFrame <- abundanceFrame %>%
                      group_by_at(firstColumn) %>%
                      summarise(across(everything(), sum)) %>%
                      ungroup() %>% 
                      as.data.frame() %>%
                      column_to_rownames(var = firstColumn)

# Sort by sum of all samples, data would be out-of-order after group_by
sortedNames <- names(sort(rowSums(abundanceFrame), decreasing = FALSE))
abundanceFrame <- abundanceFrame[sortedNames, , drop = FALSE]

## Stacked barplot
stackedBarplot(abundanceFrame)
