# Load R packages
suppressMessages(library(tidyverse))
library(Cairo)
library(argparse)
options(stringsAsFactors=F)

inputData <- read.table("infile.txt", sep='\t', header=T,quote='',check.names=F)

geneNum <- as.numeric(lapply(inputData[,4], function (x) unlist(strsplit(x, ' '))[1]))
bgGeneNum <- as.numeric(lapply(inputData[,5], function (x) unlist(strsplit(x, ' '))[1]))
ratio <- round(geneNum/bgGeneNum, 2)
inputData <- data.frame(keggLevel1 = inputData[,1],
                        keggPathway = inputData[,2],
                        number = geneNum,
                        ratio = ratio,
                        Pvalue = inputData[,6],
                        qvalue = inputData[,7])
inputData <- inputData %>%
             rename('P-value' = Pvalue,
                    'Adjusted P-value' = qvalue)
# Extract the top 20 significant KEGG pathways
## Sort by Pvalue, keggPathway
inputData <- inputData %>%
             arrange(`P-value`, keggPathway)
## Extract the top 20
if (length(inputData[,1]) > 20)
  inputData <- inputData[ c(1:20), ]
# Sort by Pvalue
inputData <- inputData %>%
             arrange(desc(`P-value`), desc(keggPathway)) %>%
             mutate(keggPathway = factor(keggPathway, levels = keggPathway))

# Customize left margin of the picture
maxPathwayLength <- max(nchar(as.character(inputData[,2])))
if (maxPathwayLength < 25) {
  leftMarginScaleSize <- 30
} else if (maxPathwayLength < 50){
  leftMarginScaleSize <- 20
} else {
  leftMarginScaleSize <- 18
}
leftMargin <- maxPathwayLength/leftMarginScaleSize

# Customize picture width
 width <- ( 5 + leftMargin )

# Customize picture height
annotItemNum <- length(inputData[,1])
if ( annotItemNum <= 10 ) {
    height <- annotItemNum/2.6
} else if ( annotItemNum <= 20 ) {
    height <- annotItemNum/4
}
if (height < 4)
  {height = 4}

# Plot
suffix = 'P-value'
p <- ggplot(inputData, aes(keggPathway, ratio, colour = `P-value`))
p <- p + geom_point(aes(size = number))
p <- p + theme_bw()
p <- p + theme(panel.grid = element_blank())
p <- p + theme(panel.border = element_rect(color='black', fill=NA))
p <- p + scale_color_gradientn(colours = c('#EE0000', '#FF8247', '#76EE00', '#00F5FF', '#B23AEE'),
                               guide = 'colourbar', limits = c(0, 1.02))
diffNumVector <- inputData$number
i <- 1
while (TRUE){
	legendKeyVector <- seq(min(diffNumVector), max(diffNumVector), by = i)
	if (length(legendKeyVector) <= 5){
		break
	}
	i <- i + 1
}
p <- p + scale_size_continuous(breaks = legendKeyVector)
p <- p + labs(x = '', y = 'Rich Factor', size = 'Count')
p <- p + guides(colour = guide_colorbar(order = 1),size = guide_legend(order = 2))
p <- p + coord_flip()
p <- p + theme(text = element_text(size = 13))
# Save png
ggsave(paste('KEGG_Enrichment_', suffix, '.png', sep=''), width = width, height = height, units = 'in')
# Save pdf
ggsave(paste('KEGG_Enrichment_', suffix, '.pdf', sep=''), width = width, height = height, units = 'in')

