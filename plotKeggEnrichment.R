# load R packages
suppressMessages(library(tidyverse))
library(Cairo)
library(argparse)

inputData <- read.table('infile.txt', sep = '\t', header = T, quote = '', check.names = F)

geneNum <- as.numeric(lapply(inputData$`Cluster_frequency`, function (x) unlist(strsplit(x, " "))[1]))
bgGeneNum <- as.numeric(lapply(inputData$`Metabolome_frequency`, function (x) unlist(strsplit(x, " "))[1]))
ratio <- round(geneNum / bgGeneNum, 2)
inputData <- data.frame(keggLevel1 = inputData$`Kegg_level_1`,
                        keggPathway = inputData$`Kegg_pathway`,
                        number = geneNum,
                        ratio = ratio,
                        Pvalue = inputData$`P-value`,
                        qvalue = inputData$`Corrected_P-value`)
# Extract the top 20 significant KEGG pathways
# Sort by Pvalue
inputData <- inputData %>%
             arrange(Pvalue, keggPathway)
# Extract the top 20
if (nrow(inputData) > 20)
  inputData <- inputData[c(1:20), ]
# diaplay pathway name Sort by Pvalue 
order<- sort(inputData$Pvalue,index.return=TRUE, decreasing = TRUE)
inputData$keggPathway <- factor(inputData$keggPathway, 
                                 levels = inputData$keggPathway[order$ix])

# Plot

p <- ggplot(inputData, aes(keggPathway, ratio, colour = Pvalue))

p <- p + geom_point(aes(size = number))+ 
         theme_bw()+ 
         theme(panel.grid = element_blank())+ 
         theme(panel.border = element_rect(color = 'black', fill = NA))+
         scale_color_gradientn(colours = c("#EE0000", "#FF8247", "#76EE00", "#00F5FF", "#B23AEE"),
                               guide = "colourbar",
                               limits = c(0, 1.02))

p <- p + scale_size_continuous(breaks = seq(min(inputData$number), max(inputData$number), by = 5))

p <- p + labs(x = "", y = "Rich Factor", colour = 'P-value', size = "Count")+
    guides(colour = guide_colorbar(order = 1),size = guide_legend(order = 2)) + 
    coord_flip() + 
    theme(text = element_text(size = 20))

# Save png
ggsave('KEGG_Enrichment.png', width = 12, height = 8, units = "in")
# Save pdf
ggsave('KEGG_Enrichment.pdf', width = 12, height = 8)
