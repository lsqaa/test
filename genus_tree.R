evolution <- function(tree, finalFrame, annotFrame, outdir){
    colors <- c('#0000AA', '#FF0000', '#006400', '#00BFFF',
                '#388E8E', '#458B00', '#FF34B3', '#FF7F00', '#6A5ACD', '#8DEEEE',
                '#B03060', '#4A708B', '#7D26CD', '#000000', '#4876FF', '#8B475D')
    phylums <- unique(finalFrame$Phylum)
	colors <- getDiscreteColor(length(phylums) + 1)
    phyCol <- data.frame(colors = colors[1:length(phylums)], Phylum = unique(finalFrame$Phylum))
    finalFrame <- inner_join(finalFrame, phyCol, by = "Phylum")
    dat1 <- finalFrame %>% 
            select(c("label", "Phylum", "colors")) %>%
            distinct() %>%
            rename(phylum = Phylum)
    dat1$phylum <- factor(dat1$phylum, levels = sort(unique(dat1$phylum)))

    dat2 <- aggregate(. ~ Phylum, annotFrame, FUN = paste, collapse = ",")
    labels <- lapply(dat2$label, function(x) {unlist(strsplit(x, split = ","))})
    names(labels) <- dat2$Phylum

    dat3 <- finalFrame %>%
            select(c("label", "Group", "Value"))
    dat3$Group <- factor(dat3$Group)

    tr <- groupOTU(tree, labels, "Phylum")
    Phylum <- NULL
    p <- ggtree(tr = tr, layout = "circular", open.angle = 15, 
                size = 0.2, aes(colour = Phylum), show.legend = F) +
         theme(legend.position = 'none')
    p1 <- p %<+% dat1 +
                 geom_tippoint(aes(colour = phylum), alpha=0) +
                 geom_tiplab(aes(colour = phylum),
                             align = T,
                             linetype = 3,
                             size = 4,
                             linesize = 0.2,
                             show.legend = F) +
                 scale_colour_manual(name = "Phylum",
                                     values = colors,
                                     breaks = unique(dat1$phylum),
                                     guide = guide_legend(
                                                keywidth = 0.5,
                                                keyheight = 0.5,
                                                order = 2,
                                                override.aes = list(
                                                                size = 2, 
                                                                alpha = 1)))
    p2 <- p1 + geom_fruit(data = dat3,
                          geom = geom_bar,
                          stat="identity", 
                          orientation = "y",  
                          position = position_stackx(), 
                          mapping = aes(
                                y = label,
                                x = Value,
                                fill = Group),
                          pwidth = 0.4, 
                          offset = 1.8) +
               scale_size_continuous(range = c(1, 3),                                
                                     guide = guide_legend(
                                         keywidth = 0.5, 
                                         keyheight = 0.5, 
                                         override.aes = list(starshape = 15), 
                                         order = 2)) +
               theme(legend.position = c(0.9, 0.85), plot.margin = margin(0, 0, 0, 0),
                     legend.title = element_text(size=14),                 
                     legend.text = element_text(size=12),                 
                     legend.spacing.y = unit(0.02, "cm"),
                     legend.margin = unit(0, "lines"))
    
    pdfFigure <- "genus_group100.tree.pdf"
    cairo_pdf(pdfFigure, height = 10, width = 10)
    print(p2)
    dev.off()
}

################################################################################
#
# Load package and initialization
#
################################################################################
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggtreeExtra))
suppressMessages(library(ggtree))
suppressMessages(library(phyloseq))
suppressMessages(library(treeio))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(ggnewscale))
suppressMessages(library(tools))

################################################################################
#
# Main code, data import and visialization 
#
################################################################################
# Load OTU table
genusFrame <- read.table("genus_quant.txt", header = T, sep = '\t', quote = '', check.names = F, fileEncoding = "UTF-8", comment.char = '')
annotFrame <- read.table("phylum_genus.annot.txt", header = T, sep = '\t', quote = '', check.names = F, fileEncoding = "UTF-8", comment.char = '')
groupFrame <- read.table("sample.txt", header = T, sep = '\t', quote = '', check.names = F, fileEncoding = "UTF-8", comment.char = '')
tree <- read.tree("genus.tre")
genusFrame <- genusFrame %>%
              group_by(Genus) %>%
              summarise(across(everything(), sum)) %>%
              ungroup() %>%
              as.data.frame() %>%
              rename(., label = Genus)
genusFrame$label <- str_replace_all(genusFrame$label, " ", "_")
genusFrame <- genusFrame[str_detect(genusFrame$label, 
                                    "uncultured|unassigned|Unassigned", 
                                    negate = T), ]
genusMeltFrame <- melt(genusFrame, 
                       id = "label", 
                       variable.name = "Sample", 
                       value.name = "Value")
mergeFrame <- inner_join(genusMeltFrame, groupFrame, by = "Sample")
mergeFrame <- subset(mergeFrame, select=c("Group","label","Value"))
meanFrame <- mergeFrame %>%
             group_by(Group, label) %>%
             summarise(Value = mean(Value)) %>%
             ungroup() %>%
             as.data.frame()

annotFrame <- annotFrame %>%
              select(Phylum, Genus)
annotFrame$Genus <- str_replace_all(annotFrame$Genus, " ", "_")
annotFrame <- annotFrame[str_detect(annotFrame$Genus, 
                                    "uncultured|unassigned|Unassigned", 
                                    negate = T), ] 
annotFrame <- annotFrame %>%
              rename(., label = Genus) %>%
              distinct()
finalFrame <- inner_join(meanFrame, annotFrame, by = "label")

evolution(tree, finalFrame, annotFrame, "./")
