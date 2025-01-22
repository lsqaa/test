library(ggplot2)
library(ggprism)
setwd("d:wangliyi/R/gene")

x <- read.table(file="KEGG2.txt",sep="\t",header=T,check.names=FALSE)
head(x)
p <-ggplot(x,aes(x= log10pvalue, y=reorder(KEGG, log10pvalue), size=genecount, fill= qvalue))+
  geom_point(color="black", shape=21)+ 
  scale_size_continuous(range = c(0.5, 10))+ 
  labs (subtitle = "KEGG enrichment ", x = "-Log10(P-value) ", y = "") +
  scale_fill_gradient(low = "#E1A123", high = "#B9402D")+
  theme_bw()+theme(panel.grid=element_blank(), panel.border=element_blank(), axis.line = element_line(), axis.text.x = element_text(angle=90))
p
