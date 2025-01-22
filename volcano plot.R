library(ggplot2)
library(openxlsx)
setwd("d:/R/Met")
x<- read.xlsx("Met.xlsx", sheet=1, rowNames = T)
logFC <- log2(x$Fold.change)
padj <- x$pvalue
data <- data.frame(logFC=logFC,padj=padj)
data$sig[(data$padj >0.05|data$padj=="NA")|(data$logFC < 1& data$logFC > -1)] <- "no"
data$sig[data$padj <= 0.05 & data$logFC >= 1] <- "up"
data$sig[data$padj <=0.05 & data$logFC <= -1] <- "down"
theme_set(theme_bw())
p <- ggplot(data,aes(logFC,-1*log10(padj),
                     color = sig))+geom_point()+
  xlim(-10,10) + labs(x="log2(FoldChange)",y="-log10(p-value)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4) +geom_point(size=2)
p <- p +theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))+ylim(0,6)
p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
p

