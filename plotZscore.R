# Load R packages
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
library(Cairo)

# read group and infile
gDF <- read_tsv('group.txt', col_types='cc', locale=locale(encoding='UTF-8'), na='') %>%
       mutate(Group = factor(Group,levels=unique(Group)))

inDF <- read_tsv('infile.txt', col_types=cols(), locale=locale(encoding='UTF-8'))

# Set colors 
colorSets <- c('#1b9e77', '#d95f02')
# Set shape
shapeSets <- c(21,22)

plotDF <- inDF %>% 
          arrange(desc(VIP)) %>% 
          head(50)

zPlotDF <- plotDF %>%
               column_to_rownames('Compounds') %>%
               select(c(gDF$Sample)) %>% 
               t %>% scale %>% t %>% 
               as.data.frame %>%
               rownames_to_column('Index') %>% 
               as_tibble %>%
               mutate(Index = factor(Index, levels=unique(Index))) %>%　# order
               select(c('Index', everything())) %>%
               gather(key='Sample', value='Z-score', -Index) %>%
               left_join(gDF, by='Sample')

g <- ggplot(zPlotDF, aes(x=`Z-score`, y=Index, colour=Group, label=Sample)) + 
     geom_point(alpha=0.7, aes(shape=Group, fill=Group)) + 
     scale_color_manual(values=colorSets) + 
     scale_fill_manual(values=colorSets) + 
     scale_shape_manual(values=shapeSets) +
     labs(y='', title='Z-score Plot') + 
     theme_bw() + 
     theme(panel.grid.major.x=element_blank(), 
           panel.grid.minor.x=element_blank(),
           legend.title=element_blank(),
           plot.title=element_text(hjust=.5)) + 
     theme(legend.position=c(1,1), 
           legend.justification=c(1,1),
           legend.background=element_rect(fill=alpha('white',0.3), 
                                          colour='black', size=0.3))
png('zscore.png', 
    width=6, height=6, res=300, units='in')
print(g)
dev.off()
cairo_pdf('zscore.pdf', 
         width=6, height=6)
print(g)
dev.off()

