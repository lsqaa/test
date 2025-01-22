# Functions -------------------------------------------------------------------
# Load package and initialization -------------------------------------
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(Cairo))
suppressMessages(library(ggrepel))

# load results of diff analysis, assign each metabolite/gene/protein to be 
#     up/down/insig according to filtering conditions ----------------------
df <- read_tsv("infile.txt", col_types = cols(), col_names = T, locale = locale(encoding = "utf-8"), quote ="", na = '')


df <- df %>% mutate(Type = case_when(Log2FC >= 0 ~ "up",
                                     Log2FC <= 0 ~ "down",
                                     TRUE ~ "insig"))


df <- df %>%
         mutate(`Log2FC_select` = ifelse(Log2FC >= 0 | 
                                         Log2FC <= 0, '-', NA))


df <- df %>%
         mutate(`Pvalue_select` = ifelse(`P-value` < 0.05, '-', NA))


df <- df %>%
          mutate(`VIP_select` = ifelse(VIP > 1, '-', NA))


filterCol <- c('Log2FC_select', 'Pvalue_select', 'VIP_select')

df[which(apply(df[filterCol], 1, function(x) any(is.na(x)))), "Type"] <- "insig"
df <- df %>% select(-one_of(filterCol))

# Prepare for plot ------------------------------------------------------
# Define text of legend
up <- filter(df, Type == 'up') %>% count() %>% unlist()
down <- filter(df, Type == 'down') %>% count() %>% unlist()
insig <- filter(df, Type == 'insig') %>% count() %>% unlist()

upM <- paste("Up: ", up, sep = "")
downM <- paste("Down: ", down, sep = "")
insigM <- paste("Insignificant: ", insig, sep = "")


df <- df %>% 
         mutate(Color = case_when(Type == 'up' ~ upM,
         Type == 'down' ~ downM, TRUE ~ insigM))


# Preprocess data frame and axis labels

df <- df %>% mutate(`-Log10pvalue` = -log(`P-value`,10))
x <- df$Log2FC
y <- df$`-Log10pvalue`
xlabel <- expression(paste(Log[2], '(Fold Change)'))
ylabel <- expression(paste(-Log[10], '(P-value)'))
xintercept <- 0
yintercept <- -log(0.05, 10)
x_max <- na.omit(df$Log2FC)
x_max <- x_max[is.finite(x_max)] %>% abs() %>% max()

# Plot ---------------------------------------------------------------------
g <- ggplot(data = df, mapping = aes(x = Log2FC, y = `-Log10pvalue`, size = VIP))

g <- g + geom_point(alpha = 0.85, aes(colour = Color))
g <- g + geom_vline(xintercept = xintercept,
                    linetype = 'longdash', colour = 'black')
g <- g + geom_hline(yintercept = yintercept,
                    linetype = 'longdash', colour = 'black')
g <- g + scale_size(range = c(0.05, 3))
g <- g + scale_colour_manual(name = "Statistics", values = c('#cc0000', '#006600', 'grey'),
                             breaks = c(upM, downM, insigM),
                             labels = c(upM, downM, insigM),
                             limits = c(upM, downM, insigM))
g <- g + scale_x_continuous(labels = c(-16, -8, -4, -2,
                                                0, 2, 4, 8, 16),
                            breaks = c(-16, -8, -4, -2,
                                       0, 2, 4, 8, 16),
                            limits = c(-x_max, x_max))
g <- g + theme_bw()
g <- g + theme(panel.grid.minor.x = element_blank())
g <- g + guides(colour =  guide_legend(order = 1),
                size = guide_legend(order = 0))
g <- g + labs(x = xlabel, y = ylabel)

g <- g + theme(text = element_text(size = 22))

png('vol.png',
    height = 6, width = 8, res = 300, units = 'in')
print(g)
dev.off()
cairo_pdf('vol.pdf',
          height = 8, width = 6)
print(g)
dev.off()