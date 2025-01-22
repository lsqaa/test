tTest <- function(dataFrame, groupFrame, groupVector, prefix, adjustMethod, top, height){
    firstColumnTitle <- colnames(dataFrame)[1]
    controlGroup <- groupVector[1]
    caseGroup <- groupVector[2]
    # Extract data of needed samples
    dataFrame <- column_to_rownames(dataFrame, var = firstColumnTitle)
    dataFrame <- dataFrame[, groupFrame$Sample, drop = FALSE]
    # Write original OTU data
    levelFile <- paste0(prefix, ".level.txt")
    tmpFrame <- rownames_to_column(dataFrame, var = firstColumnTitle)
    write.table(tmpFrame, levelFile, quote = F, sep = '\t', row.names = F)
    # Extract sample of each group
    controlSamples <- groupFrame[groupFrame$Group == controlGroup, ]$Sample
    caseSamples <- groupFrame[groupFrame$Group == caseGroup, ]$Sample
    ## juFrame if all group length bigger than 2
    result <- NULL
    if (length(controlSamples) >= 2 && length(caseSamples) >= 2){
        for(species in rownames(dataFrame)){
            g1Value <- unlist(dataFrame[species, controlSamples, drop = T])
            g2Value <- unlist(dataFrame[species, caseSamples, drop = T])
            g1Value <- as.numeric(format(g1Value, scientific = TRUE))
            g2Value <- as.numeric(format(g2Value, scientific = TRUE))
            g1Unique <- length(unique(g1Value))
            g2Unique <- length(unique(g2Value))
            if (g1Unique != 1 && g2Unique != 1){
                tResult <- t.test(g1Value, g2Value, paired = F,
                                  alternative = "two.sided", 
                                  cond.lvel=0.95)
                g1Mean <- mean(g1Value)
                g2Mean <- mean(g2Value)
                g1Sd <- sd(g1Value)
                g2Sd <- sd(g2Value)
                pValue <- tResult$p.value
                lower <- tResult$conf.int[1]
                upper <- tResult$conf.int[2]
                content <- c(species, g1Mean, g1Sd, g2Mean, g2Sd, 
                                lower, upper, pValue)
                result <- rbind(result, content)
            }
        }
    }
    if (is.null(result)){
        return()
    }
    colnames(result) <- c(firstColumnTitle, 
                          paste0("Mean(", controlGroup, ")"), 
                          paste0("SD(", controlGroup, ")"), 
                          paste0("Mean(", caseGroup, ")"), 
                          paste0("SD(", caseGroup, ")"), 
                          "Interval lower", "Interval upper", 
                          "P-value")
    rownames(result) <- NULL
    result <- as.data.frame(result)
    # Write stat results
    result$`P-value` <- as.numeric(result$`P-value`)
    sortResult <- result[order(result$`P-value`), ]
    sortResult$`Q-value` <- p.adjust(sortResult$`P-value`, method = adjustMethod)
    sortResult <- na.omit(sortResult)
    # Add quant data
    dataFrame <- dataFrame %>% rownames_to_column(var = firstColumnTitle)
    sortResult <- left_join(sortResult, dataFrame)
    testFile <- paste0(prefix, ".T-test.txt")
    write.table(sortResult, testFile, quote = F, sep = '\t', row.names = F)
    testPvalueFile <- paste0(prefix, ".T-test.psig.txt")
    sortResultPvalue <- sortResult[sortResult$`P-value` < args$pvalue, , 
                                   drop = F]
    write.table(sortResultPvalue, testPvalueFile, quote = F, sep = '\t', row.names = F)
    testQvalueFile <- paste0(prefix, ".T-test.qsig.txt")
    sortResultQvalue <- sortResult[sortResult$`Q-value` < args$qvalue, , drop = F]
    write.table(sortResultQvalue, testQvalueFile, quote = F, sep = '\t', row.names = F)
    if(nrow(sortResultPvalue) >= 1){
        filePrefix <- str_replace(testPvalueFile, '.txt', '')
        combinePlot(sortResultPvalue, controlGroup, caseGroup, filePrefix, top = top, height = height)
    }
}

formatDecimal <- function(x){
    if (abs(x) == 0){
        return(x)
    } else if (abs(x) < 0.01){
        return(sprintf('%.1e', x))
    } else {
        return(signif(x, 2))
    }
}

calcRangeStep <- function(x){
    x.scientific <- format(abs(x), scientific = TRUE)
    magnitude <- as.integer(strsplit(x.scientific, 'e')[[1]][2])
    for (index in seq(10)){
        step <- index * 10 ** magnitude
        if (abs(x)/step <= 3){
            return(index * 10 ** magnitude)
        }
    }
}

generateAxis <- function(min, max){
    if (max <= 0){
        step <- calcRangeStep(min)
        axisVector <- - rev(seq(0, abs(min) + step, step))
    } else if (min >= 0){
        step <- calcRangeStep(max)
        axisVector <- seq(0, max + step, step)
    } else if (min < 0 && max > 0) {
        leftAxisStep <- calcRangeStep(min)
        rightAxisStep <- calcRangeStep(max)
        step <- max(leftAxisStep, rightAxisStep)
        leftAxisVector <- - rev(seq(0, abs(min) + step, step))
        rightAxisVector <- seq(0, max + step, step)
        axisVector <- unique(c(leftAxisVector, rightAxisVector))
    }
    return(axisVector)
}

# Barplot and CI-dot plot combination
combinePlot <- function(tTestFrame, group1, group2, prefix, top = 30, height = NULL){
    cbbPalette <- c("#E69F00", "#56B4E9")
    names(cbbPalette) <- c(group1, group2)
    colnames(tTestFrame)[1:9] <- c("Taxa", "column2", "__internal_1", "column4", "__internal_2", 
                                   "Lower", "Upper", "P-value", "Qvalue")
    tTestFrame <- tTestFrame %>% 
                  mutate(column2 = as.numeric(column2),
                         column4 = as.numeric(column4))
    # Sort by P-value and filter
    tTestFrame <- tTestFrame %>%
                  arrange(`P-value`)
    if(nrow(tTestFrame) > top){
        tTestFrame <- tTestFrame[1:top, ]
    }
    # Adjust height of images
    if (is.null(height)){
        diffNum <- nrow(tTestFrame)
        if (diffNum <= 3){
            height <- 3
        } else {
            height <- (diffNum - 3) * 0.5 + 3
        }
    }

    data <- tTestFrame %>% 
            mutate(Diff = column2 - column4) %>%
            mutate(Group = ifelse(column2 > column4, group1, group2)) %>% 
            mutate(`P-value` = sapply(`P-value`, formatDecimal)) %>%
            select(Taxa, column2, column4, Diff, `P-value`, Lower, Upper, Group) %>%
            arrange(Diff)
    if (args$separator == 'not'){
        data$NewTaxa <- data$Taxa
    } else {
        data$NewTaxa <- guessAndSplitVector(data$Taxa, len = 2)
    } 
    data <- data %>% mutate(NewTaxa = factor(NewTaxa, levels = NewTaxa))
    data <- subset(data, select = -c(Taxa))
    data$Lower <- as.numeric(data$Lower)
    data$Upper <- as.numeric(data$Upper)
    colnames(data)[which(names(data) == "column2")] <- group1
    colnames(data)[which(names(data) == "column4")] <- group2
    min <- range(data$Diff)[1]
    max <- range(data$Diff)[2]
    if(min > 0){
        min <- min / 1.2
    } else {
        min <- min * 1.2
    }
    if(max > 0){
        max <- max * 1.2
    } else {
        max <- max / 1.2
    }
    diffAxisVector <- generateAxis(min, max)
    diffAxisLabel <- sapply(diffAxisVector, formatDecimal)

    p1 <- ggplot(data, aes(NewTaxa, Diff, fill = Group)) +
          scale_x_discrete(limits = levels(data$NewTaxa)) +
          scale_y_continuous(breaks = diffAxisVector,
                             labels = diffAxisLabel) +
          labs(x = NULL, y = "Difference between groups", 
               title = "95% confidence intervals") +
          coord_flip() + 
          theme(panel.background = element_rect(fill = 'transparent'),
                panel.grid = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y = element_blank(),
                axis.text = element_text(size = 19),
                axis.text.x = element_text(vjust = -1),
                axis.text.y = element_blank(),
                axis.title.x = element_text(size = 23, vjust = -1),
                legend.position = 'none',
                plot.title = element_text(size = 23, vjust = 0.5, hjust = 0.5))
    if (nrow(data) > 1){
        for (i in 1:(nrow(data) - 1)){
            p1 <- p1 + annotate('rect', xmin = i + 0.5, xmax = i + 1.5, 
                                ymin = -Inf, ymax = Inf,
                                fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
        }       
    }
    p1 <- p1 + geom_point(shape = 21, size = 5) + 
          geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                        size = 0.1, width = 0.1) + 
          scale_fill_manual(values = cbbPalette) + 
          geom_hline(yintercept = 0, linetype = "dashed")

    data2 <- data[, c("NewTaxa", group1, group2)]
    data2 <- gather(data2, Group, Value, -NewTaxa)   
    data2$Value <- as.numeric(data2$Value)
    data2$Group <- factor(data2$Group, levels = rev(unique(data2$Group)))
    meanAxisVector <- generateAxis(0, max(data2$Value) * 1.2)
    meanAxisLabel <- sapply(meanAxisVector, formatDecimal)

    p2 <- ggplot(data2, aes(x = NewTaxa, y = Value, fill = Group)) + 
          scale_x_discrete(limits = levels(data$NewTaxa)) +
          scale_y_continuous(breaks = meanAxisVector,
                             labels = meanAxisLabel) +
          labs(x = NULL, y = "Mean in groups") + 
          coord_flip(expand = FALSE) + 
          theme(panel.background = element_rect(fill = 'transparent'),
                panel.grid = element_blank(),
                axis.text = element_text(size = 23),
                axis.text.x = element_text(vjust = -1),
                axis.text.y = element_text(color = 'black'),
                axis.title.x = element_text(size = 23, vjust = -1),
                legend.position = 'top',
                legend.direction = "horizontal",
                legend.title = element_blank(),
                legend.text = element_text(size = 23,
                                           margin = margin(r = 20)))
    if (nrow(data) > 1){
        for (i in 1:(nrow(data) - 1)){
            p2 <- p2 + annotate('rect', xmin = i + 0.5, xmax = i + 1.5, 
                                ymin = -Inf, ymax = Inf,
                                fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
        }
    }
    p2 <- p2 + geom_bar(stat = 'identity', width = 0.3,
                        position = position_dodge(0.3)) + 
          scale_fill_manual(values = cbbPalette)
    
    p3 <- ggplot(data, aes(NewTaxa, Diff, fill = Group)) +
          geom_text(aes(y = 0, x = NewTaxa), 
                    label = data$`P-value`,
                    hjust = 0, 
                    inherit.aes = F, 
                    size = 8.5) +
          geom_text(aes(x = nrow(data) / 2 + 0.5, y = 0.85), 
                    label = "P-value",
                    srt = 90, 
                    size = 8.5) +
          coord_flip() +
          ylim(c(0, 1)) +
          theme(panel.background = element_blank(),
                panel.grid = element_blank(),
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                axis.title = element_blank())

    pdfFigure <- paste0(prefix, ".pdf")
    cairo_pdf(pdfFigure, width = 20, height = height)
    p <- p2 + p1 + p3 + plot_layout(widths = c(4, 6, 2))
    print(p)
    dev.off()
    pngFigure <- str_replace(pdfFigure, ".pdf$", ".png")
    convertCMD <- 'convert -alpha off -density 300 -quality 100 '
    system(paste0(convertCMD, pdfFigure, ' ', pngFigure))
}

# Load package and initialization ---------------------------------
suppressMessages(library(Cairo))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrepel))
suppressMessages(library(qvalue))
suppressMessages(library(readxl))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))
suppressMessages(library(tools))
suppressMessages(library(patchwork))

# Main code, data import and visialization ---------------------------------

# Load original OTU and group
totalOTUFrame <- read.table("infile.txt", header = T, sep = '\t', quote = '', check.names = F, fileEncoding = "UTF-8", comment.char = '')
totalGroupFrame <- read.table("sample.txt", header = T, sep = '\t')
# Extract sample information of diff group
diffGroupVector <- str_split("A,B", ',')[[1]]
groupFrame <- normalizeSampleFrame(totalGroupFrame, diffGroupVector)
groupFrame <- groupFrame %>%
              select(Sample, Group) %>%
              filter(Group %in% diffGroupVector)

# T-test analysis
prefix <- "./t-test"
tTest(totalOTUFrame, groupFrame, diffGroupVector, prefix, "fdr",
      30, NULL)
