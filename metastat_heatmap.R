drawComplexHeatmap <- function(heatmap, legends, locus, unit = 'inch') {
    # Draw ComplexHeatmap and get its width and height
    drawHeatmap <- draw(heatmap,
                        legend_grouping = 'adjusted',
                        annotation_legend_list = legends,
                        adjust_annotation_extension = FALSE,
                        align_heatmap_legend = locus,
                        padding = unit(c(.2, .4, .2, .4), 'inch'))  # set margins (bottom, left, top, right)
    width <- ComplexHeatmap:::width(drawHeatmap)
    height <- ComplexHeatmap:::height(drawHeatmap)
    width <- convertX(width, unit, valueOnly = TRUE)
    height <- convertY(height, unit, valueOnly = TRUE)
    return(list(width = width,
                height = height,
                heatmap = drawHeatmap))
}

# Load packages -----------------------------------------------------------------------------------
suppressMessages(library(tidyverse))
suppressMessages(library(readxl))
suppressMessages(library(writexl))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))

# Data preparation --------------------------------------------------------------------------------
# Load input data
inputDataFrame <- read.table("infile.txt", header = T, sep = '\t', quote = '', check.names = F, fileEncoding = "UTF-8", comment.char = '')
inputIndexColName <- colnames(inputDataFrame)[1]
colnames(inputDataFrame)[1] <- 'Index'  # rename index column name
# Initialize additional information
infoDataFrame <- transmute(inputDataFrame, Index = Index, Label = Index)

# Load group information and extract quantitative data
sampleDataFrame <- read.table("sample.txt", header = T, sep = '\t', quote = '', check.names = F, fileEncoding = "UTF-8", comment.char = '')
sampleVector <- sampleDataFrame[['Sample']]
quantDataFrame <- extractSampleData(inputDataFrame, sampleVector)
quantDataFrame <- column_to_rownames(quantDataFrame, var = 'Index')

# Z-score normalization
quantMatrix <- t(quantMatrix) %>% scale %>% t
colNum <- ncol(quantMatrix)
rowNum <- nrow(quantMatrix)

# Parse heatmap parameters ------------------------------------------------------------------------
# Whether to show row and column names
cellWidth <- 7
cellHeight <- 7
heatmapWidth = colNum * unit(cellWidth, 'mm')
heatmapHeight = rowNum * unit(cellHeight, 'mm')
colNameFontSize <- 12
rowNameFontSize <- 12
showRowName <- TRUE
showColName <- TRUE

# Heatmap body color, make sure middle legend correspond to 0
heatmapColorSet <- getContinuousColor(n = 9, colorPalette = "RdYlBu")
dataMin <- min(quantMatrix, na.rm = TRUE)
dataMax <- max(quantMatrix, na.rm = TRUE)
dataRange <- c(seq(dataMin, 0, length.out=5), seq(0, dataMax, length.out=5)[-1])
heatmapColorFun <- colorRamp2(dataRange, heatmapColorSet)
cellBorderColor <- '#e5e5e5'

# Set column (Group) annotation -------------------------------------------------------------------
colAnnotName <- "Group"
colAnnotOrder <- unique(sampleDataFrame[[colAnnotName]])
# Color blending
colAnnotUniqVector <- intersect(colAnnotOrder, sampleDataFrame[[colAnnotName]])
colAnnotNum <- length(colAnnotUniqVector)
colAnnotColors <- getDiscreteDarkColor(colAnnotNum)
names(colAnnotColors) <- colAnnotUniqVector
# HeatmapAnnotation
colAnnotVector <- sampleDataFrame[[colAnnotName]] %>% as.vector
names(colAnnotVector) <- sampleVector
colHmAnnot <- HeatmapAnnotation(
        colAnnot = colAnnotVector,
        annotation_label = colAnnotName,
        col = list(colAnnot = colAnnotColors),
        show_legend = FALSE,
        simple_anno_size = unit(.2, 'cm'),
        simple_anno_size_adjust = FALSE,
        annotation_name_gp = gpar(fontsize = 12),
        annotation_name_side = 'right')


# Set row label ----------------------------------------------------------------
rowNameLabel <- 'Index'
# Set row names to show
colNames <- colnames(quantMatrix)
# rowNames <- rownames(quantMatrix)
rowNames <- as.character(infoDataFrame[['Label']])
names(rowNames) <- as.character(infoDataFrame[['Index']])


# Set row annotation -----------------------------------------------------------
# Load or extract annotation columns
# Specify an annotation file
rowAnnotDataFrame <-  read.table("heatmap.signif.txt", header = T, sep = '\t', quote = '', check.names = F, fileEncoding = "UTF-8", comment.char = '')
infoDataFrame <- left_join(infoDataFrame, rowAnnotDataFrame, by = 'Index')
# Load format file of annot 
rowAnnotFormatFrame <- read.table("color_shape.txt", header = T, sep = '\t', quote = '', check.names = F, fileEncoding = "UTF-8", comment.char = '')

longerDataFrame <- pivot_longer(infoDataFrame, - c('Index', 'Label'),
                                        names_to = 'Column_name',
                                        values_to = 'Value_name')
# 绘制行注释
rowAnnotFrame <- infoDataFrame %>%
                         select(-Label) %>%
                         column_to_rownames(var = 'Index')
# Define cell function
cell_fun <- function(i, j, x, y, width, height, fill){
            for (k in seq(nrow(rowAnnotFormatFrame))){
                tmpRows <- rowAnnotFormatFrame[k, ]
                if (rowAnnotFrame[j, i] == tmpRows[, 1]){
                    if ('Shape' %in% colnames(rowAnnotFormatFrame)){
                        grid.points(x, y, pch = tmpRows$Shape,
                                    gp = gpar(col = tmpRows$Color,
                                              fill = tmpRows$Color))
                    } else {
                        grid.rect(x, y, width, height, 
                                  gp = gpar(fill = tmpRows$Color,
                                            col = cellBorderColor,
                                            lwd = .1))
                    }
                }
            }
        }
rowHmRightAnnot <- Heatmap(
            as.matrix(rowAnnotFrame),
            cell_fun = cell_fun,
            rect_gp = gpar(fill = NA, lty = 'blank'),
            row_labels = rowNames,
            row_names_gp = gpar(fontsize = rowNameFontSize),
            row_names_max_width = max_text_width(rowNames),
            column_labels = colnames(rowAnnotFrame),
            column_names_gp = gpar(fontsize = colNameFontSize),
            column_names_rot = -45,
            column_names_max_height = max_text_width(colnames(rowAnnotFrame)),
            width = ncol(rowAnnotFrame) * unit(5, "mm"),
            height = heatmapHeight,
            show_row_names = showRowName,
            show_heatmap_legend = FALSE)
        rowAnnotColors <- rowAnnotFormatFrame$Color
        names(rowAnnotColors) <- rowAnnotFormatFrame[, 1]


# Plot heatmap ------------------------------------------------------------------------------------
pdf(NULL)  # prevent generating Rplots.pdf
set.seed(666)
hm <- suppressWarnings(Heatmap(
    quantMatrix,
    col = heatmapColorFun,
    na_col = 'grey',
    cluster_rows = T,
    row_dend_width = unit(2, 'cm'),
    row_dend_gp = gpar(lwd = .5),
    cluster_columns = F,
    column_dend_height = unit(1.2, 'cm'),
    column_dend_gp = gpar(lwd = .5),
    show_row_names = F,
    row_labels = rowNames,
    row_names_gp = gpar(fontsize = rowNameFontSize),
    row_names_max_width = max_text_width(rowNames),
    show_column_names = T,
    column_labels = colNames,
    column_names_gp = gpar(fontsize = colNameFontSize),
    column_names_rot = -45,
    column_names_max_height = max_text_width(colNames),
    width = heatmapWidth,
    height = heatmapHeight,
    show_heatmap_legend = FALSE,
    top_annotation = colHmAnnot,
    left_annotation = rowHmLeftAnnot,
    use_raster = FALSE))
hmHeight <- component_height(hm, k = 'heatmap_body')
hmHeight <- convertHeight(hmHeight, 'inch', valueOnly = TRUE)  # will be 0 if heatmapHeight is NULL
# Add dataframe annotation of rows to main heatmap
hm <- hm + rowHmRightAnnot

# Add customized legends --------------------------------------------------------------------------
# Set heatmap legend
colorBarName <- 'Z-score'
heatmapLegend <- Legend(
    col_fun = heatmapColorFun,
    legend_height = unit(5, 'cm'),
    title = colorBarName,
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 11))
heatmapLegendHeight <- ComplexHeatmap:::height(heatmapLegend)
heatmapLegendHeight <- convertHeight(heatmapLegendHeight, 'inch', valueOnly = TRUE)
heatmapLegendList <- list(heatmapLegend)

# Set column (Group) annotation legend
colAnnotLegendList <- list()
colAnnotLegend <- Legend(
        labels = names(colAnnotColors),
        labels_gp = gpar(fontsize = 11),
        title = colAnnotName,
        title_gp = gpar(fontsize = 12),
        type = 'grid',
        legend_gp = gpar(fill = colAnnotColors))
# Adjust legend columns if needed
colAnnotLegendHeight <- ComplexHeatmap:::height(colAnnotLegend)
colAnnotLegendHeight <- convertHeight(colAnnotLegendHeight, 'inch', valueOnly = TRUE)
colAnnotLegendList <- list(colAnnotLegend)

# Set row annotation legend
rowAnnotLegendList <- list()
rowAnnotLegend <- Legend(
        labels = names(rowAnnotColors),
        labels_gp = gpar(fontsize = 11),
        title = rowAnnotName,
        title_gp = gpar(fontsize = 12),
        type = 'grid',
        legend_gp = gpar(fill = rowAnnotColors))
# Adjust legend columns if needed
rowAnnotLegendHeight <- ComplexHeatmap:::height(rowAnnotLegend)
rowAnnotLegendHeight <- convertHeight(rowAnnotLegendHeight, 'inch', valueOnly = TRUE)
rowAnnotLegendList <- list(rowAnnotLegend)

# Pack heatmap legend and annotation legends
legendList <- c(heatmapLegendList, colAnnotLegendList, rowAnnotLegendList)
legendPd <- packLegend(list = legendList,
                       direction = 'horizontal')


# Define the size of output device ----------------------------------------------------------------
# Adjust width
heatmapList <- drawComplexHeatmap(heatmap = hm, legends = legendPd, locus = legendLocus)
outWidth <- heatmapList$width %>% round(digits = 2)

# Adjust height
outHeight <- heatmapList$height %>% round(digits = 2)
outHeight <- max(outHeight, 4.31)  # the min height to show full heatmap legend
# Output files ------------------------------------------------------------------------------------
outPrefix <- "heatmap"
outPDF <- paste0(outPrefix, '.pdf')
outPNG <- paste0(outPrefix, '.png')

# Export PDF and PNG
graphics.off()  # close all previous devices
cairo_pdf(filename = outPDF, width = outWidth, height = outHeight)
cp <- dev.cur()  # record pdf device
png(filename = outPNG, type = 'cairo', width = outWidth, height = outHeight, res = 300, units = 'in')
dev.control('enable')
print(heatmapList$heatmap)
dev.copy(which = cp)  # copy plots from pdf to png
