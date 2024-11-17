
# 1. 处理数据 ------------------------------------------------------------------

# **1.1 feature2count -----------------------------------------------------------

#' Convert FeatureCount data into raw count matrix

#' @import tidyverse
#' @param featureCount A data frame containing FeatureCount data without the 1st row.
#' @return A matrix with raw counts where row names are gene IDs.
#' @examples
#' # Example usage
#' raw_count <- feature2count(featureCount)
feature2count <- function(featureCount) {
  rawCount <- featureCount %>% select(-c(Geneid:Length)) %>% as.matrix()
  rownames(rawCount) <- featureCount$Geneid
  return(rawCount)
}


# **1.2 Top variable -----------------------------------------------------------

#' Select top variable features
#'
#' This function identifies the top variable features from a dataset.
#' @importFrom matrixStats rowVars
#' @param data A matrix or data frame containing feature counts.
#' @param n An integer specifying the number of top variable features to select (default: 1000).
#' @return A subset of the original dataset with the top variable features.
#' @examples
#' # Example usage
#' top_features <- top_variable(data, n = 1000)
top_variable <- function(data, n = 1000) {
  rV <- matrixStats::rowVars(as.matrix(data))
  names(rV) <- rownames(data)
  rV_top <- sort(rV,decreasing = T)[1:n]
  return(data[names(rV_top),])
}


# 2. 标准化数据 ----------------------------------------------------------------

# **2.1 spike_nrom --------------------------------------------------------------------

#' Normalize data using size factors
#'
#' This function normalizes data using provided size factors.
#'
#' @param data A numeric matrix of raw counts.
#' @param sf A numeric vector of spike in reads.
#' @return A normalized numeric matrix.
#' @examples
#' # Example usage
#' normalized_data <- spike_norm(data, sf)
spike_norm <- function(data, sf) {
  sf <- sf / max(sf)
  data_norm <- t(t(data) / sf)
  return(data_norm)
}


# **2.2 ERCC loess --------------------------------------------------------------------

#' Normalize data using loess
#'
#' This function performs loess normalization on the dataset.
#' @importFrom affy normalize.loess
#' @param data A numeric matrix of raw counts.
#' @param subset A vector of indices to subset the data for normalization.
#' @return A loess-normalized numeric matrix.
#' @examples
#' # Example usage
#' normalized_data <- spike_loess(data, subset = c(1:10))
spike_loess <- function(data, subset) {
  library(affy)
  data_expr <- data[which(data <= 1)] <- 1.1
  data_loess <- normalize.loess(data_expr, subset = subset)
  return(data_loess)
}


# **2.3 RUVseq --------------------------------------------------------------------

#' Perform RUVseq normalization
#'
#' This function uses RUVseq to normalize RNA-seq data with spike-in controls.
#' @import RUVSeq
#' @param data A numeric matrix of raw counts.
#' @param groupsName A vector of group names.
#' @return A SeqExpressionSet object with normalized data.
#' @examples
#' # Example usage
#' normalized_set <- rUV_fun(data, groupsName = c("group1", "group2"))
rUV_fun <- function(data, groupsName) {
  library(RUVSeq)
  filter <- apply(data, 1, function(x) length(x[x > 5]) >= 2)
  filtered <- data[filter, ]
  spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]
  
  x <- as.factor(rep(groupsName))
  set <- newSeqExpressionSet(as.matrix(filtered),
                             phenoData = data.frame(x, row.names = colnames(filtered)))
  set <- betweenLaneNormalization(set, which = "upper")
  set1 <- RUVg(set, spikes, k = 1)
  return(set1)
}

# 3. 差异分析----------------------------------------------------------------

# **3.1 edgeR function ----------------------------------------------------------

#' Differential Expression Analysis with edgeR
#'
#' This function performs differential expression analysis using the edgeR package.
#' @import edgeR
#' @param data A matrix of raw counts.
#' @param select A vector of column indices to select from the data matrix (default: 1:4).
#' @param group A vector specifying group labels for the selected samples. Default: c("A", "A", "B", "B").
#' @param LogFC A numeric value for the log fold change threshold (default: 0.5).
#' @param spike A logical indicating if the spike in normalization has been done (default: TRUE).
#' @return A data frame containing the results of the differential expression analysis, including logFC, FDR, and group labels ("UP", "DN", "NC").
#' @examples
#' # Example usage
#' results <- edger_fun(data, select = 1:4, group = c("A", "A", "B", "B"))
edger_fun <- function(data, select = 1:4, group = NULL, LogFC = 0.5, spike = TRUE) {
  # 选择指定的列
  data <- data[, select]
  
  # 如果没有指定 group，则使用默认值
  if (is.null(group)) {
    group <- c("A", "A", "B", "B")
  }
  
  # 创建 DGEList 对象
  dge <- DGEList(counts = data, group = group)
  
  # 根据 spike 参数调整操作
  if (spike) {
    dge$samples$norm.factors <- rep(1, length(group))
    dge$samples$lib.size <- min(dge$samples$lib.size)
  }
  
  # 估计离散度并进行拟合
  dge <- estimateDisp(dge)
  fit <- glmFit(dge, design = model.matrix(~group))
  lrt <- glmLRT(fit)
  
  # 返回分析结果
  results <- topTags(lrt, n = Inf) %>% as.data.frame()
  results <- results %>% mutate(group = case_when(
    logFC >= LogFC & FDR <= 0.05 ~ "UP",
    logFC <= -LogFC & FDR <= 0.05 ~ "DN",
    TRUE ~ "NC"
  ))
  
  return(results)
}

# 4. 可视化----------------------------------------------------------------

# **4.1 火山图 ----------------------------------------------------------------

#' Generate a Volcano Plot
#'
#' This function creates a volcano plot for visualizing differential expression results.
#' @import tidyverse
#' @param x A data frame containing differential expression results.
#' @param xlim A numeric vector specifying x-axis limits (optional).
#' @param ylim A numeric vector specifying y-axis limits (optional).
#' @param sample An integer specifying the number of samples to randomly select for plotting (optional).
#' @param LogFC A numeric value for the log fold change threshold (default: 0.5).
#' @param cols A character vector of colors for the plot groups (default: c("#fc8d62", "#66c2a5")).
#' @return A ggplot object representing the volcano plot.
#' @examples
#' # Example usage
#' p <- plotV(results, xlim = c(-3, 3), ylim = c(0, 10))
#' print(p)
plotV <- function(x, xlim = NULL, ylim = NULL, sample = NULL, LogFC = 0.5, cols = c("#fc8d62", "#66c2a5")) {
  # 如果指定了 sample，则随机抽取样本行
  if (!is.null(sample)) {
    set.seed(123)
    x <- x %>% sample_n(sample)
  }
  
  data <- x %>% as_tibble() %>% mutate(group = case_when(
    logFC >= LogFC & FDR <= 0.05 ~ "UP",
    logFC <= -LogFC & FDR <= 0.05 ~ "DN",
    TRUE ~ "NC"
  ))
  data1 <- data
  data1$group <- factor(data1$group, levels = c("UP", "DN", "NC"))
  
  p <- ggplot(data1) +
    geom_vline(xintercept = c(-LogFC, LogFC), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_point(data = data1 %>% filter(group == "NC"),
               aes(logFC, -log10(FDR), size = abs(logFC)), color = "grey", alpha = 0.1) +
    geom_point(data = data1 %>% filter(group != "NC"),
               aes(logFC, -log10(FDR), size = abs(logFC), color = group), alpha = 0.5) +
    theme_bw(base_size = 8) +
    theme(legend.key.size = unit(3, "mm")) +
    scale_size_continuous(range = c(0.1, 2)) +
    scale_color_manual(values = cols)
  
  # 如果指定了 xlim 或 ylim，则添加 coord_cartesian()
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }
  
  return(p)
}

# **4.2 plotErrorBar ---------------------------------------------------------------

#' Create Error Bar Plot
#'
#' This function generates an error bar plot for summarizing and visualizing grouped data.
#' @import tidyverse
#' @param data A data frame containing the data to be plotted.
#' @param x A variable for the x-axis.
#' @param y A variable for the y-axis.
#' @param cols A character vector specifying colors for the groups (optional).
#' @return A ggplot object representing the error bar plot.
#' @examples
#' # Example usage
#' p <- plotErrorBar(data, x = "group", y = "value", cols = c("red", "blue"))
#' print(p)
plotErrorBar <- function(data, x, y, cols = NULL) {
  # 使用 !! 和 sym 动态指定列名
  x <- sym(x)
  y <- sym(y)
  
  # 计算数据摘要
  data_summary <- data %>%
    group_by(!!x) %>%
    summarise(
      Median = median(!!y, na.rm = TRUE),
      err = sd(!!y, na.rm = TRUE) / ((n() + 1)^0.5)
    )
  
  # 创建 ggplot 对象
  p <- ggplot(data_summary, aes(x = !!x, y = Median, color = !!x)) + 
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = Median - err, ymax = Median + err), width = 0.2, position = position_dodge(0.05)) +
    #geom_hline(yintercept = 0) +
    theme_classic(base_size = 8) +
    theme(legend.position = "none")
  
  # 如果提供了颜色，则添加颜色映射
  if (!is.null(cols)) {
    p <- p + scale_color_manual(values = cols)
  }
  
  return(p)
}


# 5. bed处理----------------------------------------------------------------

# **5.1 Convert to GR -------------------------------------------------------

#' Convert Data Frame to Genomic Ranges
#'
#' This function converts a data frame with genomic coordinates into a GRanges object.
#' @importFrom plyranges as_granges
#' @param data A data frame containing genomic data.
#' @param ID A character string specifying the column name for the ID field (default: "Geneid").
#' @return A GRanges object representing the genomic data.
#' @examples
#' # Example usage
#' gr <- CovertToGR(data, ID = "Geneid")
CovertToGR <- function(data, ID = "Geneid") {
  # 动态引用 ID 列
  id_col <- sym(ID)
  
  # 分割 ID 列并转换为数据框
  data_df <- data %>%
    pull(!!id_col) %>% # 提取指定列
    str_split("_", simplify = TRUE) %>%
    as_tibble()
  
  # 添加列名
  colnames(data_df) <- c("seqnames", "start", "end")
  
  # 转换为 GRanges 对象
  data_gr <- data_df %>%
    mutate(start = as.integer(start), end = as.integer(end)) %>%
    bind_cols(data) %>%
    plyranges::as_granges()
  
  return(data_gr)
}


# **5.2 FindOverlaps_gr --------------------------------------------------------------------

#' Find Overlaps Between Genomic Ranges
#'
#' This function identifies overlaps between two GRanges objects and returns a combined data frame.
#' @importFrom GenomicRanges findOverlaps
#' @param data1_gr A GRanges object representing the first genomic range.
#' @param data2_gr A GRanges object representing the second genomic range.
#' @return A data frame combining the information from both genomic ranges where overlaps occur.
#' @examples
#' # Example usage
#' overlaps <- findOverlaps_gr(gr1, gr2)
findOverlaps_gr <- function(data1_gr, data2_gr){
  overlaps <- findOverlaps(data1_gr, data2_gr)
  data1_res <- as_tibble(data1_gr)
  data2_res <- as_tibble(data2_gr)
  overlaps_res <- bind_cols(data1_res[queryHits(overlaps),], data2_res[subjectHits(overlaps),])
  return(overlaps_res)
}















