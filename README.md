
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Easyseq

Easyseq is an R package for RNA-seq data analysis and visualization. It
provides functions for data normalization, differential expression
analysis, and visualization.

## Installation

You can install the development version of Easyseq from source:

``` r
install.packages("path/to/Easyseq_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Usage

Here are examples of how to use key functions from Easyseq:

### Normalize data

``` r
library(Easyseq)

# Example usage of spike_norm
data <- matrix(runif(20, 1, 100), nrow = 5)
sf <- c(1.5, 2.0, 2.5, 3.0)
normalized_data <- spike_norm(data, sf)
print(normalized_data)
```

### Differential expression analysis

``` r
# Example usage of edger_fun
counts <- matrix(sample(1:100, 40, replace = TRUE), nrow = 10)
group <- c("A", "A", "B", "B")
results <- edger_fun(counts, group = group, LogFC = 0.5)
print(results)
```

### Generate a volcano plot

``` r
# Example usage of plotV
library(ggplot2)
results <- data.frame(
  logFC = runif(100, -3, 3),
  FDR = runif(100, 0, 0.1)
)
volcano_plot <- plotV(results, LogFC = 0.5)
print(volcano_plot)
```

## Contributing

If you’d like to contribute to Easyseq, feel free to submit issues or
pull requests on [GitHub](https://github.com/yourusername/Easyseq).

## License

This package is licensed under the GPL-3 license.
