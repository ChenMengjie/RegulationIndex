## RegulationIndex

To better distinguish between the gene regulatory states, we introduce a gene-specific statistic, ‘k-proportion’. This is the percentage of cells with gene expression less than a defined integer value, k, which is determined based on the mean gene expression count. Consider a situation where certain genes exhibit high expression levels in a limited subset of cells. This extreme expression skews the overall average upward, however the k-proportions (a measure reflecting the portion of cells with considerably lower expression) remain high. Essentially, under equivalent mean expressions, regulatory genes display significantly higher values of k-proportions compared to their homeostatic counterparts.

**RegulationIndex** is an implementation for k-proportion inflation test. 




## Installation

**RegulationIndex** can be installed from github directly as follows:

```r
devtools::install_github("ChenMengjie/RegulationIndex")
```

Read in an example data set. 
Zheng et al data was downloaded from the 10X website. Only cells that labeled as CD34+ were selected for analysis. 

```r       
library(RegulationIndex)
load(CD34)
counts <- counts[, -which(duplicated(colnames(counts)))]
testdata <- t(counts)
```

There are 277 cells with 31482 genes. We performed UMAP on log-transformed data. We used k-means to separate cells into 3 subtypes using their UMAP coordinates. The subtype information is saved in `final_cluster'. 
```r
dim(testdata)
table(final_cluster)
```
## **RegulationIndex** function 1: performing k-proportion inflation tests 

The main function is `Kprop_inflation_test`, with following arguments:
- `dat` Input data matrix for gene \times cell, preferably UMI counts. All our experiments are performed on UMI counts. The performance is not accessed for processed data.
- `grouping` The grouping IDs of cells. If the `grouping` is NULL, all cells will be assessed at once. The grouping can be cell types, individuals or any phenotype. If the `grouping` is supplied, cells will be evaluated by groups. If you have data with complex nested structures, like multiple cell types from multiple donors, you need to run `Kprop_inflation_test` multiple times. Now only one layer of grouping structure is supported. See our examples. 
- `min.cell.num.group` Minimum cell number required for each group. The default is 25. Groups with less cells will be excluded from the analysis.
- `min.depth.group` Minimum per-cell average depth (number of reads) required for each group. The default is 4000, set empirically for 10X data. Group with lower depths will be excluded from the analysis. 
- `cell.depth.ranges` A two-value vector of range of total counts per cell. The default is c(2500, 60000), set empirically for 10X data. The recommended per cell sequencing depth for 10x Genomics projects is between 30,000 and 70,000 reads. If your sequencing depth is substantially higher or lower than common practice, please adjust `min.depth.group` and `cell.depth.ranges` accordingly. 
- `kprop.cutoff` A pre-specified cut-off on K props. Only k >= k_cutoff will be used for estimating the dispersion of NB. The default is 1, meaning not using 0-prop for estimation.
- `fdr.control.method` Method to adjust for multiple comparisons. The default is "BH". See function p.adjust() for more options.
- `qvalue.cutoff` A cut-off on q-value adjusted by multiple comparisons, to call significantly inflated genes. The default is 0.05.
- `genemean.cutoff` A cut-off on mean UMI per cell, to report in the list of inflated genes. The default is 1. 
- `filter.gene` Whether reporting IncRNA, mtRNA, or ribosomal RNA in the list of inflated genes. The default is TRUE.

`Kprop_inflation_test` will return a list of summary of k-proportion test for each group.
- groupID : the ID of the group. The value is NA if there is no grouping information supplied.
- NBdispersion  :  estimated dispersion parameter
- Zscores : a matrix summarizing gene mean, zscore, p-value and q-value for all genes in `dat`. NAs are reported for genes with mean = 0.
- Inflated : a matrix of selected genes that pass `qvalue.cutoff` and `genemean.cutoff`.


```r       
cluster.zscore.summary <- Kprop_inflation_test(testdata, grouping = final_cluster, min.depth.group = 500, cell.depth.ranges = c(500, 10000), genemean.cutoff = 0.5)
```

If the `grouping` is not provided, the inflation will be evaluated treating all cells as a single group. 

```r       
cluster.zscore.summary.as.one <- Kprop_inflation_test(testdata, grouping = NULL, min.depth.group = 500, cell.depth.ranges = c(500, 10000), genemean.cutoff = 0.5)
```

## **RegulationIndex** function 2: Visualizing k-proportion inflation tests using a wave plot

Consider a situation where certain genes exhibit high expression levels in a limited subset of cells. This extreme expression skews the overall average upward, however the k-proportions (a measure reflecting the portion of cells with considerably lower expression) remain high. Essentially, under equivalent mean expressions, regulatory genes display significantly higher values of k-proportions compared to their homeostatic counterparts. This relationship between k-proportions and the mean expression can be illustrated on a 'wave plot'. In this visualization, regulatory genes surface as anomalies, appearing as distinct outliers above the general trend, akin to water droplets cast airborne from a wave. By isolating these 'droplets', we can identify potential regulatory genes, offering new avenues for characterizing gene regulation directly from single-cell data.

The function `waveplot` will visualize the k-proportion inflation test summary for one group, output from `Kprop_inflation_test`, using a wave plot. By default, only coding RNA will be plotted and top genes will be labeled. Options are available to add a 95% confidence interval based on the negative binomial distribution, and to add expected k-proportions calculated from a Poisson distribution.

The following commands draw wave plots for three different subpopulations for CD34+ cells. 
```r       
wave_one <- waveplot(cluster.zscore.summary[[1]], show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#CC8D1A")
```

```r       
wave_two <- waveplot(cluster.zscore.summary[[2]], show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#8CBEB")
```

```r       
wave_three <- waveplot(cluster.zscore.summary[[3]], show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#F06060")
```

The following commands draw a wave plot for  CD34+ cells as a whole population. 

```r       
wave_all <- waveplot(cluster.zscore.summary.as.one, show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#F06060")
```


### Author

**Mengjie Chen** (U Chicago)

Bug report, comments or questions please send to mengjiechen@uchicago.edu.
