Dunlop growth rates - habitat
================
Cassandra Wattenburger
11/10/2020

**Goal of this script:** Investigate the effects of soil habitat on
in-situ growth rates

# Import libraries

``` r
library(tidyverse)
library(cowplot)
library(data.table)
library(phyloseq)
library(lmerTest)

sessionInfo()
```

    ## R version 3.6.3 (2020-02-29)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] lmerTest_3.1-3    lme4_1.1-23       Matrix_1.2-18     phyloseq_1.30.0  
    ##  [5] data.table_1.13.0 cowplot_1.1.0     forcats_0.5.0     stringr_1.4.0    
    ##  [9] dplyr_1.0.2       purrr_0.3.4       readr_1.4.0       tidyr_1.1.2      
    ## [13] tibble_3.0.3      ggplot2_3.3.2     tidyverse_1.3.0  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-147        fs_1.5.0            lubridate_1.7.9    
    ##  [4] httr_1.4.2          numDeriv_2016.8-1.1 tools_3.6.3        
    ##  [7] backports_1.1.10    R6_2.4.1            vegan_2.5-6        
    ## [10] mgcv_1.8-31         DBI_1.1.0           BiocGenerics_0.32.0
    ## [13] colorspace_1.4-1    permute_0.9-5       ade4_1.7-15        
    ## [16] withr_2.3.0         tidyselect_1.1.0    compiler_3.6.3     
    ## [19] cli_2.0.2           rvest_0.3.6         Biobase_2.46.0     
    ## [22] xml2_1.3.2          scales_1.1.1        digest_0.6.25      
    ## [25] minqa_1.2.4         rmarkdown_2.4       XVector_0.26.0     
    ## [28] pkgconfig_2.0.3     htmltools_0.5.0     dbplyr_1.4.4       
    ## [31] rlang_0.4.8         readxl_1.3.1        rstudioapi_0.11    
    ## [34] generics_0.0.2      jsonlite_1.7.1      magrittr_1.5       
    ## [37] biomformat_1.14.0   Rcpp_1.0.5          munsell_0.5.0      
    ## [40] S4Vectors_0.24.4    Rhdf5lib_1.8.0      fansi_0.4.1        
    ## [43] ape_5.4-1           lifecycle_0.2.0     stringi_1.5.3      
    ## [46] yaml_2.2.1          MASS_7.3-51.6       zlibbioc_1.32.0    
    ## [49] rhdf5_2.30.1        plyr_1.8.6          grid_3.6.3         
    ## [52] blob_1.2.1          parallel_3.6.3      crayon_1.3.4       
    ## [55] lattice_0.20-41     Biostrings_2.54.0   haven_2.3.1        
    ## [58] splines_3.6.3       multtest_2.42.0     hms_0.5.3          
    ## [61] knitr_1.30          pillar_1.4.6        igraph_1.2.6       
    ## [64] boot_1.3-25         reshape2_1.4.4      codetools_0.2-16   
    ## [67] stats4_3.6.3        reprex_0.3.0        glue_1.4.2         
    ## [70] evaluate_0.14       modelr_0.1.8        nloptr_1.2.2.2     
    ## [73] vctrs_0.3.4         foreach_1.5.0       cellranger_1.1.0   
    ## [76] gtable_0.3.0        assertthat_0.2.1    xfun_0.18          
    ## [79] broom_0.7.1         survival_3.1-12     iterators_1.0.12   
    ## [82] IRanges_2.20.2      cluster_2.1.0       statmod_1.4.34     
    ## [85] ellipsis_0.3.1

``` r
rm(list=ls())
```

Import estimated growth rates:

``` r
gr.est = readRDS("rdata.files/gr.final.cleaned.rds")
```

# Growth properties

Average across ASVs to get community level data. ASVs within same
community are not independent of one-another, which would screw up
statistics. Unfortunately, singleton ASV estimates prevents me frrom
using a mixed-linear effects model with ASV as a random effect, which
would be ideal, so we must average instead.

``` r
gr.est.avg = gr.est %>%
  group_by(Soil, Amendment, Replicate) %>% # average across ASVs
  summarize(avg_start = mean(Start),
            avg_end = mean(End),
            avg_length = mean(Length),
            avg_k = mean(k))
```

### Start of growth

``` r
gr.est.avg %>%
  ggplot(aes(x=Soil, y=avg_start, color=Amendment)) +
  geom_boxplot() +
  labs(title="Start of growth", y="Day") +
  theme_test()
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

**Statistics**

Try
ANOVA.

``` r
strt.aov = aov(avg_start ~ Soil + Amendment + Soil:Amendment, data=gr.est.avg)
hist(resid(strt.aov))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Model residuals are not normal. Try log
transform.

``` r
strt.aov2 = aov(log(avg_start) ~ Soil + Amendment + Soil:Amendment, data=gr.est.avg)
hist(resid(strt.aov2))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
plot(predict(strt.aov2), resid(strt.aov2))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
summary(strt.aov2)
```

    ##                Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Soil            1 0.0440  0.0440   3.491   0.0987 .  
    ## Amendment       1 0.9104  0.9104  72.185 2.82e-05 ***
    ## Soil:Amendment  1 0.1338  0.1338  10.610   0.0116 *  
    ## Residuals       8 0.1009  0.0126                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Start of growth is affected by soil and depends on the amendment.

Post hoc:

Only interested in certain contrasts, ie doesn’t make sense to compare
cropped water to successional C-amended.

Using Welch’s t-tests (fewer assumptions) with Benjamini-Hochberg
correction.

``` r
# Contrasts
start.wtest1 = t.test(log(avg_start) ~ Amendment, data=gr.est.avg[gr.est.avg$Soil=="C3",], var.equal = FALSE)
start.wtest2 = t.test(log(avg_start) ~ Amendment, data=gr.est.avg[gr.est.avg$Soil=="S17",], var.equal = FALSE)
start.wtest3 = t.test(log(avg_start) ~ Soil, data=gr.est.avg[gr.est.avg$Amendment=="Y",], var.equal = FALSE)
start.wtest4 = t.test(log(avg_start) ~ Soil, data=gr.est.avg[gr.est.avg$Amendment=="N",], var.equal = FALSE)

# Adjust p-values
start.pvals = c(start.wtest1$p.value, start.wtest2$p.value, start.wtest3$p.value, start.wtest4$p.value)
start.adjpvals = p.adjust(start.pvals, method="BH", n=length(start.pvals))

# Results
start.results = data.frame(contrast = c("Cropped, C vs W", "Successional, C vs W", "C-amended, cropped vs succ.", "Water, cropped vs succ."), pvals = start.adjpvals)
start.results
```

    ##                      contrast      pvals
    ## 1             Cropped, C vs W 0.01021037
    ## 2        Successional, C vs W 0.09463981
    ## 3 C-amended, cropped vs succ. 0.09463981
    ## 4     Water, cropped vs succ. 0.10136663

### End of growth window

``` r
gr.est.avg %>%
  ggplot(aes(x=Soil, y=avg_end, color=Amendment)) +
  geom_boxplot() +
  labs(title="End of growth", y="Day") +
  theme_test()
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

**Statistics**

``` r
end.aov = aov(avg_end ~ Soil + Amendment + Soil:Amendment, data=gr.est.avg)
hist(resid(end.aov))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
plot(predict(end.aov), resid(end.aov))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
summary(end.aov)
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)  
    ## Soil            1  3.683   3.683   9.550 0.0149 *
    ## Amendment       1  1.645   1.645   4.266 0.0727 .
    ## Soil:Amendment  1  0.801   0.801   2.077 0.1875  
    ## Residuals       8  3.085   0.386                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

End of growth is affected by soil and depends on the amendment.

Post hoc:

Only interested in certain contrasts, ie doesn’t make sense to compare
cropped water to successional C-amended.

Using Welch’s t-tests (fewer assumptions) with Benjamini-Hochberg
correction.

``` r
# Contrasts
end.wtest1 = t.test(avg_end ~ Amendment, data=gr.est.avg[gr.est.avg$Soil=="C3",], var.equal = FALSE)
end.wtest2 = t.test(avg_end ~ Amendment, data=gr.est.avg[gr.est.avg$Soil=="S17",], var.equal = FALSE)
end.wtest3 = t.test(avg_end ~ Soil, data=gr.est.avg[gr.est.avg$Amendment=="Y",], var.equal = FALSE)
end.wtest4 = t.test(avg_end ~ Soil, data=gr.est.avg[gr.est.avg$Amendment=="N",], var.equal = FALSE)

# Adjust p-values
end.pvals = c(end.wtest1$p.value, end.wtest2$p.value, end.wtest3$p.value, end.wtest4$p.value)
end.adjpvals = p.adjust(end.pvals, method="BH", n=length(end.pvals))

# Results
end.results = data.frame(contrast = c("Cropped, C vs W", "Successional, C vs W", "C-amended, cropped vs succ.", "Water, cropped vs succ."), pvals = end.adjpvals)
end.results
```

    ##                      contrast     pvals
    ## 1             Cropped, C vs W 0.1765747
    ## 2        Successional, C vs W 0.7194618
    ## 3 C-amended, cropped vs succ. 0.4649738
    ## 4     Water, cropped vs succ. 0.1765747

No significance.

### Length of growth window

``` r
gr.est.avg %>%
  ggplot(aes(x=Soil, y=avg_length, color=Amendment)) +
  geom_boxplot() +
  labs(title="Length of growth", y="Day") +
  theme_test()
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

**Statistics**

``` r
length.lm = lm(avg_length ~ Soil + Amendment + Soil:Amendment, data=gr.est.avg)
hist(resid(length.lm))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Linear model fit is not normal. Try log
transform.

``` r
length.lm2 = lm(log(avg_length) ~ Soil + Amendment + Soil:Amendment, data=gr.est.avg)
hist(resid(length.lm2))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
plot(predict(length.lm2), resid(length.lm2))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

``` r
summary(length.lm2)
```

    ## 
    ## Call:
    ## lm(formula = log(avg_length) ~ Soil + Amendment + Soil:Amendment, 
    ##     data = gr.est.avg)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.16575 -0.06974 -0.01896  0.06178  0.17176 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         1.83972    0.07691  23.920 9.94e-09 ***
    ## SoilS17            -0.23148    0.10877  -2.128    0.066 .  
    ## AmendmentN         -0.12181    0.10877  -1.120    0.295    
    ## SoilS17:AmendmentN -0.04011    0.15382  -0.261    0.801    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1332 on 8 degrees of freedom
    ## Multiple R-squared:  0.6391, Adjusted R-squared:  0.5037 
    ## F-statistic: 4.722 on 3 and 8 DF,  p-value: 0.03518

Length of growth in not significant.

### Specific growth rates

``` r
gr.est.avg %>%
  ggplot(aes(x=Soil, y=log(avg_k), color=Amendment)) +
  geom_boxplot() +
  labs(title="Specific growth rate", y="Specific growth rate, ln (per Day)") +
  theme_test()
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

**Statistics**

``` r
k.aov = lm(avg_k ~ Soil + Amendment + Soil:Amendment, data=gr.est.avg)
hist(resid(k.aov))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Try log
transform.

``` r
k.aov2 = aov(log(avg_k) ~ Soil + Amendment + Soil:Amendment, data=gr.est.avg)
hist(resid(k.aov2))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
plot(predict(k.aov2), resid(k.aov2))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

``` r
summary(k.aov2)
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)
    ## Soil            1 0.0067 0.00672   0.051  0.827
    ## Amendment       1 0.1532 0.15320   1.156  0.314
    ## Soil:Amendment  1 0.0706 0.07057   0.533  0.486
    ## Residuals       8 1.0599 0.13249

No significance found.

### Create table of results:

``` r
gr.est %>%
  group_by(Soil, Amendment) %>%
  summarize(avg_k = mean(k), sd_k = sd(k), avg_start = mean(Start), sd_start = sd(Start), avg_end = mean(End), sd_end = sd(End), avg_len = mean(Length), sd_len = sd(Length))
```

    ## # A tibble: 4 x 10
    ## # Groups:   Soil [2]
    ##   Soil  Amendment avg_k  sd_k avg_start sd_start avg_end sd_end avg_len sd_len
    ##   <fct> <fct>     <dbl> <dbl>     <dbl>    <dbl>   <dbl>  <dbl>   <dbl>  <dbl>
    ## 1 C3    Y         0.268 0.159      1.75     1.08    7.31   2.12    6.57   2.35
    ## 2 C3    N         0.235 0.123      3.81     1.83    8.43   2.11    5.62   1.66
    ## 3 S17   Y         0.296 0.178      2.49     2.00    6.76   2.41    5.27   2.40
    ## 4 S17   N         0.219 0.190      3.53     3.10    7.04   2.94    4.51   1.70

### End of growth

Generation time:

``` r
gr.est %>% summarize(min_g=min(g), max_g=max(g), mean_g=mean(g))
```

    ##       min_g   max_g   mean_g
    ## 1 0.7039976 63.5031 5.043394

# Overlapping estimates between treatments

``` r
# Isolate treatment ASV names
gr.C3.asv = gr.est %>% filter(Soil=="C3") %>% distinct(ASV)
gr.S17.asv = gr.est %>% filter(Soil=="S17") %>% distinct(ASV)
gr.C3y.asv = gr.est %>% filter(Soil=="C3", Amendment=="Y") %>% distinct(ASV)
gr.C3n.asv = gr.est %>% filter(Soil=="C3", Amendment=="N") %>% distinct(ASV)
gr.S17y.asv = gr.est %>% filter(Soil=="S17", Amendment=="Y") %>% distinct(ASV)
gr.S17n.asv = gr.est %>% filter(Soil=="S17", Amendment=="N") %>% distinct(ASV)

# Overlap
inner_join(gr.C3.asv, gr.S17.asv) %>% nrow()
```

    ## [1] 88

``` r
inner_join(gr.C3y.asv, gr.C3n.asv) %>% nrow()
```

    ## [1] 77

``` r
inner_join(gr.S17y.asv, gr.S17n.asv) %>% nrow()
```

    ## [1] 12

# Growth rate distributions

Average across replicates: makes for cleaner presentation.

``` r
gr.estavg = gr.est %>% 
  group_by(Soil, Amendment, ASV, Domain, Phylum, Class, Order, Family, Genus) %>% 
  summarize(avgk = mean(k), avgStart = mean(Start), avgEnd = mean(End), avgLength = mean(Length))
```

Kernel density graphs

See:
<https://www.statsmodels.org/dev/examples/notebooks/generated/kernel_density.html>

``` r
gr.estavg %>% 
  filter(Amendment=="N") %>%
  ggplot(aes(x=log(avgk), color=Soil, fill=Soil)) +
  geom_density(alpha=0.5) +
  labs(title="Only water", x="log specific growth rate", y="Kernel Density") +
  theme_test() +
  theme(#legend.position="none",
        title=element_text(size=12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  scale_color_manual(labels = c("Cropped", "Successional"), values=c("#F8766D","#00BFC4")) +
  scale_fill_manual(labels = c("Cropped", "Successional"), values=c("#F8766D","#00BFC4")) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(-5,0))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
gr.estavg %>% 
  filter(Amendment=="Y") %>%
  ggplot(aes(x=log(avgk), color=Soil, fill=Soil)) +
  geom_density(alpha=0.5) +
  labs(title="C-amended", x="log specific growth rate", y="Kernel Density") +
  theme_test() +
  theme(#legend.position="none",
        title=element_text(size=12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  scale_color_manual(labels = c("Cropped", "Successional"), values=c("#F8766D","#00BFC4")) +
  scale_fill_manual(labels = c("Cropped", "Successional"), values=c("#F8766D","#00BFC4")) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(-5,0))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

``` r
gr.estavg %>% 
  filter(Soil=="C3") %>%
  ggplot(aes(x=log(avgk), color=Amendment, fill=Amendment)) +
  geom_density(alpha=0.5) +
  labs(title="Cropped", x="log specific growth rate", y="Kernel Density") +
  theme_test() +
  theme(#legend.position="none",
        title=element_text(size=12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  scale_color_manual(labels = c("C-added", "Water"), values=c("#E78E43", "#423EC5")) +
  scale_fill_manual(labels = c("C-added", "Water"), values=c("#E78E43", "#423EC5")) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(-5,0))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-21-3.png)<!-- -->

``` r
gr.estavg %>% 
  filter(Soil=="S17") %>%
  ggplot(aes(x=log(avgk), color=Amendment, fill=Amendment)) +
  geom_density(alpha=0.5) +
  labs(title="Successional", x="log specific growth rate", y="Kernel Density") +
  theme_test() +
  theme(#legend.position="none",
        title=element_text(size=12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  scale_color_manual(labels = c("C-added", "Water"), values=c("#E78E43", "#423EC5")) +
  scale_fill_manual(labels = c("C-added", "Water"), values=c("#E78E43", "#423EC5")) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(-5,0))
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-21-4.png)<!-- -->

# Growth rate distribution binning

Hypothesis: growth rate distributions will differ by soil and by
amendment.

In order to conduct tests on the distirbution of the data, I need to bin
it discretely. The bins will need to be equal width and have the same
start and end points across treatments so that they are comparable, and
the number of bins must be chosen to simplify the data without losing
too much of the “shape” of the distribution (see smoothed histograms
above for reference).

Ranging from 5 to 10 bins.

Overall:

``` r
for (i in 5:10) {
  print(gr.estavg %>%
    ggplot(aes(x=log(avgk))) +
    geom_histogram(bins=i, color="black", fill="gray") +
    theme_test() +
    labs(title=paste("Overall,", i, "bins", sep=" ")))
}
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-22-2.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-22-3.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-22-4.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-22-5.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-22-6.png)<!-- -->

Treatments:

``` r
for (s in c("C3", "S17")) {
  for (a in c("Y", "N")) {
    for (i in 5:10) {
      print(gr.estavg %>%
        filter(Soil==s, Amendment==a) %>%
        ggplot(aes(x=log(avgk))) +
        geom_histogram(bins=i, color="black", fill="gray") +
        theme_test() +
        labs(title=paste(s, a, i, "bins")))
    }
  }
}
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-3.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-4.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-5.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-6.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-7.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-8.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-9.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-10.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-11.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-12.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-13.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-14.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-15.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-16.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-17.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-18.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-19.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-20.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-21.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-22.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-23.png)<!-- -->![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-23-24.png)<!-- -->

It looks like 5 bins is too few, 6 or 7 seem to more adequately capture
the distributions. I’ll choose 7 bins.

There need to be the same number and width of bins in each dataset.

Caculate min, max, and range of specific growth rate values in whole
dataset:

``` r
print(paste("Min:", min(log(gr.estavg$avgk))))
```

    ## [1] "Min: -4.51760171995445"

``` r
print(paste("Max:", max(log(gr.estavg$avgk))))
```

    ## [1] "Max: -0.0155326196739766"

``` r
print(paste("Range:", max(log(gr.estavg$avgk))-min(log(gr.estavg$avgk))))
```

    ## [1] "Range: 4.50206910028048"

Calculate bin width based on range:

``` r
width = (max(log(gr.estavg$avgk))-min(log(gr.estavg$avgk)))/7
print(paste("Bin width:", width))
```

    ## [1] "Bin width: 0.643152728611497"

``` r
kmin=min(log(gr.estavg$avgk))
print(paste("Bin 1 boundaries: less than", kmin+width))
```

    ## [1] "Bin 1 boundaries: less than -3.87444899134296"

``` r
print(paste("Bin 2 boundaries:", kmin+width, "to",  kmin+width*2))
```

    ## [1] "Bin 2 boundaries: -3.87444899134296 to -3.23129626273146"

``` r
print(paste("Bin 3 boundaries:", kmin+width*2, "to",  kmin+width*3))
```

    ## [1] "Bin 3 boundaries: -3.23129626273146 to -2.58814353411996"

``` r
print(paste("Bin 4 boundaries:", kmin+width*3, "to",  kmin+width*4))
```

    ## [1] "Bin 4 boundaries: -2.58814353411996 to -1.94499080550847"

``` r
print(paste("Bin 5 boundaries:", kmin+width*4, "to",  kmin+width*5))
```

    ## [1] "Bin 5 boundaries: -1.94499080550847 to -1.30183807689697"

``` r
print(paste("Bin 6 boundaries:", kmin+width*5, "to",  kmin+width*6))
```

    ## [1] "Bin 6 boundaries: -1.30183807689697 to -0.658685348285473"

``` r
print(paste("Bin 7 boundaries:", kmin+width*6, "and greater"))
```

    ## [1] "Bin 7 boundaries: -0.658685348285473 and greater"

``` r
# Bin data based on boundaries chosen and log k
gr.est.bins = gr.est %>% mutate(bin = ifelse(log(k) < kmin+width, 1, ifelse(log(k) < kmin+width*2, 2, ifelse(log(k) < kmin+width*3, 3, ifelse(log(k) < kmin+width*4, 4, ifelse(log(k) < kmin+width*5, 5, ifelse(log(k) < kmin+width*6, 6, 7)))))))

# Calculate bin information (proportion, frequency)
gr.bins = data.frame()
for (s in as.character(unique(gr.est.bins$Soil))) { 
  for (a in as.character(unique(gr.est.bins$Amendment))) {
    for (r in as.character(unique(gr.est.bins$Replicate))) {
      total = nrow(gr.est.bins[gr.est.bins$Soil==s & gr.est.bins$Amendment==a & gr.est.bins$Replicate==r,])
      for (b in as.character(unique(gr.est.bins$bin))) {
        bin = gr.est.bins[gr.est.bins$Soil==s & gr.est.bins$Amendment==a & gr.est.bins$Replicate==r & gr.est.bins$bin==b,]
        binfreq = nrow(bin)
        binratio = binfreq/total
        thisrow = data.frame(s, a, r, b, binfreq, binratio)
        gr.bins = rbind(gr.bins, thisrow)   
      }
    }
  }
}
colnames(gr.bins) = c("Soil","Amendment","Replicate","bin","frequency","proportion")
gr.bins$bin = factor(gr.bins$bin, levels=c(1,2,3,4,5,6,7))
```

Plot bins:

``` r
ggplot(gr.bins) +
  geom_boxplot(aes(x=bin, y=proportion, color=Soil)) +
  facet_wrap(~Amendment) +
  theme_test() +
  labs(x="Slow to fast", y="Proportion of total taxa")
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
ggplot(gr.bins) +
  geom_boxplot(aes(x=bin, y=frequency, color=Soil)) +
  facet_wrap(~Amendment) +
  theme_test() +
  labs(x="Slow to fast", y="Frequency")
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-27-2.png)<!-- -->

``` r
ggplot(gr.bins) +
  geom_boxplot(aes(x=bin, y=proportion, color=Amendment)) +
  facet_wrap(~Soil) +
  scale_color_manual(values=c("#E78E43", "#423EC5")) +
  theme_test() +
  labs(x="Slow to fast", y="Proportion of total taxa")
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-27-3.png)<!-- -->

``` r
ggplot(gr.bins) +
  geom_boxplot(aes(x=bin, y=frequency, color=Amendment)) +
  facet_wrap(~Soil) +
  scale_color_manual(values=c("#E78E43", "#423EC5")) +
  theme_test() +
  labs(x="Slow to fast", y="Frequency")
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-27-4.png)<!-- -->

**Statistics**

First, I want to know if the distribution of growth rates differs by
habitat or amendment. I’ll use a contingency table and Fisher’s exact
test.

See:
<http://www.biostathandbook.com/fishers.html>

``` r
# Average replicates and round (Frischer's needs integers, not continuous numbers)
gr.bins.avg = gr.bins %>% group_by(Soil, Amendment, bin) %>% summarize(avgfreq = mean(frequency)) %>% mutate(avgfreq.rnd = round(avgfreq))

# Make contingency tables
conting.y = dcast(gr.bins.avg[gr.bins.avg$Amendment=="Y",], Soil~bin, value.var="avgfreq.rnd")
conting.y = conting.y[,-1] %>% as.matrix()

conting.n = dcast(gr.bins.avg[gr.bins.avg$Amendment=="N",], Soil~bin, value.var="avgfreq.rnd")
conting.n = conting.n[,-1] %>% as.matrix()

conting.C3 = dcast(gr.bins.avg[gr.bins.avg$Soil=="C3",], Amendment~bin, value.var="avgfreq.rnd")
conting.C3 = conting.C3[,-1] %>% as.matrix()

conting.S17 = dcast(gr.bins.avg[gr.bins.avg$Soil=="S17",], Amendment~bin, value.var="avgfreq.rnd")
conting.S17 = conting.S17[,-1] %>% as.matrix()

# Run Fisher's exact tests
set.seed(1)
fish.testy = fisher.test(conting.y, simulate.p.value = TRUE)
fish.testn = fisher.test(conting.n, simulate.p.value = TRUE)
fish.testC3 = fisher.test(conting.C3, simulate.p.value = TRUE)
fish.testS17 = fisher.test(conting.S17, simulate.p.value = TRUE)

# Multiple test correction
fish.pvals = c(fish.testn$p.value, fish.testy$p.value, fish.testC3$p.value, fish.testS17$p.value)
fish.padj = p.adjust(fish.pvals, method="holm", n=length(fish.pvals))
print(paste("Unamended cropped vs Succ., adj. p-value:", fish.padj[1]))
```

    ## [1] "Unamended cropped vs Succ., adj. p-value: 0.00399800099950025"

``` r
print(paste("Amended cropped vs Succ., adj. p-value:", fish.padj[2]))
```

    ## [1] "Amended cropped vs Succ., adj. p-value: 0.547226386806597"

``` r
print(paste("Cropped amended vs unamended, adj. p-value:", fish.padj[3]))
```

    ## [1] "Cropped amended vs unamended, adj. p-value: 0.00899550224887556"

``` r
print(paste("Succ. amended vs unamended, adj. p-value:", fish.padj[4]))
```

    ## [1] "Succ. amended vs unamended, adj. p-value: 0.254872563718141"

# Growth rate shift up/down

Among ASVs that were estimated in both cropped and successional soil
treatments, do they reliably shift up or down in growth rate?

Hypotheses:

1.  ASVs found in water and C-amended treatments will grow faster in the
    C-amended treatment in both soils because of higher resource
    availability.

2.  ASVs found in both cropped and successional soils overall will grow
    faster in the cropped soil because of greater disturbances in the
    cropped soil.

<!-- end list -->

``` r
# Isolate treatment data
gr.C3.avgk = gr.est %>% 
  filter(Soil=="C3") %>% 
  group_by(Soil, ASV) %>% 
  summarize(C3.avgk=mean(k))

gr.S17.avgk = gr.est %>%
  filter(Soil=="S17") %>% 
  group_by(Soil, ASV) %>% 
  summarize(S17.avgk=mean(k))

gr.C3y.avgk = gr.est %>% 
  filter(Soil=="C3", Amendment=="Y") %>% 
  group_by(Soil, Amendment, ASV) %>% 
  summarize(C3y.avgk=mean(k))

gr.C3n.avgk = gr.est %>% 
  filter(Soil=="C3", Amendment=="N") %>% 
  group_by(Soil, Amendment, ASV) %>% 
  summarize(C3n.avgk=mean(k))

gr.S17y.avgk = gr.est %>% 
  filter(Soil=="S17", Amendment=="Y") %>% 
  group_by(Soil, Amendment, ASV) %>% 
  summarize(S17y.avgk=mean(k))

gr.S17n.avgk = gr.est %>% 
  filter(Soil=="S17", Amendment=="N") %>% 
  group_by(Soil, Amendment, ASV) %>% 
  summarize(S17n.avgk=mean(k))
```

Visualize:

``` r
# Overall cropped vs succesional
inner_join(gr.C3.avgk, gr.S17.avgk, by="ASV") %>%
  ggplot(aes(x=C3.avgk, y=S17.avgk)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="Cropped, log specific growth rate", y="Successional, log specific growth rate") +
  theme_test()
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

``` r
# Cropped C-amended vs water
inner_join(gr.C3y.avgk, gr.C3n.avgk, by="ASV") %>% 
  ggplot(aes(x=C3y.avgk, y=C3n.avgk)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="Cropped C-amended, log specific growth rate", y="Cropped water-amended, log specific growth rate") +
  theme_test()
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-30-2.png)<!-- -->

``` r
# Succesional C-amended vs water
inner_join(gr.S17y.avgk, gr.S17n.avgk, by="ASV") %>% 
  ggplot(aes(x=S17y.avgk, y=S17n.avgk)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="Successional C-amended, log specific growth rate", y="Successional water-amended, log specific growth rate") +
  theme_test()
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-30-3.png)<!-- -->

Doesn’t look promising.

**Statistcs**

Try simple linear regression first:

``` r
# Test overall cropped vs successional using linear regression
gr.C3S17.avgk = inner_join(gr.C3.avgk, gr.S17.avgk, by="ASV")
gr.C3S17.avgk.lm = lm(C3.avgk ~ S17.avgk, data=gr.C3S17.avgk)

# Evaluate assumptions of linear regression
hist(gr.C3S17.avgk.lm$residuals) # evaluate assumption of normality
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

``` r
plot(predict(gr.C3S17.avgk.lm), resid(gr.C3S17.avgk.lm)) # evaluate assumption of homoskedasticity
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-31-2.png)<!-- -->

``` r
# P value
summary(gr.C3S17.avgk.lm)$coefficients[2,4]
```

    ## [1] 0.06631589

``` r
# Cropped C-amended vs water
gr.C3yn.avgk = inner_join(gr.C3y.avgk, gr.C3n.avgk, by="ASV")
gr.C3yn.avgk.lm = lm(C3y.avgk ~ C3n.avgk, data=gr.C3yn.avgk)

# Evaluate assumptions of linear regression
hist(gr.C3yn.avgk.lm$residuals) # evaluate assumption of normality
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

``` r
plot(predict(gr.C3yn.avgk.lm), resid(gr.C3yn.avgk.lm)) # evaluate assumption of homoskedasticity
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-32-2.png)<!-- -->

``` r
# P value
summary(gr.C3yn.avgk.lm)$coefficients[2,4]
```

    ## [1] 0.09315551

``` r
# Successional C-amended vs water
gr.S17yn.avgk = inner_join(gr.S17y.avgk, gr.S17n.avgk, by="ASV")
gr.S17yn.avgk.lm = lm(S17y.avgk ~ S17n.avgk, data=gr.S17yn.avgk)

# Evaluate assumptions of linear regression
hist(gr.S17yn.avgk.lm$residuals) # evaluate assumption of normality
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

``` r
plot(predict(gr.S17yn.avgk.lm), resid(gr.S17yn.avgk.lm)) # evaluate assumption of homoskedasticity
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-33-2.png)<!-- -->

``` r
# P value
summary(gr.S17yn.avgk.lm)$coefficients[2,4]
```

    ## [1] 0.6783426

Assumptions look fine, but the successional overlap data is just sparse.
Anyways, nothing is significant.

# Explore bimodal distributions - slow and fast groups?

The growth rate distribution graphs show distinct bimodal distributions
in the water treatment. This seems to indicate differing populations,
which I want to explore further.

Split data into “fast” and “slow” groups based on the bimodal
distributions.

**Split by the bimodal distributions in unamended soils:**

Choose thresholds for each soil:

``` r
gr.estavg %>% 
  filter(Soil=="C3") %>%
  ggplot(aes(x=log(avgk), color=Amendment, fill=Amendment)) +
  geom_density(alpha=0.5) +
  theme_test() +
  theme(#legend.position="none",
        title=element_text(size=12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  scale_color_manual(values=c("#E78E43", "#423EC5")) +
  scale_fill_manual(values=c("#E78E43", "#423EC5")) +
  labs(title="Cropped", x="ln speciifc growth rate", y="Density") +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(-5,0)) +
  geom_vline(xintercept=-2.4, linetype=2)
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

``` r
gr.estavg %>% 
  filter(Soil=="S17") %>%
  ggplot(aes(x=log(avgk), color=Amendment, fill=Amendment)) +
  geom_density(alpha=0.5) +
  theme_test() +
  theme(#legend.position="none",
        title=element_text(size=12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  scale_color_manual(values=c("#E78E43", "#423EC5")) +
  scale_fill_manual(values=c("#E78E43", "#423EC5")) +
  labs(title="Successional", x="ln speciifc growth rate", y="Density") +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(-5,0)) +
  geom_vline(xintercept=-1.9, linetype=2)
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-34-2.png)<!-- -->

Label ASVs based on avg growth rate threshold:

``` r
# Water only
gr.estavg.group = gr.estavg %>% 
  filter(Amendment=="N") %>% 
  mutate(group = ifelse(Soil=="C3" & log(avgk) < -2.4, "slow", ifelse(Soil=="S17" & log(avgk) < -1.9, "slow", "fast")))

gr.estavg.group = gr.estavg.group[,c(1,3:ncol(gr.estavg.group))]
```

``` r
# save 
saveRDS(gr.estavg.group, file="rdata.files/gr.groups.rds")
```

### Abundances of groups in field communities

Do these taxa exist in the in-field community data? If so, what are
they’re abundances? I.e. did the soil shape the growth dynamics we
observed in the microcosms?

Hypothesis: abundances of slow and fast taxa within field communities
differ between cropped and successional soils.

Import field data:

``` r
# Import data
field = readRDS("rdata.files/field.norm.cleaned.rds")

# Add group labels
field.group = field %>% 
  inner_join(gr.estavg.group, var=c("Soil", "ASV", "Domain", "Phylum", "Class", "Order", "Family", "Genus",)) %>%
  filter(!norm_abund==0) %>% # remove 0 values, ie ASV not present
  mutate(group = parse_factor(group, levels=c("slow", "fast"))) # more intuitive order for groups
```

#### Individual ASV normalized abundances in each group:

``` r
# Visualize
field.group %>% 
  ggplot(aes(x=group, y=log(norm_abund), color=as.factor(Block))) +
  geom_boxplot() +
  geom_jitter(alpha=0.2) +
  facet_wrap(~Soil) +
  labs(title="Individual ASV abundances", x="", y="natural log internal-standard normalized abundance", color="Block") +
  theme_test()
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

**Statistics**

Linear mixed effects model

  - Block as random effect
  - Repeated measures for slow and fast groups, which were each measured
    from the same sample

<!-- end list -->

``` r
# Cropped
field.group.ind.C3.lmer = lmer(log(norm_abund) ~ group + Block + (1|Replicate), data=field.group[field.group$Soil=="C3",]) # fit model
plot(resid(field.group.ind.C3.lmer)) # check normality
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

``` r
plot(predict(field.group.ind.C3.lmer), resid(field.group.ind.C3.lmer)) # check variances
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-39-2.png)<!-- -->

``` r
anova(field.group.ind.C3.lmer) # results
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##       Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## group 26.950  26.950     1 907.03  21.747 3.578e-06 ***
    ## Block 16.859  16.859     1 909.85  13.604  0.000239 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Successional
field.group.ind.S17.lmer = lmer(log(norm_abund) ~ group + Block + (1|Replicate), data=field.group[field.group$Soil=="S17",]) # fit model
plot(resid(field.group.ind.S17.lmer)) # check normality
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-39-3.png)<!-- -->

``` r
plot(predict(field.group.ind.S17.lmer), resid(field.group.ind.S17.lmer)) # check variances
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-39-4.png)<!-- -->

``` r
anova(field.group.ind.S17.lmer) # results
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##       Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
    ## group 0.3769  0.3769     1 400.04  0.2323 0.63011  
    ## Block 5.6536  5.6536     1 401.02  3.4842 0.06269 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#### Total, aggregated abundances in each group:

``` r
# Sum abundances in each group
field.group.sum = field.group %>% 
  group_by(Soil, Block, Replicate, group) %>%
  summarize(total_abund = sum(norm_abund))

# Visualize
field.group.sum %>% 
  ggplot(aes(x=group, y=log(total_abund), color=as.factor(Block))) +
  geom_boxplot() +
  geom_jitter(alpha=0.2) +
  facet_wrap(~Soil) +
  labs(title="Total, aggregated abundance", x="", y="natural log internal-standard normalized abundance", color="Block") +
  theme_test()
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

``` r
field.group.sum %>% 
  ggplot(aes(x=group, y=log(total_abund))) +
  geom_boxplot() +
  geom_jitter(alpha=0.2) +
  facet_wrap(~Soil) +
  labs(title="Total, aggregated abundance", x="", y="natural log internal-standard normalized abundance", color="Block") +
  theme_test()
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-40-2.png)<!-- -->

Strong block effect in cropped soil.

**Statistics**

Linear mixed effects model

  - Block as random effect
  - Repeated measures for slow and fast groups, which were each measured
    from the same sample

<!-- end list -->

``` r
# Cropped
field.group.sum.C3.lmer = lmer(log(total_abund) ~ group + Block + (1|Replicate), data=field.group.sum[field.group.sum$Soil=="C3",]) # fit model
plot(resid(field.group.sum.C3.lmer)) # check normality
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

``` r
plot(predict(field.group.sum.C3.lmer), resid(field.group.sum.C3.lmer)) # check variances
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-41-2.png)<!-- -->

``` r
anova(field.group.sum.C3.lmer) # results
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##        Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## group 23.7801 23.7801     1 8.4216 35.9922 0.0002642 ***
    ## Block  2.2964  2.2964     1 9.1912  3.4757 0.0944703 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Successional
field.group.sum.S17.lmer = lmer(log(total_abund) ~ group + Block + (1|Replicate), data=field.group.sum[field.group.sum$Soil=="S17",]) # fit model
plot(resid(field.group.sum.S17.lmer)) # check normality
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-41-3.png)<!-- -->

``` r
plot(predict(field.group.sum.S17.lmer), resid(field.group.sum.S17.lmer)) # check variances
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-41-4.png)<!-- -->

``` r
anova(field.group.sum.S17.lmer) # results
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##         Sum Sq  Mean Sq NumDF DenDF F value Pr(>F)
    ## group 0.244902 0.244902     1    10  0.2121 0.6550
    ## Block 0.005137 0.005137     1    10  0.0044 0.9481

Conclusion: we found evidence that, in the cropped soil, ASVs in the
“fast” group are, in-total and individually, more abundant than ASVs
in the “slow” group. No such evidence in successional soil.

# Figures

### Distribution of growth rates

``` r
p1 = gr.estavg %>%
  mutate(Soil = recode_factor(Soil, S17="Successional", C3="Cropped"),
         Amendment = recode_factor(Amendment, N="Water", Y="C-added")) %>%
  ggplot(aes(x=log(avgk), fill=Soil)) +
  geom_density(alpha=0.5) +
  facet_wrap(~Amendment, ncol=1) +
  scale_fill_manual(values=c(NA, "gray")) +
  labs(y="Kernel density", x="Specific growth rate") +
  theme_test() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_blank(),
        legend.position = c(.25,.93),
        legend.text = element_text(size=8),
        legend.key.width=unit(0.5, "cm"),
        legend.key.height=unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(margin = unit(c(1, 0, 0, 0), "mm"), size=10),
        axis.title.y = element_text(margin = unit(c(0, 1.5, 0, 0), "mm"), size=10),
        axis.text = element_text(size = 8))
p1
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

``` r
p2 = gr.bins %>%
  mutate(Soil = recode_factor(Soil, S17="Successional", C3="Cropped"),
         Amendment = recode_factor(Amendment, N="Water", Y="C-added")) %>%
  ggplot(aes(x=bin, y=frequency, fill=Soil)) +
  geom_boxplot(alpha=0.5) +
  facet_wrap(~Amendment, ncol=1) +
  scale_fill_manual(values=c(NA, "gray")) +
  labs(x="Bin number", y="Number of ASVs") +
  theme_test() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(margin = unit(c(1, 0, 0, 0), "mm"), size=10),
        axis.title.y = element_text(margin = unit(c(0, 1.5, 0, 0), "mm"), size=10),
        axis.text = element_text(size = 8))
p2
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-42-2.png)<!-- -->

``` r
habitatdist.plot = plot_grid(p1, NULL, p2, nrow=1, align="h", rel_widths=c(1,0.04,1), labels=c("a", NA, "b"), label_size=10)
habitatdist.plot
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-42-3.png)<!-- -->

``` r
ggsave(plot=habitatdist.plot, file="figures/fig_habitatdist.pdf", device=pdf, width=140, height=90, units="mm")
```

### Shift in growth rate

``` r
p3 = inner_join(gr.C3y.avgk, gr.C3n.avgk, by="ASV") %>% 
  ggplot(aes(x=log(C3y.avgk), y=log(C3n.avgk))) +
  geom_point(shape=1) +
  #geom_smooth(method="lm", color="black", linetype=2) +
  labs(x="Cropped, C-amended", y="Cropped, water") +
  scale_x_continuous(limits=c(-3,0)) +
  scale_y_continuous(limits=c(-3,0)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(margin = unit(c(1, 0, 0, 0), "mm"), size=10),
        axis.title.y = element_text(margin = unit(c(0, 1.5, 0, 0), "mm"), size=10),
        axis.text = element_text(size = 8))
p3
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

``` r
p4 = inner_join(gr.S17y.avgk, gr.S17n.avgk, by="ASV") %>% 
  ggplot(aes(x=log(S17y.avgk), y=log(S17n.avgk))) +
  geom_point(shape=1) +
  #geom_smooth(method="lm", color="black", linetype=2) +
  labs(x="Successional, C-amended", y="Succesional, water") +
  scale_x_continuous(limits=c(-3,0)) +
  scale_y_continuous(limits=c(-3,0)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(margin = unit(c(1, 0, 0, 0), "mm"), size=10),
        axis.title.y = element_text(margin = unit(c(0, 1.5, 0, 0), "mm"), size=10),
        axis.text = element_text(size = 8))
p4
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-44-2.png)<!-- -->

``` r
p5 = inner_join(gr.C3y.avgk, gr.S17y.avgk, by="ASV") %>% 
  ggplot(aes(x=log(C3y.avgk), y=log(S17y.avgk))) +
  geom_point(shape=1) +
  #geom_smooth(method="lm", color="black", linetype=2) +
  labs(x="Cropped, C-amended", y="Succesional, C-amended") +
  scale_x_continuous(limits=c(-3,0)) +
  scale_y_continuous(limits=c(-3,0)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(margin = unit(c(1, 0, 0, 0), "mm"), size=10),
        axis.title.y = element_text(margin = unit(c(0, 1.5, 0, 0), "mm"), size=10),
        axis.text = element_text(size = 8))
p5
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-44-3.png)<!-- -->

``` r
p6 = inner_join(gr.C3n.avgk, gr.S17n.avgk, by="ASV") %>% 
  ggplot(aes(x=log(C3n.avgk), y=log(S17n.avgk))) +
  geom_point(shape=1) +
  #geom_smooth(method="lm", color="black", linetype=2) +
  labs(x="Cropped, water", y="Succesional, water") +
  scale_x_continuous(limits=c(-3,0)) +
  scale_y_continuous(limits=c(-3,0)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(margin = unit(c(1, 0, 0, 0), "mm"), size=10),
        axis.title.y = element_text(margin = unit(c(0, 1.5, 0, 0), "mm"), size=10),
        axis.text = element_text(size = 8))
p6
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-44-4.png)<!-- -->

``` r
shift.plot = plot_grid(p3, p4, p5, p6)
shift.plot
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-44-5.png)<!-- -->

``` r
ggsave(shift.plot, file="figures/suppfig_shiftgr.pdf", device=pdf, width=140, height=140, units="mm")
```

### Start and end of growth

Note sure if will include.

``` r
p7 = gr.est.avg %>%
  mutate(Soil = recode_factor(Soil, C3="Cropped", S17="Successional"),
         Amendment = recode_factor(Amendment, Y="C-amended", N="Water")) %>%
  ggplot(aes(x=Soil, y=avg_start)) +
  geom_boxplot(aes(fill=Amendment)) +
  scale_fill_manual(values=c(NA, "gray")) +
  scale_y_continuous(limits=c(0,4)) +
  labs(title="Start of growth", y="Day") +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.15, 0.8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_blank(),
        title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 8))
p7
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

``` r
p8 = gr.est.avg %>%
  mutate(Soil = recode_factor(Soil, C3="Cropped", S17="Successional"),
         Amendment = recode_factor(Amendment, Y="C-amended", N="Water")) %>%
  ggplot(aes(x=Soil, y=avg_end)) +
  geom_boxplot(aes(fill=Amendment)) +
  scale_fill_manual(values=c(NA, "gray")) +
  scale_y_continuous(limits=c(0,10)) +
  labs(title="End of growth", y="Day") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_blank(),
        title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 8))
p8
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-46-2.png)<!-- -->

``` r
startend.plots = plot_grid(p7, p8, nrow=2, labels=c("a", "b"), label_size=10)
startend.plots
```

![](05_dunlopgr_habitat_files/figure-gfm/unnamed-chunk-46-3.png)<!-- -->

``` r
ggsave(startend.plots, file="figures/fig_startend.pdf", device=pdf, width=140, height=90, units="mm")
```
