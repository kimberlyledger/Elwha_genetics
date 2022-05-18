steelhead_PCA_DAPC
================
Kimberly Ledger
4/12/2022

# looking at O.mykiss population structure

load libraries

``` r
library(adegenet) #this package is used for analysis of genetic/genomic data 
library(dplyr) # data manipulation
library(tidyr) # data manipulation
library(forcats)
library(ggplot2)
```

# Part 1: create genind objects

``` r
#onmy_n_df <- read.csv("outputs/onmy_loci4genstr_outflank_hwe.csv")
onmy_n_df <- read.csv("outputs/omy_loci4genstr_outflank_hwe_oneadapt.csv")
```

here i will subset PRE and POST dam individuals

``` r
geno_pre <- onmy_n_df %>%
  filter(Time == "Pre") %>%
  dplyr::mutate(Location = fct_relevel(Location, "BD", "SBLR", "ID", "AD")) %>% 
  arrange(Location, rkm)

geno_post <- onmy_n_df %>%
  filter(Time == "Post") %>%
  dplyr::mutate(Location = fct_relevel(Location, "BD", "SBLR", "ID", "AD")) %>% 
  arrange(Location, rkm)
```

## create genind objects from pre and post dam samples

``` r
#geno_n_pre <- geno_pre[,c(2:297)]        #for "outputs/onmy_loci4genstr_outflank_hwe.csv" use 297 
geno_n_pre <- geno_pre[,c(2:265)]        #for "outputs/omy_loci4genstr_outflank_hwe_oneadapt.csv" use 265

rownames(geno_n_pre) <- geno_pre$Sample_ID
col_geno_n_pre <- gsub("\\.", "_", colnames(geno_n_pre))
colnames(geno_n_pre) <- col_geno_n_pre

#meta_pre <- geno_pre[,-c(1:297)]
meta_pre <- geno_pre[,-c(1:265)]

pre_loci <- colnames(geno_n_pre)
pre_ind <- rownames(geno_n_pre)
pre_location <- meta_pre$Location
pre_rkm <- meta_pre$rkm
pre_site <- meta_pre$Sampling_Site

onmy_genind_pre_location <- df2genind(geno_n_pre,
                         sep="/",
                         ind.names=pre_ind,
                         loc.names=pre_loci, 
                         pop = pre_location,     
                         ploidy = 2)
onmy_genind_pre_site <- df2genind(geno_n_pre,
                         sep="/",
                         ind.names=pre_ind,
                         loc.names=pre_loci, 
                         pop = pre_site,      
                         ploidy = 2)
onmy_genind_pre_location
```

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 531 individuals; 264 loci; 522 alleles; size: 1.3 Mb
    ## 
    ##  // Basic content
    ##    @tab:  531 x 522 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 1-2)
    ##    @loc.fac: locus factor for the 522 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: df2genind(X = geno_n_pre, sep = "/", ind.names = pre_ind, loc.names = pre_loci, 
    ##     pop = pre_location, ploidy = 2)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 81-166)

``` r
onmy_genind_pre_site
```

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 531 individuals; 264 loci; 522 alleles; size: 1.3 Mb
    ## 
    ##  // Basic content
    ##    @tab:  531 x 522 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 1-2)
    ##    @loc.fac: locus factor for the 522 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: df2genind(X = geno_n_pre, sep = "/", ind.names = pre_ind, loc.names = pre_loci, 
    ##     pop = pre_site, ploidy = 2)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 1-128)

``` r
#geno_n_post <- geno_post[,c(2:297)]
geno_n_post <- geno_post[,c(2:265)]

rownames(geno_n_post) <- geno_post$Sample_ID
col_geno_n_post <- gsub("\\.", "_", colnames(geno_n_post))
colnames(geno_n_post) <- col_geno_n_post

#meta_post <- geno_post[,-c(1:297)]
meta_post <- geno_post[,-c(1:265)]

post_loci <- colnames(geno_n_post)
post_ind <- rownames(geno_n_post)
post_location <- meta_post$Location
post_rkm <- meta_post$rkm
post_site <- meta_post$Sampling_Site

onmy_genind_post_location <- df2genind(geno_n_post,
                         sep="/",
                         ind.names=post_ind,
                         loc.names=post_loci, 
                         pop = post_location,
                         ploidy = 2)

onmy_genind_post_site <- df2genind(geno_n_post,
                         sep="/",
                         ind.names=post_ind,
                         loc.names=post_loci, 
                         pop = post_site,
                         ploidy = 2)
onmy_genind_post_location
```

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 934 individuals; 264 loci; 525 alleles; size: 2.1 Mb
    ## 
    ##  // Basic content
    ##    @tab:  934 x 525 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 1-2)
    ##    @loc.fac: locus factor for the 525 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: df2genind(X = geno_n_post, sep = "/", ind.names = post_ind, loc.names = post_loci, 
    ##     pop = post_location, ploidy = 2)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 25-445)

``` r
onmy_genind_post_site
```

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 934 individuals; 264 loci; 525 alleles; size: 2.1 Mb
    ## 
    ##  // Basic content
    ##    @tab:  934 x 525 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 1-2)
    ##    @loc.fac: locus factor for the 525 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: df2genind(X = geno_n_post, sep = "/", ind.names = post_ind, loc.names = post_loci, 
    ##     pop = post_site, ploidy = 2)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 1-280)

# Part 2: Hardy-Weinberg equilibrium - AS OF 17 May 2022 HWE was checkin in “steeloutliertest.Rmd”

# Part 4: convet genind to genlight objects

``` r
#install_github("green-striped-gecko/dartR")
library(dartR)

genlit_pre_site <- gi2gl(onmy_genind_pre_site)
```

    ## Starting gi2gl 
    ## Starting gl.compliance.check 
    ##   Processing genlight object with SNP data
    ##   Checking coding of SNPs
    ##     SNP data scored NA, 0, 1 or 2 confirmed
    ##   Checking locus metrics and flags
    ##   Recalculating locus metrics
    ##   Checking for monomorphic loci
    ##     Dataset contains monomorphic loci
    ##   Checking whether individual names are unique.
    ##   Checking for individual metrics
    ##   Warning: Creating a slot for individual metrics
    ##   Checking for population assignments
    ##     Population assignments confirmed
    ##   Spelling of coordinates checked and changed if necessary to lat/lon
    ## Completed: gl.compliance.check 
    ## Completed: gi2gl

``` r
genlit_pre_location <- gi2gl(onmy_genind_pre_location)
```

    ## Starting gi2gl 
    ## Starting gl.compliance.check 
    ##   Processing genlight object with SNP data
    ##   Checking coding of SNPs
    ##     SNP data scored NA, 0, 1 or 2 confirmed
    ##   Checking locus metrics and flags
    ##   Recalculating locus metrics
    ##   Checking for monomorphic loci
    ##     Dataset contains monomorphic loci
    ##   Checking whether individual names are unique.
    ##   Checking for individual metrics
    ##   Warning: Creating a slot for individual metrics
    ##   Checking for population assignments
    ##     Population assignments confirmed
    ##   Spelling of coordinates checked and changed if necessary to lat/lon
    ## Completed: gl.compliance.check 
    ## Completed: gi2gl

``` r
genlit_post_site <- gi2gl(onmy_genind_post_site)
```

    ## Starting gi2gl 
    ## Starting gl.compliance.check 
    ##   Processing genlight object with SNP data
    ##   Checking coding of SNPs
    ##     SNP data scored NA, 0, 1 or 2 confirmed
    ##   Checking locus metrics and flags
    ##   Recalculating locus metrics
    ##   Checking for monomorphic loci
    ##     Dataset contains monomorphic loci
    ##   Checking whether individual names are unique.
    ##   Checking for individual metrics
    ##   Warning: Creating a slot for individual metrics
    ##   Checking for population assignments
    ##     Population assignments confirmed
    ##   Spelling of coordinates checked and changed if necessary to lat/lon
    ## Completed: gl.compliance.check 
    ## Completed: gi2gl

``` r
genlit_post_location <- gi2gl(onmy_genind_post_location)
```

    ## Starting gi2gl 
    ## Starting gl.compliance.check 
    ##   Processing genlight object with SNP data
    ##   Checking coding of SNPs
    ##     SNP data scored NA, 0, 1 or 2 confirmed
    ##   Checking locus metrics and flags
    ##   Recalculating locus metrics
    ##   Checking for monomorphic loci
    ##     Dataset contains monomorphic loci
    ##   Checking whether individual names are unique.
    ##   Checking for individual metrics
    ##   Warning: Creating a slot for individual metrics
    ##   Checking for population assignments
    ##     Population assignments confirmed
    ##   Spelling of coordinates checked and changed if necessary to lat/lon
    ## Completed: gl.compliance.check 
    ## Completed: gi2gl

# consider filtering genlight objects prior to running a PCA or PCoA?

prep for pca (a) Filter stringently on call rate, using a threshold of
at least 95% loci called. (b) Remove individuals for which call rate is
exceptionally low, say \<80%. (c) Impute the remaining missing values on
a population‐by‐population basis, where populations can be considered
panmictic.

``` r
#genlit_filter <- gl.filter.callrate(genlit_pre_site, method="loc", threshold=0.90)
#genlit_filter <- gl.filter.callrate(genlit_filter, method="ind", threshold=0.80)
#genlit_filter <- gl.impute(genlit_filter, method="random")
#pcoa <- gl.pcoa(genlit_filter)
```

this is for pre-dam samples only

``` r
#library(directlabels)
#gl.pcoa.plot(pcoa, genlit_filter, xaxis = 1, yaxis =2, ellipse = TRUE, plevel = 0.9)
```

##PCA using filtered genlight object set color palette

``` r
library(viridisLite)
mycol = viridisLite::viridis(4, alpha = 0.6)
```

``` r
pca_pre <- glPca(genlit_pre_location, center = T, scale = F, nf = 50)
s.class(pca_pre$scores, pop(genlit_pre_location),
        xax=1, yax=2, col=mycol,
        axesel=FALSE, cstar=0, cpoint=3)
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
pca_post <- glPca(genlit_post_location, center = T, scale = F, nf = 50)
s.class(pca_post$scores, pop(genlit_post_location),
        xax=1, yax=2, col=mycol,
        axesel=FALSE, cstar=0, cpoint=3)
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

# Part 5: DAPC

i am following this tutorial:
<https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf>

find clusters

``` r
grp_pre <- find.clusters(genlit_pre_location, max.n.clust=40, n.pca=100)
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

    ## Choose the number of clusters (>=2):

check Kstat

``` r
grp_pre$Kstat
```

    ##      K=1      K=2      K=3      K=4      K=5      K=6      K=7      K=8 
    ## 1427.851 1413.072 1410.153 1409.635 1409.718 1410.957 1412.912 1413.764 
    ##      K=9     K=10     K=11     K=12     K=13     K=14     K=15     K=16 
    ## 1415.700 1418.400 1420.697 1423.682 1427.035 1428.476 1430.436 1435.951 
    ##     K=17     K=18     K=19     K=20     K=21     K=22     K=23     K=24 
    ## 1436.980 1441.007 1442.865 1447.353 1450.575 1454.071 1457.614 1460.196 
    ##     K=25     K=26     K=27     K=28     K=29     K=30     K=31     K=32 
    ## 1464.500 1467.343 1471.504 1474.434 1478.670 1482.271 1486.044 1488.451 
    ##     K=33     K=34     K=35     K=36     K=37     K=38     K=39     K=40 
    ## 1492.445 1497.222 1500.135 1502.794 1507.434 1512.342 1513.836 1518.069

will use K=2 and K=3 since BIC stops decreasing all that much beyond
there…

``` r
grp_pre_2 <- find.clusters(genlit_pre_location, max.n.clust=40, n.clust = 2)
```

    ## Choose the number PCs to retain (>=1):

``` r
grp_pre_3 <- find.clusters(genlit_pre_location, max.n.clust=40, n.clust = 3)
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

    ## Choose the number PCs to retain (>=1):

run dapc

``` r
dapc_pre_2 <- dapc(genlit_pre_location, grp_pre_2$grp, n.pca = 40, n.da = 10) #keep first 40 PCs; and 10 n.da
dapc_pre_3 <- dapc(genlit_pre_location, grp_pre_3$grp, n.pca = 40, n.da = 10) #keep first 40 PCs; and 10 n.da
```

dapc output

``` r
dapc_pre_2
```

    ##  #################################################
    ##  # Discriminant Analysis of Principal Components #
    ##  #################################################
    ## class: dapc
    ## $call: dapc.genlight(x = genlit_pre_location, pop = grp_pre_2$grp, n.pca = 40, 
    ##     n.da = 10)
    ## 
    ## $n.pca: 40 first PCs of PCA used
    ## $n.da: 1 discriminant functions saved
    ## $var (proportion of conserved variance): 0.453
    ## 
    ## $eig (eigenvalues): 2922  vector    length content                   
    ## 1 $eig      1      eigenvalues               
    ## 2 $grp      531    prior group assignment    
    ## 3 $prior    2      prior group probabilities 
    ## 4 $assign   531    posterior group assignment
    ## 5 $pca.cent 264    centring vector of PCA    
    ## 6 $pca.norm 264    scaling vector of PCA     
    ## 7 $pca.eig  258    eigenvalues of PCA        
    ## 
    ##   data.frame    nrow ncol content                                          
    ## 1 $tab          531  40   retained PCs of PCA                              
    ## 2 $means        2    40   group means                                      
    ## 3 $loadings     40   1    loadings of variables                            
    ## 4 $ind.coord    531  1    coordinates of individuals (principal components)
    ## 5 $grp.coord    2    1    coordinates of groups                            
    ## 6 $posterior    531  2    posterior membership probabilities               
    ## 7 $pca.loadings 264  40   PCA loadings of original variables               
    ## 8 $var.contr    264  1    contribution of original variables

``` r
dapc_pre_3
```

    ##  #################################################
    ##  # Discriminant Analysis of Principal Components #
    ##  #################################################
    ## class: dapc
    ## $call: dapc.genlight(x = genlit_pre_location, pop = grp_pre_3$grp, n.pca = 40, 
    ##     n.da = 10)
    ## 
    ## $n.pca: 40 first PCs of PCA used
    ## $n.da: 2 discriminant functions saved
    ## $var (proportion of conserved variance): 0.453
    ## 
    ## $eig (eigenvalues): 1540 588.1  vector    length content                   
    ## 1 $eig      2      eigenvalues               
    ## 2 $grp      531    prior group assignment    
    ## 3 $prior    3      prior group probabilities 
    ## 4 $assign   531    posterior group assignment
    ## 5 $pca.cent 264    centring vector of PCA    
    ## 6 $pca.norm 264    scaling vector of PCA     
    ## 7 $pca.eig  258    eigenvalues of PCA        
    ## 
    ##   data.frame    nrow ncol content                                          
    ## 1 $tab          531  40   retained PCs of PCA                              
    ## 2 $means        3    40   group means                                      
    ## 3 $loadings     40   2    loadings of variables                            
    ## 4 $ind.coord    531  2    coordinates of individuals (principal components)
    ## 5 $grp.coord    3    2    coordinates of groups                            
    ## 6 $posterior    531  3    posterior membership probabilities               
    ## 7 $pca.loadings 264  40   PCA loadings of original variables               
    ## 8 $var.contr    264  2    contribution of original variables

plot DAPC using scatter()

``` r
scatter(dapc_pre_2, col = viridis(2, alpha = 0.6))
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
scatter(dapc_pre_3, col = viridis(3, alpha = 0.6))
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

plot DAPC using ggplot

``` r
dapc_pre_2_scores <- as.data.frame(dapc_pre_2$tab)
dapc_pre_2_scores$Assigned_Pop <- dapc_pre_2$grp
dapc_pre_2_scores$Original_Pop <- meta_pre$Location

set.seed(9)
dapc_p2 <- ggplot(dapc_pre_2_scores, aes(x=PC1, y=PC2, colour=Assigned_Pop)) 
dapc_p2 <- dapc_p2 + geom_point(size=2, aes(shape = Original_Pop)) 
dapc_p2 <- dapc_p2 + stat_ellipse(level = 0.95, size = 1)
dapc_p2 <- dapc_p2 + scale_color_manual(values = viridis(3)) 
dapc_p2 <- dapc_p2 + geom_hline(yintercept = 0) 
dapc_p2 <- dapc_p2 + geom_vline(xintercept = 0) 
dapc_p2 <- dapc_p2 + theme_bw()

dapc_p2
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

plot DAPC using ggplot

``` r
dapc_pre_3_scores <- as.data.frame(dapc_pre_3$tab)
dapc_pre_3_scores$Assigned_Pop <- dapc_pre_3$grp
dapc_pre_3_scores$Original_Pop <- meta_pre$Location

set.seed(9)
dapc_p3 <- ggplot(dapc_pre_3_scores, aes(x=PC1, y=PC2, colour=Assigned_Pop)) 
dapc_p3 <- dapc_p3 + geom_point(size=2, aes(shape = Original_Pop)) 
dapc_p3 <- dapc_p3 + stat_ellipse(level = 0.95, size = 1)
dapc_p3 <- dapc_p3 + scale_color_manual(values = viridis(3)) 
dapc_p3 <- dapc_p3 + geom_hline(yintercept = 0) 
dapc_p3 <- dapc_p3 + geom_vline(xintercept = 0) 
dapc_p3 <- dapc_p3 + theme_bw()
dapc_p3
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

which loci have the greatest contribution?

``` r
set.seed(4)
loadingplot(dapc_pre_2$var.contr, axis=1,
                       thres=.07, lab.jitter=1)
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

    ## NULL

``` r
loadingplot(dapc_pre_3$var.contr, axis=1,
                       thres=.07, lab.jitter=1)
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

    ## NULL

``` r
loadingplot(dapc_pre_3$var.contr, axis=2,
                       thres=.07, lab.jitter=1)
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->

    ## NULL

look at loci with \>0.05 contribution to pc axis

``` r
data.frame(dapc_pre_3$var.contr) %>%
  filter(LD2 > 0.05)
```

    ##             LD1        LD2
    ## 146 0.005303561 0.05985673

``` r
pre_loci[146]
```

    ## [1] "NC_035094_1_50984953"

## this is Omy_RAD3209-10 on omy 18 (Adaptive. Basin-wide, top-outlier)

``` r
#summary(dapc_pre_2)
#summary(dapc_pre_3)
```

genotype compostion plot - individuals are ordered from BD to AD

``` r
compoplot(dapc_pre_2, posi="bottomright",
          txt.leg=paste("Cluster", 1:2), xlab="individuals", col=viridis(2))
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
compoplot(dapc_pre_3, posi="bottomright",
          txt.leg=paste("Cluster", 1:3), xlab="individuals", col=viridis(3))
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->

make easier to interpret DAPC plot using ggplot

``` r
dapc.results <- as.data.frame(dapc_pre_2$posterior)
dapc.results$pop <- meta_pre$Location
dapc.results$indNames <- rownames(dapc.results)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))  #pivot table to be able to use ggplot 
```

rename columns and plot

``` r
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

p2 <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p2 <- p2 + geom_bar(stat='identity') 
p2 <- p2 + scale_fill_manual(values = viridis(3)) 
p2 <- p2 + facet_grid(~Original_Pop, scales = "free")
p2 <- p2 + theme(axis.text.x = element_blank())
p2
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

make easier to interpret DAPC plot using ggplot

``` r
dapc.results <- as.data.frame(dapc_pre_3$posterior)
dapc.results$pop <- meta_pre$Location
dapc.results$indNames <- rownames(dapc.results)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))  #pivot table to be able to use ggplot 
```

rename columns and plot

``` r
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

p3 <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p3 <- p3 + geom_bar(stat='identity') 
p3 <- p3 + scale_fill_manual(values = viridis(3)) 
p3 <- p3 + facet_grid(~Original_Pop, scales = "free")
p3 <- p3 + theme(axis.text.x = element_blank())
p3
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

make combined figure for predam individuals…

``` r
library(cowplot)

pre_plot <- plot_grid(dapc_p2, p2, dapc_p3, p3, labels = "AUTO", ncol =2)
pre_plot
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
ggsave2("outputs/pre_plot.png", width = 14, height = 8)
```

# Now analyses post-removal dam individuals

find clusters

``` r
grp_post <- find.clusters(genlit_post_location, max.n.clust=40, n.pca=100)
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

    ## Choose the number of clusters (>=2):

check Kstat

``` r
grp_post$Kstat
```

    ##      K=1      K=2      K=3      K=4      K=5      K=6      K=7      K=8 
    ## 2453.333 2442.355 2436.224 2433.806 2429.950 2428.272 2427.024 2427.334 
    ##      K=9     K=10     K=11     K=12     K=13     K=14     K=15     K=16 
    ## 2427.743 2429.750 2431.376 2432.859 2434.350 2435.933 2438.228 2440.638 
    ##     K=17     K=18     K=19     K=20     K=21     K=22     K=23     K=24 
    ## 2442.950 2445.880 2448.181 2451.054 2454.132 2455.391 2458.771 2462.276 
    ##     K=25     K=26     K=27     K=28     K=29     K=30     K=31     K=32 
    ## 2465.358 2468.679 2471.162 2472.366 2477.947 2480.711 2482.787 2485.909 
    ##     K=33     K=34     K=35     K=36     K=37     K=38     K=39     K=40 
    ## 2490.395 2493.894 2497.229 2500.861 2504.294 2506.271 2512.049 2515.124

will use K=2 to K=5 since BIC stops decreasing all that much beyond
there…

``` r
grp_post_2 <- find.clusters(genlit_post_location, max.n.clust=40, n.clust = 2)
```

    ## Choose the number PCs to retain (>=1):

``` r
grp_post_3 <- find.clusters(genlit_post_location, max.n.clust=40, n.clust = 3)
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

    ## Choose the number PCs to retain (>=1):

``` r
grp_post_4 <- find.clusters(genlit_post_location, max.n.clust=40, n.clust = 4)
```

    ## Choose the number PCs to retain (>=1):

``` r
grp_post_5 <- find.clusters(genlit_post_location, max.n.clust=40, n.clust = 5)
```

    ## Choose the number PCs to retain (>=1):

run dapc

``` r
dapc_post_2 <- dapc(genlit_post_location, grp_post_2$grp, n.pca = 40, n.da = 10) #keep first 40 PCs; and 10 n.da
dapc_post_3 <- dapc(genlit_post_location, grp_post_3$grp, n.pca = 40, n.da = 10) #keep first 40 PCs; and 10 n.da
dapc_post_4 <- dapc(genlit_post_location, grp_post_4$grp, n.pca = 40, n.da = 10) #keep first 40 PCs; and 10 n.da
dapc_post_5 <- dapc(genlit_post_location, grp_post_5$grp, n.pca = 40, n.da = 10) #keep first 40 PCs; and 10 n.da
```

dapc outputs

``` r
dapc_post_2
```

    ##  #################################################
    ##  # Discriminant Analysis of Principal Components #
    ##  #################################################
    ## class: dapc
    ## $call: dapc.genlight(x = genlit_post_location, pop = grp_post_2$grp, 
    ##     n.pca = 40, n.da = 10)
    ## 
    ## $n.pca: 40 first PCs of PCA used
    ## $n.da: 1 discriminant functions saved
    ## $var (proportion of conserved variance): 0.402
    ## 
    ## $eig (eigenvalues): 2273  vector    length content                   
    ## 1 $eig      1      eigenvalues               
    ## 2 $grp      934    prior group assignment    
    ## 3 $prior    2      prior group probabilities 
    ## 4 $assign   934    posterior group assignment
    ## 5 $pca.cent 264    centring vector of PCA    
    ## 6 $pca.norm 264    scaling vector of PCA     
    ## 7 $pca.eig  261    eigenvalues of PCA        
    ## 
    ##   data.frame    nrow ncol content                                          
    ## 1 $tab          934  40   retained PCs of PCA                              
    ## 2 $means        2    40   group means                                      
    ## 3 $loadings     40   1    loadings of variables                            
    ## 4 $ind.coord    934  1    coordinates of individuals (principal components)
    ## 5 $grp.coord    2    1    coordinates of groups                            
    ## 6 $posterior    934  2    posterior membership probabilities               
    ## 7 $pca.loadings 264  40   PCA loadings of original variables               
    ## 8 $var.contr    264  1    contribution of original variables

``` r
dapc_post_3
```

    ##  #################################################
    ##  # Discriminant Analysis of Principal Components #
    ##  #################################################
    ## class: dapc
    ## $call: dapc.genlight(x = genlit_post_location, pop = grp_post_3$grp, 
    ##     n.pca = 40, n.da = 10)
    ## 
    ## $n.pca: 40 first PCs of PCA used
    ## $n.da: 2 discriminant functions saved
    ## $var (proportion of conserved variance): 0.402
    ## 
    ## $eig (eigenvalues): 1057 714.7  vector    length content                   
    ## 1 $eig      2      eigenvalues               
    ## 2 $grp      934    prior group assignment    
    ## 3 $prior    3      prior group probabilities 
    ## 4 $assign   934    posterior group assignment
    ## 5 $pca.cent 264    centring vector of PCA    
    ## 6 $pca.norm 264    scaling vector of PCA     
    ## 7 $pca.eig  261    eigenvalues of PCA        
    ## 
    ##   data.frame    nrow ncol content                                          
    ## 1 $tab          934  40   retained PCs of PCA                              
    ## 2 $means        3    40   group means                                      
    ## 3 $loadings     40   2    loadings of variables                            
    ## 4 $ind.coord    934  2    coordinates of individuals (principal components)
    ## 5 $grp.coord    3    2    coordinates of groups                            
    ## 6 $posterior    934  3    posterior membership probabilities               
    ## 7 $pca.loadings 264  40   PCA loadings of original variables               
    ## 8 $var.contr    264  2    contribution of original variables

``` r
dapc_post_4
```

    ##  #################################################
    ##  # Discriminant Analysis of Principal Components #
    ##  #################################################
    ## class: dapc
    ## $call: dapc.genlight(x = genlit_post_location, pop = grp_post_4$grp, 
    ##     n.pca = 40, n.da = 10)
    ## 
    ## $n.pca: 40 first PCs of PCA used
    ## $n.da: 3 discriminant functions saved
    ## $var (proportion of conserved variance): 0.402
    ## 
    ## $eig (eigenvalues): 767.8 530 410.3  vector    length content                   
    ## 1 $eig      3      eigenvalues               
    ## 2 $grp      934    prior group assignment    
    ## 3 $prior    4      prior group probabilities 
    ## 4 $assign   934    posterior group assignment
    ## 5 $pca.cent 264    centring vector of PCA    
    ## 6 $pca.norm 264    scaling vector of PCA     
    ## 7 $pca.eig  261    eigenvalues of PCA        
    ## 
    ##   data.frame    nrow ncol content                                          
    ## 1 $tab          934  40   retained PCs of PCA                              
    ## 2 $means        4    40   group means                                      
    ## 3 $loadings     40   3    loadings of variables                            
    ## 4 $ind.coord    934  3    coordinates of individuals (principal components)
    ## 5 $grp.coord    4    3    coordinates of groups                            
    ## 6 $posterior    934  4    posterior membership probabilities               
    ## 7 $pca.loadings 264  40   PCA loadings of original variables               
    ## 8 $var.contr    264  3    contribution of original variables

``` r
dapc_post_5
```

    ##  #################################################
    ##  # Discriminant Analysis of Principal Components #
    ##  #################################################
    ## class: dapc
    ## $call: dapc.genlight(x = genlit_post_location, pop = grp_post_5$grp, 
    ##     n.pca = 40, n.da = 10)
    ## 
    ## $n.pca: 40 first PCs of PCA used
    ## $n.da: 4 discriminant functions saved
    ## $var (proportion of conserved variance): 0.402
    ## 
    ## $eig (eigenvalues): 589.7 413.8 284.6 226.7  vector    length content                   
    ## 1 $eig      4      eigenvalues               
    ## 2 $grp      934    prior group assignment    
    ## 3 $prior    5      prior group probabilities 
    ## 4 $assign   934    posterior group assignment
    ## 5 $pca.cent 264    centring vector of PCA    
    ## 6 $pca.norm 264    scaling vector of PCA     
    ## 7 $pca.eig  261    eigenvalues of PCA        
    ## 
    ##   data.frame    nrow ncol content                                          
    ## 1 $tab          934  40   retained PCs of PCA                              
    ## 2 $means        5    40   group means                                      
    ## 3 $loadings     40   4    loadings of variables                            
    ## 4 $ind.coord    934  4    coordinates of individuals (principal components)
    ## 5 $grp.coord    5    4    coordinates of groups                            
    ## 6 $posterior    934  5    posterior membership probabilities               
    ## 7 $pca.loadings 264  40   PCA loadings of original variables               
    ## 8 $var.contr    264  4    contribution of original variables

plot DAPC using scatter()

``` r
scatter(dapc_post_2, col = viridis(2, alpha = 0.6))
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

``` r
scatter(dapc_post_3, col = viridis(3, alpha = 0.6))
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-34-2.png)<!-- -->

``` r
scatter(dapc_post_4, col = viridis(4, alpha = 0.6))
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-34-3.png)<!-- -->

``` r
scatter(dapc_post_5, col = viridis(5, alpha = 0.6))
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-34-4.png)<!-- -->

plot DAPC using ggplot

``` r
dapc_post_3_scores <- as.data.frame(dapc_post_3$tab)
dapc_post_3_scores$Assigned_Pop <- dapc_post_3$grp
dapc_post_3_scores$Original_Pop <- meta_post$Location

set.seed(9)
dapc_d3 <- ggplot(dapc_post_3_scores, aes(x=PC1, y=PC2, colour=Assigned_Pop)) 
dapc_d3 <- dapc_d3 + geom_point(size=2, aes(shape = Original_Pop)) 
dapc_d3 <- dapc_d3 + stat_ellipse(level = 0.95, size = 1)
dapc_d3 <- dapc_d3 + scale_color_manual(values = viridis(4)) 
dapc_d3 <- dapc_d3 + geom_hline(yintercept = 0) 
dapc_d3 <- dapc_d3 + geom_vline(xintercept = 0) 
dapc_d3 <- dapc_d3 + theme_bw()
dapc_d3
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

plot DAPC using ggplot

``` r
dapc_post_4_scores <- as.data.frame(dapc_post_4$tab)
dapc_post_4_scores$Assigned_Pop <- dapc_post_4$grp
dapc_post_4_scores$Original_Pop <- meta_post$Location

set.seed(9)
dapc_d4 <- ggplot(dapc_post_4_scores, aes(x=PC1, y=PC2, colour=Assigned_Pop)) 
dapc_d4 <- dapc_d4 + geom_point(size=2, aes(shape = Original_Pop)) 
dapc_d4 <- dapc_d4 + stat_ellipse(level = 0.95, size = 1)
dapc_d4 <- dapc_d4 + scale_color_manual(values = viridis(4)) 
dapc_d4 <- dapc_d4 + geom_hline(yintercept = 0) 
dapc_d4 <- dapc_d4 + geom_vline(xintercept = 0) 
dapc_d4 <- dapc_d4 + theme_bw()
dapc_d4
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

genotype compostion plot for 2 to 5 clusters

``` r
compoplot(dapc_post_2, posi="bottomright",
          txt.leg=paste("Cluster", 1:2), xlab="individuals", col=viridis(2))
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

``` r
compoplot(dapc_post_3, posi="bottomright",
          txt.leg=paste("Cluster", 1:3), xlab="individuals", col=viridis(3))
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-37-2.png)<!-- -->

``` r
compoplot(dapc_post_4, posi="bottomright",
          txt.leg=paste("Cluster", 1:4), xlab="individuals", col=viridis(4))
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-37-3.png)<!-- -->

``` r
compoplot(dapc_post_5, posi="bottomright",
          txt.leg=paste("Cluster", 1:5), xlab="individuals", col=viridis(5))
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-37-4.png)<!-- -->

make easier to interpret membership probability plot using ggplot

``` r
dapc.results.d3 <- as.data.frame(dapc_post_3$posterior)
dapc.results.d3$pop <- meta_post$Location
dapc.results.d3$indNames <- rownames(dapc.results.d3)
dapc.results.d3  <- pivot_longer(dapc.results.d3, -c(pop, indNames))  #pivot table to be able to use ggplot 
```

rename columns and plot

``` r
colnames(dapc.results.d3) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

d3 <- ggplot(dapc.results.d3, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
d3 <- d3 + geom_bar(stat='identity') 
d3 <- d3 + scale_fill_manual(values = viridis(4)) 
d3 <- d3 + facet_grid(~Original_Pop, scales = "free")
d3 <- d3 + theme(axis.text.x = element_blank())
d3
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

make easier to interpret membership probability plot using ggplot

``` r
dapc.results.d4 <- as.data.frame(dapc_post_4$posterior)
dapc.results.d4$pop <- meta_post$Location
dapc.results.d4$indNames <- rownames(dapc.results.d4)
dapc.results.d4  <- pivot_longer(dapc.results.d4, -c(pop, indNames))  #pivot table to be able to use ggplot 
```

rename columns and plot

``` r
colnames(dapc.results.d4) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

d4 <- ggplot(dapc.results.d4, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
d4 <- d4 + geom_bar(stat='identity') 
d4 <- d4 + scale_fill_manual(values = viridis(4)) 
d4 <- d4 + facet_grid(~Original_Pop, scales = "free")
d4 <- d4 + theme(axis.text.x = element_blank())
d4
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

make combined figure for post dam samples …

``` r
library(cowplot)

post_plot <- plot_grid(dapc_d3, d3, dapc_d4, d4, labels = "AUTO", ncol =2)
post_plot
```

![](steelhead_PCA_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

``` r
ggsave2("outputs/post_plot.png", width = 14, height = 8)
```
