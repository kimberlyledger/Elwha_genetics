steelhead\_adaptive\_analyses
================
Jong Yoon Jeon
Mar 14 2022

\#\#Load packages

``` r
library(vcfR)
```

    ## 
    ##    *****       ***   vcfR   ***       *****
    ##    This is vcfR 1.12.0 
    ##      browseVignettes('vcfR') # Documentation
    ##      citation('vcfR') # Citation
    ##    *****       *****      *****       *****

``` r
library(adegenet)
```

    ## Loading required package: ade4

    ## 
    ##    /// adegenet 2.1.5 is loaded ////////////
    ## 
    ##    > overview: '?adegenet'
    ##    > tutorials/doc/questions: 'adegenetWeb()' 
    ##    > bug reports/feature requests: adegenetIssues()

``` r
library(pegas)
```

    ## Loading required package: ape

    ## Registered S3 method overwritten by 'pegas':
    ##   method      from
    ##   print.amova ade4

    ## 
    ## Attaching package: 'pegas'

    ## The following object is masked from 'package:ape':
    ## 
    ##     mst

    ## The following object is masked from 'package:ade4':
    ## 
    ##     amova

    ## The following objects are masked from 'package:vcfR':
    ## 
    ##     getINFO, write.vcf

``` r
library(poppr)
```

    ## Warning: multiple methods tables found for 'direction'

    ## Warning: multiple methods tables found for 'gridDistance'

    ## This is poppr version 2.9.3. To get started, type package?poppr
    ## OMP parallel support: available

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(readr)
library(readxl)
library(tibble)
library(dartR)
```

    ## Loading required package: ggplot2

    ## Registered S3 method overwritten by 'GGally':
    ##   method from   
    ##   +.gg   ggplot2

    ## Registered S3 method overwritten by 'genetics':
    ##   method      from 
    ##   [.haplotype pegas

    ## **** Welcome to dartR ****

    ## Be aware that owing to CRAN requirements and compatibility reasons not all functions of the packages may run yet, as some dependencies could be missing. Hence for a most enjoyable experience we recommend to run the function

    ## gl.install.vanilla.dartR()

    ## This installs all missing and required packages for your version of dartR. 
    ## For citation information please use:

    ## citation('dartR')

    ## 
    ## **** Have fun using dartR! ****

\#\#Load data

``` r
steelhead.metadata <- read_csv("C:/Users/jyj55/Desktop/LGC-DGS_2022/Team_project/Steelhead_genotype/Elwha_Steelhead_Formatted.csv")
```

    ## Rows: 1693 Columns: 15
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: ","
    ## chr (10): Sample_ID, Smolt, NvH_Origin, Sex, Date, Time, Location, Run_Timin...
    ## dbl  (5): Year, Fork_Length, Lat, Long, rkm
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
steelhead.vcf <- read.vcfR("C:/Users/jyj55/Desktop/LGC-DGS_2022/Team_project/Steelhead_genotype/Elwha_GTSeq_Sans_CCT.vcf.gz")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 35
    ##   header_line: 36
    ##   variant count: 336
    ##   column count: 1178
    ## Meta line 35 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 336
    ##   Character matrix gt cols: 1178
    ##   skip: 0
    ##   nrows: 336
    ##   row_num: 0
    ## Processed variant: 336
    ## All variants processed

``` r
steelhead.genind <- vcfR2genind(steelhead.vcf, sep = "/")
steelhead.gt.vcf <- as.data.frame(t(extract.gt(steelhead.vcf, return.alleles = TRUE))) #t for transpose rows and columns
colnames(steelhead.gt.vcf) <- gsub("\\.", "_", colnames(steelhead.gt.vcf)) # Replace ".' to "_" in loci names
#steelhead.gt.genind <- df2genind(steelhead.gt.vcf, sep = "/")
steelhead.snp <- read_csv("C:/Users/jyj55/Desktop/LGC-DGS_2022/Team_project/Steelhead_genotype/SNP_Coordinates_CRITFC.csv")
```

    ## Rows: 367 Columns: 7
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: ","
    ## chr (4): Locus, chromosome, Scaffold, SNP
    ## dbl (3): snp coordinate in genome, Other, Physical_Position
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Steelhead.Locus_Key <- read_excel("C:/Users/jyj55/Desktop/LGC-DGS_2022/Team_project/Steelhead_genotype/Steelhead_Locus_Key.xlsx")
```

\#\#Make a dataframe to store population info

``` r
steelhead_pop <- matrix(NA, nrow=nrow(steelhead.genind@tab), ncol=8)
steelhead_pop <- as.data.frame(steelhead_pop)
names(steelhead_pop) <- c("Sample_ID", "NvH", "Time", "Location", "Run_Timing", "Life_History_Trait", "Lat", "Long")
```

\#\#Store pop info for each individual

``` r
for (i in 1:nrow(steelhead.genind@tab)){
  steelhead_pop$Sample_ID[i] <- rownames(steelhead.genind@tab)[i]
  steelhead_pop$NvH[i] <- steelhead.metadata %>% filter(Sample_ID == rownames(steelhead.genind@tab)[i]) %>% select(NvH_Origin)
  steelhead_pop$Time[i] <- steelhead.metadata %>% filter(Sample_ID == rownames(steelhead.genind@tab)[i]) %>% select(Time)
  steelhead_pop$Location[i] <- steelhead.metadata %>% filter(Sample_ID == rownames(steelhead.genind@tab)[i]) %>% select(Location)
  steelhead_pop$Run_Timing[i] <- steelhead.metadata %>% filter(Sample_ID == rownames(steelhead.genind@tab)[i]) %>% select(Run_Timing)
  steelhead_pop$Life_History_Type[i] <- steelhead.metadata %>% filter(Sample_ID == rownames(steelhead.genind@tab)[i]) %>% select(Life_History_Type)
  steelhead_pop$X[i] <- steelhead.metadata %>% filter(Sample_ID == rownames(steelhead.genind@tab)[i]) %>% select(Long)
  steelhead_pop$Y[i] <- steelhead.metadata %>% filter(Sample_ID == rownames(steelhead.genind@tab)[i]) %>% select(Lat)
}
invisible(steelhead_pop$NvH[lengths(steelhead_pop$NvH) == 0] <- NA_character_)
invisible(steelhead_pop$Time[lengths(steelhead_pop$Time) == 0] <- NA_character_)
invisible(steelhead_pop$Location[lengths(steelhead_pop$Location) == 0] <- NA_character_)
invisible(steelhead_pop$Run_Timing[lengths(steelhead_pop$Run_Timing) == 0] <- NA_character_)
invisible(steelhead_pop$Life_History_Type[lengths(steelhead_pop$Life_History_Type) == 0] <- NA_character_)
invisible(steelhead_pop$X[lengths(steelhead_pop$X) == 0] <- NA_character_)
invisible(steelhead_pop$Y[lengths(steelhead_pop$Y) == 0] <- NA_character_)
steelhead_coord <- cbind(steelhead_pop$Sample_ID, steelhead_pop$X, steelhead_pop$Y)
steelhead.genind@other$xy <- steelhead_coord
#steelhead.gt.genind@other$xy <- steelhead_coord
```

\#\#Divide pop by “Time” and “Location”

``` r
strata(steelhead.genind) <- steelhead_pop
setPop(steelhead.genind) <- ~Time/Location
steelhead.genind <- steelhead.genind[!is.na(steelhead.genind@strata$Time) & !is.na(steelhead.genind@strata$Location)] #Remove individuals without pop info
#strata(steelhead.gt.genind) <- steelhead_pop
#setPop(steelhead.gt.genind) <- ~Time/Location
#steelhead.gt.genind <- steelhead.gt.genind[!is.na(steelhead.gt.genind@strata$Time) & !is.na(steelhead.gt.genind@strata$Location)] #Remove individuals without pop info
```

\#\#Filter adaptive loci only, removing “character(0)” loci

``` r
steelhead_adaptive <- Steelhead.Locus_Key %>% filter(grepl('Adaptive', Steelhead.Locus_Key$`SNPeff Annotation output`))
for (i in 1:nrow(steelhead_adaptive)){
  steelhead_adaptive$SnpPos[i] <- steelhead.snp %>% filter(Locus == steelhead_adaptive[i,]$`SNPPIT or Alias`) %>% select(SNP)
}
```

    ## Warning: Unknown or uninitialised column: `SnpPos`.

``` r
toRetain <- steelhead_adaptive$SnpPos[lengths(steelhead_adaptive$SnpPos) != 0] #Remove loci without name (character(0))
toRetain <- gsub("\\.", "_", unlist(toRetain)) # Replace ".' to "_"
steelhead_adaptive.genind <- steelhead.genind[loc=unlist(toRetain)]
```

    ## Warning: the following specified loci do not exist: NC_035086_1_10773803,
    ## NC_035098_1_40502475

``` r
steelhead_adaptive.gt.vcf <- steelhead.gt.vcf[, names(steelhead.gt.vcf) %in% toRetain] #Retain adaptive loci only
#steelhead_adaptive.genind locus name change from SnpPos to SnpName? - manually after allele frequency calculation
```

\#\#HWE test for each pop (by “Time” or by “Location”)

``` r
steelhead_pop <- steelhead_pop %>% filter(!is.na(steelhead_pop$Time) & !is.na(steelhead_pop$Location)) #Remove individuals without pop info
strata(steelhead_adaptive.genind) <- steelhead_pop

setPop(steelhead_adaptive.genind) <- ~Time #by "Time"
steelhead_time_hwt <- seppop(steelhead_adaptive.genind) %>% lapply(hw.test, B = 0)
write.table(steelhead_time_hwt, file = "steelhead_time_hwt.txt", sep = "\t")

setPop(steelhead_adaptive.genind) <- ~Location #by "Location"
steelhead_location_hwt <- seppop(steelhead_adaptive.genind) %>% lapply(hw.test, B = 0)
write.table(steelhead_location_hwt, file = "steelhead_location_hwt.txt", sep = "\t")
```

\#\#Transform to genpop object and calculate allele frequencies

``` r
#Assigne population by "Time" + "Location"
setPop(steelhead_adaptive.genind) <- ~Time/Location
steelhead_adaptive.genpop <- genind2genpop(steelhead_adaptive.genind)
```

    ## 
    ##  Converting data from a genind to a genpop object... 
    ## 
    ## ...done.

``` r
steelhead_af <- makefreq(steelhead_adaptive.genpop, missing = NA)
```

    ## 
    ##  Finding allelic frequencies from a genpop object... 
    ## 
    ## ...done.

\#\#PCA using only adaptive loci

``` r
library(ggplot2)
steelhead_adaptive.genlight <- gi2gl(steelhead_adaptive.genind)
```

    ## Starting gi2gl 
    ## Completed: gi2gl

``` r
steelhead_adaptive.pca <- glPca(steelhead_adaptive.genlight, center = F, scale = F, nf =50)
steelhead_adaptive.PCA <- as.data.frame(steelhead_adaptive.pca$scores)
steelhead_pop$Location <- factor(steelhead_pop$Location, levels = c("AD", "ID", "BD", "SBLR")) 
steelhead_adaptive.PCA$Color <- steelhead_pop$Location
steelhead_adaptive.PCA$Color <- unlist(steelhead_adaptive.PCA$Color)
steelhead_pop$Time <- factor(steelhead_pop$Time, levels = c("Pre", "During", "Post")) 
steelhead_adaptive.PCA$Shape <- steelhead_pop$Time
steelhead_adaptive.PCA$Shape <- unlist(steelhead_adaptive.PCA$Shape)
steelhead_adaptive.PCA %>% ggplot(aes(x=PC1, y = PC2)) + geom_point(aes(col = Color, shape = Shape), size = 3, alpha=0.5) + 
  labs(color = "Location", shape = "Time")+
  theme_classic() + 
  theme(axis.text = element_text(size = 10, face = "bold"), 
        axis.title = element_text(size = 13, face = "bold"), 
        axis.line = element_line(size = 1.2),
        axis.title.y = element_text(angle = 90),
        legend.position = "bottom", 
        legend.box = "horizontal", 
        legend.title = element_text(face = "bold", size = 16), 
        legend.text = element_text(size = 13, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"))
```

![](steelhead_adaptive_analyses_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

\#\#Plot allele frequency heatmaps (\*check summer-run/winter-run
alleles and modify)

``` r
library(RColorBrewer)
library(gplots)
```

    ## Registered S3 method overwritten by 'gplots':
    ##   method         from 
    ##   reorder.factor gdata

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

``` r
#Change row order as following: Pre_AD, Pre_ID, Pre_BD, During_BD, Post_AD, Post_ID, Post_BD, Post_SBLR
steelhead_af_ordered <- steelhead_af[c(2, 7, 1, 6, 3, 4, 5, 8),]
write.table(steelhead_af_ordered, "steelhead_af.txt", sep = "\t")
steelhead_af_pruned <- read.delim("C:/Users/jyj55/Desktop/LGC-DGS_2022/Team_project/Analysis/Steelhead analysis/steelhead_af_pruned.txt", row.names=1) #Reread manually pruned table to have major allele frequency only for each locus 
steelhead_af_pruned <- as.matrix(steelhead_af_pruned)
steelhead_af_pruned <- steelhead_af_pruned[,order(colnames(steelhead_af_pruned))]
heatmap_color <- colorRampPalette(brewer.pal(10,"RdYlBu"))
steelhead.heatmap_all <- heatmap.2(steelhead_af_pruned, col=heatmap_color, Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace = "none", density.info = "none", xlab = "Locus", labCol = "", margins = c(2.5,10))
```

![](steelhead_adaptive_analyses_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
steelhead.heatmap_all_dendro <- heatmap.2(steelhead_af_pruned, col=heatmap_color, trace = "none", density.info = "none", xlab = "Locus", labCol = "", margins = c(2.5,10)) #Heatmap with dendrogram option
```

![](steelhead_adaptive_analyses_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
steelhead_af_005 <- read.delim("C:/Users/jyj55/Desktop/LGC-DGS_2022/Team_project/Analysis/Steelhead analysis/steelhead_af_pruned_005.txt", row.names=1) #Manually filtered loci to include higher 5% of them in terms of absolute value of allele frequency change between PreDam and PostDam for any of AD, ID, BD
steelhead_af_005 <- as.matrix(steelhead_af_005)
steelhead.heatmap_005 <- heatmap.2(steelhead_af_005, col=heatmap_color, Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace = "none", density.info = "none", srtCol = 45, margins = c(10,10))
```

![](steelhead_adaptive_analyses_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->

``` r
steelhead.heatmap_005_dendro <- heatmap.2(steelhead_af_005, col=heatmap_color, trace = "none", density.info = "none", srtCol = 45, margins = c(10,10)) #Heatmap with dendrogram option
```

![](steelhead_adaptive_analyses_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->

\#\#Calculate (observed) genotype frequencies for each pop
(Time/Location) and save them to manually re-arrange them

``` r
library(mixIndependR)
steelhead_pop$Time_Location <- paste0(steelhead_pop$Time,"_",steelhead_pop$Location)
steelhead_adaptive.gt.vcf <- as_tibble(steelhead_adaptive.gt.vcf, rownames = "Sample_ID")
steelhead_adaptive.gt.vcf <- left_join(steelhead_adaptive.gt.vcf, steelhead_pop %>% select(Sample_ID, Time_Location), by = c("Sample_ID" = "Sample_ID"))
steelhead_adaptive.gt.vcf <- steelhead_adaptive.gt.vcf %>% filter(!is.na(steelhead_adaptive.gt.vcf$Time_Location)) #Remove individuals without pop info
steelhead_adaptive$SnpPos <- gsub("\\.", "_", steelhead_adaptive$SnpPos)
steelhead_adaptiveSNP.gt.vcf <- steelhead_adaptive.gt.vcf

#Change loci names from SnpPos to SnpName (="SNPPIT or Alias")
for (i in 1:ncol(steelhead_adaptiveSNP.gt.vcf[-c(1,length(steelhead_adaptiveSNP.gt.vcf))])){
  colnames(steelhead_adaptiveSNP.gt.vcf)[i+1] <- as.character(steelhead_adaptive[(which(steelhead_adaptive$SnpPos == colnames(steelhead_adaptiveSNP.gt.vcf)[i+1])),3]) #3rd column of "steelhead_adaptive" includes SnpName
}

steelhead_adaptiveSNP.gt.vcf_PreAD <- steelhead_adaptiveSNP.gt.vcf[steelhead_adaptiveSNP.gt.vcf$Time_Location == "Pre_AD",]
steelhead_adaptiveSNP.gt.vcf_PostAD <- steelhead_adaptiveSNP.gt.vcf[steelhead_adaptiveSNP.gt.vcf$Time_Location == "Post_AD",]
steelhead_adaptiveSNP.gt.vcf_PreID <- steelhead_adaptiveSNP.gt.vcf[steelhead_adaptiveSNP.gt.vcf$Time_Location == "Pre_ID",]
steelhead_adaptiveSNP.gt.vcf_PostID <- steelhead_adaptiveSNP.gt.vcf[steelhead_adaptiveSNP.gt.vcf$Time_Location == "Post_ID",]
steelhead_adaptiveSNP.gt.vcf_PreBD <- steelhead_adaptiveSNP.gt.vcf[steelhead_adaptiveSNP.gt.vcf$Time_Location == "Pre_BD",]
steelhead_adaptiveSNP.gt.vcf_DuringBD <- steelhead_adaptiveSNP.gt.vcf[steelhead_adaptiveSNP.gt.vcf$Time_Location == "During_BD",]
steelhead_adaptiveSNP.gt.vcf_PostBD <- steelhead_adaptiveSNP.gt.vcf[steelhead_adaptiveSNP.gt.vcf$Time_Location == "Post_BD",]
steelhead_adaptiveSNP.gt.vcf_PostSBLR <- steelhead_adaptiveSNP.gt.vcf[steelhead_adaptiveSNP.gt.vcf$Time_Location == "Post_SBLR",]

steelhead_gf_PreAD <- GenotypeFreq(steelhead_adaptiveSNP.gt.vcf_PreAD[,2:(ncol(steelhead_adaptiveSNP.gt.vcf_PreAD)-1)],"/",expect=FALSE) #expect=FALSE to calculate "observed" genotype frequencies
steelhead_gf_PreAD[] <- sapply(steelhead_gf_PreAD, function(x){x/sum(x)}) #Make genotype frequency as "frequency"; default is observed numbers"
steelhead_gf_PreAD$Time_Location <- "Pre_AD"
steelhead_gf_PostAD <- GenotypeFreq(steelhead_adaptiveSNP.gt.vcf_PostAD[,2:(ncol(steelhead_adaptiveSNP.gt.vcf_PostAD)-1)],"/",expect=FALSE)
steelhead_gf_PostAD[] <- sapply(steelhead_gf_PostAD, function(x){x/sum(x)})
steelhead_gf_PostAD$Time_Location <- "Post_AD"
write.table(steelhead_gf_PreAD, file = "steelhead_gf_PreAD.txt", sep = "\t")
write.table(steelhead_gf_PostAD, file = "steelhead_gf_PostAD.txt", sep = "\t")

steelhead_gf_PreID <- GenotypeFreq(steelhead_adaptiveSNP.gt.vcf_PreID[,2:(ncol(steelhead_adaptiveSNP.gt.vcf_PreID)-1)],"/",expect=FALSE)
steelhead_gf_PreID[] <- sapply(steelhead_gf_PreID, function(x){x/sum(x)})
steelhead_gf_PreID$Time_Location <- "Pre_ID"
steelhead_gf_PostID <- GenotypeFreq(steelhead_adaptiveSNP.gt.vcf_PostID[,2:(ncol(steelhead_adaptiveSNP.gt.vcf_PostID)-1)],"/",expect=FALSE)
steelhead_gf_PostID[] <- sapply(steelhead_gf_PostID, function(x){x/sum(x)})
steelhead_gf_PostID$Time_Location <- "Post_ID"
write.table(steelhead_gf_PreID, file = "steelhead_gf_PreID.txt", sep = "\t")
write.table(steelhead_gf_PostID, file = "steelhead_gf_PostID.txt", sep = "\t")

steelhead_gf_PreBD <- GenotypeFreq(steelhead_adaptiveSNP.gt.vcf_PreBD[,2:(ncol(steelhead_adaptiveSNP.gt.vcf_PreBD)-1)],"/",expect=FALSE)
steelhead_gf_PreBD[] <- sapply(steelhead_gf_PreBD, function(x){x/sum(x)})
steelhead_gf_PreBD$Time_Location <- "Pre_BD"
steelhead_gf_DuringBD <- GenotypeFreq(steelhead_adaptiveSNP.gt.vcf_DuringBD[,2:(ncol(steelhead_adaptiveSNP.gt.vcf_DuringBD)-1)],"/",expect=FALSE)
steelhead_gf_DuringBD[] <- sapply(steelhead_gf_DuringBD, function(x){x/sum(x)})
steelhead_gf_DuringBD$Time_Location <- "During_BD"
steelhead_gf_PostBD <- GenotypeFreq(steelhead_adaptiveSNP.gt.vcf_PostBD[,2:(ncol(steelhead_adaptiveSNP.gt.vcf_PostBD)-1)],"/",expect=FALSE)
steelhead_gf_PostBD[] <- sapply(steelhead_gf_PostBD, function(x){x/sum(x)})
steelhead_gf_PostBD$Time_Location <- "Post_BD"
write.table(steelhead_gf_PreBD, file = "steelhead_gf_PreBD.txt", sep = "\t")
write.table(steelhead_gf_DuringBD, file = "steelhead_gf_DuringBD.txt", sep = "\t")
write.table(steelhead_gf_PostBD, file = "steelhead_gf_PostBD.txt", sep = "\t")

steelhead_gf_PostSBLR <- GenotypeFreq(steelhead_adaptiveSNP.gt.vcf_PostSBLR[,2:(ncol(steelhead_adaptiveSNP.gt.vcf_PostSBLR)-1)],"/",expect=FALSE)
steelhead_gf_PostSBLR[] <- sapply(steelhead_gf_PostSBLR, function(x){x/sum(x)})
steelhead_gf_PostSBLR$Time_Location <- "Post_SBLR"
write.table(steelhead_gf_PostSBLR, file = "steelhead_gf_PostSBLR.txt", sep = "\t")

#After this, in Excel, manually separate loci according to positions on Greb1 (Premature/Mature and Omy5 (Ancestral/Recombined) and combine files of same Location (e.g. steelhead_gf_BD; combined PreBD + DuringBD + PostBD), because "GenotypeFreq" function removes Sample ID. But you can just make combined files in advacne here too, as they retain Pop info.
```

\#\#Plot allele frequency heatmaps (\*check summer-run/winter-run
alleles and modify)

``` r
#Read files of GREB1L loci genotype frequencies for each run timing (summer-run, heterozygotes, winter-run)
steelhead_gf_summerrun <- read_excel("C:/Users/jyj55/Desktop/LGC-DGS_2022/Team_project/Analysis/Steelhead analysis/Elwha_genetics/steelhead_gf_summerrun.xlsx")
```

    ## New names:
    ## * `` -> ...1

``` r
steelhead_gf_summerrun <- as.matrix(steelhead_gf_summerrun)
Popinfo <- steelhead_gf_summerrun[,1]
steelhead_gf_summerrun <- steelhead_gf_summerrun[,-1]
steelhead_gf_summerrun <- apply(steelhead_gf_summerrun, 2, as.numeric)
rownames(steelhead_gf_summerrun) <- Popinfo

heatmap_color <- colorRampPalette(brewer.pal(10,"RdYlBu"))

steelhead.heatmap_summerrun <- heatmap.2(steelhead_gf_summerrun, col=heatmap_color, Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace = "none", srtCol = 45, density.info = "none", margins = c(10,10))
```

![](steelhead_adaptive_analyses_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
steelhead.heatmap_summerrun_dendro <- heatmap.2(steelhead_gf_summerrun, col=heatmap_color, trace = "none", density.info = "none", srtCol = 45, margins = c(10,10)) #Heatmap with dendrogram option
```

![](steelhead_adaptive_analyses_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

``` r
steelhead_gf_heterozygote <- read_excel("C:/Users/jyj55/Desktop/LGC-DGS_2022/Team_project/Analysis/Steelhead analysis/Elwha_genetics/steelhead_gf_heterozygote.xlsx")
```

    ## New names:
    ## * `` -> ...1

``` r
steelhead_gf_heterozygote <- as.matrix(steelhead_gf_heterozygote)
steelhead_gf_heterozygote <- steelhead_gf_heterozygote[,-1]
steelhead_gf_heterozygote <- apply(steelhead_gf_heterozygote, 2, as.numeric)
rownames(steelhead_gf_heterozygote) <- Popinfo

steelhead.heatmap_heterozygote <- heatmap.2(steelhead_gf_heterozygote, col=heatmap_color, Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace = "none", srtCol = 45, density.info = "none", margins = c(10,10))
```

![](steelhead_adaptive_analyses_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->

``` r
steelhead.heatmap_heterozygote_dendro <- heatmap.2(steelhead_gf_heterozygote, col=heatmap_color, trace = "none", density.info = "none", srtCol = 45, margins = c(10,10)) #Heatmap with dendrogram option
```

![](steelhead_adaptive_analyses_files/figure-gfm/unnamed-chunk-12-4.png)<!-- -->

``` r
steelhead_gf_winterrun <- read_excel("C:/Users/jyj55/Desktop/LGC-DGS_2022/Team_project/Analysis/Steelhead analysis/Elwha_genetics/steelhead_gf_winterrun.xlsx")
```

    ## New names:
    ## * `` -> ...1

``` r
steelhead_gf_winterrun <- as.matrix(steelhead_gf_winterrun)
steelhead_gf_winterrun <- steelhead_gf_winterrun[,-1]
steelhead_gf_winterrun <- apply(steelhead_gf_winterrun, 2, as.numeric)
rownames(steelhead_gf_winterrun) <- Popinfo

steelhead.heatmap_winterrun <- heatmap.2(steelhead_gf_winterrun, col=heatmap_color, Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace = "none", srtCol = 45, density.info = "none", margins = c(10,10))
```

![](steelhead_adaptive_analyses_files/figure-gfm/unnamed-chunk-12-5.png)<!-- -->

``` r
steelhead.heatmap_winterrun_dendro <- heatmap.2(steelhead_gf_winterrun, col=heatmap_color, trace = "none", density.info = "none", srtCol = 45, margins = c(10,10)) #Heatmap with dendrogram option
```

![](steelhead_adaptive_analyses_files/figure-gfm/unnamed-chunk-12-6.png)<!-- -->