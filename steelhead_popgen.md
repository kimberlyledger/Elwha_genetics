Elwha Steelhead Population Genetics
================
Kimberly Ledger
8 March 2022

This code will look at neutral population genetics of steelhead/rainbow
trout in the Elwha River, Washington - this dataset includes individuals
sampled before and after two dams were removed (lower dam removal
completed in 2012; upper dam removal completed in 2015)

load libraries

``` r
library(vcfR) #this package is used to visualize and manipulate VCF files
library(adegenet) #this package is used for analysis of genetic/genomic data 
library(dplyr) # data manipulation
library(tidyr) # data manipulation
library(pegas) #a package for pop gen data analysis
library(poppr) #a package for pop gen data analysis
```

# Part 1: Prepare the data

## read in the steelhead (Oncorhynchus mykiss) metadata

``` r
onmy_metadata <- read.csv("~/Desktop/LG_Proj4/Elwha_datafiles/Elwha_Steelhead_Formatted.csv")
head(onmy_metadata)
```

    ##   Sample_ID Year   Smolt Fork_Length NvH_Origin  Sex    Date Time Location
    ## 1  33649_17 2004 Unknown          NA          N <NA> 7/14/04  Pre       ID
    ## 2  33649_18 2004 Unknown          NA          N <NA> 7/14/04  Pre       ID
    ## 3  33649_19 2004 Unknown          NA          N <NA> 7/14/04  Pre       ID
    ## 4  33649_20 2004 Unknown          NA          N <NA> 7/14/04  Pre       ID
    ## 5  33649_23 2004 Unknown          NA          N <NA> 7/14/04  Pre       ID
    ## 6  33649_26 2004 Unknown          NA          N <NA> 7/14/04  Pre       ID
    ##   Run_Timing Life_Stage Life_History_Type      Lat      Long rkm Sampling_Site
    ## 1    Unknown   Juvenile       Land_Locked 48.11933 -123.5535  NA  little_river
    ## 2    Unknown   Juvenile       Land_Locked 48.06303 -123.5770  NA  little_river
    ## 3    Unknown   Juvenile       Land_Locked 48.06303 -123.5770  NA  little_river
    ## 4    Unknown   Juvenile       Land_Locked 48.06303 -123.5770  NA  little_river
    ## 5    Unknown   Juvenile       Land_Locked 48.06303 -123.5770  NA  little_river
    ## 6    Unknown   Juvenile       Land_Locked 48.06303 -123.5770  NA  little_river

### check out a quick summary tables of the metadata

``` r
onmy_metadata %>%
  group_by(Time, Location, Life_Stage) %>%
  summarize(total = n())
```

    ## # A tibble: 11 × 4
    ## # Groups:   Time, Location [9]
    ##    Time   Location Life_Stage total
    ##    <chr>  <chr>    <chr>      <int>
    ##  1 During BD       Adult         44
    ##  2 Post   AD       Adult         61
    ##  3 Post   AD       Juvenile     420
    ##  4 Post   BD       Adult        512
    ##  5 Post   BD       Juvenile       1
    ##  6 Post   ID       Adult         61
    ##  7 Post   SBLR     Juvenile      27
    ##  8 Pre    AD       Juvenile     208
    ##  9 Pre    BD       Juvenile     104
    ## 10 Pre    ID       Juvenile     169
    ## 11 Pre    SBLR     Juvenile      86

``` r
onmy_metadata %>%
  group_by(Time, Sampling_Site) %>%
  summarize(total = n()) %>%
  pivot_wider(names_from = "Time", values_from = "total")
```

    ## # A tibble: 31 × 4
    ##    Sampling_Site      During  Post   Pre
    ##    <chr>               <int> <int> <int>
    ##  1 elwha_river_lower      40    71   104
    ##  2 elwha_river_mouth       3     4    NA
    ##  3 lekt_outlet             1   320    NA
    ##  4 altaire                NA    14    32
    ##  5 boulder                NA     6    NA
    ##  6 chicago_camp           NA     6     7
    ##  7 ds_fishermans_bend     NA    10    NA
    ##  8 ds_ranger_station      NA     7    NA
    ##  9 elkhorn                NA   117    37
    ## 10 elwha_river_middle     NA    11    NA
    ## # … with 21 more rows

## read in the VCF file and convert to genind object

``` r
onmy_vcf <- read.vcfR("~/Desktop/LG_Proj4/Elwha_datafiles/Elwha_GTSeq_Sans_CCT.vcf")
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
onmy_genind <- vcfR2genind(onmy_vcf, sep = "/")
onmy_genind
```

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 1,169 individuals; 336 loci; 664 alleles; size: 3.3 Mb
    ## 
    ##  // Basic content
    ##    @tab:  1169 x 664 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 1-2)
    ##    @loc.fac: locus factor for the 664 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: adegenet::df2genind(X = t(x), sep = sep)
    ## 
    ##  // Optional content
    ##    - empty -

## add population information into the genind object - Location or Sampling_Site????

``` r
# create empty data frame 
onmy_pop <- matrix(NA, nrow=nrow(onmy_genind@tab), ncol=2)
onmy_pop <- as.data.frame(onmy_pop)
names(onmy_pop) <- c("Sample_ID", "pop")

# add population info for each individual to data frame
for (i in 1:nrow(onmy_genind@tab)){
  onmy_pop$Sample_ID[i] <- rownames(onmy_genind@tab)[i]
  onmy_pop$pop[i] <- onmy_metadata %>% filter(Sample_ID == rownames(onmy_genind@tab)[i]) %>% select(Location)
}

# add population info to the genind obj
strata(onmy_genind) <- onmy_pop
setPop(onmy_genind) <- ~pop

onmy_genind
```

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 1,169 individuals; 336 loci; 664 alleles; size: 3.4 Mb
    ## 
    ##  // Basic content
    ##    @tab:  1169 x 664 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 1-2)
    ##    @loc.fac: locus factor for the 664 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: adegenet::df2genind(X = t(x), sep = sep)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 8-568)
    ##    @strata: a data frame with 2 columns ( Sample_ID, pop )

# Part 2: Genetic Diversity

## global hardy-weinberg equilibrium

``` r
onmy_hw <- pegas::hw.test(onmy_genind, B=0)
head(onmy_hw)
```

    ##                             chi^2 df  Pr(chi^2 >)
    ## NC_035086_1_7508221    3.61939683  1 5.710957e-02
    ## NC_035086_1_10852282  27.31379350  1 1.729731e-07
    ## NC_035086_1_30372084 897.68130228  1 0.000000e+00
    ## NC_035086_1_32872704   2.27159271  1 1.317644e-01
    ## NC_035086_1_38704654   0.01658003  1 8.975448e-01
    ## NC_035086_1_40486618   0.05964253  1 8.070618e-01

## LD test using “poppr” package and “genclone” object

``` r
#onmy_genclone <- as.genclone(onmy_genind)
#LDpair <- onmy_genclone %>% pair.ia
#LDpair <- tibble::rownames_to_column(as.data.frame(LDpair), var = "pairID")
```
