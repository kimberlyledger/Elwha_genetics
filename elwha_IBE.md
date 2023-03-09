elwha_IBE
================
Kimberly Ledger
2022-11-21

## last update: 8 March 2023

## starting by just working with post-dam removal scores

load libraries

``` r
library(tidyverse)
library(ggplot2)

## check later what i'm actually using...
#library(nlme)
library(usdm)
#library(corMLPE)
library(lme4)
library(stargazer)
#library(MuMIn)
```

## steelhead

load data notes: added back in OUT_DIST manually

``` r
omy_neutral_scores <- read.csv("~/Desktop/LG_Proj4/Elwha_genetics/outputs/MAR2023/onmy_n_post_dapc.csv") %>%
  rename(Sample_ID = X) %>%
  dplyr::select(Sample_ID, PC1)

omy_adaptive_scores <- readRDS("~/Desktop/LG_Proj4/Elwha_datafiles/2023/IBE_input_Steelhead.rds") %>%
  #rename(Sample_ID = sampleID) %>%
  dplyr::select(Sample_ID, alleleS_count1)

enviro <- read.csv("~/Desktop/LG_Proj4/Elwha_environmentaldata/outputs/final/enviro_summary_reduced_22Feb2023.csv") %>%
  rename(Sampling_Site = X)

omy_metadata <- read.csv("~/Desktop/LG_Proj4/Elwha_datafiles/2023/Elwha_Steelhead_Formatted_Post_022023_sitesfixed_kjl.csv") %>%
  filter(Time == "Post")
```

join the environmental data to the list of omy collection sites

``` r
omy_sites <- omy_metadata %>%
  distinct(Sampling_Site)
#omy_sites

omy_sites_enviro <- omy_sites %>%
  left_join(enviro, by = "Sampling_Site")
```

create a data frame

``` r
omy_df <- omy_metadata %>%
  dplyr::select(Sample_ID, Life_Stage, Life_History_Type, Sampling_Site) %>%
  left_join(omy_sites_enviro, by = "Sampling_Site") %>%
  left_join(omy_neutral_scores, by = "Sample_ID") %>%
  left_join(omy_adaptive_scores, by = "Sample_ID") %>%
  arrange(OUT_DIST) %>%
  drop_na(PC1)
```

## run some mixed effect models

### sort of following worked examples:

<https://bookdown.org/hhwagner1/LandGenCourse_book/WE_6.html#WE_6>
<https://bookdown.org/hhwagner1/LandGenCourse_book/WE_12.html#WE_12>

### center and scale explainatory variables

``` r
sample <- omy_df$Sample_ID
PC1 <- omy_df$PC1
PropA <- omy_df$alleleS_count1
site <- omy_df$Sampling_Site
dist <- omy_df$OUT_DIST

for (i in 1:length(colnames(omy_df))){
  if (is.numeric(omy_df[, i])==TRUE)
    omy_df[, i] <- as.numeric(scale(omy_df[, i]))
  else
    omy_df[, i] <- omy_df[, i]
}

## replace values i don't want scaled
omy_df$PC1 <- PC1
omy_df$alleleS_count1 <- PropA
```

check data frame structure

``` r
omy_df$Sampling_Site <- as.factor(omy_df$Sampling_Site)
str(omy_df)
```

    ## 'data.frame':    802 obs. of  20 variables:
    ##  $ Sample_ID        : chr  "51761_17_005" "51761_17_008" "51686_15_92" "51686_15_93" ...
    ##  $ Life_Stage       : chr  "Adult" "Adult" "Adult" "Adult" ...
    ##  $ Life_History_Type: chr  "Steelhead" "Steelhead" "Steelhead" "Steelhead" ...
    ##  $ Sampling_Site    : Factor w/ 19 levels "aldwell","altaire",..: 9 9 15 15 15 15 15 15 15 15 ...
    ##  $ CANOPY           : num  -0.303 -0.303 -0.442 -0.442 -0.442 ...
    ##  $ SLOPE            : num  -1.4 -1.4 -1.4 -1.4 -1.4 ...
    ##  $ PRECIP           : num  -1.27 -1.27 -1.27 -1.27 -1.27 ...
    ##  $ FlowVel          : num  -0.97 -0.97 -0.938 -0.938 -0.938 ...
    ##  $ AWAT             : num  0.874 0.874 0.874 0.874 0.874 ...
    ##  $ IP_Steelhd       : num  -1.32 -1.32 -1.32 -1.32 -1.32 ...
    ##  $ IP_Chinook       : num  1.22 1.22 1.22 1.22 1.22 ...
    ##  $ Pool_freq        : num  -0.843 -0.843 -0.843 -0.843 -0.843 ...
    ##  $ Logjams          : num  0.774 0.774 0.774 0.774 0.774 ...
    ##  $ Spawnable        : num  -0.08 -0.08 -0.226 -0.226 -0.226 ...
    ##  $ fines            : num  -1.554 -1.554 0.946 0.946 0.946 ...
    ##  $ gravels          : num  -2.227 -2.227 0.823 0.823 0.823 ...
    ##  $ MWMT             : num  0.726 0.726 0.726 0.726 0.726 ...
    ##  $ OUT_DIST         : num  -0.973 -0.973 -0.948 -0.948 -0.948 ...
    ##  $ PC1              : num  -0.186 0.637 0.62 0.967 0.725 ...
    ##  $ alleleS_count1   : num  0.0556 0 0.4444 0 NA ...

response = PC1 (continuous), the loading score of PC1 from the DAPC
analysis using neutral loci OR alleleS_count1 (continuous), the
proportion of early run alleles explanatory = environmental vars
(centered and scaled) random effect = sampling site

## look at the distribution of the response variable

``` r
hist(omy_df$PC1)
```

![](elwha_IBE_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
hist(omy_df$alleleS_count1)
```

![](elwha_IBE_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

PC1 okay, but will need to use a different distribution for alleleS

## subset dataframe to removed environmental data missing from some sites.

``` r
omy_df_reduced <- omy_df %>%
  dplyr::select(!gravels) %>%
  dplyr::select(!fines)
```

## fit candidate models for datasets

1.  Null model: OUT_DIST (aka. Rkm)
2.  Climate model: MWMT (max weekly max temp) and PRECIP
3.  Flow model: FlowVel, SLOPE
4.  Habitat type: CANOPY, IP_Steelhd, Pool_freq, Logjams, spawnable
5.  Full model: CANOPY, SLOPE, PRECIP, FlowVel, MWMT, IP_Steelhd,
    Pool_freq, Logjams, Spawnable

## check for multicollinearity - climate model

``` r
df <- with(omy_df_reduced, data.frame(MWMT, PRECIP))
usdm::vif(df) 
```

    ##   Variables      VIF
    ## 1      MWMT 5.302078
    ## 2    PRECIP 5.302078

will proceed only with MWMT

## check for multicollinearity - flow model

``` r
df <- with(omy_df_reduced, data.frame(FlowVel, SLOPE))
usdm::vif(df) 
```

    ##   Variables      VIF
    ## 1   FlowVel 1.540041
    ## 2     SLOPE 1.540041

will keep both

## check for multicollinearity - habitat model

``` r
df <- with(omy_df_reduced, data.frame(CANOPY, IP_Steelhd, Logjams, Pool_freq, Spawnable))
usdm::vif(df) 
```

    ##    Variables       VIF
    ## 1     CANOPY  4.191717
    ## 2 IP_Steelhd  2.141861
    ## 3    Logjams 10.256301
    ## 4  Pool_freq  7.084923
    ## 5  Spawnable  1.346930

drop Logjams

``` r
df <- with(omy_df_reduced, data.frame(CANOPY, IP_Steelhd, Pool_freq, Spawnable))
usdm::vif(df) 
```

    ##    Variables      VIF
    ## 1     CANOPY 2.213114
    ## 2 IP_Steelhd 2.110795
    ## 3  Pool_freq 1.155398
    ## 4  Spawnable 1.143541

drop Logjams

## check for multicollinearity - full

``` r
df <- with(omy_df_reduced, data.frame(MWMT, FlowVel, CANOPY, IP_Steelhd, Pool_freq, Spawnable))
usdm::vif(df) 
```

    ##    Variables      VIF
    ## 1       MWMT 3.456305
    ## 2    FlowVel 8.794318
    ## 3     CANOPY 2.715157
    ## 4 IP_Steelhd 9.087993
    ## 5  Pool_freq 2.357833
    ## 6  Spawnable 1.298098

drop IP_Steelhd

``` r
df <- with(omy_df_reduced, data.frame(MWMT, FlowVel, CANOPY, Pool_freq, Spawnable))
usdm::vif(df) 
```

    ##   Variables      VIF
    ## 1      MWMT 2.454570
    ## 2   FlowVel 2.046027
    ## 3    CANOPY 2.671420
    ## 4 Pool_freq 2.178273
    ## 5 Spawnable 1.137516

good to go.

create subsets for (1) jv all, (2) adult all, (3) adults known steelhead

``` r
omy_jv <- omy_df_reduced %>%
  filter(Life_Stage == "Juvenile")

omy_ad <- omy_df_reduced %>%
  filter(Life_Stage == "Adult")
  
omy_ad_stlhd <- omy_ad %>%
  filter(Life_History_Type == "Steelhead")
```

# juveniles

## set up models - no random effects

``` r
mod1 <- lm(PC1 ~ OUT_DIST, data = omy_jv)
mod2 <- lm(PC1 ~ MWMT, data = omy_jv)
mod3 <- lm(PC1 ~ FlowVel + SLOPE, data = omy_jv)
mod4 <- lm(PC1 ~ CANOPY + IP_Steelhd +Pool_freq + Spawnable, data = omy_jv)
mod5 <- lm(PC1 ~ MWMT + FlowVel + CANOPY + Pool_freq + Spawnable, data = omy_jv)
```

## compare models

``` r
stargazer(mod1, mod2, mod3, mod4, mod5, 
          type = "text",
          digits = 3, 
          ci = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
```

    ## 
    ## =============================================================================================================================================
    ##                                                                        Dependent variable:                                                   
    ##                     -------------------------------------------------------------------------------------------------------------------------
    ##                                                                                PC1                                                           
    ##                               (1)                      (2)                      (3)                     (4)                     (5)          
    ## ---------------------------------------------------------------------------------------------------------------------------------------------
    ## OUT_DIST                   -0.517***                                                                                                         
    ##                         (-0.605, -0.430)                                                                                                     
    ##                                                                                                                                              
    ## MWMT                                                 0.270***                                                                  0.253         
    ##                                                   (0.218, 0.322)                                                          (-0.018, 0.524)    
    ##                                                                                                                                              
    ## FlowVel                                                                       -0.126                                          -0.365         
    ##                                                                           (-0.355, 0.104)                                 (-0.830, 0.101)    
    ##                                                                                                                                              
    ## SLOPE                                                                         -0.201*                                                        
    ##                                                                          (-0.358, -0.044)                                                    
    ##                                                                                                                                              
    ## CANOPY                                                                                               -0.341***                 0.049         
    ##                                                                                                  (-0.473, -0.209)         (-0.337, 0.434)    
    ##                                                                                                                                              
    ## IP_Steelhd                                                                                             0.088                                 
    ##                                                                                                   (-0.093, 0.269)                            
    ##                                                                                                                                              
    ## Pool_freq                                                                                              0.061                                 
    ##                                                                                                   (-0.018, 0.140)                            
    ##                                                                                                                                              
    ## Spawnable                                                                                                                                    
    ##                                                                                                                                              
    ##                                                                                                                                              
    ## Constant                     0.113                  -0.246***                -0.424***               -0.649***               -0.378***       
    ##                         (-0.007, 0.233)          (-0.324, -0.169)        (-0.611, -0.236)        (-0.759, -0.539)        (-0.559, -0.196)    
    ##                                                                                                                                              
    ## ---------------------------------------------------------------------------------------------------------------------------------------------
    ## Observations                  273                      273                      273                     273                     273          
    ## R2                           0.332                    0.276                    0.254                   0.406                   0.406         
    ## Adjusted R2                  0.330                    0.273                    0.248                   0.400                   0.400         
    ## Residual Std. Error     0.463 (df = 271)         0.482 (df = 271)        0.490 (df = 270)        0.438 (df = 269)        0.438 (df = 269)    
    ## F Statistic         134.957*** (df = 1; 271) 103.058*** (df = 1; 271) 45.916*** (df = 2; 270) 61.377*** (df = 3; 269) 61.377*** (df = 3; 269)
    ## =============================================================================================================================================
    ## Note:                                                                                                           *p<0.05; **p<0.01; ***p<0.001

week 6 tutorial suggests using the ML fit, as AIC is not valid for REML…
“The function extractAIC refits the models with ‘REML=FALSE’ to obtain
AIC values that are comparable between models with different fixed
effects or between models fitted with functions lm and lmer.”

``` r
aic_vals <- c(extractAIC(mod1)[2], extractAIC(mod2)[2], extractAIC(mod3)[2], 
              extractAIC(mod4)[2], extractAIC(mod5)[2])
names(aic_vals) <- c("mod1","mod2","mod3", "mod4", "mod5")
aic_vals
```

    ##      mod1      mod2      mod3      mod4      mod5 
    ## -418.9334 -396.5918 -386.5292 -446.9679 -446.9679

**hmm here mod4 (habitat model) and mod5 (full model) shows lowest AIC**

## check residuals of habitat model

``` r
plot(mod4, abline = c(0,0))
```

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-19-4.png)<!-- -->

``` r
#par(mfrow=c(1,2))
#hist(residuals(mod4)) 
#qqnorm(residuals(mod4))
```

looks okay.

## check residuals of full model

``` r
plot(mod5, abline = c(0,0))
```

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-20-3.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-20-4.png)<!-- -->

``` r
#par(mfrow=c(1,2))
#hist(residuals(mod5)) 
#qqnorm(residuals(mod5))
```

# adults

## set up models - no random effects

``` r
mod1 <- lm(PC1 ~ OUT_DIST, data = omy_ad)
mod2 <- lm(PC1 ~ MWMT, data = omy_ad)
mod3 <- lm(PC1 ~ FlowVel + SLOPE, data = omy_ad)
mod4 <- lm(PC1 ~ CANOPY + IP_Steelhd +Pool_freq + Spawnable, data = omy_ad)
mod5 <- lm(PC1 ~ MWMT + FlowVel + CANOPY + Pool_freq + Spawnable, data = omy_ad)
```

## compare models

``` r
stargazer(mod1, mod2, mod3, mod4, mod5, 
          type = "text",
          digits = 3, 
          ci = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
```

    ## 
    ## ===========================================================================================================================================
    ##                                                                       Dependent variable:                                                  
    ##                     -----------------------------------------------------------------------------------------------------------------------
    ##                                                                               PC1                                                          
    ##                               (1)                     (2)                     (3)                     (4)                     (5)          
    ## -------------------------------------------------------------------------------------------------------------------------------------------
    ## OUT_DIST                   -0.398***                                                                                                       
    ##                        (-0.508, -0.288)                                                                                                    
    ##                                                                                                                                            
    ## MWMT                                               0.355***                                                                 0.309**        
    ##                                                 (0.220, 0.490)                                                          (0.083, 0.534)     
    ##                                                                                                                                            
    ## FlowVel                                                                    0.175***                                         0.054*         
    ##                                                                         (0.100, 0.250)                                  (0.002, 0.106)     
    ##                                                                                                                                            
    ## SLOPE                                                                       -0.099                                                         
    ##                                                                         (-0.201, 0.003)                                                    
    ##                                                                                                                                            
    ## CANOPY                                                                                              0.108**                 0.111**        
    ##                                                                                                 (0.041, 0.176)          (0.042, 0.181)     
    ##                                                                                                                                            
    ## IP_Steelhd                                                                                           0.022                                 
    ##                                                                                                 (-0.033, 0.077)                            
    ##                                                                                                                                            
    ## Pool_freq                                                                                          -0.178***                -0.047         
    ##                                                                                                (-0.248, -0.108)         (-0.165, 0.070)    
    ##                                                                                                                                            
    ## Spawnable                                                                                            0.031                   0.018         
    ##                                                                                                 (-0.008, 0.069)         (-0.021, 0.057)    
    ##                                                                                                                                            
    ## Constant                     0.016                   0.082                 0.202***                0.170***                  0.065         
    ##                         (-0.065, 0.097)        (-0.0003, 0.164)         (0.139, 0.266)          (0.112, 0.227)          (-0.024, 0.154)    
    ##                                                                                                                                            
    ## -------------------------------------------------------------------------------------------------------------------------------------------
    ## Observations                  528                     528                     528                     528                     528          
    ## R2                           0.087                   0.048                   0.063                   0.120                   0.134         
    ## Adjusted R2                  0.085                   0.046                   0.060                   0.113                   0.125         
    ## Residual Std. Error    0.496 (df = 526)        0.506 (df = 526)        0.503 (df = 525)        0.488 (df = 523)        0.485 (df = 522)    
    ## F Statistic         50.195*** (df = 1; 526) 26.613*** (df = 1; 526) 17.670*** (df = 2; 525) 17.761*** (df = 4; 523) 16.121*** (df = 5; 522)
    ## ===========================================================================================================================================
    ## Note:                                                                                                         *p<0.05; **p<0.01; ***p<0.001

``` r
aic_vals <- c(extractAIC(mod1)[2], extractAIC(mod2)[2], extractAIC(mod3)[2], 
              extractAIC(mod4)[2], extractAIC(mod5)[2])
names(aic_vals) <- c("mod1","mod2","mod3", "mod4", "mod5")
aic_vals
```

    ##      mod1      mod2      mod3      mod4      mod5 
    ## -738.8682 -716.8039 -723.1407 -751.9965 -758.5630

\*\* full model (mod5) has lowest AIC \*\*

## check residuals of full model

``` r
plot(mod5, abline = c(0,0))
```

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-24-2.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-24-3.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-24-4.png)<!-- -->

``` r
#par(mfrow=c(1,2))
#hist(residuals(mod5)) 
#qqnorm(residuals(mod5))
```

# adult steelhead

## set up models - no random effects

``` r
mod1 <- lm(PC1 ~ OUT_DIST, data = omy_ad_stlhd)
mod2 <- lm(PC1 ~ MWMT, data = omy_ad_stlhd)
mod3 <- lm(PC1 ~ FlowVel + SLOPE, data = omy_ad_stlhd)
mod4 <- lm(PC1 ~ CANOPY + IP_Steelhd +Pool_freq + Spawnable, data = omy_ad_stlhd)
mod5 <- lm(PC1 ~ MWMT + FlowVel + CANOPY + Pool_freq + Spawnable, data = omy_ad_stlhd)
```

## compare models

``` r
stargazer(mod1, mod2, mod3, mod4, mod5, 
          type = "text",
          digits = 3, 
          ci = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
```

    ## 
    ## ===========================================================================================================================================
    ##                                                                       Dependent variable:                                                  
    ##                     -----------------------------------------------------------------------------------------------------------------------
    ##                                                                               PC1                                                          
    ##                               (1)                     (2)                     (3)                     (4)                     (5)          
    ## -------------------------------------------------------------------------------------------------------------------------------------------
    ## OUT_DIST                   -0.390***                                                                                                       
    ##                        (-0.508, -0.272)                                                                                                    
    ##                                                                                                                                            
    ## MWMT                                               0.378***                                                                 0.362**        
    ##                                                 (0.232, 0.525)                                                          (0.115, 0.610)     
    ##                                                                                                                                            
    ## FlowVel                                                                    0.163***                                         0.060*         
    ##                                                                         (0.084, 0.242)                                  (0.007, 0.113)     
    ##                                                                                                                                            
    ## SLOPE                                                                       -0.079                                                         
    ##                                                                         (-0.188, 0.031)                                                    
    ##                                                                                                                                            
    ## CANOPY                                                                                              0.096**                 0.111**        
    ##                                                                                                 (0.024, 0.169)          (0.037, 0.186)     
    ##                                                                                                                                            
    ## IP_Steelhd                                                                                           0.034                                 
    ##                                                                                                 (-0.023, 0.090)                            
    ##                                                                                                                                            
    ## Pool_freq                                                                                          -0.193***                -0.040         
    ##                                                                                                (-0.270, -0.117)         (-0.169, 0.089)    
    ##                                                                                                                                            
    ## Spawnable                                                                                            0.048                   0.032         
    ##                                                                                                (-0.00001, 0.096)        (-0.016, 0.080)    
    ##                                                                                                                                            
    ## Constant                     0.022                   0.074                 0.218***                0.168***                  0.043         
    ##                         (-0.066, 0.111)         (-0.016, 0.164)         (0.149, 0.286)          (0.106, 0.229)          (-0.054, 0.140)    
    ##                                                                                                                                            
    ## -------------------------------------------------------------------------------------------------------------------------------------------
    ## Observations                  495                     495                     495                     495                     495          
    ## R2                           0.079                   0.049                   0.063                   0.125                   0.142         
    ## Adjusted R2                  0.077                   0.047                   0.059                   0.118                   0.133         
    ## Residual Std. Error    0.496 (df = 493)        0.503 (df = 493)        0.500 (df = 492)        0.484 (df = 490)        0.480 (df = 489)    
    ## F Statistic         42.183*** (df = 1; 493) 25.584*** (df = 1; 493) 16.485*** (df = 2; 492) 17.550*** (df = 4; 490) 16.128*** (df = 5; 489)
    ## ===========================================================================================================================================
    ## Note:                                                                                                         *p<0.05; **p<0.01; ***p<0.001

``` r
aic_vals <- c(extractAIC(mod1)[2], extractAIC(mod2)[2], extractAIC(mod3)[2], 
              extractAIC(mod4)[2], extractAIC(mod5)[2])
names(aic_vals) <- c("mod1","mod2","mod3", "mod4", "mod5")
aic_vals
```

    ##      mod1      mod2      mod3      mod4      mod5 
    ## -693.0982 -677.5022 -682.5658 -712.7343 -720.0171

\*\* full model (mod5) has lowest AIC \*\*

## check residuals of full model

``` r
plot(mod5, abline = c(0,0))
```

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-28-2.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-28-3.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-28-4.png)<!-- -->

``` r
#par(mfrow=c(1,2))
#hist(residuals(mod5)) 
#qqnorm(residuals(mod5))
```

## now using the adaptive response variable… this is a proportion so go with “family = binomial” in glmer

``` r
library(bbmle) # for AICtab
```

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'bbmle'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     slice

# juvenile

GLMs

``` r
mod1_glm <- glm(alleleS_count1  ~ OUT_DIST, family = "binomial", data = omy_jv)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
mod2_glm <- glm(alleleS_count1  ~ MWMT, family = "binomial", data = omy_jv)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
mod3_glm <- glm(alleleS_count1  ~ FlowVel + SLOPE, family = "binomial", data = omy_jv)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
mod4_glm <- glm(alleleS_count1  ~ CANOPY + IP_Steelhd +Pool_freq + Spawnable, family = "binomial", data = omy_jv)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
mod5_glm <- glm(alleleS_count1  ~ MWMT + FlowVel + CANOPY + Pool_freq + Spawnable, family = "binomial", data = omy_jv)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
AICtab(mod1_glm, mod2_glm, mod3_glm, mod4_glm, mod5_glm)
```

    ##          dAIC df
    ## mod4_glm  0.0 4 
    ## mod5_glm  0.0 4 
    ## mod1_glm  6.0 2 
    ## mod2_glm 14.6 2 
    ## mod3_glm 62.1 3

mod4 (habitat) and mod5 (full) have lowest AIC

``` r
stargazer(mod1_glm, mod2_glm, mod3_glm, mod4_glm, mod5_glm, 
          type = "text",
          digits = 3, 
          ci = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
```

    ## 
    ## ====================================================================================================
    ##                                                  Dependent variable:                                
    ##                   ----------------------------------------------------------------------------------
    ##                                                     alleleS_count1                                  
    ##                         (1)              (2)              (3)             (4)              (5)      
    ## ----------------------------------------------------------------------------------------------------
    ## OUT_DIST              1.487***                                                                      
    ##                    (1.051, 1.923)                                                                   
    ##                                                                                                     
    ## MWMT                                  -0.819***                                          -0.343     
    ##                                    (-1.063, -0.576)                                  (-1.777, 1.092)
    ##                                                                                                     
    ## FlowVel                                                 -0.454                           -0.138     
    ##                                                     (-1.440, 0.532)                  (-2.655, 2.379)
    ##                                                                                                     
    ## SLOPE                                                   0.947**                                     
    ##                                                     (0.265, 1.628)                                  
    ##                                                                                                     
    ## CANOPY                                                                  1.456***          0.589     
    ##                                                                      (0.789, 2.123)  (-1.472, 2.650)
    ##                                                                                                     
    ## IP_Steelhd                                                              -1.012*                     
    ##                                                                     (-1.910, -0.115)                
    ##                                                                                                     
    ## Pool_freq                                                                0.034                      
    ##                                                                     (-0.395, 0.464)                 
    ##                                                                                                     
    ## Spawnable                                                                                           
    ##                                                                                                     
    ##                                                                                                     
    ## Constant             -1.320***         -0.340*          -0.346          0.926**           0.132     
    ##                   (-1.881, -0.760) (-0.679, -0.001) (-1.157, 0.465)  (0.337, 1.514)  (-0.771, 1.035)
    ##                                                                                                     
    ## ----------------------------------------------------------------------------------------------------
    ## Observations            269              269              269             269              269      
    ## Log Likelihood        -120.004         -124.308        -147.075         -115.005        -115.005    
    ## Akaike Inf. Crit.     244.009          252.616          300.151         238.010          238.010    
    ## ====================================================================================================
    ## Note:                                                                  *p<0.05; **p<0.01; ***p<0.001

model checking and diagnostics

## check residuals of model

``` r
plot(mod4_glm, abline = c(0,0))
```

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-33-2.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-33-3.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-33-4.png)<!-- -->

``` r
plot(mod5_glm, abline = c(0,0))
```

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-33-5.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-33-6.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-33-7.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-33-8.png)<!-- -->

# adults

GLMs

``` r
mod1_glm <- glm(alleleS_count1  ~ OUT_DIST, family = "binomial", data = omy_ad)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
mod2_glm <- glm(alleleS_count1  ~ MWMT, family = "binomial", data = omy_ad)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
mod3_glm <- glm(alleleS_count1  ~ FlowVel + SLOPE, family = "binomial", data = omy_ad)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
mod4_glm <- glm(alleleS_count1  ~ CANOPY + IP_Steelhd +Pool_freq + Spawnable, family = "binomial", data = omy_ad)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
mod5_glm <- glm(alleleS_count1  ~ MWMT + FlowVel + CANOPY + Pool_freq + Spawnable, family = "binomial", data = omy_ad)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
AICtab(mod1_glm, mod2_glm, mod3_glm, mod4_glm, mod5_glm)
```

    ##          dAIC df
    ## mod1_glm  0.0 2 
    ## mod4_glm  0.3 5 
    ## mod5_glm  0.8 6 
    ## mod3_glm 13.7 3 
    ## mod2_glm 18.0 2

mod1 (rKM), mod4 (habitat) and mod5 (full) have lowest AIC

``` r
stargazer(mod1_glm, mod2_glm, mod3_glm, mod4_glm, mod5_glm, 
          type = "text",
          digits = 3, 
          ci = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
```

    ## 
    ## ======================================================================================================
    ##                                                   Dependent variable:                                 
    ##                   ------------------------------------------------------------------------------------
    ##                                                      alleleS_count1                                   
    ##                         (1)              (2)              (3)              (4)              (5)       
    ## ------------------------------------------------------------------------------------------------------
    ## OUT_DIST              1.558***                                                                        
    ##                    (0.851, 2.266)                                                                     
    ##                                                                                                       
    ## MWMT                                   -1.419**                                            -0.625     
    ##                                    (-2.353, -0.486)                                   (-2.298, 1.048) 
    ##                                                                                                       
    ## FlowVel                                                 -0.653*                            -0.153     
    ##                                                     (-1.184, -0.123)                  (-0.550, 0.244) 
    ##                                                                                                       
    ## SLOPE                                                    0.427                                        
    ##                                                     (-0.247, 1.101)                                   
    ##                                                                                                       
    ## CANOPY                                                                    -0.310           -0.289     
    ##                                                                      (-0.795, 0.175)  (-0.830, 0.252) 
    ##                                                                                                       
    ## IP_Steelhd                                                                -0.023                      
    ##                                                                      (-0.436, 0.389)                  
    ##                                                                                                       
    ## Pool_freq                                                                0.849***          0.610      
    ##                                                                       (0.404, 1.294)  (-0.134, 1.354) 
    ##                                                                                                       
    ## Spawnable                                                                 -0.086           -0.069     
    ##                                                                      (-0.357, 0.185)  (-0.350, 0.212) 
    ##                                                                                                       
    ## Constant              -0.550*          -0.758**        -1.235***        -1.092***         -0.859*     
    ##                   (-1.047, -0.052) (-1.309, -0.208) (-1.680, -0.789) (-1.477, -0.708) (-1.551, -0.167)
    ##                                                                                                       
    ## ------------------------------------------------------------------------------------------------------
    ## Observations            291              291              291              291              291       
    ## Log Likelihood        -97.776          -106.772         -103.624         -94.940          -94.170     
    ## Akaike Inf. Crit.     199.553          217.544          213.248          199.881          200.340     
    ## ======================================================================================================
    ## Note:                                                                    *p<0.05; **p<0.01; ***p<0.001

model checking and diagnostics

## check residuals of model

``` r
plot(mod1_glm, abline = c(0,0))
```

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-37-2.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-37-3.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-37-4.png)<!-- -->

``` r
plot(mod4_glm, abline = c(0,0))
```

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-37-5.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-37-6.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-37-7.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-37-8.png)<!-- -->

``` r
plot(mod5_glm, abline = c(0,0))
```

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-37-9.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-37-10.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-37-11.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-37-12.png)<!-- -->

# steelhead

GLMs

``` r
mod1_glm <- glm(alleleS_count1  ~ OUT_DIST, family = "binomial", data = omy_ad_stlhd)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
mod2_glm <- glm(alleleS_count1  ~ MWMT, family = "binomial", data = omy_ad_stlhd)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
mod3_glm <- glm(alleleS_count1  ~ FlowVel + SLOPE, family = "binomial", data = omy_ad_stlhd)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
mod4_glm <- glm(alleleS_count1  ~ CANOPY + IP_Steelhd +Pool_freq + Spawnable, family = "binomial", data = omy_ad_stlhd)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
mod5_glm <- glm(alleleS_count1  ~ MWMT + FlowVel + CANOPY + Pool_freq + Spawnable, family = "binomial", data = omy_ad_stlhd)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
AICtab(mod1_glm, mod2_glm, mod3_glm, mod4_glm, mod5_glm)
```

    ##          dAIC df
    ## mod1_glm  0.0 2 
    ## mod5_glm  3.2 6 
    ## mod4_glm  3.4 5 
    ## mod3_glm 15.3 3 
    ## mod2_glm 19.1 2

mod1 (rKM) has lowest AIC

``` r
stargazer(mod1_glm, mod2_glm, mod3_glm, mod4_glm, mod5_glm, 
          type = "text",
          digits = 3, 
          ci = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
```

    ## 
    ## ====================================================================================================
    ##                                                  Dependent variable:                                
    ##                   ----------------------------------------------------------------------------------
    ##                                                     alleleS_count1                                  
    ##                         (1)             (2)              (3)              (4)              (5)      
    ## ----------------------------------------------------------------------------------------------------
    ## OUT_DIST             1.700***                                                                       
    ##                   (0.953, 2.446)                                                                    
    ##                                                                                                     
    ## MWMT                                  -1.608**                                           -0.945     
    ##                                   (-2.590, -0.625)                                   (-2.775, 0.886)
    ##                                                                                                     
    ## FlowVel                                                -0.747**                          -0.169     
    ##                                                    (-1.296, -0.198)                  (-0.576, 0.237)
    ##                                                                                                     
    ## SLOPE                                                   0.582                                       
    ##                                                    (-0.125, 1.289)                                  
    ##                                                                                                     
    ## CANOPY                                                                   -0.297          -0.318     
    ##                                                                     (-0.809, 0.214)  (-0.899, 0.263)
    ##                                                                                                     
    ## IP_Steelhd                                                               -0.036                     
    ##                                                                     (-0.460, 0.388)                 
    ##                                                                                                     
    ## Pool_freq                                                               0.881***          0.523     
    ##                                                                      (0.413, 1.350)  (-0.295, 1.340)
    ##                                                                                                     
    ## Spawnable                                                                -0.048          -0.018     
    ##                                                                     (-0.330, 0.234)  (-0.310, 0.275)
    ##                                                                                                     
    ## Constant              -0.461          -0.673*         -1.177***        -1.081***         -0.733     
    ##                   (-0.995, 0.073) (-1.258, -0.088) (-1.643, -0.711) (-1.484, -0.678) (-1.479, 0.014)
    ##                                                                                                     
    ## ----------------------------------------------------------------------------------------------------
    ## Observations            281             281              281              281              281      
    ## Log Likelihood        -89.586         -99.161          -96.216          -88.302          -87.167    
    ## Akaike Inf. Crit.     183.173         202.322          198.431          186.604          186.335    
    ## ====================================================================================================
    ## Note:                                                                  *p<0.05; **p<0.01; ***p<0.001

model checking and diagnostics

## check residuals of model

``` r
plot(mod1_glm, abline = c(0,0))
```

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-41-2.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-41-3.png)<!-- -->

    ## Warning in plot.window(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "abline" is not a graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in axis(side = side, at = at, labels = labels, ...): "abline" is not a
    ## graphical parameter

    ## Warning in box(...): "abline" is not a graphical parameter

    ## Warning in title(...): "abline" is not a graphical parameter

    ## Warning in plot.xy(xy.coords(x, y), type = type, ...): "abline" is not a
    ## graphical parameter

    ## Warning in title(sub = sub.caption, ...): "abline" is not a graphical parameter

![](elwha_IBE_files/figure-gfm/unnamed-chunk-41-4.png)<!-- -->
