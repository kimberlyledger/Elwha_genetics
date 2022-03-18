library(vcfR)
library(adegenet, dep=TRUE)
library(pegas)
library(poppr)
library(dplyr)

library(readr)
steelhead.metadata <- read_csv("C:/Users/jyj55/Desktop/LGC-DGS_2022/Team_project/Steelhead_genotype/Elwha_Steelhead_Formatted.csv")

steelhead.vcf <- read.vcfR("C:/Users/jyj55/Desktop/LGC-DGS_2022/Team_project/Steelhead_genotype/Elwha_GTSeq_Sans_CCT.vcf.gz")
steelhead.genind <- vcfR2genind(steelhead.vcf, sep = "/")

# Make a dataframe to store population info
steelhead_pop <- matrix(NA, nrow=nrow(steelhead.genind@tab), ncol=2)
steelhead_pop <- as.data.frame(steelhead_pop)
names(steelhead_pop) <- c("Sample_ID", "pop")

# Store population info for each individual
for (i in 1:nrow(steelhead.genind@tab)){
  steelhead_pop$Sample_ID[i] <- rownames(steelhead.genind@tab)[i]
  steelhead_pop$pop[i] <- steelhead.metadata %>% filter(Sample_ID == rownames(steelhead.genind@tab)[i]) %>% select(Location)
}

# Add population info into genind object
strata(steelhead.genind) <- steelhead_pop
setPop(steelhead.genind) <- ~pop

# HWE test
steelhead_hwt <- hw.test(steelhead.genind, B=0) # Because the number of individuals (= sample size) is large enough, parametric test (B=0) would be ok?

# Store loci with p < 0.05 into new object
steelhead_HWDloci <- matrix(NA, nrow=nrow(steelhead_hwt), ncol=2)
colnames(steelhead_HWDloci) <- c("name_loci", "p_value")
for (i in 1:nrow(steelhead_hwt)){
  if (steelhead_hwt[i,3] < 0.05){
    steelhead_HWDloci[i,1] <- rownames(steelhead_hwt)[i]
    steelhead_HWDloci[i,2] <- steelhead_hwt[i,3]
  }
}
write.csv(cbind(steelhead_HWDloci, steelhead_hwt), "loci_HWD.csv")

# LD test using "poppr" package and "genclone" object
steelhead.genclone <- as.genclone(steelhead.genind)
LDpair <- steelhead.genclone %>% pair.ia
LDpair <- tibble::rownames_to_column(as.data.frame(LDpair), var = "pairID")

# Store loci pairs in LD into new object
LD_over07 <- LDpair %>% filter(rbarD > 0.7) %>% select('pairID', 'rbarD') # I set up rbarD threshold for LD as 0.7 but I am not sure if it is the right value.
write.csv(LD_over07, "locipair_LD.csv")