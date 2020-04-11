# script for analysis of copepod methylation data for Homework 3
### using New data Reid uploaded

library(methylKit)
library(tidyverse)
library(pheatmap)
library(magrittr)
library(ggthemes)
library(ggbiplot)
library(gt)

# first, we want to read in the raw methylation calls with methylkit

# set directory with absolute path (why is this necessary? I have no idea, but gz files wont work with relative paths)

dir <- "/Users/baxterworthing/Downloads/Epigenetics_data" 

# read in the sample ids

samples <- read.table("/Users/baxterworthing/Downloads/Epigenetics_data/sample_id.txt", header = F)

# now point to coverage files

files <- file.path(dir, samples$V1)
all(file.exists(files))

# convert to list

file.list <- as.list(files)

# get the names only for naming our samples

nmlist <- as.list(gsub("_1_bismark_bt2_pe.bismark.cov.gz","",samples$V1))

# use methRead to read in the coverage files

myobj <- methRead(location= file.list,
                  sample.id =   nmlist,
                  assembly = "atonsa", # this is just a string. no actual database
                  dbtype = "tabix", 
                  context = "CpG",
                  resolution = "base", # SNP data
                  mincov = 20, # filetr low cov bases (high bc we have pooled data, generally would use 5-10)
                  treatment = 
                    c(0,0,0,0, #these correspond to the nmlist (first 4 are AA_F0, next 4 are AA_F25 etc, this is just giving single number for each rep ID)
                      1,1,1,1,
                      2,2,2,2,
                      3,3,3,3,
                      4,4,4,4),
                  pipeline = "bismarkCoverage",
                  dbdir = "~/Storage/methylkit")


######
# visualize coverage and filter
######

# We can look at the coverage for individual samples with getCoverageStats()
getCoverageStats(myobj[[1]], plot = T) # makes histo of per base coverage, high values are probably error 
# so this is AA_F00_1

getCoverageStats(myobj[[2]])




# I absolutely hate this fuction because it auto plots and auto prints, so you cant get the stats in an object without capture.output() and cant save the plot as an object without recordPlot()

# I hate myslef for doing this
test <- capture.output(getCoverageStats(myobj[[3]])) 
. <- unlist(str_split(test[4], "   "))

testOut <- data.frame(Median=.[3], Mean=.[4])

# make function

medMen <- function(x){
  
  test <- capture.output(getCoverageStats(x)) 
  . <- unlist(str_split(test[4], "   "))
  # this is sooooo annoying bc length can be either 5 or 6
  ifelse(length(.)==5, 
  return(data.frame(Median=as.numeric(.[3]), Mean=as.numeric(.[4]))), return(data.frame(Median=as.numeric(.[4]), Mean=as.numeric(.[5]))))

    
}

# get all median and mean coverage stats 
covList <- myobj %>%
  map(medMen)

# would rather do DF
covDF <- myobj %>%
  map_df(medMen)

# add names 
covDF$Sample <- unlist(nmlist)


# filter out my samples 

myCovDF <- covDF %>%
filter(str_detect(Sample, 'AA_F25|HH|AH')) 

#### fig 1 shoud just be coverage ####

covGT <- gt(myCovDF) %>% 
  tab_header(title = " Coverage Across Samples") 

gtsave(covGT, "covGT.png")

# filter samples by depth with filterByCoverage() 

filtered.myobj <- filterByCoverage(myobj, lo.count = 20, lo.perc = NULL,
                                   hi.count = NULL, high.perc=97.5, db.dir="~/Storage/methylkit")



# read merged samples 
meth <- methylKit:::readMethylBaseDB(
  dbpath = "/Users/baxterworthing/Downloads/Epigenetics_data/methylBase_united.txt.bgz",
  dbtype = "tabix",
  sample.id =   unlist(nmlist),
  assembly = "atonsa", # this is just a string. no actual database
  context = "CpG",
  resolution = "base",
  treatment = c(0,0,0,0,
                1,1,1,1,
                2,2,2,2,
                3,3,3,3,
                4,4,4,4),
  destrand = FALSE)



# percMethylation() calculates the percent methylation for each site and sample

pm <- percMethylation(meth) # makes huge mat of meth% at each site in each sample

#plot methylation histograms
ggplot(gather(as.data.frame(pm)), aes(value)) + 
  geom_histogram(bins = 10, color="black", fill="grey") + 
  facet_wrap(~key)

# calculate and plot mean methylation

sp.means <- colMeans(pm)

p.df <- data.frame(sample=names(sp.means),
                   group = substr(names(sp.means), 1,6),
                   methylation = sp.means)

ggplot(p.df, aes(x=group, y=methylation, color=group)) + 
  stat_summary(color="black") + geom_jitter(width=0.1, size=3)

#### I want to subset my samples for fig 2 ####

groups <- c("AA_F25", "HH_F25", "AH_F25")

myP.df <- filter(p.df, group %in% groups)

names(myP.df)[1] <- "Sample"

# make custom summary

stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data = fun, colour = "grey", geom = geom, width = 0.3)
}

# custom pal 

Pal <- c("#5DB1DDFF","#F0E685FF", "#D60047FF")


gmP <- ggplot(myP.df, aes(x=group, y=methylation, color=group)) + 
  stat_sum_df("mean_se", mapping = aes(group = group)) +
  geom_jitter(width=0.1, size=2) +
  labs(y= "Mean % Methylation", x= "Treatment", title = "Mean Methylation Across All Samples") +
  theme_solarized(light = F) +
  theme(legend.position ="none") +
  theme(axis.text =element_text(color="white")) +
  theme(axis.title =element_text(color="white")) +
  theme(title =element_text(color="white")) +
  scale_x_discrete(labels=c("AA", "AH",
                            "HH")) +
  scale_color_manual(values = Pal)

ggsave("gmP.jpg", gmP, height = 5, width = 5)

# each point is replicate, point is mean % of reads that had methylated C across all SNPs box shows mean +/- stand error

####  is there a relationship between coverage and methylation
####

covMethDf <- full_join(myP.df, myCovDF, by=)

cor.test(covMethDf$Mean, covMethDf$methylation, method = "spearman") # no 






# PCA
PCASamples(meth, screeplot=T)
PCASamples(meth)




#### better PCA ####

pm.df <- as.data.frame(t(pm))
pm.df <- pm.df[,colSums(pm.df != 0)]
pm.df$Sample <- row.names(pm.df)


pm.df %<>%
  filter(str_detect(Sample, 'AA_F25|HH|AH')) 

sampleFullNames <- pm.df$Sample

pm.df %<>%
  separate(Sample, into= c("Treatment", "Gen", "Rep"), sep = "_") 
  
# hate that I'm doing it this way bc I know theres a better way with magrittr
pcaGroups <- pm.df$Treatment

pm.df %<>%
  select(-c("Treatment", "Gen", "Rep"))

pm.pca <- prcomp(pm.df, center = TRUE,scale. = T)

base::summary(pm.pca)

pcaMeth <- ggbiplot(pm.pca, groups = pcaGroups, var.axes = F, alpha = .8) +
  geom_point(aes(color=pcaGroups), size= 4) +
  scale_color_manual( values = Pal, aes(title="Treatment")) +
  labs(x ="PC1 (85.7%)", y="PC2 (6.6%)", title = "PCA of Methylation Across All SNPs") +
  geom_vline(xintercept = 0, lty = 2, color="grey") +
  geom_hline(yintercept = 0, lty = 2, color="grey") +
  theme_solarized(light = F) +
  theme(legend.key = element_blank()) +
  theme(axis.text =element_text(color="white")) +
  theme(axis.title =element_text(color="white")) +
  theme(title =element_text(color="white")) +
  theme(legend.text  =element_text(color="white")) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) 
  
# I guess this is my fig 2

ggsave("methPCAsmall.jpg", pcaMeth, height = 4.5, width = 4.5)
  


### find differentially methylated sites between two groups ####

# this is where I break from the in class script becasue I want AA_F25 and HH_F25 and AA_F25 and AH_F25
Meth_sub_HH <- reorganize(meth,
                       sample.ids =c("AA_F25_1","AA_F25_2","AA_F25_3", "AA_F25_4",
                                     "HH_F25_1","HH_F25_2","HH_F25_3","HH_F25_4"),
                       treatment = c(0,0,0,0,1,1,1,1),
                       save.db=FALSE) 


Meth_sub_AH <- reorganize(meth,
                          sample.ids =c("AA_F25_1","AA_F25_2","AA_F25_3", "AA_F25_4",
                                        "AH_F25_1","AH_F25_2","AH_F25_3","AH_F25_4"),
                          treatment = c(0,0,0,0,1,1,1,1),
                          save.db=FALSE) 

myDiff_HH <- calculateDiffMeth(Meth_sub_HH, overdispersion = "MN", mc.cores = 1, suffix= "AA_HH", adjust = "qvalue", test="Chisq")

myDiff_AH <- calculateDiffMeth(Meth_sub_AH, overdispersion = "MN", mc.cores = 1, suffix= "AA_AH", adjust = "qvalue", test="Chisq")


hhDiffbase <- getMethylDiff(myDiff_HH, qvalue = .05, difference = 10)

ahDiffbase <- getMethylDiff(myDiff_AH, qvalue = .05, difference = 10)



#### plot differenital methylation btwn treatments across all differentially methylated SNPs that have diff> 10% and adusted P < .05  ####

hist(getData(hhDiffbase)$meth.diff)

#### fig 3 ####

hDatHH <- data.frame(diff= getData(hhDiffbase)$meth.diff)

hDatAH <- data.frame(diff= getData(ahDiffbase)$meth.diff) 

pal1 <- Pal[3]

hist1 <- ggplot(hDatHH, aes(x=diff)) +
  geom_histogram( alpha =.7, bins=40, fill=pal1, color="grey")+
  labs(x= "Methylation Difference", y="Count", title = "HH vs. AA Differential Methylation") +
  theme_solarized(light=F) +
  theme(axis.text =element_text(color="white")) +
  theme(axis.title =element_text(color="white")) +
  theme(title =element_text(color="white")) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) 

pal2 <- Pal[2]

hist2 <- ggplot(hDatAH, aes(x=diff)) +
  geom_histogram(alpha =.85, bins=40, fill=pal2, color="grey94")+
  labs(x= "Methylation Difference", y="Count", title = "AH vs. AA Differential Methylation") +
  theme_solarized(light=F) +
  theme(axis.text =element_text(color="white")) +
  theme(axis.title =element_text(color="white")) +
  theme(title =element_text(color="white")) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) 


ggsave("hist1.jpg", hist1, device = "jpg", height = 3.5, width = 4)

ggsave("hist2.jpg", hist2, device = "jpg", height = 3.5, width = 4)

#hyper

hyperHH <- getMethylDiff(myDiff_HH,difference=10,qvalue=0.05,type="hyper")
hyperAH <- getMethylDiff(myDiff_AH,difference=10,qvalue=0.05,type="hyper")

# get hypo methylated bases

hypoHH <- getMethylDiff(myDiff_HH,difference=10,qvalue=0.05,type="hypo")
hypoAH <- getMethylDiff(myDiff_AH,difference=10,qvalue=0.05,type="hypo")

# get percent methylation matrix
pmHH <- percMethylation(Meth_sub_HH)
pmAH <- percMethylation(Meth_sub_AH)
# make a dataframe with snp id's, methylation, etc.


pm.sigHH <- pmHH[as.numeric(row.names(hhDiffbase)),]
pm.sigAH <- pmAH[as.numeric(row.names(ahDiffbase)),]


# add snp, chr, start, stop


dinHH <- getData(hhDiffbase)[,1:3]
df.outHH <- cbind(paste(getData(hhDiffbase)$chr, getData(hhDiffbase)$start, sep=":"), dinHH, pm.sigHH)
colnames(df.outHH) <- c("snp", colnames(dinHH), colnames(df.outHH[5:ncol(df.outHH)]))
df.outHH <- (cbind(df.outHH,getData(hhDiffbase)[,5:7]))


dinAH <- getData(ahDiffbase)[,1:3]
df.outAH <- cbind(paste(getData(ahDiffbase)$chr, getData(ahDiffbase)$start, sep=":"), dinAH, pm.sigAH)
colnames(df.outAH) <- c("snp", colnames(dinAH), colnames(df.outAH[5:ncol(df.outAH)]))
df.outAH <- (cbind(df.outAH,getData(ahDiffbase)[,5:7]))


# which SNPs the same between the 2?

which(df.outAH$snp %in% df.outHH$snp) # only 4/30 in common



# specific gene or snp

df.plotHH <- df.outHH[,c(1,5:12)] %>% pivot_longer(-snp, values_to = "methylation")
df.plotHH$group <- substr(df.plotHH$name,1,2)
head(df.plotHH)

df.plotAH <- df.outAH[,c(1,5:12)] %>% pivot_longer(-snp, values_to = "methylation")
df.plotAH$group <- substr(df.plotHH$name,1,2)
head(df.plotAH)

# looking at snp LS049205.1:248

df.plotHH %>% filter(snp=="LS049205.1:248") %>% 
  ggplot(., aes(x=group, y=methylation, color=group, fill=group)) +
  stat_summary(fun.data = "mean_se", size = 2) +
  geom_jitter(width = 0.1, size=3, pch=21, color="black")

# same, SNP, differnt treatment group

df.plotAH %>% filter(snp=="LS049205.1:248") %>% 
  ggplot(., aes(x=group, y=methylation, color=group, fill=group)) +
  stat_summary(fun.data = "mean_se", size = 2) +
  geom_jitter(width = 0.1, size=3, pch=21, color="black")


# oringinal contrast
mP <- df.plot %>% filter(snp=="LS049205.1:248") %>% 
  ggplot(., aes(x=group, y=methylation, color=group, fill=group)) +
  stat_summary(fun.data = "mean_se", size = 2) +
  geom_jitter(width = 0.1, size=3, pch=21, color="black")

ggsave("mP.jpg", mP, device = "jpg", width = 5, height =5)


# points are SNPs in each of the 4 reps fat things are mean +- SD 

# I think what's going on here is that we are identifying C/T SNPs and asking, for each sample, what % of reads are C (methylated) vs T (unmethylated). Reference has had all Ts converted to Ts that's why this works 

write.table(file = "~/Documents/EcoGenom/diffmethHH.bed",
            data.frame(chr= df.outHH$chr, start = df.outHH$start, end = df.outHH$end),
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

write.table(file = "~/Documents/EcoGenom/diffmethAH.bed",
            data.frame(chr= df.outAH$chr, start = df.outAH$start, end = df.outAH$end),
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


# now realizing SNPs and chrs are factors  


df.outHH$chr <- as.character(df.outHH$chr)
df.outAH$chr <- as.character(df.outAH$chr)

df.outHH$snp <- as.character(df.outHH$snp)
df.outAH$snp <- as.character(df.outAH$snp)


# read in filtered bed files with gannotation info

annoHH <- read_delim("HH_uniqueHitsNames.bed", delim = "\t", col_names = F) %>%
  unite("snp", X1:X2, sep = ":")

annoAH  <- read_delim("AH_uniqueHitsNames.bed", delim = "\t", col_names = F)%>%
  unite("snp", X1:X2, sep = ":")

# need to add col to this table that says if the given gene was hyper or hypo methylated 

hypoSNPhh <- df.outHH %>%
  filter(meth.diff<0)

hypoSNPah <- df.outAH %>%
  filter(meth.diff<0)


annoHH %<>%
  mutate(Methylation=ifelse(snp %in% hypoSNPhh$snp, "hypo", "hyper")) 

annoAH %<>%
  mutate(Methylation=ifelse(snp %in% hypoSNPah$snp, "hypo", "hyper")) 

# there are too many damn columns in these 

annoHHtab <- annoHH[,-c(1,4,5,9)] 

names(annoHHtab) <- c("Gene Name 1", "Gene Name 2", "Cell Component", "Mol. Function", "Bol. Process", "Methylation") 

annoAHtab <- annoAH[,-c(1,4,5,9)] 

names(annoAHtab) <- c("Gene Name 1", "Gene Name 2", "Cell Component", "Mol. Function", "Bol. Process", "Methylation") 

write.csv(annoAHtab,"annoAH.tab")

write.csv(annoHHtab,"annoHH.tab")

# jk going to need to sacrifice some stuff


annoAHlite <- annoAHtab[,-c(2,3)]
names(annoAHlite)[1] <- "Gene Name"

annoHHlite <- annoHHtab[,-c(2,3)]
names(annoHHlite)[1] <- "Gene Name"

write.csv(annoAHlite,"annoAH.tab")

write.csv(annoHHlite,"annoHH.tab")

AHgt <-  gt(annoAHlite) %>%
  tab_options(table.width = pct(50))

HHgt <- gt(annoHHlite)#%>%
  tab_options(table.width = pct(80))


gtsave(AHgt, "AHgt.png")
gtsave(HHgt, "HHgt.png")
