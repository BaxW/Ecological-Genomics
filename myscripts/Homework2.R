####Questions:
# Do trees from different source climates differ in gene expression?
# Is there a gene expression response to heat stress?
        # Does response change when drought stress is added?
# Is there an interaction between a tree's source climate and the stress treatment in terms of expression
# What genes show the greatest biological response 

#!#!#!#!#!#!# note that I'm only using data from day 10 of each treatment

library(tidyverse)
library(magrittr)
library(ggthemes)
library(DESeq2)
library(patchwork)


# Import the counts matrix
countsTable <- read.table("RS_counts_samples/RS_cds2kb_countsMatrix.txt", header=TRUE, row.names=1)

#DEseq wants integers
countsTable %<>% round()


# only want day 10 
ten_countsTable <- countsTable %>%
  select(contains("_10_")) 

dim(ten_countsTable)

# need conditions for day 10

conds <- read.delim("RS_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1, colClasses=c('factor', 'factor', 'factor', 'factor'))

conds10  <- conds %>%
  filter(day=="10")


# average counts per gene 
mean(rowSums(ten_countsTable)) #1296.096
median(rowSums(ten_countsTable)) #10 (very high dispersion)
sd(rowSums(ten_countsTable)) # 154537.8 ("")

# DEseq object for model of gene expression as respoonse to source clim, treatment, and interaction between the two

dds <- DESeqDataSetFromMatrix(countData = ten_countsTable, colData = conds10, 
design = ~ climate + treatment + climate:treatment)

# without interaction

ddsNint <- DESeqDataSetFromMatrix(countData = ten_countsTable, colData = conds10, 
                                  design = ~ climate + treatment)


# swap source clim for the individual populations (which are nested wihin source clim)
# becomes singular with interaction 

ddsPop <- DESeqDataSetFromMatrix(countData = ten_countsTable, colData = conds10, 
                                 design = ~ pop + treatment )


#filter genes with low reads 

dds <- dds[rowSums(counts(dds)) > 30]
ddsNint <- ddsNint[rowSums(counts(ddsNint)) > 30]
ddsPop <- ddsPop[rowSums(counts(ddsPop)) > 30]

## Run the DESeq model to test for differential gene expression: 
# 1) estimate size factors (per sample), 2) estimate dispersion (per gene), 
# 3) run negative binomial glm
dds %<>% DESeq()
ddsNint %<>% DESeq()
ddsPop %<>% DESeq()

conts <- resultsNames(dds)[-1]
#[1] "Intercept"    "climate_HD_vs_CW"     "treatment_D_vs_C"    #[4] "[4] treatment_H_vs_C"     "climateHD.treatmentD" "climateHD.treatmentH"

nintConts <- resultsNames(ddsNint)[-1]
#[1] "Intercept"        "climate_HD_vs_CW" "treatment_D_vs_C"
#[4] "treatment_H_vs_C"

popConts <-  resultsNames(ddsPop)[-1]
#[1] "Intercept"            "pop_BRU_05_vs_ASC_06" "pop_CAM_02_vs_ASC_06"
#[4] "pop_ESC_01_vs_ASC_06" "pop_JAY_02_vs_ASC_06" "pop_KAN_04_vs_ASC_06"
#[7] "pop_LOL_02_vs_ASC_06" "pop_MMF_13_vs_ASC_06" "pop_NOR_02_vs_ASC_06"
#[10] "pop_XBM_07_vs_ASC_06" "treatment_D_vs_C"     "treatment_H_vs_C"

# build function to generate results for each contrast in dds 
# this will require building a big list 
# maybe not enough time/ no reason to do that for ddspop, but should compare the two contrasts in common because I'm not sure if a differnt model will change the same contrast 

getRes <- function(contrast, dds=dds) {
  
  #!#!#! not sure if I should use name or contrast
  res <-results(dds, name = contrast
                  ,alpha = .05) # set signif cutoff
 
   # only want the good stuff, and want tibble output
  res@listData$IDs <- res@rownames
  
  resOut <- as_tibble(res@listData) %>%
    arrange(padj)

  
  return(resOut)
}


# then conts is input to map() to make list of results for every contrast
# just trying to save some lines of code here 

ddsRes <- map(conts, getRes)
names(ddsRes) <- conts #need to add names so I know what's what

map(ddsRes, nrow) # 24,300 genes in each contrast

resDC <-results(dds, name = "treatment_D_vs_C" 
              ,alpha = .05)
summary(resDC)

resClim <-results(dds, name = "climate_HD_vs_CW"
                ,alpha = .05)
summary(resClim)


######## got a bit ahead of myself up there

# question 1

# need to subset for just control 

cont_ten_countsTable <- ten_countsTable %>%
  select(contains("_C_")) 

conds10Cont  <- conds10 %>%
  filter(treatment=="C")

ddsCont <- DESeqDataSetFromMatrix(countData = cont_ten_countsTable, colData = conds10Cont, 
                                 design = ~ climate )

ddsCont <- ddsCont[rowSums(counts(ddsCont)) > 30]

ddsCont %<>% DESeq()

resultsNames(ddsCont)

resCont <-  results(ddsCont, name = "climate_HD_vs_CW") 

summary(resCont) # no signif DE

# what if there are differences between source clim across all treatments?

resDCcon <-results(dds, contrast = c("climate", "HD","CW")
                 ,alpha = .05)

summary(resDCcon) # one gene significantly down regulated, kinda cool

# which gene?

ddsRes$climate_HD_vs_CW$IDs[1] #"MA_129323g0010"



# if I do the pop-based model can I pull out contrasts between pops with biggest difference 


ddsPopRes <- map(popConts, getRes, dds= ddsPop)
names(ddsPopRes) <- popConts


# now I need to map a function to count DE genes in each contrast


popDEcount<- map_df(ddsPopRes,function(x) {filter(x, padj<.05) %>%
    nrow()}) %>%
  pivot_longer(everything(),names_to = "contrast", values_to = "DEgenes")

# Question 2 

# hot vs control 
resHC <-results(dds, contrast = c("treatment", "H","C")
                   ,alpha = .05)
summary(resHC) # one up-regulated, one down-regualted

#with no interaction

resNintHC <-results(ddsNint, contrast = c("treatment", "H","C")
                ,alpha = .05)
summary(resNintHC) # WAY MORE DE

# I guess when I add the interaction term, varaicne in expression is better explained by the interaction, so DE genes get sucked away from the treatment contrast

#pop level

resPopHC <-results(ddsPop, contrast = c("treatment", "H","C")
                    ,alpha = .05)
summary(resPopHC) # Damn

summary((results(ddsPop, contrast = c("pop", "BRU_05","ASC_06")
        ,alpha = .05)))

# I like this pop model because I'm accounting for within-treatment differences at the finest scale I can. IMO that is a better way to ask if there is an overall expression change in response to the treatments

# I already have the results for this in ddsPopRes, and I just need the last 2 DFs in that list for this question 

# based on popDEcount, I already know that I have 122 genes DE btwn H and C and 2430 DE btwn D and C

# how many are in common?

deD <- filter(ddsPopRes$treatment_D_vs_C, padj<.05) %>%
  select(IDs)



deH <- filter(ddsPopRes$treatment_H_vs_C, padj<.05) %>%
  select(IDs)



dedeH %>%
  filter(IDs %in% deD$IDs) %>%
  nrow() # 90 in common, makes sense 

# I just want up-regulated genes for enrichment analysis 

upD <- filter(ddsPopRes$treatment_D_vs_C, padj<.05) %>%
  filter(log2FoldChange>0) %>%
  select(IDs)

upH <- filter(ddsPopRes$treatment_H_vs_C, padj<.05) %>%
  filter(log2FoldChange>0) %>%
  select(IDs)

write.table(upH,"upH.txt",quote = F, row.names = F)

write.table(upD,"upD.txt", quote = F, row.names = F)

# heat map of Go enrichment results 

hGO <- readxl::read_xlsx("SpruceGOenrich.xlsx", sheet = 3) %>% 
  mutate(Treatment= "Heat")

dGO <- readxl::read_xlsx("SpruceGOenrich.xlsx", sheet = 4) %>%
  mutate(Treatment= "Heat+Drought")

goDat <- full_join(dGO,hGO)

goDat <- base::rbind(dGO,hGO)

which(hGO$GO %in% dGO$GO)


goP <- ggplot(goDat, aes(y=Description, fill=-log10(pval), x=Treatment)) +
  geom_tile() +
  scale_fill_gradient(high = "#CE3D32FF", low = "#F4E8CFFF", na.value = "grey")  +
  labs(y= "GO Term", fill= "-log10(P-value)", title="GO Enrichment of Up-regulated Genes") + theme(plot.title = element_text(hjust=-4)) +
   theme_classic() 
  
ggsave("goP.jpg", goP, device = "jpg", height  =7, width = 7)

# volcano plots for both of these 
pal <- c("#5A655EFF","#CE3D32FF" )
pal2 <- c("#5A655EFF","#FFC20AFF" )

# pull out data and add a significance column based on -10log of adjusted p value 

D_volcDat <- ddsPopRes$treatment_D_vs_C %>% 
 mutate(signif=case_when(padj <.05~ "Significant", padj >.05 ~ "Insignificant")) %>%
  mutate(alpha=case_when(padj <.05~ .5, padj >.05 ~ .3) )

volD <- ggplot(D_volcDat, aes(x=log2FoldChange, y=-log10(padj), color=signif)) +
  geom_point(aes(alpha=alpha), size= 1.5) + 
  labs(x= "log2(Fold Change)", y="-log10(P-Value)", color="", title = "Differential Expression: Heat + Drought vs. Control") +
  theme_solarized(light = F) +
  scale_alpha(range=c(.2, .5)) + # I clearly don't know how to scale alpha right
  scale_color_manual(values=pal) +
  theme(axis.text = element_text(size = 19)) +
  theme(axis.title =element_text(size = 19)) +
  theme(legend.position = "none")



# do same for the H treatment

H_volcDat <- ddsPopRes$treatment_H_vs_C %>% 
  mutate(signif=case_when(padj <.05~ "Significant", padj >.05 ~ "Insignificant")) %>%
  mutate(alpha=case_when(padj <.05~ .5, padj >.05 ~ .3) )

volH <- ggplot(H_volcDat, aes(x=log2FoldChange, y=-log10(padj), color=signif)) +
  geom_point(aes(alpha=alpha), size= 1.5) + 
  labs(x= "log2(Fold Change)", y="-log10(P-Value)", color="", title = "Differential Expression: Heat vs. Control") +
  theme_solarized(light = F) +
  scale_alpha(range=c(.2, .5)) + 
  scale_color_manual(values=pal2) +
  theme(axis.text = element_text(size = 19)) +
  theme(axis.title =element_text(size = 19)) +
  theme(legend.position = "none")

#ggsave("./volH.jpg", volH, device = "jpg", height = 7, width = 7)

volComb <-  volH + volD

ggsave("./volComb.jpg", volComb, device = "jpg", height = 7, width = 12)

# question 3
# interaction term will be crucial for question 3 because the interaction results ask if the effect of a given treatment depends on climate 


resIntD <-results(dds, name = "climateHD.treatmentD"
                ,alpha = .05)

summary(resIntD) # nada

resIntH <-results(dds, name = "climateHD.treatmentH"
                   ,alpha = .05)

summary(resIntH) # also nada

# there is another way to test for this (subset down to one treatment and look for DE between source clim)

D_ten_countsTable <- ten_countsTable %>%
  select(contains("_D_H_")) 
conds10D  <- conds10 %>%
  filter(treatment=="D")

ddsD <- DESeqDataSetFromMatrix(countData = D_ten_countsTable, colData = conds10D, 
                                  design = ~ climate )

ddsD <- ddsD[rowSums(counts(ddsD)) > 30]

ddsD %<>% DESeq()
resultsNames(ddsD)
resD <-  results(ddsD, name = "climate_HD_vs_CW") 
summary(resD) # just 2 genes 


# at this point, I'm increadibly confused bc Melissa's script reports over 1000 genes being differentially expressed in climateHD.treatmentD


#attempt pca nonetheless 

vsd <- vst(dds, blind=FALSE)


data <- plotPCA(vsd,intgroup=c("climate","treatment"),returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

data$treatment <- factor(data$treatment, levels=c("C","H","D"), labels = c("C","H","D"))


pcaPal <- c("#466983FF", "#BA6338FF")

pcaClim <- ggplot(data, aes(PC1, PC2, color=climate, shape=treatment)) +
  geom_point(size=3.9, alpha=0.8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(color= "Climate", shape= "Treatment", title = "PCA of Gene Expression")+
  geom_vline(xintercept = 0, lty = 2, color="grey") +
  geom_hline(yintercept = 0, lty = 2, color="grey") +
  theme_few() +
  scale_color_manual(values = pcaPal)

#ggsave("pcaClim.jpg", pcaClim, device="jpg", height = 7, width = 7 )
