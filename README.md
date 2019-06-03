# ReserveBenefit--RDA_adaptation

# Stacks Workflow

RADseq workflow using [STACKS](http://creskolab.uoregon.edu/stacks/)

Developed by [Laura Benestan](https://github.com/laurabenestan) in
[Stephanie Manel](https://sites.google.com/site/stephaniemanel/home)'s
laboratory.

## Concept and definition
This README.txt is widely inspired by the vignette created by Brenna Forester (see https://popgen.nescent.org/2018-03-27_RDA_GEA.html).
We used of Redundancy Analysis (RDA) as a genotype-environment association (GEA) method to simultaneously assess the percent of genomic variation explained by environmental variables and to detect loci under selection (see the relevant paper of Forester et al., 2018). 
RDA is a two-step analysis in which genetic and environmental data are analyzed using multivariate linear regression. 
Then PCA of the fitted values is used to produce canonical axes, which are linear combinations of the predictors (Legendre & Legendre, 2012). 
here, we performed the RDA on a individual-based sampling design.

## Application to our dataset
Here, we applied RDA to genomic data from 276 and 335 indviduals of the white seabream (Diplodus sargus) and the stripped red mullet (Mullus surmuletus) sampled across the Mediterranean sea. 
Results of the RDA at the full set of 18,512 and 14,318 single nucleotide polymorphism (SNP) markers will be soon available . 
We are interested in understanding which environmental factor may drive the genomic divergence of several species. 
In this case, the data are individual-based, and are input as allele counts (i.e. 0/1/2) for each locus for each individual fish. 

# Using R to perform the analysis
## Install packages
First, we need to install the necessary packages and then download the corresponding libraries.
```
library(psych)    
library(vegan)
library(adegenet)
library(dplyr)
library(fmsb)
library(gsl)
```
## Read genetic data

```
#plink_diplodus <- read.PLINK("18512snps-276ind-diplodus.raw", map.file = "18512snps-276ind-diplodus.map", quiet = FALSE)
plink_diplodus <- read.table("18512snps-276ind-diplodus.raw", header=TRUE, sep=" ", row.names=1)
plink_mullus <- read.table("14318snps_312ind.raw", header=TRUE, sep=" ", row.names=1)
dim(plink_diplodus)
```

## Prepare genetic data

We need to remove the names of the samples as we won't consider it to fill the missing values.
Then missing values are filled using the overall average of the allele frequency variation.

```
### Remove names
gen <- plink_diplodus[,6:37029]
gen <- plink_mullus[,6:28641]
### Calculate the NA
sum(is.na(gen))# 553,232 NAs in the matrix (~3% missing data)
### Fill the NAs
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp)) # No NAs
```

## Read the environmental data
```
### Download environmental data 
env <- read.table("276ind-24env-variables.txt",sep="\t",header=T) # for diplodus
env <- read.table("24env_variables_mullus.txt",sep="\t",header=T) # for mullus
str(env)
### Add habitats variable
habitat <- read.table("habitat_diplodus_276ind.txt",sep="\t",header=T)
env <- cbind(env, habitat)
```

## Prepare the environmental data
```
### Make individual names characters (not factors)
env$labels <- as.character(env$labels)
```

## Remove one of the highly correlated predictors

We visually and quantitatively checked the correlation among predictors using the function called `pairs.panels`.
### Visualize correlations among predictors
pairs.panels(env[,2:25], scale=T)

! Screen Shot 2019-04-16 at 19.58.38.png

```
### Create a correlation matrix
matrix_cor <- cor(env[,2:25])

### Selection for the non-correlated variables
env_selected <- select(env, salinity_surface_repro_min,chlo_surface_max,chlo_benthic_min) # for diplodus
env_selected <- select(env, salinity_water_column_whole_yr_max,chlo_surface_max,chlo_benthic_max,chlo_benthic_min)# for mullus
table(env_selected$salinity_water_column_whole_yr_max)
```

## Run the Redundancy Analysis

### Perform the RDA calculation
This step can take a while depending on your dataset.
The R2 inform about the percent of genomic variation that can be explained by one of the predictor.
Here for instance, we found a adjusted R2 of 0.0002, which means that our constrained ordination explains about 0.02% of the variation.
This low explanatory power is not surprising given that we expect that most of the SNPs in our dataset will not show a relationship with the environmental predictors (e.g., most SNPs will be neutral).
```
diplodus.rda <- rda(gen.imp ~ ., data=pred, scale=T)
diplodus.rda 
# O
### Calculate the adjusted R2
RsquareAdj(diplodus.rda)

### The eigenvalues for the constrained axes reflect the variance explained by each canonical axis
summary(eigenvals(diplodus.rda, model = "constrained")
screeplot(diplodus.rda)
```
Then we tested the significance of the RDA.
```
load.rda <- scores(diplodus.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes
##  we can see that the first three constrained axes explain all the variance.
signif.full <- anova.cca(diplodus.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full
# 0.01 ** 
#0.49 for mullus
```
## Represent the RDA outcomes in a nice graph
```
### Plot the RDA
plot(diplodus.rda, scaling=3) 
plot(diplodus.rda, choices = c(1, 3), scaling=3)
```
We add the grouping information, here the Marine Protected Areas info.
```
### Add the MPA information
mpa_group <- read.table("population-map-276ind-diplodus-mpa.txt",header=TRUE,sep="\t")
mpa_group <- read.table("population-map-335ind-mullus-mpa.txt",header=TRUE,sep="\t")
mpa_env <- cbind(env, mpa_group)
eco <- mpa_group$STRATA
```
We also add a second grouping information, here the "inside" or "outside" an MPA category.

```
### Add outside inside information
inside_outside <- read.table("Inside_outside_diplodus.txt",header=FALSE,sep="\t")
mpa_env_inside_outside <- merge(x=mpa_env,y=inside_outside, by.x=c("INDIVIDUALS"), by.y=c("V1"))
reserve <- mpa_env_inside_outside$V5
### 6 nice colors for the MPA
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c")
bg <-  c("#ff7f00","darkred","#1f78b4","#ffff33","pink","#a6cee3","#33a02c","#e31a1c","blueviolet")

### Make a nice RDA plot
pdf("RDA_outside_inside.pdf",width=10,height=10)
plot(diplodus.rda, type="n", scaling=3)
points(diplodus.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(diplodus.rda, display="sites", cex=1.3, scaling=3, col = bg[as.numeric(eco)],pch =c(16, 17)[as.numeric(reserve)]) # the wolves
text(diplodus.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()

### Extract individuals information
rda_indv <- as.data.frame(scores(diplodus.rda, display=c("sites")))
write.table(rda_indv, "rda_diplodus_276ind.txt", quote=FALSE, row.names=TRUE)
```

#################### Identify candidate SNPs involved in local adaptation #########
### Extract the SNP loadings from the three significant constrained axes
load.rda <- scores(diplodus.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes

#histograms of the loadings on each RDA axis
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

### Identifying as many potential candidate loci as possible 
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

### Apply it to each significant constrained axis
cand1 <- outliers(load.rda[,1],3) # 82 for diplodus # 95 for mullus
cand2 <- outliers(load.rda[,2],3) # 69 for diplodus # 131 for mullus
cand3 <- outliers(load.rda[,3],3) # 131 for diplodus # 64 for mullus

### Calculate the total of candidates loci
ncand <- length(cand1) + length(cand2) + length(cand3)
ncand
write.table(ncand, "290outliers_mullus.txt", quote=FALSE, sep= )
### Organize our results by making one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
# Give colnames
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
# Paste all the candidates
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

### Add in the correlations of each candidate SNP with the three environmental predictors:
foo <- matrix(nrow=(ncand), ncol=3)  # 3 columns for 3 predictors
colnames(foo) <- c("salinity_surface_repro_min","chlo_surface_max","chlo_benthic_min")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

### Now we have a data frame of 282 candidate SNPs and their correlation with our 3 environmental predictors.
cand <- cbind.data.frame(cand,foo)  
head(cand)

### Investigate the duplicated selection
length(cand$snp[duplicated(cand$snp)]) # 0 duplicated loci
foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) # no duplicates on axis 2
table(foo[foo[,1]==3,2]) # no duplicates on axis 3
cand <- cand[!duplicated(cand$snp),] # remove duplicate detections

### See which of the predictors each candidate SNP is most strongly correlated with
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation
}

colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

### Most SNPs are associated with the chlorophylle variable
table(cand$predictor) 
#chlo_benthic_min           chlo_surface_max 
#127                         75 
#salinity_surface_repro_min 
#80 

### Plot the SNPs
sel <- cand$snp
env <- cand$predictor
env[env=="chlo_benthic_min"] <- '#1f78b4'
env[env=="salinity_surface_repro_min" ] <- '#a6cee3'
env[env=="chlo_surface_max"] <- '#6a3d9a'

### Color by predictor:
col.pred <- rownames(diplodus.rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("chr",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a')

# axes 1 & 2
plot(diplodus.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(diplodus.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(diplodus.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(diplodus.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("chlo_benthic_min","salinity_surface_repro_min","chlo_surface_max"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(wolf.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(wolf.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(wolf.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(wolf.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=c("AP","cvP","MDR","AMT","NDVI","Elev","sdT","Tree"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)