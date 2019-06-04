# ReserveBenefit--RDA_adaptation

# Stacks Workflow

RADseq workflow using [STACKS](http://creskolab.uoregon.edu/stacks/)

Developed by [Laura Benestan](https://github.com/laurabenestan) in
[Stephanie Manel](https://sites.google.com/site/stephaniemanel/home)'s
laboratory.

## Concept and definition
This README.txt is widely inspired by the vignette created by Brenna Forester at [popgennescent](https://popgen.nescent.org/2018-03-27_RDA_GEA.html).
We used of Redundancy Analysis (RDA) as a genotype-environment association (GEA) method to simultaneously assess the percent of genomic variation explained by environmental variables and to detect loci under selection (see the relevant paper of Forester et al., 2018). 
RDA is a two-step analysis in which genetic and environmental data are analyzed using multivariate linear regression. 
Then PCA of the fitted values is used to produce canonical axes, which are linear combinations of the predictors (Legendre & Legendre, 2012). 
here, we performed the RDA on a individual-based sampling design.

## Application to our dataset
Here, we apply RDA to genomic data from 276 individuals of the white seabream (Diplodus sargus) sampled across the Mediterranean sea. 
Results of the RDA at the full set of 18,512 single nucleotide polymorphism (SNP) markers will be soon available. 
The goal of this study is to understand which environmental factor may drive the genomic divergence observed in this species. 
In this case, the data is individual-based, and inputs are allele counts (i.e. 0/1/2) for each locus and individual fish. 

# Using R to perform the analysis
## Install packages
First, we install the necessary packages and then download the corresponding libraries.

``` {r}
library(psych)    
library(vegan)
library(adegenet)
library(dplyr)
library(fmsb)
library(gsl)
```
## Read genetic data

```{r}
plink_diplodus <- read.table("18512snps-276ind-diplodus.raw", header=TRUE, sep=" ", row.names=1)
dim(plink_diplodus)
```

## Prepare genetic data

We remove the names of the samples as we don't consider these values when we will fill the missing values.
Yet other way to handle the missing values and before doing this step visualize patterns of missingness and perform map-independent imputations of missing genotypes using [grur package](https://github.com/thierrygosselin/grur)

```{r}
gen <- plink_diplodus[,6:37029]
```

We calculate the percent of missing data.
Multivariate analysis are very sensitive to missing data and these missing needs to be filled.
```{r}
sum(is.na(gen))# 553,232 NAs in the matrix (~3% missing data)
```

Missing values are then filled using the overall average of the allele frequency variation.
```{r}
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp)) # No NAs
```

## Read the environmental data

We download environmental data in R.
```{r}
env <- read.table("276ind-24env-variables.txt",sep="\t",header=T) # for diplodus
str(env)
```

## Prepare the environmental data

We change the class of the labels.
```{r}
env$labels <- as.character(env$labels)
```

We remove one of the highly correlated predictors by visually and quantitatively checking the correlation among predictors using the function called `pairs.panels`.

```{r}
pairs.panels(env[,2:25], scale=T)
```

![Correlations among our predictors.](Correlation_predictor.png)

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
# Run the RDA
diplodus.rda <- rda(gen.imp ~ ., data=pred, scale=T)
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
### 6 and 9 nice colors for the MPA
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
```

![RDA considering two categories.](RDA_both_species.pdf)

```
### Extract individuals information
rda_indv <- as.data.frame(scores(diplodus.rda, display=c("sites")))
write.table(rda_indv, "rda_diplodus_276ind.txt", quote=FALSE, row.names=TRUE)
```
