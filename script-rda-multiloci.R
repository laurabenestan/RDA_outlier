# Download libraries
library(psych)    
library(vegan)
library(adegenet)
library(dplyr)
library(fmsb)
library(gsl)

### Download data in plink format
#plink_diplodus <- read.PLINK("18512snps-276ind-diplodus.raw", map.file = "18512snps-276ind-diplodus.map", quiet = FALSE)
plink_diplodus <- read.table("18512snps-276ind-diplodus.raw", header=TRUE, sep=" ", row.names=1)
plink_mullus <- read.table("14318snps_312ind.raw", header=TRUE, sep=" ", row.names=1)
dim(plink_diplodus)

### Remove names
gen <- plink_diplodus[,6:37029]
gen <- plink_mullus[,6:28641]

### Calculate the NA
sum(is.na(gen))# 553,232 NAs in the matrix (~3% missing data)

### Fill the NAs
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp)) # No NAs

### Download environmental data 
env <- read.table("276ind-24env-variables.txt",sep="\t",header=T) # for diplodus
env <- read.table("24env_variables_mullus.txt",sep="\t",header=T) # for mullus
str(env)

### Add habitats variable
habitat <- read.table("habitat_diplodus_276ind.txt",sep="\t",header=T)
env <- cbind(env, habitat)

### Make individual names characters (not factors)
env$labels <- as.character(env$labels)

### Visualize correlations among predictors
pairs.panels(env[,2:25], scale=T)

### Create a correlation matrix
matrix_cor <- cor(env[,2:25])

### Correction for VIF
# Call the vif function
vif_func<-function(in_frame,thresh=10,trace=T,...){
  library(fmsb)
  if(any(!'data.frame' %in% class(in_frame))) in_frame<-data.frame(in_frame)
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    in_dat<-in_frame
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      vif_vals<-NULL
      var_names <- names(in_dat)
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      vif_max<-as.numeric(vif_vals[max_row,2])
      if(vif_max<thresh) break
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
    }
    return(names(in_dat))
  }
}
vif_func(in_frame=matrix_cor,thresh=5,trace=T)

### Selection for the non-correlated variables
env_selected <- select(env, salinity_surface_repro_min,chlo_surface_max,chlo_benthic_min) # for diplodus
env_selected <- select(env, salinity_water_column_whole_yr_max,chlo_surface_max,chlo_benthic_max,chlo_benthic_min)# for mullus
table(env_selected$salinity_water_column_whole_yr_max)

#The distribution of factor levels is highly skewed towards classes around 36
#leaving the remaining classes with small numbers of individuals. These characteristics make it unlikely to be a highly informative predictor...

### Visualize correlations among selected predictors
pairs.panels(env_selected[,1:4], scale=T)
pred <- env_selected[,1:3]

################# RUN THE RDA ###############
diplodus.rda <- rda(gen.imp ~ ., data=pred, scale=T)
diplodus.rda 
# Our constrained ordination explains about 0.02% of the variation; this low explanatory power is not surprising given that we expect that most of the SNPs in our dataset will not show a relationship with the environmental predictors (e.g., most SNPs will be neutral).

### Calculate the adjusted R2
RsquareAdj(diplodus.rda)

### The eigenvalues for the constrained axes reflect the variance explained by each canonical axis
summary(eigenvals(diplodus.rda, model = "constrained")
screeplot(diplodus.rda)

load.rda <- scores(diplodus.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes
##  we can see that the first three constrained axes explain all the variance.
signif.full <- anova.cca(diplodus.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full
# 0.01 ** 
#0.49 for mullus

### Plot the RDA
plot(diplodus.rda, scaling=3) 
plot(diplodus.rda, choices = c(1, 3), scaling=3)

### Add the MPA information
mpa_group <- read.table("population-map-276ind-diplodus-mpa.txt",header=TRUE,sep="\t")
mpa_group <- read.table("population-map-335ind-mullus-mpa.txt",header=TRUE,sep="\t")
mpa_env <- cbind(env, mpa_group)
eco <- mpa_group$STRATA

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