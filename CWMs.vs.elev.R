#<<TABLE OF CONTENTS>>
# 0 - Load & Prep Data
# 1 - Calculate CWM of each plot in 1970 and 2010 
# 2 - Traits vs elevation
#   2A - HERBIVORY LAYER - LMA & Leaf thickness weakly correlated with elevation

#   2B - CANOPY LAYER  

#<<WORKSPACES>>
wrk.dir<-("C:/Users/Julie/Desktop/Postdoc/Megantic Traits/Workspaces/") # Workspaces
data.dir<-(("C:/Users/Julie/Desktop/Postdoc/Megantic Traits/Data/")) # data
res.dir<-("C:/Users/Julie/Desktop/Postdoc/Megantic Traits/Results/")  # Results
grp.dir<-("C:/Users/Julie/Desktop/Postdoc/Megantic Traits/Graphs/")   # Graphs
fct.dir<-("C:/Users/Julie/Desktop/Postdoc/Megantic Traits/Functions/") # Functions

#<<LIBRARIES>>
library(car) # for powerTransform
library(vegan) # for decostand (standardizing data)
library(tree) # for tree()
library(MASS) # for glm.nb()  
library(mgcv) # for gam()
library(FD) # calculate CWMs and other Functional Diversity indices

#<<LOAD FILES>>

load(file=paste0(wrk.dir,'list.trait.names.herbivory.layer.Rdata')) # H.trait.names
load(file=paste0(wrk.dir,"Herbaceous.Layer.Traits.Standardized.RData")) # H.traits
load(file=paste0(wrk.dir,'species-level.traits.Herbaceous.layer.Rdata')) # sp.H.traits

load(file=paste0(wrk.dir,'list.trait.names.canopy.layer.Rdata')) # C.trait.names
load(file=paste0(wrk.dir,"Canopy.Layer.Traits.Standardized.RData")) # C.traits
load(file=paste0(wrk.dir,'species-level.traits.Canopy.layer.Rdata')) # sp.C.traits

# 0 - Load and prep data ####
#===========================#

# load elevation data
plot.elev<-read.csv(paste0(data.dir,'Plot_elevation_megantic.csv'),header=T)
rownames(plot.elev)<-paste0('X',plot.elev$plot)

# Create plot x sp matrices
#----

plot70.dat<-read.csv(paste0(data.dir,'Megantic_1970_vasc_2.csv'),header=T)
plot12.dat<-read.csv(paste0(data.dir,"Megantic_2012_vasc_2.csv"),header=T)

str(plot70.dat)
str(plot12.dat)
all(rownames(plot70.dat)==rownames(plot12.dat)) # TRUE
all(colnames(plot70.dat)==colnames(plot12.dat)) # TRUE

# Make sp names of sp.H.traits and sp.C.traits dataframes match sp names in plot70.dat and plot12.dat
#----

head(rownames(sp.H.traits)) # [1] "ARNU" "ARTR" "ATFI" "CAIN" "CASC" "CIAL"
head(plot70.dat$Accepted_name_Novembre_2016)# Abies_balsamea      Acer_pensylvanicum  Acer_rubrum

plot70.dat$Accepted_name_Novembre_2016==sp.names

# Split name to create species code
for(x in 1:length(plot70.dat$Accepted_name_Novembre_2016)){
  plot70.dat[x,'Genus']<-c(strsplit(as.character(plot70.dat$Accepted_name_Novembre_2016)[x],"_")[1][[1]][1])
  plot70.dat[x,'epi']<-c(strsplit(as.character(plot70.dat$Accepted_name_Novembre_2016)[x],"_")[1][[1]][2])
}

for(x in 1:length(plot12.dat$Accepted_name_Novembre_2016)){
  plot12.dat[x,'Genus']<-c(strsplit(as.character(plot12.dat$Accepted_name_Novembre_2016)[x],"_")[1][[1]][1])
  plot12.dat[x,'epi']<-c(strsplit(as.character(plot12.dat$Accepted_name_Novembre_2016)[x],"_")[1][[1]][2])
}

# Create species code to match trait data
for(x in 1:length(plot70.dat$Accepted_name_Novembre_2016)){
  plot70.dat[x,'sp.code']<-paste0(toupper(strtrim(plot70.dat$Genus[x],2)),toupper(strtrim(plot70.dat$epi[x],2)))
}

for(x in 1:length(plot12.dat$Accepted_name_Novembre_2016)){
  plot12.dat[x,'sp.code']<-paste0(toupper(strtrim(plot12.dat$Genus[x],2)),toupper(strtrim(plot12.dat$epi[x],2)))
}

# Recorder colummns
str(plot70.dat)
plot70.dat<-plot70.dat[c(1:2,54:56,3:5,6:53)]
rownames(plot70.dat)<-plot70.dat$sp.code
# non-unique values when setting 'row.names': ‘CADI’, ‘EUMA’, ‘OSCL’, ‘TACA’

# clean Duplicates
plot70.dat[plot70.dat$sp.code%in%c('CADI', 'EUMA', 'OSCL', 'TACA'),1:5]
#     Accepted_name_Novembre_2016     Nom_ABS_uniforme      Genus         epi sp.code
# 18               Carex_disperma       Carex disperma      Carex    disperma    CADI
# 36           Cardamine_diphylla    Dentaria diphylla  Cardamine    diphylla    CADI
# 46         Eupatorium_maculatum Eupatorium maculatum Eupatorium   maculatum    EUMA
# 47          Eurybia_macrophylla  Eurybia macrophylla    Eurybia macrophylla    EUMA
# 78          Osmorhiza_claytonii  Osmorhiza claytonii  Osmorhiza   claytonii    OSCL
# 80          Osmunda_claytoniana  Osmunda claytoniana    Osmunda claytoniana    OSCL
# 111        Taraxacum_campylodes Taraxacum officinale  Taraxacum  campylodes    TACA
# 112            Taxus_canadensis     Taxus canadensis      Taxus  canadensis    TACA

# don't keep either CADI - not in traits measured
plot70.dat[plot70.dat$Accepted_name_Novembre_2016=='Carex_disperma','sp.code']<-'CADIS'
plot70.dat[plot70.dat$Accepted_name_Novembre_2016=='Cardamine_diphylla','sp.code']<-'CADIP'
# don't keep either EUMA - - not in traits measured
plot70.dat[plot70.dat$Accepted_name_Novembre_2016=='Eupatorium_maculatum','sp.code']<-'EUPMA'
plot70.dat[plot70.dat$Accepted_name_Novembre_2016=='Eurybia_macrophyll','sp.code']<-'EURMA'
# Keep Osmorhiza claytonii
plot70.dat[plot70.dat$Accepted_name_Novembre_2016=='Osmunda_claytoniana','sp.code']<-'OSCLAY'
# don't Keep either TACA -- not in traits measured
plot70.dat[plot70.dat$Accepted_name_Novembre_2016=='Taraxacum_campylodes','sp.code']<-'TARCA'
plot70.dat[plot70.dat$Accepted_name_Novembre_2016=='Taxus_canadensis','sp.code']<-'TAXCA'

rownames(plot70.dat)<-plot70.dat$sp.code

str(plot12.dat)
plot12.dat<-plot12.dat[c(1:2,54:56,3:5,6:53)]
rownames(plot12.dat)<-plot12.dat$sp.code
# non-unique values when setting 'row.names': ‘CADI’, ‘EUMA’, ‘OSCL’, ‘TACA’ 

# don't keep either CADI
plot12.dat[plot12.dat$Accepted_name_Novembre_2016=='Carex_disperma','sp.code']<-'CADIS'
plot12.dat[plot12.dat$Accepted_name_Novembre_2016=='Cardamine_diphylla','sp.code']<-'CADIP'
# don't keep either EUMA
plot12.dat[plot12.dat$Accepted_name_Novembre_2016=='Eupatorium_maculatum','sp.code']<-'EUPMA'
plot12.dat[plot12.dat$Accepted_name_Novembre_2016=='Eurybia_macrophyll','sp.code']<-'EURMA'
# Keep Osmorhiza claytonii
plot12.dat[plot12.dat$Accepted_name_Novembre_2016=='Osmunda_claytoniana','sp.code']<-'OSCLAY'
# don't Keep either TACA
plot12.dat[plot12.dat$Accepted_name_Novembre_2016=='Taraxacum_campylodes','sp.code']<-'TARCA'
plot12.dat[plot12.dat$Accepted_name_Novembre_2016=='Taxus_canadensis','sp.code']<-'TAXCA'

rownames(plot12.dat)<-plot12.dat$sp.code

# Create empty plot x trait matrices
#----

H.CWT.70<-data.frame(matrix(nrow=48,ncol=16))
colnames(H.CWT.70)<-c('plot','elev',paste0('CWM.',names(sp.H.traits)),paste0('CWV.',names(sp.H.traits)))
rownames(H.CWT.70)<-colnames(plot70.dat[6:length(colnames(plot70.dat))])

H.CWT.12<-data.frame(matrix(nrow=48,ncol=16))
colnames(H.CWT.12)<-c('plot','elev',paste0('CWM.',names(sp.H.traits)),paste0('CWV.',names(sp.H.traits)))
rownames(H.CWT.12)<-colnames(plot12.dat[6:length(colnames(plot70.dat))])

C.CWT.70<-data.frame(matrix(nrow=48,ncol=16))
colnames(C.CWT.70)<-c('plot','elev',paste0('CWM.',names(sp.H.traits)),paste0('CWV.',names(sp.H.traits)))
rownames(C.CWT.70)<-colnames(plot70.dat[6:length(colnames(plot70.dat))])

C.CWT.12<-data.frame(matrix(nrow=48,ncol=16))
colnames(C.CWT.12)<-c('plot','elev',paste0('CWM.',names(sp.H.traits)),paste0('CWV.',names(sp.H.traits)))
rownames(C.CWT.12)<-colnames(plot12.dat[6:length(colnames(plot70.dat))])

# 1 - Calculate CWM of each plot in 1970 and 2010 ####
#====================================================#

## Super. important - missing trait data for coniferous species, which are dominant = problem. 
# Calculate CWMs only for herbaceous layer

# need to transpose plot70.dat because rows need to be plots
# need to subset species in plot70.dat to match list of species with traits

H.plot.by.sp.70<-t(plot70.dat[rownames(sp.H.traits),9:ncol(plot70.dat)])
dim(H.plot.by.sp.70) # 48 45

# Subset plot-by-species matrix to remove species with 0 abundance in 1970
H.plot.by.sp.70<-H.plot.by.sp.70[,colnames(H.plot.by.sp.70)%in%rownames(plot70.dat)]
dim(H.plot.by.sp.70) # 48 41

# Subset species-by-trait matrix to remove species with 0 abundance in 1970
H.sp.by.tr.70<-sp.H.traits[rownames(sp.H.traits)%in%rownames(plot70.dat),]
dim(H.sp.by.tr.70) # 41 7

# test compatibility of the two matrices
all(colnames(H.plot.by.sp.70)==rownames(H.sp.by.tr.70)) # TRUE

H.cwm.70<-functcomp(H.sp.by.tr.70,
                   H.plot.by.sp.70,
                   CWM.type='all')

cwm70.dat<-merge(plot.elev,H.cwm.70,by.y=0,by.x=0)
cwm70.dat$Row.names<-NULL # remove useless column entitled 'Row.names'
dim(cwm70.dat) # 48 9

save(H.cwm.70,file=paste0(wrk.dir,'1970.plot.cwm.herbaceous.layer.Rdata'))
save(H.plot.by.sp.70,file=paste0(wrk.dir,'plot.by.species.1970.herbaceous.layer.Rdata'))
save(H.sp.by.tr.70,file=paste0(wrk.dir,'sp.by.traits.1970.herbaceous.layer.Rdata'))
save(cwm70.dat,file=paste0(wrk.dir,'1970.plot.cwm.and.elevation.herbaceous.layer.Rdata'))

# 2 - Trait vs Elevation ####
#===========================#

  # 2.1 - Scatterplots - pairwise interactions? #### 
  
  pairs(cwm70.dat[,2:9])
  # LMA, Lamina.thck, myc.frac seem to be correlated with elevation. 
  
  plot(density(cwm70.dat$MeanElev)) #bimodal with peaks at 600m and 1000m
  
  # 2.2 - GAM - non-linearities?  ####
  
  par(mfrow=c(2,2))
  plot(gam(cwm70.dat$MeanElev~s(Ht.veg)+s(myc.frac)+
             s(Min.Root.Loca),data=cwm70.dat))
  
  #myc.frac appears and Min.Root.Loca appear correlated
  
  plot(gam(cwm70.dat$MeanElev~s(Ht.veg)+s(Lamina.thck)+
             s(LMA),data=cwm70.dat))
  # Lamina thickness appears non-linear
  
  plot(gam(cwm70.dat$MeanElev~s(Leaf.Area),data=cwm70.dat))
  # No relationships
  
  par(mfrow=c(1,1))
  
#=====================================================
  
  # Diagnostic plots
  
  # 1) Residuals vs Fitted - detects residual non-linear relationships between x & y variables 
  #    and heteroscedasticity (wedge)
  # 2) QQ plot - shows whether residuals are normally distributed (will follow straight line)
  # 3) Scale-Location - shows whether variance (residuals) increase with incrasing mean (also tests
  # for homoscedasticity). Another version of plot 1
  # 4) Residuals vs Leverage - shows whether individual datapoints have a lot of 'weight' on 
  # the regression
  
  # 2.3 - Regressions ####
  #===========================================================

  # Ht.veg
  
  summary(m1<-lm(cwm70.dat$Ht.veg~MeanElev+I(MeanElev^2),data=cwm70.dat)) # N.S.
  
  # Min.Root.Loca
  
  summary(m2<-lm(cwm70.dat$Min.Root.Loca~MeanElev+I(MeanElev^2),data=cwm70.dat))# Significant
  
  # Coefficients:
  #                 Estimate Std. Error t value Pr(>|t|)   
  # (Intercept)   -6.440e-01  1.492e+00  -0.432  0.66814   
  # MeanElev       1.076e-02  4.091e-03   2.629  0.01168 * 
  # I(MeanElev^2) -7.602e-06  2.683e-06  -2.833  0.00687 **
  #   ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # Residual standard error: 0.5368 on 45 degrees of freedom
  # Multiple R-squared:  0.1936,	Adjusted R-squared:  0.1578 
  # F-statistic: 5.402 on 2 and 45 DF,  p-value: 0.007893
  
  plot(cwm70.dat$Min.Root.Loca~MeanElev,data=cwm70.dat,
       xlab='Plot elevation (m)', ylab='Rooting Depth CWM',family='serif',pch=20)
  points(x=cwm70.dat$MeanElev, 
         y=fitted(m2),
         col='red', pch=20)
  text(900,4.0,labels= 'y= 0.01x -7e-06 x^2 \n AdjR2=0.16, p=0.008',family='serif',cex=0.8)
  
  # Lamina.thck
  
  summary(m3<-lm(cwm70.dat$Lamina.thck~MeanElev, data=cwm70.dat)) # Significant
  # Significant, AdjR2 = 0.15
  
  # Coefficients:
  #               Estimate Std. Error t value Pr(>|t|)   
  # (Intercept) -0.9309745  0.3343314  -2.785  0.00775 **
  #  MeanElev    0.0014116  0.0004623   3.053  0.00376 **
  #  ---
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # Residual standard error: 0.542 on 46 degrees of freedom
  # Multiple R-squared:  0.1685,	Adjusted R-squared:  0.1504 
  # F-statistic: 9.323 on 1 and 46 DF,  p-value: 0.003755
  
  plot(cwm70.dat$Lamina.thck~MeanElev,data=cwm70.dat,
       xlab='Plot elevation (m)', ylab='Leaf thickness',family='serif',pch=20)
  points(x=cwm70.dat$MeanElev, 
         y=fitted(m3),
         col='red', pch=20)
  text(650,1.3,labels= 'y= -0.93 + 0.001x \n AdjR2=0.15, p=0.003',family='serif',cex=0.9)
  
  # LMA
  
  summary(m4<-lm(cwm70.dat$LMA~MeanElev,data=cwm70.dat))
  # significant, AdjR2=0.08
  
  # Coefficients:
  #              Estimate Std. Error t value Pr(>|t|)  
  # (Intercept) -0.3332923  0.1617211  -2.061   0.0450 *
  #  MeanElev    0.0005099  0.0002236   2.280   0.0273 *
  # 
  # Residual standard error: 0.2622 on 46 degrees of freedom
  # Multiple R-squared:  0.1015,	Adjusted R-squared:  0.08201 
  # F-statistic: 5.199 on 1 and 46 DF,  p-value: 0.02728
  
  plot(cwm70.dat$LMA~MeanElev,data=cwm70.dat,
       xlab='Plot elevation (m)', ylab='LMA',family='serif',pch=20)
  points(x=cwm70.dat$MeanElev, 
         y=fitted(m4),
         col='red', pch=20)
  text(650,1.3,labels= 'y= -0.93 + 0.001x \n AdjR2=0.15, p=0.003',family='serif',cex=0.9)
  
  # find outlier point
  text(x=cwm70.dat$MeanElev, y=cwm70.dat$LMA, labels = cwm70.dat$plot) # plot 02 has a very low CWM-LMA
  
  # run model without that outlier
  summary(lm(cwm70.dat[cwm70.dat$plot!='62',]$LMA~MeanElev,data=cwm70.dat[cwm70.dat$plot!='62',]))
  # model not significnat without outlier
  
  # LDMC
  
  summary(m5<-lm(cwm70.dat$LDMC~MeanElev,data=cwm70.dat)) # N.S.
  
  # Leaf.Area
  
  summary(m6<-lm(cwm70.dat$Leaf.Area~MeanElev,data=cwm70.dat)) # p=val=0.008
  
  # Coefficients:
  # Estimate Std. Error t value Pr(>|t|)   
  # (Intercept) -0.7371092  0.2327403  -3.167  0.00273 **
  # MeanElev     0.0008945  0.0003218   2.780  0.00786 **
  # ---
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # Residual standard error: 0.3773 on 46 degrees of freedom
  # Multiple R-squared:  0.1438,	Adjusted R-squared:  0.1252 
  # F-statistic: 7.726 on 1 and 46 DF,  p-value: 0.007859
  
  plot(cwm70.dat$Leaf.Area~MeanElev,data=cwm70.dat,
       xlab='Plot elevation (m)', ylab='Leaf Area',family='serif',pch=20)
  points(x=cwm70.dat$MeanElev, 
         y=fitted(m6),
         col='red', pch=20)
  text(700,0.6,labels= 'y= -0.73 + 0.001x  \n AdjR2=0.13, p=0.008',family='serif',cex=0.9)
    
  # myc.frac
  
  summary(m7<-lm(cwm70.dat$myc.frac~MeanElev+I(MeanElev^2),data=cwm70.dat)) # NS
  plot(cwm70.dat$myc.frac~,data=cwm70.dat)
  text(cwm70.dat$MeanElev,cwm70.dat$myc.frac,labels=cwm70.dat$plot)
