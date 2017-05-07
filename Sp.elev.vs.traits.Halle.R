#<<TABLE OF CONTENTS>>
# 0 - Bin plots into elevation categories
# 1 - Calculate mean elevation of each species
# 2 - Create dataframes
# A - HERBIVORY LAYER
# A3- Data Exploration - Scatterplots, Tree models and GAMs
# A4 - Model Selection
#   A4.1 - lm - elev ~ traits. Traits do not predict herbaceous species elevation
#   A4.2 - glm - elev ~ traits. Traits do not predict herbaceous species elevation
# B - CANOPY LAYER
# B3 - Data Exploration - Scatterplots, Tree models and GAMs
# B4 - Model selection ####
#   B4.1 lm - elev ~ traits - Traits do not predict canopy species elevation !
#   B4.2 glm, gamma link fct - elev ~ traits - Traits do not predict canopy species elevation!

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

#<<LOAD FILES>>

load(file=paste0(wrk.dir,'list.trait.names.herbivory.layer.Rdata')) # H.trait.names
load(file=paste0(wrk.dir,"Herbaceous.Layer.Traits.Standardized.RData")) # H.traits
load(file=paste0(wrk.dir,'species-level.traits.Herbaceous.layer.Rdata')) # sp.H.traits

load(file=paste0(wrk.dir,'list.trait.names.canopy.layer.Rdata')) # C.trait.names
load(file=paste0(wrk.dir,"Canopy.Layer.Traits.Standardized.RData")) # C.traits
load(file=paste0(wrk.dir,'species-level.traits.Canopy.layer.Rdata')) # sp.C.traits

# load(file=paste0(wrk.dir,'Species.abundances.full.data.Rdata')) # object = abund
load(file=paste0(wrk.dir,'Species.abundances.Inf.removed.Rdata')) # object = abund.c
#=====================================================================================

# 0 - Bin plots into elevation categories ####
#============================================#

  #0.1 - Presence/Absence at each elevation ####
  #============================================#

  # Hall = presence/absence of species along small plots spread along elevational gradient
  Hall<-read.csv(paste0(data.dir,'G.Hall_sp_by_site_corrected.csv'),header=T)
  sp.names<-colnames(Hall)[,12:ncol(Hall)]
  
  # Create vector of intervals
  intervals<-seq(from=400,to=1060,by=30)
  length(intervals) # 23
  
  # Create vector of lower and upper intervals
  lower.int<-intervals[1:22]+1
  upper.int<-intervals[2:23]
  mid.int<-(lower.int+upper.int)/2
  
  # Create empty matrix to fill - currently using 14 bins of 50 m elevation.
  
  Hall.bins.pa<-as.data.frame(matrix(NA,nrow=22,ncol=ncol(Hall)-11))
  colnames(Hall.bins.pa)<-sp.names
  elev.range<-paste0(lower.int,'-',upper.int)
  no.plots<-rep(NA,times=22)
  Hall.bins.pa<-cbind(elev.range,lower.int,upper.int,mid.int,no.plots,Hall.bins.pa)
  
  # fill empty dataframe
  
  for (i in 1:nrow(Hall.bins.pa)) { # loop through each elevation intervals (currently 14), i
    l<-Hall.bins.pa[i,'lower.int']
    u<-Hall.bins.pa[i,'upper.int']
    subset.Hall<-Hall[which(Hall$Elev_m >= l & Hall$Elev_m <= u),] # create a data subset with only plots in elev. bin
    
    for(s in sp.names) { # loop through each sp in the colnames (459 spp.), s
      
      ifelse(
        test= any(subset.Hall[,s]==1)==TRUE,               # if this test is true - if sp. present in any plot in subset
        yes = Hall.bins.pa[Hall.bins.pa$elev.range==elev.range[i],s]<-1,   # set the value for that elevation bin for that species to 1
        no= Hall.bins.pa[Hall.bins.pa$elev.range==elev.range[i],s]<-0      # else, set it to 0
      )
      
      Hall.bins.pa[elev.range==elev.range[i],'no.plots']<-nrow(subset.Hall) # count no. of plots into each elevation bin
      
    }
  }
  
  head(Hall.bins.pa)[1:10] # check that all looks good.
  save(Hall.bins.pa, file=paste0(wrk.dir,'Hall.data.plots.binned.by.elevation.slice.Rdata'))
  
  # 0.2 - Relative abundance at each elevation ####
  #===============================================#
  
  # Create empty matrix to fill - currently using 22 bins of 30 m elevation.
  
  intervals.by.sp<-as.data.frame(matrix(NA,nrow=22,ncol=ncol(Hall)-11))
  colnames(intervals.by.sp)<-sp.names
  elev.range<-paste0(lower.int,'-',upper.int)
  no.plots<-rep(NA,times=22)
  Hall.bins.prop.presence<-cbind(elev.range,lower.int,upper.int,mid.int,no.plots,intervals.by.sp)
  
  # fill empty dataframe
  
  for (i in 1:nrow(Hall.bins.prop.presence)) { # loop through each elevation intervals (currently 22), i
    l<-Hall.bins.prop.presence[i,'lower.int']
    u<-Hall.bins.prop.presence[i,'upper.int']
    subset.Hall<-Hall[which(Hall$Elev_m >= l & Hall$Elev_m <= u),] # create a data subset with only plots in elev. bin
    
    for(s in sp.names) { # loop through each sp in the colnames (459 spp.), s
      
      Hall.bins.prop.presence[elev.range==elev.range[i],'no.plots']<-nrow(subset.Hall) # count no. of plots in each elevation bin
      
      elev.presence<-sum(subset.Hall[,s])   # Number of plots where species s occurs within elevation bin
      p<-elev.presence/Hall.bins.prop.presence$no.plots[i] # proportion of plots with sp. s present of a given elevation bin
      Hall.bins.prop.presence[Hall.bins.prop.presence$elev.range==elev.range[i],s]<-p 
      # assign proportion of species occurence for each elevation bin 
      
    }
  }
  
  head(Hall.bins.prop.presence)[1:10]
  
# 1 - Calculate mean elevation of each species ####
#=================================================#
  
  # s- species s
  # i- elevation bin i
  # prop.pres(s.i)*elev(i)/sum(prop.pres(i))

# Calculate mean elevation for each species 

mean.elev<-data.frame(matrix(ncol=5,nrow=459))
names(mean.elev)<-c('Species','avg.elev.pa','mid.elev','wghtd.elev')

# Calculate mean species elevation
for(s in sp.names) {
  mean.elev[which(sp.names==s),'Species']<-s                  # Assign species name
  
  # Mean elevation based on presence/absence data
  mean.elev[which(mean.elev$Species==s),'avg.elev.pa']<-
    round(sum(Hall.bins.pa$mid.int * Hall.bins.pa[,s]) / 
            sum(Hall.bins.pa[,s]),digits=0)           
  
  # mid-range of elevation as median between min and max
  max.e<-max(Hall.bins.pa[Hall.bins.pa[,s]==1,'mid.int'])
  min.e<-min(Hall.bins.pa[Hall.bins.pa[,s]==1,'mid.int'])
  mean.elev[which(mean.elev$Species==s),'mid.elev']<-min.e+((max.e-min.e)/2)
                                                              
  # mean elevation weighted by relative abundace at each elevation 
  mean.elev[which(mean.elev$Species==s),'wghtd.elev']<-sum(
    Hall.bins.prop.presence$mid.int *
      Hall.bins.prop.presence[,s]
    ) / sum(Hall.bins.prop.presence[,s])

# Split name to create species code
for(x in 1:length(sp.names)){
  mean.elev[x,'Genus']<-c(strsplit(sp.names[x],"_")[1][[1]][1])
  mean.elev[x,'epi']<-c(strsplit(sp.names[x],"_")[1][[1]][2])
}
# Create species code to match trait data
for(x in 1:length(sp.names)){
mean.elev[x,'sp.code']<-paste0(toupper(strtrim(mean.elev$Genus[x],2)),toupper(strtrim(mean.elev$epi[x],2)))
}

str(mean.elev)
# re-order columns
mean.elev<-mean.elev[c('Species','Genus','epi','sp.code','avg.elev.pa','mid.elev','wghtd.elev','wghtd.elev2')]
str(mean.elev)
  # data.frame':	459 obs. of  5 variables:
  # $ Species: chr  "Abies_balsamea" "Acer_pensylvanicum" "Acer_rubrum" "Acer_saccharum" ...
  # $ Genus  : chr  "Abies" "Acer" "Acer" "Acer" ...
  # $ epi    : chr  "balsamea" "pensylvanicum" "rubrum" "saccharum" ...
  # $ sp.code: chr  "ABBA" "ACPE" "ACRU" "ACSA" ...
  # $ avg.elev.pa   : num  625 582 523 572 622 ...
  # $ mid.elev: num  744 676 660 680 742 ...
  # $ wghtd.elev2: num  17792 9783 7295 9216 16081 ...

save(mean.elev,file=paste0(wrk.dir,'species.mean.elevation.based.on.Hall.plot.surveys.RData'))

# 2- Create dataframes ####
#=========================#

# Herbaceous dataframe - merge trait data and elevation based on species code
#----

H.dat<-merge(sp.H.traits.c,mean.elev,by.y='sp.code',by.x=0) # '0' indicates rownames
str(H.dat)
H.dat[,c('Species','Genus','epi')]<-NULL
row.names(H.dat)<-H.dat$Row.names
# Error - non-unique values of ‘CASC’, ‘COCA’, ‘COTR’, ‘POPU’ 

# clean Duplicates
mean.elev[mean.elev$sp.code%in%c('CASC', 'COCA', 'COTR', 'POPU'),]
  #                   Species        Genus        epi sp.code avg.elev.pa
  # 84         Carex_scabrata        Carex   scabrata    CASC  599
  # 85         Carex_scoparia        Carex   scoparia    CASC  550
  # 112     Conyza_canadensis       Conyza canadensis    COCA  458
  # 113       Coptis_trifolia       Coptis   trifolia    COTR  592
  # 115  Corallorhiza_trifida Corallorhiza    trifida    COTR  530
  # 117     Cornus_canadensis       Cornus canadensis    COCA  669
  # 317 Polygonatum_pubescens  Polygonatum  pubescens    POPU  583
  # 330  Potamogeton_pusillus  Potamogeton   pusillus    POPU  490

# Keep Carex scabrata
mean.elev[mean.elev$Species=='Carex_scoparia','sp.code']<-'CASCO'
# Keep Cornus canadensis
mean.elev[mean.elev$Species=='Conyza_canadensis','sp.code']<-'COCAN'
# Keep Polygonatum pubescens
mean.elev[mean.elev$Species=='Potamogeton_pusillus','sp.code']<-'POPUS'
# Keep Coptis trifolia
mean.elev[mean.elev$Species=='Corallorhiza_trifida','sp.code']<-'COTRI'

# Overwrite object with new changes
save(mean.elev,file=paste0(wrk.dir,'species.mean.elevation.based.on.Hall.plot.surveys.RData'))

# Recreate database
H.dat<-merge(sp.H.traits.c,mean.elev,by.y='sp.code',by.x=0) # '0' indicates rownames
str(H.dat)
H.dat[,c('Species','Genus','epi')]<-NULL

row.names(H.dat)<-H.dat$Row.names
H.dat$Row.names<-NULL

plot(density(H.dat$avg.elev.pa)) # left-skewed
shapiro.test((H.dat$avg.elev.pa)) # p-val = 0.01 ! Not bad!
plot(density(H.dat$wghtd.elev2)) # This carries same problem of unevely sampled no. of plots along elevation 
plot(H.dat$avg.elev.pa,H.dat$wghtd.elev2)
     #xlim=c(500,800),ylim=c(500,800))
# points(x=H.dat$mid.elev,
#        y=predict(lm(H.dat$avg.elev.pa~H.dat$mid.elev)),
#        col='red', pch=20) 
# abline(a=0,b=1)

summary(lm(H.dat$avg.elev.pa~H.dat$wghtd.elev2)) #AdjR2 = 0.30
# These two elevations are almost the same!


save(H.dat,file=paste0(wrk.dir,'Herbaceous.layer.dataframe.with.species.mean.elevation.and.traits.RData'))

# Canopy dataframe - merge trait data and elevation based on species code
#----

C.dat<-merge(sp.C.traits.c,mean.elev,by.y='sp.code',by.x=0) # '0' indicates rownames
str(C.dat)
C.dat[,c('Species','Genus','epi')]<-NULL

row.names(C.dat)<-C.dat$Row.names
# Error - non-unique values when setting 'row.names': ‘ACSP’, ‘POTR’

# Clean Duplicates
mean.elev[mean.elev$sp.code%in%c('ACSP', 'POTR'),]
  #                 Species   Genus         epi sp.code avg.elev.pa
  # 5         Acer_spicatum    Acer    spicatum    ACSP  622
  # 7             Actaea_sp  Actaea          sp    ACSP  593
  # 316       Poa_trivialis     Poa   trivialis    POTR  549
  # 328 Populus_tremuloides Populus tremuloides    POTR  528

# Keep Acer spicatum
mean.elev[mean.elev$Species=='Actaea_sp','sp.code']<-'ACTSP'
# Keep Populus tremuloides
mean.elev[mean.elev$Species=='Poa_trivialis','sp.code']<-'POTRI'

# Overwrite object with new changes
save(mean.elev,file=paste0(wrk.dir,'species.mean.elevation.based.on.Hall.plot.surveys.RData'))

# Recreate database
C.dat<-merge(sp.C.traits.c,mean.elev,by.y='sp.code',by.x=0) # '0' indicates rownames
str(C.dat)
C.dat[,c('Species','Genus','epi')]<-NULL

row.names(C.dat)<-C.dat$Row.names
C.dat$Row.names<-NULL

plot(density(C.dat$avg.elev.pa)) # looks quite normal 
shapiro.test((C.dat$avg.elev.pa)) # p-val = 0.10 ! Woo-hoo! Don't transform.
plot(density(C.dat$mid.elev)) # looks quite normal except for a little right-skewed
shapiro.test((C.dat$avg.elev.pa)) # p-val = 0.20 ! Woo-hoo! Don't transform.

save(C.dat,file=paste0(wrk.dir,'Canopy.layer.dataframe.with.species.mean.elevation.and.traits.RData'))


# A - HERBIVORY LAYER ####
#==========================#

  # A3 - Data Exploration ####
  #==========================#

        # A3.1 Scatterplots - pairwise interactions? ####
        #----

        names(H.dat) 
        # [1] "Ht.veg"        "Min.Root.Loca" "Lamina.thck"   "LMA"           "LDMC"          "Leaf.Area"     "myc.frac"     
        # [8] "avg.elev.pa"   "mid.elev"      "wghtd.elev"    "wghtd.elev2"

        pairs(H.dat[,c(1:7,11)],panel=panel.smooth)
        # Variance seems homoscedastic
        # correlations among Leaf.Area-Ht.veg; Lamina.thck-LMA; LDMc-LMA

        # A3.2 - Tree models - non-linearities ####
        #----
        
        elev.tree.model<-tree(H.dat$avg.elev.pa~Ht.veg+Min.Root.Loca+Lamina.thck+LMA+LDMC+Leaf.Area+myc.frac,
                              data=H.dat)
        plot(elev.tree.model); text(elev.tree.model) ; title('herbaceous layer, avg.elev.pa vs traits')
        # LDMC is most important variable, 
        # for those with high LDMC, LDMC matters
        # For those with low LDMC min.root.loca matters
        elev.tree.model #  null deviance = 60840
        1-(deviance(elev.tree.model)/60840) # 0.54
        
        # A3.3 GAM - non-linearities?  ####
        #----
        
        par(mfrow=c(2,2))
        plot(gam(H.dat$avg.elev.pa~s(Ht.veg)+s(Min.Root.Loca)
                 ,data=H.dat))
        # no curvature
        
        plot(gam(H.dat$avg.elev.pa~s(Lamina.thck)+s(LMA)
                 ,data=H.dat))
        # no curvature
        
        plot(gam(H.dat$avg.elev.pa~s(LDMC)+s(Leaf.Area),
                 data=H.dat))
        # no curvature     
        
        plot(gam(H.dat$avg.elev.pa~s(myc.frac),
                 data=H.dat))
        # no curvature      
        
    #=====================================================
    
    # Diagnostic plots
    
    # 1) Residuals vs Fitted - detects residual non-linear relationships between x & y variables 
    #    and heteroscedasticity (wedge)
    # 2) QQ plot - shows whether residuals are normally distributed (will follow straight line)
    # 3) Scale-Location - shows whether variance (residuals) increase with incrasing mean (also tests
    # for homoscedasticity). Another version of plot 1
    # 4) Residuals vs Leverage - shows whether individual datapoints have a lot of 'weight' on 
    # the regression
    
    # A4 - Model selection ####
    #===========================================================
        
        #  A4.1 lm - elev ~ traits - nothing significant ! ####
        
        # A) Check trait-trait interactions ####
        #----
        
        # short on degrees of freedom. Do it in two models.
        summary(lm(H.dat$avg.elev.pa~(Ht.veg+Min.Root.Loca+Lamina.thck+LMA+LDMC+Leaf.Area)^2,
                   data=H.dat))
        #LMA:Leaf.Area  marginally significant
        
        summary(lm(H.dat$avg.elev.pa~(myc.frac+Min.Root.Loca+Lamina.thck+LMA+LDMC+Leaf.Area)^2,
                   data=H.dat))
        
        # nothing significant
        
        # B) Check non-linear relationships ####
        #----
        
        summary(lm(H.dat$avg.elev.pa~.+I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
                     I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
                data=H.dat))
        # 
        # I(LMA^2)  significant
        
        # C) Fit best minimal model - not significant####
        #----
        
        H1<-lm(H.dat$avg.elev.pa~Ht.veg+Min.Root.Loca+Lamina.thck+LMA+LDMC+I(LMA^2),data=H.dat)
        H0<-lm(H.dat$avg.elev.pa~1,data=H.dat)
        
        elev.best.lm.forw<-step(H0,scope=list(lower=H0,upper=H1),trace=T,direction='both')
        summary(elev.best.lm.forw) # Min.Root.Loca retained, but not significant
        AIC(elev.best.lm.forw) #307.97
        
        # Null model just as good

        elev.best.lm.back<-step(H1,scope=list(lower=H0,upper=H1),trace=T,direction='both')
        summary(elev.best.lm.back) # Min.Root.Loca

        # Coefficients:
        #               Estimate Std. Error t value Pr(>|t|)    
        # (Intercept)    710.341     21.896  32.441   <2e-16 ***
        # Min.Root.Loca  -11.319      7.869  -1.438    0.162 
        # model not significant
        
        # Plot results 
        par(mfrow=c(1,1))
        
        plot(H.dat$avg.elev.pa,H.dat$Ht.veg,
             pch=19, mgp=c(2,1,0),
             xlab='',ylab= 'Plant Height (cm)',family='serif')
         
        plot(H.dat$avg.elev.pa,H.dat$LMA,
             pch=19, mgp=c(2,1,0),
             xlab='',ylab= 'LMA (mg/cm2)',family='serif')
        
        plot(H.dat$avg.elev.pa,H.dat$Leaf.Area,
             pch=19, mgp=c(2,1,0),
             xlab='',ylab= 'Leaf Area (cm2)',family='serif')
        
        plot(H.dat$avg.elev.pa,H.dat$Min.Root.Loca,
             pch=19, mgp=c(2,1,0),
             xlab='Species mean elevation (m)',ylab= 'Rooting depth',family='serif')
        
        plot(H.dat$avg.elev.pa,H.dat$myc.frac,
             pch=19, mgp=c(2,1,0),
             xlab= 'Species mean elevation (m)', ylab= 'fraction of roots with \n mycorhizal associations',family='serif')
        
        # A4.2 glm, gamma link fct - 
        
        # A) Check trait-trait interactions ####
        # ----
        
        summary(glm(H.dat$avg.elev.pa~.+(Min.Root.Loca+Ht.veg+Min.Root.Loca+Lamina.thck+LMA+
                                                  LDMC+Leaf.Area+myc.frac)^2,
                    data=H.dat,family=Gamma(link='log')),dispersion=1)
        
        # nothing significant
        
        # B) Check higher-order effects ####
        #----
          
        summary(glm(H.dat$avg.elev.pa~I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
                      I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
                    data=H.dat),family=Gamma(link='log'))  
        
        # nothing significant
        
        # C) Fit best model ####
        #----
        
        H1<-glm(H.dat$avg.elev.pa~.,data=H.dat,family= Gamma(link='log'))
        H0<-glm(H.dat$avg.elev.pa~1,data=H.dat,family=Gamma(link='log'))
        
        summary(elev.best.glm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=F),
               dispersion=1)
        
        # nothing significant
        
      # CANOPY LAYER
        
        # B3.1 Scatterplots - pairwise interactions? ####
        #----
        
        names(C.dat) 
        # [1] "Lamina.thck" "LMA"         "LDMC"        "Leaf.Area"   "elev"   
        
        pairs(C.dat,panel=panel.smooth)
        # Variance seems homoscedastic
        # yes, Lamina.thck-LMA; LDMc-LMA
        
        # B3.2 - Tree models - non-linearities ####
        #----
        
        C.elev.tree.model<-tree(C.dat$avg.elev.pa~.,data=C.dat)
        plot(C.elev.tree.model)
        text(C.elev.tree.model)
        title('canopy layer, elev vs traits')
        # LDMC is most important variable, 
        # for those with high LDMC, high LDMC matters
        C.elev.tree.model #  null deviance = 38840
        1-(deviance(C.elev.tree.model)/38840) # 0.55
        
        # B3.3 GAM - non-linearities?  ####
        #----
        
        plot(gam(C.dat$avg.elev.pa~s(Lamina.thck),data=C.dat))
        # curvy 
        
        plot(gam(H.dat$avg.elev.pa~s(LDMC),data=H.dat))
        # no relationship     
        
        plot(gam(H.dat$avg.elev.pa~s(LMA),data=H.dat))
        # negative ?
        
        plot(gam(H.dat$avg.elev.pa~s(Leaf.Area),data=H.dat))
        # positive?
        
        # B4 - Model selection ####
        #===========================================================
        
        #  B4.1 lm - elev ~ traits - nothing significant ! ####
        
        # A) Check trait-trait interactions ####
        #----
        
        # short on degrees of freedom. Do it in two models.
        summary(lm(C.dat$avg.elev.pa~(Lamina.thck+LMA+LDMC+Leaf.Area)^2,
                   data=C.dat))
        
        # LMA, LDMC & all interactions significant!!!
        
        vif(lm(C.dat$avg.elev.pa~(Lamina.thck+LMA+LDMC+Leaf.Area)^2,
               data=C.dat)) 
        
        # LMA, LDMC and Lamina.thck:LDMC are redundant
        
        # B) Check non-linear relationships ####
        #----
        
        summary(lm(C.dat$avg.elev.pa~Lamina.thck+I(Lamina.thck^2)+LMA+I(LMA^2)+LDMC+I(LDMC^2)+Leaf.Area+I(Leaf.Area^2),
                   data=C.dat))
        
        # nothing significant
        
        # C) Fit  minimal model ####
        #----
        
        H1<-lm(C.dat$avg.elev.pa~.^2,data=C.dat)
        H0<-lm(C.dat$avg.elev.pa~1,data=C.dat)
        
        C.elev.best.lm<-step(H0,scope=list(upper=H1,lower=H0),direction="both")
        summary(C.elev.best.lm)
        
        #             Estimate Std. Error t value Pr(>|t|)    
        # (Intercept)   566.71      11.68  48.529   <2e-16 ***
        # LDMC          -26.05      15.22  -1.712    0.108    
        # ---

        # Residual standard error: 46.54 on 15 degrees of freedom
        # Multiple R-squared:  0.1634,	Adjusted R-squared:  0.1076 
        # F-statistic:  2.93 on 1 and 15 DF,  p-value: 0.1075
        
        # starting with empty model, nothing is significant!
        # Null model is retained
        AIC(H0) 
        # 185
        
        C.elev.best.lm<-step(H1,scope=list(upper=H1,lower=H0),direction="both")
        summary(C.elev.best.lm)
        
        # Coefficients:
        #                       Estimate Std. Error t value Pr(>|t|)    
        # (Intercept)             549.80      11.32  48.551 5.12e-09 ***
        # Lamina.thck              24.69      27.49   0.898  0.40379    
        # LMA                    -268.09      70.39  -3.809  0.00888 ** 
        # LDMC                    393.06      93.44   4.207  0.00564 ** 
        # Leaf.Area                30.66      16.50   1.858  0.11248    
        # Lamina.thck:LMA        -139.57      53.46  -2.611  0.04008 *  
        # Lamina.thck:LDMC        516.52     102.64   5.032  0.00237 ** 
        # Lamina.thck:Leaf.Area   135.53      34.15   3.969  0.00738 ** 
        # LMA:LDMC                -44.63      17.86  -2.499  0.04659 *  
        # LMA:Leaf.Area          -159.35      46.54  -3.424  0.01407 *  
        # LDMC:Leaf.Area          104.86      27.45   3.820  0.00876 ** 
        # 
        # Residual standard error: 26.07 on 6 degrees of freedom
        # Multiple R-squared:  0.895,	Adjusted R-squared:   0.72 
        # F-statistic: 5.115 on 10 and 6 DF,  p-value: 0.02926
        
        # Starting with full model, everything is retained BUT model has non-significant terms...
        AIC(H1) # 165
        
        plot(C.elev.best.lm)
        # variance homoscedastic
        # normal residuals
        # PRPE & FRAM Cook's = 0.5
        
        a<-update(C.elev.best.lm,~.-Lamina.thck)
        a<-update(a,~.-LMA:Lamina.thck)
        a<-update(a,~.-Leaf.Area:Lamina.thck)
        a<-update(a,~.-LDMC:Lamina.thck )
        
        # D) Redo test without outliers
        dim(C.dat) #17 5
        C.dat.no.outliers<-C.dat[!rownames(C.dat)%in%c('PRPE','FRAM'),]
        dim(C.dat) #33 4
        
        H0<-lm(C.dat.no.outliers$avg.elev.pa~1,data=C.dat.no.outliers)
        H1<-lm(C.dat.no.outliers$avg.elev.pa~.^2,data=C.dat.no.outliers)
        
        C.elev.best.lm.no.outliers<-step(H1,scope=list(upper=H1,lower=H0),direction="both" )
        summary(C.elev.best.lm.no.outliers)
        AIC(C.elev.best.lm.no.outliers) # 148.7      
        
        C.elev.best.lm.no.outliers<-update(C.elev.best.lm,
                                           subset=(C.dat[!rownames(C.dat)%in%c('PRPE','FRAM'),]))
        
        # Redo model without LMA bc it is highly correlated wtih both Lamina.thck & LDMC
        
        H1<-lm(C.dat$avg.elev.pa~(LDMC*Lamina.thck*Leaf.Area),data=C.dat)
        H0<-lm(C.dat$avg.elev.pa~1,data=C.dat)
        
        C.elev.best.lm<-step(H1,scope=list(upper=H1,lower=H0),direction="both")
        summary(C.elev.best.lm)
        C.elev.best.lm<-update(C.elev.best.lm,~.-Lamina.thck:Leaf.Area)
        summary(C.elev.best.lm)
        C.elev.best.lm<-update(C.elev.best.lm,~.-Leaf.Area)
        summary(C.elev.best.lm)
        
        # Coefficients:
        #   Estimate Std. Error t value Pr(>|t|)    
        # (Intercept)        549.53      15.38  35.724  2.3e-14 ***
        # LDMC                59.87      41.47   1.444    0.172    
        # Lamina.thck        -28.82      23.76  -1.213    0.247    
        # LDMC:Lamina.thck   114.30      52.10   2.194    0.047 *  
        #   ---
        #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
        # 
        # Residual standard error: 42.68 on 13 degrees of freedom
        # Multiple R-squared:  0.3905,	Adjusted R-squared:  0.2498 
        # F-statistic: 2.776 on 3 and 13 DF,  p-value: 0.08334
        
        #$$$$
        #Non-significant model
        
        # B4.2 glm, gamma link fct - elev ~ traits - nothing significant ####
        
        # A) Check trait-trait interactions ####
        # ----
        
        summary(glm(C.dat$avg.elev.pa~.+(Lamina.thck+LMA+LDMC+Leaf.Area)^2,
                    data=C.dat,family=Gamma(link='log')),dispersion=1)
        
        # nothing significant
        
        # B) Check higher-order effects ####
        #----
        
        summary(glm(C.dat$avg.elev.pa~.+I(Lamina.thck^2)+I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2),
                    data=C.dat),family=Gamma(link='log'))  
        
        # nothing significant
        
        # C) Fit best model ####
        #----
        
        H1<-glm(C.dat$avg.elev.pa~.^2,data=C.dat,family= Gamma(link='log'))
        H0<-glm(C.dat$avg.elev.pa~1,data=C.dat,family=Gamma(link='log'))
        
        summary(C.elev.best.glm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=F),
                dispersion=1)
        
        # nothing significant
        
        summary(C.elev.best.glm<-step(H1,scope=list(lower=H0,upper=H1),direction='both',trace=F),
                dispersion=1)
        
        # nothing significant
        
        
        
# summary(lm) - reports t-values
# The t tests are the marginal impact of the variables in question, 
# given the presence of all the other variables. i.e. the beta coefficients are partial regression coefficients, 
# and we are testing their significance

# anova(lm) - reports F-values
# the reductions in the residual sum of squares as each term of the formula is added in turn
# anova(lm(formula)) = aov(formula). Order greatly affects results

  # The F tests are sequential - so they test for the importance of var1 in the presence of nothing 
  # but the intercept, of var2 in the presence of nothing but the intercept and V1, and of the interaction 
  # in the presence of all the above

# The t-test statistic (and its p-value) is a test of whether coefficient=0. The F-test on the anova() printout
# is whether the added variable significantly reduces the residual sum of squares.
        