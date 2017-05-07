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
      
      Hall.bins.prop.presence[Hall.bins.prop.presence$elev.range==elev.range[i],s]<-p # assign proportion of species occurence for each elevation bin 
      
    }
  }
  
  head(Hall.bins.prop.presence)[1:10]
  
# 1 - Calculate mean elevation of each species ####
#=================================================#
  
  # s- species s
  # i- elevation bin i
  # prop.pres(s.i)*elev(i)/sum(prop.pres(i))

# Calculate mean elevation for each species 

mean.elev<-data.frame(matrix(ncol=4,nrow=459))
names(mean.elev)<-c('Species','wghtd.elev','mid.elev','wghtd.elev')

# Calculate mean species elevation
for(s in sp.names) {
  mean.elev[which(sp.names==s),'Species']<-s                  # Assign species name
  
  # Mean elevation based on presence/absence data
  mean.elev[which(mean.elev$Species==s),'wghtd.elev']<-
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
}

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
mean.elev<-mean.elev[c('Species','Genus','epi','sp.code','avg.elev.pa','mid.elev','wghtd.elev')]
str(mean.elev)
  # data.frame':	459 obs. of  5 variables:
  # $ Species: chr  "Abies_balsamea" "Acer_pensylvanicum" "Acer_rubrum" "Acer_saccharum" ...
  # $ Genus  : chr  "Abies" "Acer" "Acer" "Acer" ...
  # $ epi    : chr  "balsamea" "pensylvanicum" "rubrum" "saccharum" ...
  # $ sp.code: chr  "ABBA" "ACPE" "ACRU" "ACSA" ...
  # $ avg.elev.pa   : num  625 582 523 572 622 ...
  # $ mid.elev: num  744 676 660 680 742 ...

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
mean.elev[mean.elev$sp.code%in%c('CASC', 'COCA', 'COTR', 'POPU'),] # what are the duplicate rows?

  #                   Species        Genus        epi sp.code avg.elev.pa
  # 84         Carex_scabrata        Carex   scabrata    CASC  599
  # 85         Carex_scoparia        Carex   scoparia    CASC  550
  # 112     Conyza_canadensis       Conyza canadensis    COCA  458
  # 113       Coptis_trifolia       Coptis   trifolia    COTR  592
  # 115  Corallorhiza_trifida Corallorhiza    trifida    COTR  530
  # 117     Cornus_canadensis       Cornus canadensis    COCA  669
  # 317 Polygonatum_pubescens  Polygonatum  pubescens    POPU  583
  # 330  Potamogeton_pusillus  Potamogeton   pusillus    POPU  490

# Keep Carex scabrata original sp.code
mean.elev[mean.elev$Species=='Carex_scoparia','sp.code']<-'CASCO' # rename Carex scoparia
# Keep Cornus canadensis original sp.code
mean.elev[mean.elev$Species=='Conyza_canadensis','sp.code']<-'COCAN' # rename Conyza canadensis
# Keep Polygonatum pubescens original sp.code
mean.elev[mean.elev$Species=='Potamogeton_pusillus','sp.code']<-'POPUS' # rename Potamogeton pusillus
# Keep Coptis trifolia original sp.code
mean.elev[mean.elev$Species=='Corallorhiza_trifida','sp.code']<-'COTRI'# rename Corallorhiza trifida

# Overwrite object with new changes
save(mean.elev,file=paste0(wrk.dir,'species.mean.elevation.based.on.Hall.plot.surveys.RData'))

# Recreate database
H.dat<-merge(sp.H.traits.c,mean.elev,by.y='sp.code',by.x=0) # '0' indicates rownames
str(H.dat)
H.dat[,c('Species','Genus','epi')]<-NULL

row.names(H.dat)<-H.dat$Row.names
H.dat$Row.names<-NULL

plot(density(H.dat$avg.elev.pa)) # left-skewed
shapiro.test((H.dat$avg.elev.pa)) # p-val = 0.005 ! Not bad!

summary(lm(H.dat$avg.elev.pa~H.dat$wghtd.elev)) #AdjR2 = 0.75
# These two elevations are 3/4 similar


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

# Keep Acer spicatum original sp.code
mean.elev[mean.elev$Species=='Actaea_sp','sp.code']<-'ACTSP' # rename Actaea_sp
# Keep Populus tremuloides original sp.code
mean.elev[mean.elev$Species=='Poa_trivialis','sp.code']<-'POTRI' # rename Poa trivialis

# Overwrite object with new changes
save(mean.elev,file=paste0(wrk.dir,'species.mean.elevation.based.on.Hall.plot.surveys.RData'))

# Recreate database
C.dat<-merge(sp.C.traits.c,mean.elev,by.y='sp.code',by.x=0) # '0' indicates rownames
str(C.dat)
C.dat[,c('Species','Genus','epi')]<-NULL

row.names(C.dat)<-C.dat$Row.names
C.dat$Row.names<-NULL

plot(density(C.dat$avg.elev.pa)) # bimodal - no suprise: conifers vs deciduous 
shapiro.test((C.dat$avg.elev.pa)) # p-val = 0.03 ! normal enough! Don't transform.
plot(density(C.dat$mid.elev)) # bimodal
shapiro.test((C.dat$mid.elev)) # p-val = 0.009 - Close enough

save(C.dat,file=paste0(wrk.dir,'Canopy.layer.dataframe.with.species.mean.elevation.and.traits.RData'))


# A - HERBACEOUS LAYER ####
#==========================#

  # A3 - Data Exploration ####
  #==========================#

        # A3.1 Scatterplots - pairwise interactions? ####
        #----

        names(H.dat) 
        # [1] "Ht.veg"        "Min.Root.Loca" "Lamina.thck"   "LMA"           "LDMC"          "Leaf.Area"     "myc.frac"     
        # [8] "avg.elev.pa"   "mid.elev"      "wghtd.elev"

        pairs(H.dat,panel=panel.smooth)
        # Variance seems homoscedastic
        # correlations among Leaf.Area-Ht.veg; Lamina.thck-LMA; LDMc-LMA
        # myc.frac and min.root.loca seem related to weighted elevation

        # A3.2 - Tree models - non-linearities ####
        #----
        
        elev.tree.model<-tree(H.dat$wghtd.elev~Ht.veg+Min.Root.Loca+Lamina.thck+LMA+LDMC+Leaf.Area+myc.frac,
                              data=H.dat)
        plot(elev.tree.model); text(elev.tree.model) ; title('herbaceous layer, avg.elev.pa vs traits')
        # Min.Root.Loca is most important variable, 
        # for those with high Min.Root.Loca, LDMC matters
        # For those with low min.root.loca, myc.frac matters
        elev.tree.model #  null deviance = 133000
        1-(deviance(elev.tree.model)/133000) # 0.53
        
        # A3.3 GAM - non-linearities?  ####
        #----
        
        par(mfrow=c(2,2))
        plot(gam(H.dat$wghtd.elev~s(Ht.veg)+s(Min.Root.Loca)
                 ,data=H.dat))
        # relationship w Min.Root.Loca, but no curvature
        
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
        summary(lm(H.dat$wghtd.elev~(Ht.veg+Min.Root.Loca+Lamina.thck+LMA+LDMC+Leaf.Area)^2,
                   data=H.dat))
        #LMA:Leaf.Area  marginally significant
        
        summary(lm(H.dat$wghtd.elev~(myc.frac+Min.Root.Loca+Lamina.thck+LMA+LDMC+Leaf.Area)^2,
                   data=H.dat))
        
        # nothing significant
        
        # B) Check non-linear relationships ####
        #----
        
        summary(lm(H.dat$wghtd.elev~I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
                     I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
                data=H.dat))
        # 
        # nothing significant
        
        # C) Fit best minimal model - not significant####
        #----
        
        H1<-lm(H.dat$wghtd.elev~Ht.veg+Min.Root.Loca+Lamina.thck+LMA+LDMC+LMA:Leaf.Area,data=H.dat)
        H0<-lm(H.dat$wghtd.elev~1,data=H.dat)
        
        elev.best.lm.forw<-step(H0,scope=list(lower=H0,upper=H1),trace=T,direction='both')
        summary(elev.best.lm.forw) # Min.Root.Loca retained, but marginally significant
        # AdjR2=0.08
        AIC(elev.best.lm.forw) #329.17
        
        # Null model just as good

        elev.best.lm.back<-step(H1,scope=list(lower=H0,upper=H1),trace=T,direction='both')
        summary(elev.best.lm.back) 

        # Coefficients:
        #               Estimate Std. Error t value Pr(>|t|)    
        # (Intercept)     728.87      30.51  23.887   <2e-16 ***
        # Min.Root.Loca   -21.55      11.01  -1.957   0.0620 .  
        # Lamina.thck      53.06      30.03   1.767   0.0899 .  
        # LMA             -89.57      46.49  -1.927   0.0659 .  
        # LDMC             76.73      34.78   2.206   0.0372 * 
        
        # Residual standard error: 63.68 on 24 degrees of freedom
        # Multiple R-squared:  0.2685,	Adjusted R-squared:  0.1466 
        # F-statistic: 2.202 on 4 and 24 DF,  p-value: 0.09914
        
        AIC(elev.best.lm.back) # 329.73
        
        # Continue model simplification
        
        m4<-update(elev.best.lm.back,~.-Lamina.thck); AIC(m4) #33128
        summary(m4)
        m5<-update(m4,~.-LMA); AIC(m4) #331
        summary(m5)
        m6<-update(m5,~.-LDMC); AIC(m6) # 329.17
        summary(m6)
        
        #----
        # back to minimal model from forward selection. 
        # Two models with similar AIC: one parameter explains 8% of variance, or 4 parameters explain 15%. 
        
        # Plot results 
        par(mfrow=c(1,1))
        
        plot(H.dat$wghtd.elev,H.dat$Ht.veg,
             pch=19, mgp=c(2,1,0),
             xlab='',ylab= 'Plant Height (cm)',family='serif')
         
        plot(H.dat$wghtd.elev,H.dat$LMA,
             pch=19, mgp=c(2,1,0),
             xlab='',ylab= 'LMA (mg/cm2)',family='serif')
        
        plot(H.dat$wghtd.elev,H.dat$Leaf.Area,
             pch=19, mgp=c(2,1,0),
             xlab='',ylab= 'Leaf Area (cm2)',family='serif')
        
        plot(H.dat$wghtd.elev,H.dat$Min.Root.Loca,
             pch=19, mgp=c(2,1,0),
             xlab='Species mean elevation (m)',ylab= 'Rooting depth',family='serif')
        
        plot(H.dat$wghtd.elev,H.dat$myc.frac,
             pch=19, mgp=c(2,1,0),
             xlab= 'Species mean elevation (m)', ylab= 'fraction of roots with \n mycorhizal associations',family='serif')
        
        # A4.2 glm, gamma link fct - 
        
        # A) Check trait-trait interactions ####
        # ----
        
        summary(glm(H.dat$wghtd.elev~(Min.Root.Loca+Ht.veg+Min.Root.Loca+Lamina.thck+LMA+
                                                  LDMC+Leaf.Area+myc.frac)^2,
                    data=H.dat,family=Gamma(link='log')),dispersion=1)
        
        # nothing significant
        
        # B) Check higher-order effects ####
        #----
          
        summary(glm(H.dat$wghtd.elev~I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
                      I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
                    data=H.dat),family=Gamma(link='log'))  
        
        # nothing significant
        
        # C) Fit best model ####
        #----
        
        H1<-glm(H.dat$wghtd.elev~.,data=H.dat,family= Gamma(link='log'))
        H0<-glm(H.dat$wghtd.elev~1,data=H.dat,family=Gamma(link='log'))
        
        summary(elev.best.glm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=T),
               dispersion=1)
        
        # Best model has non-significant terms
        
        # Coefficients:
        #              Estimate Std. Error z value Pr(>|z|)  
        # (Intercept)  5.199121   2.813546   1.848   0.0646 .
        # avg.elev.pa  0.001920   0.004117   0.466   0.6409  
        # myc.frac    -0.019395   0.172150  -0.113   0.9103  
        # Ht.veg      -0.025589   0.317374  -0.081   0.9357 
        
        AIC(elev.best.glm) # 280
        AIC(H0) # 330
        
        elev.best.glm #  null deviance = 0.29
        1-(deviance(elev.best.glm)/0.29) # 0.85
        
        #---
        # A model where none of the terms are significant explains 85% of the variance ?
        
      # B - CANOPY LAYER ####
      #=====================#
        
        # B3.1 Scatterplots - pairwise interactions? ####
        #----
        
        names(C.dat) 
        # [1] "Lamina.thck" "LMA"         "LDMC"        "Leaf.Area"   "avg.elev.pa"   "mid.elev"    "wghtd.elev"   
        
        pairs(C.dat,panel=panel.smooth)
        # Variance seems homoscedastic
        # yes, Lamina.thck-LMA; LDMc-LMA
        # maybe weighted elevation and Leaf Area?
        
        # B3.2 - Tree models - non-linearities ####
        #----
        
        C.elev.tree.model<-tree(C.dat$wghtd.elev~Lamina.thck+LMA+LDMC+Leaf.Area,data=C.dat)
        plot(C.elev.tree.model) ; text(C.elev.tree.model) ; title('canopy layer, elev vs traits')
        # LDMC is most important variable - non-linear relationships appear 
        # for those with high LDMC, high LDMC matters
        
        C.dat[C.dat$LDMC<=-2.0,] # Acer Spicatum. weird. 
        
        C.elev.tree.model #  null deviance = 143400
        1-(deviance(C.elev.tree.model)/143400) # 0.53
        
        # B3.3 GAM - non-linearities?  ####
        #----
        
        plot(gam(C.dat$wghtd.elev~s(Lamina.thck),data=C.dat))
        # curvy 
        
        plot(gam(H.dat$wghtd.elev~s(LDMC),data=H.dat))
        # no relationship     
        
        plot(gam(H.dat$wghtd.elev~s(LMA),data=H.dat))
        # no relationships
        
        plot(gam(H.dat$wghtd.elev~s(Leaf.Area),data=H.dat))
        # no relationship
        
        # B4 - Model selection ####
        #===========================================================
        
        #  B4.1 lm - elev ~ traits - nothing significant ! ####
        
        # A) Check trait-trait interactions ####
        #----
        
        # short on degrees of freedom. Do it in two models.
        summary(lm(C.dat$wghtd.elev~(Lamina.thck+LMA+LDMC+Leaf.Area)^2,
                   data=C.dat))
        
        # LMA, LDMC & almost all interactions significant!!!
        
        vif(lm(C.dat$wghtd.elev~(Lamina.thck+LMA+LDMC+Leaf.Area)^2,
               data=C.dat)) 
        
        # LMA, LDMC and Lamina.thck:LDMC are redundant, so remove LDMC from analyses
        
        # B) Check non-linear relationships ####
        #----
        
        summary(lm(C.dat$wghtd.elev~Lamina.thck+I(Lamina.thck^2)+LMA+I(LMA^2)+LDMC+I(LDMC^2)+Leaf.Area+I(Leaf.Area^2),
                   data=C.dat))
        
        # LDMC^2 and Lamina.thck^2 marginal
        
        # C) Fit  minimal model ####
        #----
        
        H1<-lm(C.dat$wghtd.elev~Lamina.thck+LMA+Leaf.Area+
                 I(Lamina.thck^2)+
                 +LMA:Leaf.Area+Lamina.thck:Leaf.Area+Lamina.thck:LMA,
               data=C.dat)
        H0<-lm(C.dat$wghtd.elev~1,data=C.dat)
        
        C.elev.best.lm.back<-step(H1,scope=list(upper=H1,lower=H0),direction="both")
        summary(C.elev.best.lm.back)
        
        #Coefficients:
        # Estimate Std. Error t value Pr(>|t|)    
        # (Intercept)       562.04      37.82  14.860 1.26e-08 ***
        #   Lamina.thck        16.29      68.95   0.236    0.818    
        # LMA               -10.23      72.49  -0.141    0.890    
        # Leaf.Area         -63.81      38.82  -1.644    0.128    
        # LMA:Leaf.Area    -130.40      88.11  -1.480    0.167    
        # Lamina.thck:LMA    87.67      73.88   1.187    0.260    
        # ---
        #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
        # 
        # Residual standard error: 92.8 on 11 degrees of freedom
        # Multiple R-squared:  0.3393,	Adjusted R-squared:  0.03901 
        # F-statistic:  1.13 on 5 and 11 DF,  p-value: 0.4002
        AIC(C.elev.best.lm.back) # 208
        AIC(H0) # 205 
        # Why didn't it keep the null model?
        
        C.elev.best.lm.forw<-step(H0,scope=list(upper=H1,lower=H0),direction="both")
        summary(C.elev.best.lm.forw)
        
        # Null model retained
        
        
        # B4.2 glm, gamma link fct - elev ~ traits - nothing significant ####
        
        # A) Check trait-trait interactions ####
        # ----
        
        summary(glm(C.dat$wghtd.elev~(Lamina.thck+LMA+LDMC+Leaf.Area)^2,
                    data=C.dat,family=Gamma(link='log')),dispersion=1)
        
        # nothing significant
        
        # B) Check higher-order effects ####
        #----
        
        summary(glm(C.dat$wghtd.elev~+I(Lamina.thck^2)+I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2),
                    data=C.dat),family=Gamma(link='log'))  
        
        # nothing significant
        
        # C) Fit best model ####
        #----
        
        H1<-glm(C.dat$wghtd.elev~(Lamina.thck+LMA+LDMC+Leaf.Area)^2,data=C.dat,family= Gamma(link='log'))
        H0<-glm(C.dat$wghtd.elev~1,data=C.dat,family=Gamma(link='log'))
        
        summary(C.elev.best.glm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=F),
                dispersion=1)
        
        # Null model retained
        
        summary(C.elev.best.glm<-step(H1,scope=list(lower=H0,upper=H1),direction='both',trace=F),
                dispersion=1)
        
        # everything retained, but model not significant. 
        # Redo of analyses AQUI
        
        
        
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
        