#<<TABLE OF CONTENTS>>
# 1 - Calculate mean elevation of each species
# 2- Create dataframes
# A - HERBIVORY LAYER
# A3- Data Exploration - Scatterplots, Tree models and GAMs
# A4 - Model Selection
#   A4.1 - lm - elev ~ traits. Traits do not predict herbaceous species elevation
#   A4.2 - glm - elev ~ traits. Traits do not predict herbaceous species elevation
# B - CANOPY LAYER
# B3 - Data Exploration - Scatterplots, Tree models and GAMs

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

# 1 - Calculate mean elevation of each species
# Hall = presence/absence of species along small plots spread along elevational gradient
Hall<-read.csv(paste0(data.dir,'G.Hall_sp_by_site_corrected.csv'),header=T)

# Calculate mean elevation for each species 
sp.names<-names(Hall[,12:length(names(Hall))])

mean.elev<-data.frame(matrix(ncol=2,nrow=459))
names(mean.elev)<-c('Species','elev')

# Calculate mean species elevation
for(x in sp.names) {
  mean.elev[which(sp.names==x),'Species']<-x
  mean.elev[which(mean.elev$Species==x),'elev']<-round(sum(Hall$Elev_m * Hall[,x]) / sum(Hall[,x]),digits=0) 
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
mean.elev<-mean.elev[c('Species','Genus','epi','sp.code','elev')]
str(mean.elev)
  # data.frame':	459 obs. of  5 variables:
  # $ Species: chr  "Abies_balsamea" "Acer_pensylvanicum" "Acer_rubrum" "Acer_saccharum" ...
  # $ Genus  : chr  "Abies" "Acer" "Acer" "Acer" ...
  # $ epi    : chr  "balsamea" "pensylvanicum" "rubrum" "saccharum" ...
  # $ sp.code: chr  "ABBA" "ACPE" "ACRU" "ACSA" ...
  # $ elev   : num  625 582 523 572 622 ...

# 2- Create dataframes

# Herbaceous dataframe - merge trait data and elevation based on species code
H.dat<-merge(sp.H.traits.c,mean.elev,by.y='sp.code',by.x=0) # '0' indicates rownames
str(H.dat)
H.dat[,c('Species','Genus','epi')]<-NULL
row.names(H.dat)<-H.dat$Row.names
# Error - non-unique values of ‘CASC’, ‘COCA’, ‘COTR’, ‘POPU’ 


# clean Duplicates
mean.elev[mean.elev$sp.code%in%c('CASC', 'COCA', 'COTR', 'POPU'),]
  #                   Species        Genus        epi sp.code elev
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

# Recreate database
H.dat<-merge(sp.H.traits.c,mean.elev,by.y='sp.code',by.x=0) # '0' indicates rownames
str(H.dat)
H.dat[,c('Species','Genus','epi')]<-NULL

row.names(H.dat)<-H.dat$Row.names
H.dat$Row.names<-NULL

plot(density(H.dat$elev)) # looks quite normal except for bums in tails and a little right-skewed
shapiro.test((H.dat$elev)) # p-val = 0.053 ! Woo-hoo!
shapiro.test(log(H.dat$elev)) #p-val=0.11 
plot(density(log(H.dat$elev))) # not any different. don't transform.


# Canopy dataframe - merge trait data and elevation based on species code
C.dat<-merge(sp.C.traits.c,mean.elev,by.y='sp.code',by.x=0) # '0' indicates rownames
str(C.dat)
C.dat[,c('Species','Genus','epi')]<-NULL

row.names(C.dat)<-C.dat$Row.names
# Error - non-unique values when setting 'row.names': ‘ACSP’, ‘POTR’

# Clean Duplicates
mean.elev[mean.elev$sp.code%in%c('ACSP', 'POTR'),]
  #                 Species   Genus         epi sp.code elev
  # 5         Acer_spicatum    Acer    spicatum    ACSP  622
  # 7             Actaea_sp  Actaea          sp    ACSP  593
  # 316       Poa_trivialis     Poa   trivialis    POTR  549
  # 328 Populus_tremuloides Populus tremuloides    POTR  528

# Keep Acer spicatum
mean.elev[mean.elev$Species=='Actaea_sp','sp.code']<-'ACT.SP'
# Keep Populus tremuloides
mean.elev[mean.elev$Species=='Poa_trivialis','sp.code']<-'POTRI'

# Recreate database
C.dat<-merge(sp.C.traits.c,mean.elev,by.y='sp.code',by.x=0) # '0' indicates rownames
str(C.dat)
C.dat[,c('Species','Genus','epi')]<-NULL

row.names(C.dat)<-C.dat$Row.names
C.dat$Row.names<-NULL

plot(density(C.dat$elev)) # looks quite normal except for a little right-skewed
shapiro.test((C.dat$elev)) # p-val = 0.21 ! Woo-hoo! Don't transform.


# A - HERBIVORY LAYER ####
#==========================#

  # A3 - Data Exploration ####
  #==========================#

        # 3.1 Scatterplots - pairwise interactions? ####
        #----

        names(H.dat) 
        # [1] "Row.names"     "Ht.veg"        "Min.Root.Loca" "Lamina.thck"   "LMA"           "LDMC"          "Leaf.Area"    
        # [8] "myc.frac"      "elev" 

        pairs(H.dat,panel=panel.smooth)
        # Variance seems homoscedastic
        # yes, Leaf.Area-Ht.veg; Lamina.thck-LMA; LDMc-LMA

        # 3.2 - Tree models - non-linearities ####
        #----
        
        elev.tree.model<-tree(H.dat$elev~.,data=H.dat)
        plot(elev.tree.model)
        text(elev.tree.model)
        title('herbaceous layer, elev vs traits')
        # LMA is most important variable, 
        # for those with high LMA, myc.fra matters
        # For those with with LMA and low myc.fra, ht.veg matters
        elev.tree.model #  null deviance = 33790.0
        1-(deviance(elev.tree.model)/33790.0) # 0.43
        
        # 3.3 GAM - non-linearities?  ####
        #----
        
        plot(gam(H.dat$elev~s(Ht.veg)+s(Min.Root.Loca)
                 ,data=H.dat))
        # no relationship
        
        plot(gam(H.dat$elev~s(Lamina.thck)+s(LMA)
                 ,data=H.dat))
        # no relationship
        
        plot(gam(H.dat$elev~s(LDMC)+s(Leaf.Area),
                 data=H.dat))
        # no relationship     
        
        plot(gam(H.dat$elev~s(myc.frac),
                 data=H.dat))
        # negative curved relationship appears      
        
    #=====================================================
    
    # Diagnostic plots
    
    # 1) Residuals vs Fitted - detects residual non-linear relationships between x & y variables 
    #    and heteroscedasticity (wedge)
    # 2) QQ plot - shows whether residuals are normally distributed (will follow straight line)
    # 3) Scale-Location - shows whether variance (residuals) increase with incrasing mean (also tests
    # for homoscedasticity). Another version of plot 1
    # 4) Residuals vs Leverage - shows whether individual datapoints have a lot of 'weight' on 
    # the regression
    
    # 4 - Model selection ####
    #===========================================================
        
        #  A4.1 lm - elev ~ traits - nothing significant ! ####
        
        # A) Check trait-trait interactions ####
        #----
        
        # short on degrees of freedom. Do it in two models.
                summary(lm(H.dat$elev~(Ht.veg+Min.Root.Loca+Lamina.thck+LMA+LDMC+Leaf.Area)^2,
                   data=H.dat))
        # LMA-Leaf.Area & LDMc-Leaf.Area marginally significant
        
        summary(lm(H.dat$elev~(myc.frac+Min.Root.Loca+Lamina.thck+LMA+LDMC+Leaf.Area)^2,
                   data=H.dat))
        # LMA-Leaf.Area significant
        # LDMC-Leaf.Area, LMA-LDMC marginal
        
        
        summary(lm(H.dat$elev~(myc.frac+Ht.veg)^2,
                   data=H.dat))
        # myc.frac:Ht.veg significant
        
        # B) Check non-linear relationships ####
        #----
        
        summary(lm(H.dat$elev~.+I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
                     I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
                data=H.dat))
        # nothing significant
        
        # C) Fit best minimal model ####
        #----
        
        H1<-lm(H.dat$elev~.+LMA:Leaf.Area+LDMC:Leaf.Area+LMA:LDMC+myc.frac:Ht.veg,data=H.dat)
        H0<-lm(H.dat$elev~1,data=H.dat)
        
        elev.best.lm<-step(H0,scope=list(lower=H0,upper=H1),trace=T,direction='both')
        summary(elev.best.lm)
        
        # nothing significant! 
        
        # A4.2 glm, gamma link fct - 
        
        # A) Check trait-trait interactions ####
        # ----
        
        summary(glm(H.dat$elev~.+(Min.Root.Loca+Ht.veg+Min.Root.Loca+Lamina.thck+LMA+
                                                  LDMC+Leaf.Area+myc.frac)^2,
                    data=H.dat,family=Gamma(link='log')),dispersion=1)
        
        # nothing significant
        
        # B) Check higher-order effects ####
        #----
          
        summary(glm(H.dat$elev~I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
                      I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
                    data=H.dat),family=Gamma(link='log'))  
        
        # nothing significant
        
        # C) Fit best model ####
        #----
        
        H1<-glm(H.dat$elev~.,data=H.dat,family= Gamma(link='log'))
        H0<-glm(H.dat$elev~1,data=H.dat,family=Gamma(link='log'))
        
        summary(elev.best.glm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=F),
               dispersion=1)
        
        # nothing significant
        
      # CANOPY LAYER
        
        # 4.1 Scatterplots - pairwise interactions? ####
        #----
        
        names(C.dat) 
        # [1] "Lamina.thck" "LMA"         "LDMC"        "Leaf.Area"   "elev"   
        
        pairs(C.dat,panel=panel.smooth)
        # Variance seems homoscedastic
        # yes, Lamina.thck-LMA; LDMc-LMA
        
        # 4.2 - Tree models - non-linearities ####
        #----
        
        C.elev.tree.model<-tree(C.dat$elev~.,data=C.dat)
        plot(C.elev.tree.model)
        text(C.elev.tree.model)
        title('canopy layer, elev vs traits')
        # LDMC is most important variable, 
        # for those with high LDMC, high LDMC matters
        C.elev.tree.model #  null deviance = 38840
        1-(deviance(C.elev.tree.model)/38840) # 0.55
        
        # 4.3 GAM - non-linearities?  ####
        #----
        
        plot(gam(C.dat$elev~s(Lamina.thck),data=C.dat))
        # curvy 
        
        plot(gam(H.dat$elev~s(LDMC),data=H.dat))
        # no relationship     
        
        plot(gam(H.dat$elev~s(LMA),data=H.dat))
        # negative ?
        
        plot(gam(H.dat$elev~s(Leaf.Area),data=H.dat))
        # positive?