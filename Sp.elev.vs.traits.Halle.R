#<<TABLE OF CONTENTS>>
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

# 1 - Calculate mean elevation of each species
# Hall = presence/absence of species along small plots spread along elevational gradient
Hall<-read.csv(paste0(data.dir,'G.Hall_sp_by_site_corrected.csv'),header=T)

# Calculate mean elevation for each species 
sp.names<-names(Hall[,12:length(names(Hall))])

mean.elev<-data.frame(matrix(ncol=3,nrow=459))
names(mean.elev)<-c('Species','avg.elev','mid.elev')

# Calculate mean species elevation
for(s in sp.names) {
  mean.elev[which(sp.names==s),'Species']<-s
  mean.elev[which(mean.elev$Species==s),'avg.elev']<-round(sum(Hall$Elev_m * Hall[,s]) / sum(Hall[,s]),digits=0) 
  max.e<-max(Hall[Hall[,s]==1,'Elev_m'])
  min.e<-min(Hall[Hall[,s]==1,'Elev_m'])
  mean.elev[which(mean.elev$Species==s),'mid.elev']<-min.e+((max.e-min.e)/2)
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
mean.elev<-mean.elev[c('Species','Genus','epi','sp.code','avg.elev','mid.elev')]
str(mean.elev)
  # data.frame':	459 obs. of  5 variables:
  # $ Species: chr  "Abies_balsamea" "Acer_pensylvanicum" "Acer_rubrum" "Acer_saccharum" ...
  # $ Genus  : chr  "Abies" "Acer" "Acer" "Acer" ...
  # $ epi    : chr  "balsamea" "pensylvanicum" "rubrum" "saccharum" ...
  # $ sp.code: chr  "ABBA" "ACPE" "ACRU" "ACSA" ...
  # $ avg.elev   : num  625 582 523 572 622 ...
  # $ mid.elev: num  744 676 660 680 742 ...

save(mean.elev,file=paste0(wrk.dir,'species.mean.elevation.based.on.Hall.plot.surveys.RData'))

# 2- Create dataframes

# Herbaceous dataframe - merge trait data and elevation based on species code
#----

H.dat<-merge(sp.H.traits.c,mean.elev,by.y='sp.code',by.x=0) # '0' indicates rownames
str(H.dat)
H.dat[,c('Species','Genus','epi')]<-NULL
row.names(H.dat)<-H.dat$Row.names
# Error - non-unique values of ‘CASC’, ‘COCA’, ‘COTR’, ‘POPU’ 

# clean Duplicates
mean.elev[mean.elev$sp.code%in%c('CASC', 'COCA', 'COTR', 'POPU'),]
  #                   Species        Genus        epi sp.code avg.elev
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

plot(density(H.dat$avg.elev)) # looks quite normal except for bums in tails and a little right-skewed
shapiro.test((H.dat$avg.elev)) # p-val = 0.053 ! Woo-hoo!
plot(density(H.dat$mid.avg.elev))
plot(H.dat$mid.elev,H.dat$avg.elev,
     xlim=c(500,800),ylim=c(500,800))
points(x=H.dat$mid.elev,
       y=predict(lm(H.dat$avg.elev~H.dat$mid.elev)),
       col='red', pch=20) 
abline(a=0,b=1)

summary(lm(H.dat$avg.elev~H.dat$mid.elev)) #AdjR2 = 0.30
# These two elevations are quite different!


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
  #                 Species   Genus         epi sp.code avg.elev
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

plot(density(C.dat$avg.elev)) # looks quite normal except for a little right-skewed
shapiro.test((C.dat$avg.elev)) # p-val = 0.21 ! Woo-hoo! Don't transform.
plot(density(C.dat$mid.elev)) # looks quite normal except for a little right-skewed
shapiro.test((C.dat$avg.elev)) # p-val = 0.21 ! Woo-hoo! Don't transform.

save(C.dat,file=paste0(wrk.dir,'Canopy.layer.dataframe.with.species.mean.elevation.and.traits.RData'))


# A - HERBIVORY LAYER ####
#==========================#

  # A3 - Data Exploration ####
  #==========================#

        # A3.1 Scatterplots - pairwise interactions? ####
        #----

        names(H.dat) 
        # [1] "Row.names"     "Ht.veg"        "Min.Root.Loca" "Lamina.thck"   "LMA"           "LDMC"          "Leaf.Area"    
        # [8] "myc.frac"      "avg.elev"      "mid.elev"

        pairs(H.dat,panel=panel.smooth)
        # Variance seems homoscedastic
        # correlations among Leaf.Area-Ht.veg; Lamina.thck-LMA; LDMc-LMA

        # A3.2 - Tree models - non-linearities ####
        #----
        
        elev.tree.model<-tree(H.dat$avg.elev~Ht.veg+Min.Root.Loca+Lamina.thck+LMA+LDMC+Leaf.Area+myc.frac,
                              data=H.dat)
        plot(elev.tree.model); text(elev.tree.model) ; title('herbaceous layer, avg.elev vs traits')
        # LMA is most important variable, 
        # for those with high LMA, myc.fra matters
        # For those with with LMA and low myc.fra, ht.veg matters
        elev.tree.model #  null deviance = 33790.0
        1-(deviance(elev.tree.model)/33790.0) # 0.43
        
        mid.elev.tree.model<-tree(H.dat$mid.elev~Ht.veg+Min.Root.Loca+Lamina.thck+LMA+LDMC+Leaf.Area+myc.frac,
                                  data=H.dat)
        plot(mid.elev.tree.model); text(mid.elev.tree.model) ; title('herbaceous layer, mid.elev vs traits')
        
        # mid.elevation doesn't work as well as average
        
        # A3.3 GAM - non-linearities?  ####
        #----
        
        plot(gam(H.dat$avg.elev~s(Ht.veg)+s(Min.Root.Loca)
                 ,data=H.dat))
        # no relationship
        
        plot(gam(H.dat$avg.elev~s(Lamina.thck)+s(LMA)
                 ,data=H.dat))
        # no relationship
        
        plot(gam(H.dat$avg.elev~s(LDMC)+s(Leaf.Area),
                 data=H.dat))
        # no relationship     
        
        plot(gam(H.dat$avg.elev~s(myc.frac),
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
    
    # A4 - Model selection ####
    #===========================================================
        
        #  A4.1 lm - elev ~ traits - nothing significant ! ####
        
        # A) Check trait-trait interactions ####
        #----
        
        # short on degrees of freedom. Do it in two models.
        summary(lm(H.dat$avg.elev~(Ht.veg+Min.Root.Loca+Lamina.thck+LMA+LDMC+Leaf.Area)^2,
                   data=H.dat))
        # LMA-Leaf.Area & LDMc-Leaf.Area marginally significant
        
        summary(lm(H.dat$avg.elev~(myc.frac+Min.Root.Loca+Lamina.thck+LMA+LDMC+Leaf.Area)^2,
                   data=H.dat))
        # LMA-Leaf.Area significant
        # LDMC-Leaf.Area, LMA-LDMC marginal
        
        
        summary(lm(H.dat$avg.elev~(myc.frac+Ht.veg)^2,
                   data=H.dat))
        # myc.frac:Ht.veg significant
        
        # B) Check non-linear relationships ####
        #----
        
        summary(lm(H.dat$avg.elev~.+I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
                     I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
                data=H.dat))
        # nothing significant
        
        # C) Fit best minimal model ####
        #----
        
        H1<-lm(H.dat$avg.elev~.+LMA:Leaf.Area+LDMC:Leaf.Area+LMA:LDMC+myc.frac:Ht.veg,data=H.dat)
        H0<-lm(H.dat$avg.elev~1,data=H.dat)
        
        elev.best.lm<-step(H0,scope=list(lower=H0,upper=H1),trace=T,direction='both')
        summary(elev.best.lm)
        
        # nothing significant! 
        
        # A4.2 glm, gamma link fct - 
        
        # A) Check trait-trait interactions ####
        # ----
        
        summary(glm(H.dat$avg.elev~.+(Min.Root.Loca+Ht.veg+Min.Root.Loca+Lamina.thck+LMA+
                                                  LDMC+Leaf.Area+myc.frac)^2,
                    data=H.dat,family=Gamma(link='log')),dispersion=1)
        
        # nothing significant
        
        # B) Check higher-order effects ####
        #----
          
        summary(glm(H.dat$avg.elev~I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
                      I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
                    data=H.dat),family=Gamma(link='log'))  
        
        # nothing significant
        
        # C) Fit best model ####
        #----
        
        H1<-glm(H.dat$avg.elev~.,data=H.dat,family= Gamma(link='log'))
        H0<-glm(H.dat$avg.elev~1,data=H.dat,family=Gamma(link='log'))
        
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
        
        C.elev.tree.model<-tree(C.dat$avg.elev~.,data=C.dat)
        plot(C.elev.tree.model)
        text(C.elev.tree.model)
        title('canopy layer, elev vs traits')
        # LDMC is most important variable, 
        # for those with high LDMC, high LDMC matters
        C.elev.tree.model #  null deviance = 38840
        1-(deviance(C.elev.tree.model)/38840) # 0.55
        
        # B3.3 GAM - non-linearities?  ####
        #----
        
        plot(gam(C.dat$avg.elev~s(Lamina.thck),data=C.dat))
        # curvy 
        
        plot(gam(H.dat$avg.elev~s(LDMC),data=H.dat))
        # no relationship     
        
        plot(gam(H.dat$avg.elev~s(LMA),data=H.dat))
        # negative ?
        
        plot(gam(H.dat$avg.elev~s(Leaf.Area),data=H.dat))
        # positive?
        
        # B4 - Model selection ####
        #===========================================================
        
        #  B4.1 lm - elev ~ traits - nothing significant ! ####
        
        # A) Check trait-trait interactions ####
        #----
        
        # short on degrees of freedom. Do it in two models.
        summary(lm(C.dat$avg.elev~(Lamina.thck+LMA+LDMC+Leaf.Area)^2,
                   data=C.dat))
        
        # LMA, LDMC & all interactions significant!!!
        
        vif(lm(C.dat$avg.elev~(Lamina.thck+LMA+LDMC+Leaf.Area)^2,
               data=C.dat)) 
        
        # LMA, LDMC and Lamina.thck:LDMC are redundant
        
        # B) Check non-linear relationships ####
        #----
        
        summary(lm(C.dat$avg.elev~Lamina.thck+I(Lamina.thck^2)+LMA+I(LMA^2)+LDMC+I(LDMC^2)+Leaf.Area+I(Leaf.Area^2),
                   data=C.dat))
        
        # nothing significant
        
        # C) Fit  minimal model ####
        #----
        
        H1<-lm(C.dat$avg.elev~.^2,data=C.dat)
        H0<-lm(C.dat$avg.elev~1,data=C.dat)
        
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
        
        H0<-lm(C.dat.no.outliers$avg.elev~1,data=C.dat.no.outliers)
        H1<-lm(C.dat.no.outliers$avg.elev~.^2,data=C.dat.no.outliers)
        
        C.elev.best.lm.no.outliers<-step(H1,scope=list(upper=H1,lower=H0),direction="both" )
        summary(C.elev.best.lm.no.outliers)
        AIC(C.elev.best.lm.no.outliers) # 148.7      
        
        C.elev.best.lm.no.outliers<-update(C.elev.best.lm,
                                           subset=(C.dat[!rownames(C.dat)%in%c('PRPE','FRAM'),]))
        
        # Redo model without LMA bc it is highly correlated wtih both Lamina.thck & LDMC
        
        H1<-lm(C.dat$avg.elev~(LDMC*Lamina.thck*Leaf.Area),data=C.dat)
        H0<-lm(C.dat$avg.elev~1,data=C.dat)
        
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
        
        summary(glm(C.dat$avg.elev~.+(Lamina.thck+LMA+LDMC+Leaf.Area)^2,
                    data=C.dat,family=Gamma(link='log')),dispersion=1)
        
        # nothing significant
        
        # B) Check higher-order effects ####
        #----
        
        summary(glm(C.dat$avg.elev~.+I(Lamina.thck^2)+I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2),
                    data=C.dat),family=Gamma(link='log'))  
        
        # nothing significant
        
        # C) Fit best model ####
        #----
        
        H1<-glm(C.dat$avg.elev~.^2,data=C.dat,family= Gamma(link='log'))
        H0<-glm(C.dat$avg.elev~1,data=C.dat,family=Gamma(link='log'))
        
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
        