#<<TABLE OF CONTENTS>>
# 0-  Make trait and abundance dataframes correspond based on shared species
# 1- Data Exploration 
# 2- Abund~Traits, Individual variables - look at diagnostic plots
# 3- log-transform response, lm

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

load(file=paste0(wrk.dir,'list.trait.names.herbivory.layer.Rdata'))
load(file=paste0(wrk.dir,"Herbaceous.Layer.Traits.Standardized.RData"))
load(file=paste0(wrk.dir,'species-level.traits.Herbaceous.layer.Rdata'))

load(file=paste0(wrk.dir,'list.trait.names.canopy.layer.Rdata'))
load(file=paste0(wrk.dir,"Canopy.Layer.Traits.Standardized.RData"))
load(file=paste0(wrk.dir,'species-level.traits.Canopy.layer.Rdata'))

# load(file=paste0(wrk.dir,'Species.abundances.full.data.Rdata')) # object = abund
load(file=paste0(wrk.dir,'Species.abundances.Inf.removed.Rdata')) # object = abund.c
#============================================================================================

# A - HERBIVORY LAYER, ABUNDANCE

# 0-  Make trait and abundance dataframes correspond based on shared species ####
#================================================================================

  # 0A - Herbivory layer
    names(H.abund)
    H.abund<-abund.c[which(rownames(abund.c)%in%rownames(sp.H.traits)),c('abund.ratio','log.abund.ratio',
                                                                         'ratio.log.abund',"PT.abund.ratio" )]
    H.abund<-H.abund[order(rownames(H.abund)),]
    View(H.abund)
    dim(H.abund) #45 4
    
    dim(sp.H.traits)
    # 53 8
    sp.H.traits<-sp.H.traits[which(rownames(sp.H.traits)%in%rownames(H.abund)),]
    dim(sp.H.traits) 
    #45 8
    all(rownames(H.abund)==rownames(sp.H.traits)) # TRUE

# 1 - Data Exploration - Scatterplots, Tree models and GAMs #### 
#===========================================================

    # Scatterplots - pairwise interactions?
    
    names(H.abund) # "abund.ratio"     "log.abund.ratio" "ratio.log.abund" "pt.abund.ratio" 
    pairs(cbind(sp.H.traits,H.abund[,'abund.ratio']),panel=panel.smooth)
    #yes, LMA-LDMC ; Root.Loca-myc.frac ; min.root.Loca - max.root.loca
    # variance appears homoscedastic
    
    pairs(cbind(sp.H.traits,H.abund[,'log.abund.ratio']),panel=panel.smooth)
    # variance appears homoscedastic except for one outlier
    
    pairs(cbind(sp.H.traits,H.abund[,'ratio.log.abund']),panel=panel.smooth)
    # can't tell variance homogeneity - one outlier too large
    
    pairs(cbind(sp.H.traits,H.abund[,'pt.abund.ratio']),panel=panel.smooth)
    # variance appears homoscedastic
    
    # Tree models - non-linearities
  
    tree.model<-tree(H.abund$abund.ratio~.,data=sp.H.traits)
    plot(tree.model)
    text(tree.model)
    title('herbaceous layer, abund.ratio vs traits')
    # Myc.frac is most important variable, 
    # for those with high myc.frac, LDMC matters
    tree.model #  null deviance = 167.00
    1-(deviance(tree.model)/167) # 0.29 
    
    tree.model.l<-tree(H.abund$log.abund.ratio~.,data=sp.H.traits)
    plot(tree.model.l)
    text(tree.model.l)
    title('herbaceous layer, log.abund.ratio vs traits')
    # complicated - LDMc, myc.frac, leaf.area, LMA...
    tree.model.l # null deviance = 45.4800
    1-(deviance(tree.model.l)/45.48) # 0.43
    
    tree.model.rl<-tree(H.abund$ratio.log.abund~.,data=sp.H.traits)
    plot(tree.model.rl)
    text(tree.model.rl)
    title('herbaceous layer, ratio.log.abund vs traits')
    # only LDMc retained
    tree.model.rl # null deviance = 592.50
    1-(deviance(tree.model.rl)/592.50) # 0.20
    
    tree.model.pt<-tree(H.abund$PT.abund.ratio~.,data=sp.H.traits)
    plot(tree.model.pt)
    text(tree.model.pt)
    title('herbaceous layer, PT.ratio.abund vs traits')
    # LDMC & myc.fraction retained
    tree.model.pt # null deviance = 1.281
    1-(deviance(tree.model.pt)/1.281) # 0.45
    
    # GAM - non-linearities?  
    
    par(mfrow=c(2,2))
    plot(gam(H.abund$abund.ratio~s(Ht.veg)+s(myc.frac)+
               s(Min.Root.Loca),data=sp.H.traits))
    
    #myc.frac appears non-linear
    plot(gam(H.abund$abund.ratio~s(Ht.veg)+s(Lamina.thck)+
             s(LMA),data=sp.H.traits))
    
    # No relationships
    plot(gam(H.abund$abund.ratio~s(Leaf.Area),data=sp.H.traits))
    # No relationships
    
    par(mfrow=c(1,1))
    
    #------------------------------------#
    # myc.frac and maybe LDMC appear non-linear
        
#=====================================================
     
  # Diagnostic plots
      
  # 1) Residuals vs Fitted - detects residual non-linear relationships between x & y variables 
  #    and heteroscedasticity (wedge)
  # 2) QQ plot - shows whether residuals are normally distributed (will follow straight line)
  # 3) Scale-Location - shows whether variance (residuals) increase with incrasing mean (also tests
      # for homoscedasticity). Another version of plot 1
  # 4) Residuals vs Leverage - shows whether individual datapoints have a lot of 'weight' on 
      # the regression
               
# 2 - Model selection ####
#===========================================================
 
    # Data Prep
  
    names(sp.H.traits)
    #[1] "Ht.veg"        "Min.Root.Loca" "Max.Root.Loca" "Lamina.thck"   "LMA"           "LDMC"          "Leaf.Area"    
    #[8] "myc.frac"
    sp.H.traits$Max.Root.Loca<-NULL
    dim(sp.H.traits)
    # 45 7
    
    # Remove species with NAs in myc.frac or Leaf.Area, to be able to run stepwise selection
    sp.H.traits.c<-sp.H.traits[!is.na(sp.H.traits$myc.frac),]
    sp.H.traits.c<-sp.H.traits.c[!is.na(sp.H.traits.c$Leaf.Area),]
    dim(sp.H.traits.c)
    # 35 7
    
    # Adjust species list for abundance data to have matching datasets
    H.abund.c<-H.abund[which(rownames(H.abund)%in%rownames(sp.H.traits.c)),]
    dim(H.abund.c)
    #35 4
    
  # 2.1 lm - abundance Ratio - doesn't meet assumptions & only significant bc of outliers ####
  #========================================================================#

  # A) Check trait-trait interactions ####
    ----
   summary(lm(H.abund.c$abund.ratio~(Min.Root.Loca+Ht.veg+Min.Root.Loca+Lamina.thck+LMA+
                                        LDMC+Leaf.Area+myc.frac)^2,
              data=sp.H.traits.c))
      
      # $$$$
      # No significant interactions
    
  
  # B) Check non-linear relationships ####
    ----
    summary(lm(H.abund.c$abund.ratio~I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
           I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
         data=sp.H.traits.c))
    
      # $$$$
      # myc.frac^2 significant

  # C) Fit best minimal model ####
    ----
    
  H1<-lm(H.abund.c$abund.ratio~.+I(myc.frac^2),data=sp.H.traits.c)
  H0<-lm(H.abund.c$abund.ratio~1,data=sp.H.traits.c)
  
  best.lm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace = F)
  
  summary(best.lm)
  AIC(best.lm) # 141.8
  best.lm<-update(best.lm,~.-Min.Root.Loca)
  summary(best.lm)
  AIC(best.lm) # 141.9
  best.lm<-update(best.lm,~.-myc.frac)
  AIC(best.lm) # 142.3
  summary(best.lm)
  # Coefficients:
  #               Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)     0.6054     0.3619   1.673    0.104    
  # I(myc.frac^2)   0.9402     0.2021   4.651 5.13e-05 ***
  #   ---
  # 
  # Residual standard error: 1.748 on 33 degrees of freedom
  # Multiple R-squared:  0.396,	Adjusted R-squared:  0.3777 
  # F-statistic: 21.64 on 1 and 33 DF,  p-value: 5.134e-05
  
  # D) Diagnostic plots ####
  ----
  
  plot(best.lm)
  # Variance in residuals decreases with mean - homoscedastic, but an x-y relationship remain unexplained by model
  # QQ-plot -  error not normal - J-shaped - bad
  # Leverage - CASC & LYAN have a lot of influence on model - try model without that datapoint.
  
  # E) Model without outliers ####
  ----
  dim(H.abund.c) #35 4
  H.abund.c.no.outliers<-H.abund.c[!rownames(H.abund.c)%in%c('CASC','LYAN'),]
  dim(H.abund.c.no.outliers) #33 4
  
  sp.H.traits.c.no.outliers<-sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('CASC','LYAN'),]
  dim(sp.H.traits.c.no.outliers) #33 8
  
  RsquareAdj(lm(
    H.abund.c.no.outliers$abund.ratio~I(myc.frac^2)+myc.frac,data=sp.H.traits.c.no.outliers
                          ))
  # -0.06
  
    #$$$$
    # model was significant only because of CASC and LYAN
  
  
  # 2.2 lm - log abundance ratio - MEETS ASSUMPTIONS, but only significnat bc of outliers ####
  #=================================================================#

  # A) Check trait-trait interactions ####
  ----
    
    summary(lm(H.abund.c$log.abund.ratio~(Min.Root.Loca+Ht.veg+Min.Root.Loca+Lamina.thck+LMA+
                                        LDMC+Leaf.Area+myc.frac)^2,
               data=sp.H.traits.c))
  
    # $$$$
    # no significant interactions
  
  # B) Check higher-order effects ####
  ----
  summary(lm(H.abund.c$log.abund.ratio~I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
           I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
         data=sp.H.traits.c))

  # myc.frac^2 is significant.
  
  # C) Fit best model ####
  ----
  H1<-lm(H.abund.c$log.abund.ratio~.+I(myc.frac^2),data=sp.H.traits.c)
  H0<-lm(H.abund.c$log.abund.ratio~1,data=sp.H.traits.c)
  
  best.log.lm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=F)
  
  summary(best.log.lm)
  AIC(best.log.lm) #103.4
  best.log.lm<-update(best.log.lm,~.-Ht.veg)
  summary(best.log.lm)
  AIC(best.log.lm) #104.21
  best.log.lm<-update(best.log.lm,~.-LDMC)
  summary(best.log.lm)
  AIC(best.log.lm) #105.4
  
  # Coefficients:
  #               Estimate Std. Error t value Pr(>|t|)   
  # (Intercept)    -0.4659     0.2134  -2.183  0.03625 * 
  # I(myc.frac^2)   0.3728     0.1192   3.127  0.00367 **

  # 
  # Residual standard error: 1.031 on 33 degrees of freedom
  # Multiple R-squared:  0.2286,	Adjusted R-squared:  0.2052 
  # F-statistic:  9.78 on 1 and 33 DF,  p-value: 0.003672
  
  # D) diagnostic plots ####
  ----
  
  plot(best.log.lm)
  # Variance homoscedastic
  # Residuals - normal (except OXMO) !
  # Leverage - CASC, LYAN & OXMO pulling regression
  
 # E) fit model without outliers ####
  ----
  
  dim(H.abund.c) #35 4
  H.abund.c.no.outliers<-H.abund.c[!rownames(H.abund.c)%in%c('CASC','LYAN','OXMO'),]
  dim(H.abund.c.no.outliers) #32 4
  
  sp.H.traits.c.no.outliers<-sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('CASC','LYAN','OXMO'),]
  dim(sp.H.traits.c.no.outliers) #32 7
  
  RsquareAdj(lm(
    H.abund.c.no.outliers$log.abund.ratio~I(myc.frac^2),data=sp.H.traits.c.no.outliers
  ))  
      
  # -0.02
  
  # $$$$
  # LDMC was only significant because of outliers 
  
  
  # 2.3 lm - PowerTransform response variable - MEETS ASSUMPTIONS, but only significant bc of outliers ####
  #=================================================================================
  
  # A) test trait-trait interactions ####
  ----
  
  summary(lm(H.abund.c$PT.abund.ratio~(Min.Root.Loca+Ht.veg+Min.Root.Loca+Lamina.thck+LMA+
                                    LDMC+Leaf.Area+myc.frac)^2,
        data=sp.H.traits.c))

  # nothing significant
  
  # B) Test higher-level effects ####
  ----
  
  summary(lm(H.abund.c$PT.abund.ratio~I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
          I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
        data=sp.H.traits.c))
  
  # myc.frac^2 is significant. 
  
  # C) fit best minimal model ####
  ----
  H1<-lm(H.abund.c$PT.abund.ratio~.+I(myc.frac^2),data=sp.H.traits.c)  
  H0<-lm(H.abund.c$PT.abund.ratio~1,data=sp.H.traits.c)
  
  best.pt.lm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=F)
  AIC(best.pt.lm) #-23.97
  summary(best.pt.lm)
  best.pt.lm<-update(best.pt.lm,~.-LDMC)
  AIC(best.pt.lm) # -23.3
  summary(best.pt.lm)
  
  # Coefficients:
  #               Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)    0.93031    0.03396  27.393  < 2e-16 ***
  # I(myc.frac^2)  0.07241    0.01897   3.817 0.000563 ***
  # 
  # Residual standard error: 0.1641 on 33 degrees of freedom
  # Multiple R-squared:  0.3063,	Adjusted R-squared:  0.2853 
  # F-statistic: 14.57 on 1 and 33 DF,  p-value: 0.0005631
  
  
  # D) Diagnostic plots ####
  ----
    
  plot(best.pt.lm)
  # Variance homoscedastic !
  # Residuals - almost normal
  # Leverage - CASC, LYAN & OXMO pulling regression
  
  # E) fit model without the outliers ####
  ----
  
  # use same data as in 2.2  
  
  RsquareAdj(lm(
    H.abund.c.no.outliers$PT.abund.ratio~I(myc.frac^2),data=sp.H.traits.c.no.outliers
  ))
  
  # - 0.03
  
  # model only fits the outliers !!! :-(
  
  # F) plot regression ####
  ----
  
  # plot data    
  plot(H.abund.c$PT.abund.ratio~myc.frac,data=sp.H.traits.c,type='n')
  text(sp.H.traits.c$myc.frac,H.abund.c$PT.abund.ratio,label=rownames(sp.H.traits.c),cex=0.7)
  
  # add curve
  points(x=sp.H.traits.c$myc.frac, y=fitted(best.pt.lm), col='red', pch=20)
  
  
  # G) plot plot regression without outliers ####
  ----
    
  plot(H.abund.c.no.outliers$PT.abund.ratio~sp.H.traits.c.no.outliers$myc.frac,type='n')
  
  text(sp.H.traits.c.no.outliers$myc.frac,H.abund.c.no.outliers$PT.abund.ratio,
       label=rownames(sp.H.traits.c),cex=0.7)
  
  # add curve
  points(x=sp.H.traits.c.no.outliers$myc.frac, 
         y=fitted(lm(H.abund.c.no.outliers$PT.abund.ratio~I(myc.frac^2),data=sp.H.traits.c.no.outliers)),
         col='red', pch=20)

  
  # 2.4 - Ratio of log abundance - Not sure model meets assumptions & only significant bc of outliers ####
                                   # model without outlier does not retain any variable!
  
  #================================================
  
  # A) Test trait-trait interactions ####
  ----
   
  summary(lm(H.abund.c$ratio.log.abund~(Min.Root.Loca+Ht.veg+Min.Root.Loca+Lamina.thck+LMA+
                                         LDMC+Leaf.Area+myc.frac)^2,
             data=sp.H.traits.c))
  
  # Ht.veg:LMA, Min.Root.Loca:LMA significant
  # LMA:Leaf.Area, Ht.veg:LDMC, Min.Root.Loca:Leaf.Area, Min.Root.Loca:Leaf.Area marginal
  
  # B) test higher-order effects ####
  ----
  summary(lm(H.abund.c$ratio.log.abund~I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
            I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
          data=sp.H.traits.c))
          
  #LMA^2, Leaf.Area^2 significant
  
  # C) Fit best model ####
  ----
  
  H1<-lm(H.abund.c$ratio.log.abund~.+I(LMA^2)+I(Leaf.Area^2)+I(myc.frac^2)+
           Ht.veg:LMA+Min.Root.Loca:LMA+LMA:Leaf.Area+Ht.veg:LDMC+Min.Root.Loca:Leaf.Area+Min.Root.Loca:Leaf.Area,
         data=sp.H.traits.c)  
  H0<-lm(H.abund.c$ratio.log.abund~1,data=sp.H.traits.c)
  
  best.rl.lm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=F)
  AIC(best.rl.lm) #184.9
  summary(best.rl.lm)
  best.rl.lm<-update(best.rl.lm,~.-LMA:Ht.veg)
  AIC(best.rl.lm) #  185.0
  summary(best.rl.lm)
  best.rl.lm<-update(best.rl.lm,~.-Ht.veg)
  AIC(best.rl.lm) # 187.2 
  # Actually better to keep Ht.veg
  best.rl.lm<-update(best.rl.lm,~.+Ht.veg)
  summary(best.rl.lm)
  
    # Coefficients:
    #             Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)   2.2252     0.6591   3.376  0.00199 ** 
    # I(LMA^2)     -2.5256     0.5412  -4.667 5.57e-05 ***
    # LMA           1.6825     0.6481   2.596  0.01429 *  
    # Ht.veg       -1.9898     1.0003  -1.989  0.05557 .  
    
    # 
    # Residual standard error: 3.13 on 31 degrees of freedom
    # Multiple R-squared:  0.4873,	Adjusted R-squared:  0.4377 
    # F-statistic: 9.823 on 3 and 31 DF,  p-value: 0.0001041
  
 
  # D) Diagnostic plots ####
  ----
    
  plot(best.rl.lm)
  # Variance- homoscedastic, very negative relationship between residuals and fitted values 
  # - non-linear relationship + one huge outlier
  # Residuals - normal, except IMCA
  # Leverage - LYOB & IMCA pulling regression
  
  # E) fit model without outliers ####
  ----
  best.rl.lm$call  # H.abund.c$ratio.log.abund ~ I(LMA^2) + LMA + Ht.veg, data = sp.H.traits.c
  
  dim(H.abund.c) #35 4
  H.abund.c.no.outliers<-H.abund.c[!rownames(H.abund.c)%in%c('LYOB','IMCA'),]
  dim(H.abund.c.no.outliers) #33 4
  
  dim(sp.H.traits.c) #35 7
  sp.H.traits.c.no.outliers<-sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('LYOB','IMCA'),]
  dim(sp.H.traits.c.no.outliers) #33 7
  
  RsquareAdj(lm(
    H.abund.c.no.outliers$ratio.log.abund ~ I(LMA^2) + LMA + Ht.veg, data = sp.H.traits.c.no.outliers
  )) # -0.05
  
  #$$$$$
  # model entirely driven by outliers
  
  # F) Re-run model without outlier
  
  #find outlier to remove
  plot(density(H.abund.c$ratio.log.abund))
  rownames(H.abund.c[H.abund.c$ratio.log.abund<=-20,]) #IMCA
  
  # find significant trait-trait interactions
  
  summary(lm(H.abund.c[rownames(H.abund.c)!='IMCA',]$ratio.log.abund~(Min.Root.Loca+Ht.veg+Min.Root.Loca+Lamina.thck+LMA+
                                          LDMC+Leaf.Area+myc.frac)^2,
             data=sp.H.traits.c[rownames(H.abund.c)!='IMCA',]))
  
  # Lamina.thck:Leaf.Area marginal
  
  # find significant higher-order effects
  summary(lm(H.abund.c[rownames(H.abund.c)!='IMCA',]$ratio.log.abund~I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
               I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
             data=sp.H.traits.c[rownames(H.abund.c)!='IMCA',]))
  
  # no significant terms
  
  # fit best minimal model
  H1<-lm(H.abund.c[rownames(H.abund.c)!='IMCA',]$ratio.log.abund~.+Lamina.thck:Leaf.Area,
         data=sp.H.traits.c[rownames(H.abund.c)!='IMCA',])  
  H0<-lm(H.abund.c[rownames(H.abund.c)!='IMCA',]$ratio.log.abund~1,data=sp.H.traits.c[rownames(H.abund.c)!='IMCA',])
  
  step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=T)
  
  #$$$$
  # Does not select anything !! 
  

  # 2.5 glm, gamma link fct - Meets assumption, but model not significant wihtout outliers ####
  #=========================================#
  # glm.nb - doesn't work bc non-integer values
  # glm gamma doesn't work bc NaNs produced
  #
   
  # A) Check Trait-trait interactions #####
  ----
  
  summary(glm(H.abund.c$abund.ratio+0.01~(Min.Root.Loca+Ht.veg+Min.Root.Loca+Lamina.thck+LMA+
                                          LDMC+Leaf.Area+myc.frac)^2,
                 data=sp.H.traits.c,family=Gamma(link='log')),dispersion=1)
  
  # No significant terms
  
  # B) Check higher-order effects ####
  ----
  
  summary(glm(H.abund.c$abund.ratio~I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
               I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
             data=sp.H.traits.c),family=Gamma(link='log'))
  
  # I(myc.frac^2) significant
      
    
  # C) Fit best model ####
  ----

  H1<-glm(H.abund.c$abund.ratio+0.01~.+I(myc.frac^2),data=sp.H.traits.c,family= Gamma(link='log'))
  # Error: no valid set of coefficients has been found: please supply starting values
  # In addition: Warning message:
  #   In log(ifelse(y == 0, 1, y/mu)) : NaNs produced
  H0<-glm(H.abund.c$abund.ratio+0.01~1,data=sp.H.traits.c,family=Gamma(link='log'))
  
  summary(best.glm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=F),
          dispersion=1)
  
  # lm(formula = H.abund.c$abund.ratio + 0.01 ~ I(myc.frac^2), family = Gamma(link = "log"), 
  # data = sp.H.traits.c)
  
  # Coefficients:
  # Estimate Std. Error z value Pr(>|z|)   
  # (Intercept)   -0.02954    0.20700  -0.143   0.8865   
  # I(myc.frac^2)  0.31347    0.11562   2.711   0.0067 **
  #   ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # (Dispersion parameter for Gamma family taken to be 1)
  # 
  # Null deviance: 37.981  on 34  degrees of freedom
  # Residual deviance: 26.255  on 33  degrees of freedom
  # AIC: 93.921
  
  # fitting Gamma with dispersion=1 is supposed to be equivalent to fitting lognormal. 
  # see: http://stats.stackexchange.com/questions/21447/how-to-specify-a-lognormal-distribution-in-the-glm-family-argument-in-r
  # see: https://stats.stackexchange.com/questions/67547/when-to-use-gamma-glms
  # see: http://seananderson.ca/2014/04/08/gamma-glms.html
  
  # the Gamma family the links inverse, identity and log
  gamma.dispersion(H1) # 0.5984368
  gamma.shape(H1)
  # Alpha: 1.671020
  # SE:    0.366515
  
  # Null deviance = total variance & model deviance = variance explained by model. hence, 
  # D2 = (model$null.deviance - model$deviance) / model$null.deviance
  (best.glm$null.deviance-best.glm$deviance)/best.glm$deviance
  # 0.44
  
  # D) Diagnostic plot ####
  ----
  
  plot(best.glm)
  # Variance in residuals is homoscedastic !
  # QQ-plot - distribution is pretty linear !!
  # Leverage - CASC & LYAN have a lot of influence on model - try model without that datapoint.
  
  # E) fit model without outliers
  ----
  
  best.glm$call  
  # glm(formula = H.abund.c$abund.ratio + 0.01 ~ I(myc.frac^2), family = Gamma(link = "log"), 
  # data = sp.H.traits.c)
  
  dim(H.abund.c) #35 4
  H.abund.c.no.outliers<-H.abund.c[!rownames(H.abund.c)%in%c('CASC','LYAN'),]
  dim(H.abund.c.no.outliers) #33 4
  
  dim(sp.H.traits.c) #35 7
  sp.H.traits.c.no.outliers<-sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('CASC','LYAN'),]
  dim(sp.H.traits.c.no.outliers) #33 7
  
  glm.no.out<-glm(
    H.abund.c.no.outliers$abund.ratio + 0.01 ~ I(myc.frac^2),
    data = sp.H.traits.c.no.outliers,
    family = Gamma(link = "log")
    )
  
  (glm.no.out$null.deviance-glm.no.out$deviance)/glm.no.out$deviance
  # 0.0006
  
  
  # model isn't significant without high-leverage datapoints

  
