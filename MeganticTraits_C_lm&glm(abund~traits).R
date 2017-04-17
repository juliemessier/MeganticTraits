#<<TABLE OF CONTENTS>>
# 0-  Make trait and abundance dataframes correspond based on shared species
# 1- Data Exploration, scatterplot, tree models & GAM
# 2- Abund~Traits models
#   2.1 lm - abundance Ratio - model ns.  Does not meet assumptions
#   2.2 lm - log abundance ratio - Model significant and MEETS ASSUMPTIONS
#   2.3 lm - PowerTransform response variable - Model significiant & MEETS ASSUMPTIONS 
#   2.4 lm - Ratio of log abundance - Model marginally significant & does not meets assumptions 
#   2.5g lm (gamma) - abunance ratio - meets assumption when 3 outliers are removed


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

load(file=paste0(wrk.dir,'list.trait.names.canopy.layer.Rdata'))# C.trait.names
load(file=paste0(wrk.dir,"Canopy.Layer.Traits.Standardized.RData")) # C.traits
load(file=paste0(wrk.dir,'species-level.traits.Canopy.layer.Rdata')) # sp.C.traits

# load(file=paste0(wrk.dir,'Species.abundances.full.data.Rdata')) # object = abund
load(file=paste0(wrk.dir,'Species.abundances.Inf.removed.Rdata')) # object = abund.c
#============================================================================================

# B - CANOPY LAYER, ABUNDANCE

C.abund<-abund.c[abund.c$Layer=='C',c('abund.ratio','log.abund.ratio','ratio.log.abund',"PT.abund.ratio" )]
dim(C.abund)

# 21 4

# 0-  Make trait and abundance dataframes correspond based on shared species ####
#================================================================================

  # 0A - Canopy layer

    C.abund<-abund.c[which(rownames(abund.c)%in%rownames(sp.C.traits)),c('abund.ratio','log.abund.ratio',
                                                                         'ratio.log.abund',"PT.abund.ratio" )]
    
    C.abund<-C.abund[order(rownames(C.abund)),]
    View(C.abund)
    dim(C.abund) #21 4
    
    dim(sp.C.traits)
    # 24 4
    sp.C.traits<-sp.C.traits[which(rownames(sp.C.traits)%in%rownames(C.abund)),]
    dim(sp.C.traits) 
    # 21 4
    all(rownames(C.abund)==rownames(sp.C.traits)) # TRUE

# 1 - Data Exploration - Scatterplots, Tree models and GAMs #### 
#===========================================================

    # Scatterplots - pairwise interactions?
    
    names(C.abund) # "abund.ratio"     "log.abund.ratio" "ratio.log.abund" "pt.abund.ratio" 
    pairs(cbind(sp.C.traits,C.abund[,'abund.ratio']),panel=panel.smooth)
    #yes, LMA-LDMC ; Lamina.thck-LMA
    # variance appears homoscedastic
    
    pairs(cbind(sp.C.traits,C.abund[,'log.abund.ratio']),panel=panel.smooth)
    # variance appears homoscedastic except for one outlier
    
    pairs(cbind(sp.C.traits,C.abund[,'ratio.log.abund']),panel=panel.smooth)
    # can't tell variance homogeneity - two outliers too large
    
    pairs(cbind(sp.C.traits,C.abund[,'PT.abund.ratio']),panel=panel.smooth)
    # variance appears homoscedastic
    
    # Tree models - non-linearities
  
    C.tree.model<-tree(C.abund$abund.ratio~.,data=sp.C.traits)
    plot(C.tree.model)
    text(C.tree.model)
    title('canopy layer, abund.ratio vs traits')
    # Lamina.thck is most important variable, 
    # for those with low Lamina.thck, Leaf.Area matters
    C.tree.model #  null deviance = 184.9
    1-(deviance(C.tree.model)/184.9) # 0.32
    
    C.tree.model.l<-tree(C.abund$log.abund.ratio~.,data=sp.C.traits)
    plot(C.tree.model.l)
    text(C.tree.model.l)
    title('canopy layer, log.abund.ratio vs traits')
    # Lamina.thck is most important variable, 
    # for those with low Lamina.thck, Leaf.Area matters
    C.tree.model.l # null deviance = 48.60
    1-(deviance(C.tree.model.l)/48.6) # 0.55
    
    C.tree.model.rl<-tree(C.abund[!rownames(C.abund)%in%c('FRNI','TSCA'),]$ratio.log.abund~.,
                          data=sp.C.traits[!rownames(sp.C.traits)%in%c('FRNI','TSCA'),])
    plot(C.tree.model.rl)
    text(C.tree.model.rl)
    title('canopy layer, ratio.log.abund vs traits')
    # only LDMc retained, non-linear relationship
    C.tree.model.rl # null deviance = 2647
    1-(deviance(C.tree.model.rl)/2647) # 0.19
    
    C.tree.model.pt<-tree(C.abund$PT.abund.ratio~.,data=sp.C.traits)
    plot(C.tree.model.pt)
    text(C.tree.model.pt)
    title('canopy layer, PT.ratio.abund vs traits')
    # Lamina.thck is most important variable, 
    # for those with low Lamina.thck, Leaf.Area matters
    C.tree.model.pt # null deviance = 1.41
    1-(deviance(C.tree.model.pt)/1.41) # 0.57
    
    # GAM - non-linearities?  
    
    names(sp.C.traits) 
    # "Lamina.thck" "LMA"         "LDMC"        "Leaf.Area"  
    
    par(mfrow=c(2,1))
    plot(gam(C.abund$abund.ratio~s(Lamina.thck)+s(LMA)
               ,data=sp.C.traits))
    
    #Lamina.thck appears non-linear
    plot(gam(C.abund$abund.ratio~s(LDMC),data=sp.C.traits))
    # No relationships
    
    plot(gam(C.abund$abund.ratio~s(Leaf.Area),data=sp.C.traits))
    # Leaf.Area appears non-linear
    
    par(mfrow=c(1,1))
    
    #------------------------------------#
    # Lamina.thck and Leaf.Area appear non-linear
        
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
  
    names(sp.C.traits)
    ## "Lamina.thck" "LMA"         "LDMC"        "Leaf.Area"  
    dim(sp.C.traits)
    # 21 4
    
    # Remove species with NAs in myc.frac or Leaf.Area, to be able to run stepwise selection
    sp.C.traits.c<-sp.C.traits[!is.na(sp.C.traits$Leaf.Area),]
    dim(sp.C.traits.c)
    # 14 4
    
    # Adjust species list for abundance data to have matching datasets
    C.abund.c<-C.abund[which(rownames(C.abund)%in%rownames(sp.C.traits.c)),]
    dim(C.abund.c)
    # 17 4
    
  # 2.1 lm - abundance Ratio - model ns.  Does not meet assumptions ####
  #========================================================================#

  # A) Check trait-trait interactions ####
    ----
   summary(lm(C.abund.c$abund.ratio~(Lamina.thck+LMA+LDMC+Leaf.Area)^2,
              data=sp.C.traits.c))
      
      # $$$$
      # LDMC:Leaf.Area, LMA:LDMC significant 
      # LMA:Leaf.Area, Lamina.thck:LMA marginally significant
    
  
  # B) Check non-linear relationships ####
    ----
    summary(lm(C.abund.c$abund.ratio~I(Lamina.thck^2)+I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2),
         data=sp.C.traits.c))
    
      # $$$$
      #I(Lamina.thck^2) significant
      #I(Lamina.thck^2) marginally significant 

  # C) Fit best minimal model ####
    ----
    
  H1<-lm(C.abund.c$abund.ratio~.+I(Lamina.thck^2)+I(Lamina.thck^2)+
           (LDMC:Leaf.Area)+(LMA:LDMC)+(LMA:Leaf.Area)+(Lamina.thck:LMA)
           ,data=sp.C.traits.c)
  H0<-lm(C.abund.c$abund.ratio~1,data=sp.C.traits.c)
  
  C.best.lm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace = T)
  
  summary(C.best.lm)
  AIC(C.best.lm) # 90.6
  C.best.lm<-update(C.best.lm,~.-LMA)
  summary(C.best.lm)
  AIC(C.best.lm) # 91.4
  
  # Coefficients:
  #             Estimate Std. Error t value Pr(>|t|)  
  # (Intercept)    1.246      1.013   1.229   0.2379  
  # Lamina.thck   -2.716      1.494  -1.817   0.0892 .
  # ---
  # 
  # Residual standard error: 3.179 on 15 degrees of freedom
  # Multiple R-squared:  0.1805,	Adjusted R-squared:  0.1258 
  # F-statistic: 3.303 on 1 and 15 DF,  p-value: 0.08919
  
  #$$$$
  # Model is n.s.
  
  # D) Diagnostic plots ####
  ----
  
  plot(C.best.lm)
  # Variance - homoscedastic
  # QQ-plot -  error not normal - J-shaped - bad
  # Leverage - PRPE almost at 0.5 on Cook's distance
  
  # 2.2 lm - log abundance ratio - Model significant and MEETS ASSUMPTIONS  ####
                                  # Lamina.thickness explains 32% of var
                                  # robust to outliers 
  #=================================================================#

  # A) Check trait-trait interactions ####
  ----
    
    summary(lm(C.abund.c$log.abund.ratio~(Lamina.thck+LMA+LDMC+Leaf.Area)^2,
               data=sp.C.traits.c))
  
  # $$$$
  # LDMC:Leaf.Area,  significant 
  # LMA:Leaf.Area, LMA:LDMC marginally significant
  
    # $$$$
    # no significant interactions
  
  # B) Check higher-order effects ####
  ----
    summary(lm(C.abund.c$log.abund.ratio~I(Lamina.thck^2)+I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2),
               data=sp.C.traits.c))
  
  # $$$$
  #I(Lamina.thck^2) significant
  
  # C) Fit best model ####
  ----
  H1<-lm(C.abund.c$log.abund.ratio~.+I(Lamina.thck^2),data=sp.C.traits.c)
  H0<-lm(C.abund.c$log.abund.ratio~1,data=sp.C.traits.c)
  
  C.best.log.lm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=F)
  
  summary(C.best.log.lm)
  AIC(C.best.log.lm) # 64.54
  
  # Coefficients:
  #             Estimate Std. Error t value Pr(>|t|)  
  # (Intercept)  -0.9990     0.4595  -2.174   0.0461 *
  # Lamina.thck  -1.9653     0.6776  -2.901   0.0110 *
  # 
  # Residual standard error: 1.441 on 15 degrees of freedom
  # Multiple R-squared:  0.3593,	Adjusted R-squared:  0.3166 
  # F-statistic: 8.413 on 1 and 15 DF,  p-value: 0.01098
  
  # D) diagnostic plots ####
  ----
  
  plot(C.best.log.lm)
  # Variance homoscedastic
  # Residuals - normal !
  # Leverage - FRNI almost outlier
  
 # E) fit model without outliers ####
  ----
  
  dim(C.abund.c) #17 4
  C.abund.c.no.outliers<-C.abund.c[!rownames(C.abund.c)%in%c('FRNI'),]
  dim(C.abund.c.no.outliers) #16 4
  
  sp.C.traits.c.no.outliers<-sp.C.traits.c[!rownames(sp.C.traits.c)%in%c('FRNI'),]
  dim(sp.C.traits.c.no.outliers) #16 4
  
  RsquareAdj(lm(
    C.abund.c.no.outliers$log.abund.ratio~Lamina.thck,data=sp.C.traits.c.no.outliers
  ))  # 0.35
  # $$$$
  # model still significant despitethat datapoint
  
      
  # 0.35
  
  # F) plot regression ####
  ----
    
  # plot data    
  plot(C.abund.c$log.abund.ratio~Lamina.thck,data=sp.C.traits.c,type='n',
       main='canopy layer - log.abund.ratio vs Lamina thickness')
  text(sp.C.traits.c$Lamina.thck,C.abund.c$log.abund.ratio,label=rownames(sp.C.traits.c),cex=0.7)
  
  # add curve
  points(x=sp.C.traits.c$Lamina.thck, y=fitted(C.best.log.lm), col='red', pch=20)
  abline(C.best.log.lm,col='red')
  text(0.2,2,labels='AdjR2 = 0.32 \n y= -1.0 -2.0x',family='serif' )
  
  
  # 2.3 lm - PowerTransform response variable - Model significiant & MEETS ASSUMPTIONS ####
            # Lamina thickness explains 0.33 of variance   
            # model robust to outlier
  #=================================================================================
  
  # A) Check trait-trait interactions ####
  ----
    
    summary(lm(C.abund.c$PT.abund.ratio~(Lamina.thck+LMA+LDMC+Leaf.Area)^2,
               data=sp.C.traits.c))
  
  # $$$$
  # LDMC:Leaf.Area, LMA:LDMC  significant 
  # LMA:Leaf.Area, marginally significant
  
  # $$$$
  # no significant interactions
  
  # B) Check higher-order effects ####
  ----
    summary(lm(C.abund.c$PT.abund.ratio~I(Lamina.thck^2)+I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2),
               data=sp.C.traits.c))
  
  # $$$$
  #I(Lamina.thck^2) significant
  
  # C) fit best minimal model ####
  ----
  H1<-lm(C.abund.c$PT.abund.ratio~.+I(Lamina.thck^2)+LDMC:Leaf.Area+LMA:LDMC+LMA:Leaf.Area,
         data=sp.C.traits.c)
  H0<-lm(C.abund.c$PT.abund.ratio~1,data=sp.C.traits.c)
  
  C.best.pt.lm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=T)
  
  summary(C.best.pt.lm)
  AIC(C.best.pt.lm) #  4.13
  
  #Coefficients:
  #             Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)  0.87021    0.07776  11.191 1.11e-08 ***
  # Lamina.thck -0.33977    0.11466  -2.963  0.00967 ** 
  # ---
  # Residual standard error: 0.2439 on 15 degrees of freedom
  # Multiple R-squared:  0.3692,	Adjusted R-squared:  0.3272 
  # F-statistic:  8.78 on 1 and 15 DF,  p-value: 0.009672
  # 
  
  # D) Diagnostic plots ####
  ----
    
  plot(C.best.pt.lm)
  # Variance homoscedastic !
  # Residuals - completely normal ! :-D
  # Leverage - no outliers
  
  # E) fit model without the outliers ####
  ----
  
  # use same data as in 2.2  
  
  RsquareAdj(lm(
    C.abund.c.no.outliers$PT.abund.ratio~Lamina.thck,data=sp.C.traits.c.no.outliers
  ))
  
  # 0.31
  
  # model robust to outlier !
  
  # F) plot regression ####
  ----
  
  # plot data    
  plot(C.abund.c$PT.abund.ratio~Lamina.thck,data=sp.C.traits.c,type='n',
       main='canopy layer - PT.abund.ratio vs Lamina thickness')
  text(sp.C.traits.c$Lamina.thck,C.abund.c$PT.abund.ratio,label=rownames(sp.C.traits.c),cex=0.7)
  
  # add curve
  points(x=sp.C.traits.c$Lamina.thck, y=fitted(C.best.pt.lm), col='red', pch=20)
  abline(C.best.pt.lm,col='red')
  text(0.2,1.4,labels='AdjR2 = 0.33 \n y= 0.87 -0.34x',family='serif' )
  
  
  # 2.4 lm - Ratio of log abundance - Model marginally significant & does not meets assumptions ####
  #================================================
  
  # A) Test trait-trait interactions ####
  ----
    
  C.abund.cc<-C.abund.c[C.abund.c$ratio.log.abund!=Inf,]  
  dim(C.abund.cc) # 16 4
  sp.C.traits.cc<-sp.C.traits.c[C.abund.c$ratio.log.abund!=Inf,]
  dim(sp.C.traits.cc) # 16 4 
   
  summary(lm(C.abund.cc$ratio.log.abund~(Lamina.thck+LMA+LDMC+Leaf.Area)^2,
             data=sp.C.traits.cc))
  
  # no significant interactions
  
  # B) test higher-order effects ####
  ----
  summary(lm(C.abund.cc$ratio.log.abund~I(Lamina.thck^2)+I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2),
          data=sp.C.traits.cc))
          
  # I(Lamina.thck^2) marginally significant
  
  # C) Fit best model ####
  ----
  
  H1<-lm(C.abund.cc$ratio.log.abund~.+I(Lamina.thck^2),
         data=sp.C.traits.cc)  
  H0<-lm(C.abund.cc$ratio.log.abund~1,data=sp.C.traits.cc)
  
  C.best.rl.lm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=F)
  AIC(C.best.rl.lm) #129.7
  summary(C.best.rl.lm) 

  
    # Coefficients:
    #             Estimate Std. Error t value Pr(>|t|)  
    # (Intercept)    1.948      3.085   0.631   0.5380  
    # Leaf.Area      6.058      3.295   1.839   0.0873 .
    # 
    # Residual standard error: 12.34 on 14 degrees of freedom
    # Multiple R-squared:  0.1945,	Adjusted R-squared:  0.137 
    # F-statistic:  3.38 on 1 and 14 DF,  p-value: 0.08728
  
  #$$$$
  # marginally significant
  
 
  # D) Diagnostic plots ####
  ----
    
  plot(C.best.rl.lm)
  # Variance- homoscedastic, negative relationship between residuals and fitted values 
  # - non-linear relationship
  # Residuals - normal, except 3 outliers
  # Leverage - FRAm & ACPE pulling regression
  

  # 2.5 glm (gamma) - abunance ratio - meets assumption when 3 outliers are removed ####
          # significant  Model without outliers (3 species), Lamina.thck^2, Pseudo AdjR2 = 0.74
  #=========================================#
   
  # A) Check Trait-trait interactions #####
  ----
  
  summary(glm(C.abund.c$abund.ratio+0.01~(LMA+LDMC+Lamina.thck+Leaf.Area)^2,
                 data=sp.C.traits.c,family=Gamma(link='log')),dispersion=1)
  
  # LDMC:Leaf.Area, LMA:LDMC,  significant terms
  # LDMC:Lamina.thck marginally significant 
  
  # B) Check higher-order effects ####
  ----
  
  summary(glm(C.abund.c$abund.ratio~I(Lamina.thck^2)+ I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2),
             data=sp.C.traits.c),family=Gamma(link='log'))
  
  # I(Lamina.thck^2) significant
  # I(Leaf.Area^2) marginally significant
      
    
  # C) Fit best model ####
  ----

  H1<-glm(C.abund.c$abund.ratio+0.01~.+I(Lamina.thck^2)+I(Leaf.Area^2)+
            LDMC:Leaf.Area+LMA:LDMC+LDMC:Lamina.thck,
          data=sp.C.traits.c,family= Gamma(link='log'))

  H0<-glm(C.abund.c$abund.ratio+0.01~1,data=sp.C.traits.c,family=Gamma(link='log'),maxit=100)
  
  summary(C.best.glm<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=F),
          dispersion=1)
  
  # Call:
  # glm(formula = C.abund.c$abund.ratio + 0.01 ~ I(Leaf.Area^2) + 
  #       Leaf.Area + I(Lamina.thck^2), family = Gamma(link = "log"), 
  #       data = sp.C.traits.c)
  # 
  # Deviance Residuals: 
  #     Min       1Q   Median       3Q      Max  
  # -2.1057  -0.4911  -0.1556   0.2178   1.3767  
  # 
  # Coefficients:
  #                  Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)       -1.8448     0.4638  -3.977 6.96e-05 ***
  # I(Leaf.Area^2)     0.8109     0.2462   3.294 0.000986 ***
  # Leaf.Area         -0.9517     0.3483  -2.732 0.006295 ** 
  # I(Lamina.thck^2)   3.2533     0.6363   5.113 3.17e-07 ***
  # ---
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # (Dispersion parameter for Gamma family taken to be 1)
  # 
  #     Null deviance: 35.035  on 16  degrees of freedom
  # Residual deviance: 14.849  on 13  degrees of freedom
  # AIC: 53.814
  # 
  # Number of Fisher Scoring iterations: 9
  
  # Warning messages:
  #   1: glm.fit 1-4: algorithm did not converge 
  
  # fitting Gamma with dispersion=1 is supposed to be equivalent to fitting lognormal. 
  # see: http://stats.stackexchange.com/questions/21447/how-to-specify-a-lognormal-distribution-in-the-glm-family-argument-in-r
  # see: https://stats.stackexchange.com/questions/67547/when-to-use-gamma-glms
  # see: http://seananderson.ca/2014/04/08/gamma-glms.html
  
  # the Gamma family the links inverse, identity and log
  gamma.dispersion(C.best.glm) # 0.78
  gamma.shape(C.best.glm)
  # Alpha: 1.2860094
  # SE:    0.3967661
  
  # Null deviance = total variance & model deviance = variance explained by model. hence, 
  # D2 = (model$null.deviance - model$deviance) / model$null.deviance
  (C.best.glm$null.deviance-C.best.glm$deviance)/C.best.glm$null.deviance
  # 0.58
   C.best.glm # null deviance =35.03
   1-(deviance(C.best.glm)/35.03)
  # 0.58
  
  # D) Diagnostic plot ####
  ----
  
  plot(C.best.glm)
  # Variance in residuals is homoscedastic !
  # QQ-plot - residual distribution not linear on LEFT hand-side !! (osvi, frni, potr)
  # Leverage - PRPE & OSVI have a lot of influence on model - try model without that datapoint.
  
  # E) fit model without outliers ####
  ----
  
  C.best.glm$call  
  # glm(formula = C.abund.c$abund.ratio + 0.01 ~ I(Leaf.Area^2) + 
  # Leaf.Area + I(Lamina.thck^2), family = Gamma(link = "log"), 
  # data = sp.C.traits.c, maxit = 100)
  
  dim(C.abund.cc) #16 4
  C.abund.cc.no.outliers<-C.abund.cc[!rownames(C.abund.cc)%in%c('PRPE','OSVI'),]
  dim(C.abund.cc.no.outliers) #14 4
  
  dim(sp.C.traits.cc) #16 4
  sp.C.traits.cc.no.outliers<-sp.C.traits.cc[!rownames(sp.C.traits.cc)%in%c('PRPE','OSVI'),]
  dim(sp.C.traits.cc.no.outliers) #14 4
  
  C.glm.no.out<-glm(
    C.abund.cc.no.outliers$abund.ratio + 0.01 ~ I(Leaf.Area^2) + 
    Leaf.Area + I(Lamina.thck^2), 
    data = sp.C.traits.cc.no.outliers,family = Gamma(link = "log")
    )
  
  (C.glm.no.out$null.deviance-C.glm.no.out$deviance)/C.glm.no.out$null.deviance
  # 0.78
  AIC(C.glm.no.out) # 34.52
  
  summary(C.glm.no.out)
  
  # Coefficients:
  #                  Estimate Std. Error t value Pr(>|t|)   
  # (Intercept)       -1.4165     0.4572  -3.098  0.01129 * 
  # I(Leaf.Area^2)     0.3817     0.2313   1.650  0.12993   
  # Leaf.Area         -0.2349     0.3832  -0.613  0.55356   
  # I(Lamina.thck^2)   2.8878     0.7179   4.023  0.00243 **
  #
  # (Dispersion parameter for Gamma family taken to be 0.5000209)
  # 
  # Null deviance: 20.773  on 13  degrees of freedom
  # Residual deviance:  4.556  on 10  degrees of freedom
  # AIC: 34.523
  
  # two terms n.s. 
  
  # F) update model without outliers ####
  ----
  
  C.glm.no.out<-update(C.glm.no.out,~.-Leaf.Area)
  summary(C.glm.no.out) # Leaf.Area^2 marginal
  AIC(C.glm.no.out) # 32.94
  C.glm.no.out<-update(C.glm.no.out,~.-I(Leaf.Area^2))
  summary(C.glm.no.out)  
  AIC(C.glm.no.out) #37.90
  # put Leaf.Area back in. 
  C.glm.no.out<-update(C.glm.no.out,~.+I(Leaf.Area^2))
  
  # G) check assumptions of model without outliers ####
  ----
  
  plot(C.glm.no.out)
  # variance homoscedastic
  # residuals normal !
  # FRAM has large leverage
  
  # H) remove new outliers ####
  ----
    C.glm.no.out$call  
  # glm(formula = C.abund.cc.no.outliers$abund.ratio + 0.01 ~ I(Lamina.thck^2) + 
  # I(Leaf.Area^2), family = Gamma(link = "log"), data = sp.C.traits.cc.no.outliers)

  C.abund.cc.no.outliers2<-C.abund.cc.no.outliers[!rownames(C.abund.cc.no.outliers)%in%c('FRAM'),]
  dim(C.abund.cc.no.outliers2) #13 4
  
  sp.C.traits.cc.no.outliers2<-sp.C.traits.cc.no.outliers[!rownames(sp.C.traits.cc.no.outliers)%in%c('FRAM'),]
  dim(sp.C.traits.cc.no.outliers2) #13 4
  
  C.glm.no.out2<-glm(
    C.abund.cc.no.outliers2$abund.ratio + 0.01 ~ I(Leaf.Area^2) + 
      I(Lamina.thck^2), 
    data = sp.C.traits.cc.no.outliers2,family = Gamma(link = "log")
  )
  
  (C.glm.no.out2$null.deviance-C.glm.no.out2$deviance)/C.glm.no.out2$null.deviance
  # 0.74
  AIC(C.glm.no.out2) # 26 - better than 37
  
  summary(C.glm.no.out2)
  # leaf area not significant so update
  
  C.glm.no.out2<-update(C.glm.no.out2,~.-I(Leaf.Area^2))
  summary(C.glm.no.out2)
  #Coefficients:
  #                 Estimate Std. Error t value Pr(>|t|)    
  #(Intercept)       -1.0464     0.2588  -4.043 0.001939 ** 
  #I(Lamina.thck^2)   2.4548     0.4608   5.327 0.000242 ***
  #AIC(C.glm.no.out2) # 24.37 - better than 26
  
  # I) check assumptions of new model without outliers ####
  #----
  
  plot(C.glm.no.out2)
  ## variance homoscedastic
  # errors normal
  # no more outliers, but PRSE is close to cook's value of 0.5
  
  C.glm.no.out2$call
  # glm(formula = C.abund.cc.no.outliers2$abund.ratio + 0.01 ~ I(Lamina.thck^2), 
  # family = Gamma(link = "log"), data = sp.C.traits.cc.no.outliers2)
  