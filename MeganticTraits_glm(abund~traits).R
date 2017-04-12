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
    plot(gam(H.abund$abund.ratio~s(myc.frac)+s(Min.Root.Loca)+
             s(LDMC)+s(LDMC),data=sp.H.traits))
    
    #myc.frac appears non-linear
    #Min.Root.Loca appears non-linear
    m<-gam(H.abund.c$abund.ratio~s(sp.H.traits.c$Ht.veg)+s(sp.H.traits.c$Lamina.thck)+
             s(sp.H.traits.c$LMA),data=sp.H.traits.c)
    plot(m)
    #LMA appears correlated
    m<-gam(H.abund.c$abund.ratio~s(sp.H.traits.c$Leaf.Area),data=sp.H.traits.c)
    plot(m)
    
    par(mfrow=c(1,1))
    
    #------------------------------------#
    # myc.frac. and LDMC appear non-linear
        
# 2 - Abund~Traits, Individual variables & Diagnostic plots ####
#=====================================================
     
  # Diagnostic plots
      
  # 1) Residuals vs Fitted - detects non-linear relationships between x & y variables + heteroscedasticity (wedge)
  # 2) QQ plot - shows whether residuals are normally distributed (will follow straight line)
  # 3) Scale-Location - shows whether variance (residuals) increase with incrasing mean (tests
      # for homoscedasticity)
  # 4) Residuals vs Leverage - shows whether individual datapoints have a lot of 'weight' on 
      # the regression
               
  # 2A - Abundance.ratio 
    # pick a few predictor variables
      
    names(sp.H.traits)
    #[1] "Ht.veg"        "Min.Root.Loca" "Max.Root.Loca" "Lamina.thck"   "LMA"           "LDMC"          "Leaf.Area"    
    #[8] "myc.frac"
    
    H.ht<-lm(H.abund$abund.ratio~Ht.veg,data=sp.H.traits)
    summary(H.ht) #n.s.
    plot(H.ht) 
    # Diagnostic plots - residuals not normal
    
    H.min.root<-lm(H.abund$abund.ratio~Min.Root.Loca,data=sp.H.traits)
    summary(H.min.root)#n.s.
    plot(H.min.root)
    # Diagnostic plots - residuals are right skewed; mini heteroscedasticity; EPHE large weight
    
    
    H.myc<-lm(H.abund$abund.ratio~myc.frac,data=sp.H.traits)
    summary(H.myc) # AdjR2 = 0.23, p-val=0.001
    plot(H.myc)
    # Diagnostic plot - residual-fitted relationship; non-normal variance ; CASC large leverage
    
    
    #2B - log.abundance.ratio
    
    H.ht.l<-lm(H.abund$log.abund.ratio~Ht.veg,data=sp.H.traits)
    summary(H.ht.l) # ns
    plot(H.ht.l)
    # Diagnostic plots - homoscedastic, flatter, but non-yet normal residuals, OXMO has large leverage
    
    H.min.root.l<-lm(H.abund$log.abund.ratio~Min.Root.Loca,data=sp.H.traits)
    summary(H.min.root.l)#n.s.
    plot(H.min.root.l)
    # Diagnostic plots - homoscedastic, flatter, but non-yet normal residuals, OXMO has large leverage
    
    H.myc.l<-lm(H.abund$log.abund.ratio~myc.frac,data=sp.H.traits)
    summary(H.myc.l) # ns
    plot(H.myc.l)
    # Diagnostic plot - homoscedastic, flatish, 
    
    #2C - ratio.log.abundances
    
    H.ht.rl<-lm(H.abund$ratio.log.abund~Ht.veg,data=sp.H.traits)
    summary(H.ht.rl) # ns
    plot(H.ht.rl)
    # Diagnostic plots - homoscedastic, flatter, but non-yet normal residuals, IMCA has large leverage
    
    H.min.root.rl<-lm(H.abund$ratio.log.abund~Min.Root.Loca,data=sp.H.traits)
    summary(H.min.root.rl)#n.s.
    plot(H.min.root.rl)
    # Diagnostic plots - homoscedastic, flat, but IMCA has large leverage
    
    H.myc.rl<-lm(H.abund$ratio.log.abund~myc.frac,data=sp.H.traits)
    summary(H.myc.rl) # ns
    plot(H.myc.rl)
    # Diagnostic plot -homoscedastic, flat, but IMCA has large leverage
    
    
    ##=========================
    # Residuals are non-normally distributed & non-linearities with predictor variables. 
    

# 3 - Model selection ####
#===========================================================
 
    # Data Prep
  
    # Remove species with NAs in myc.frac or Leaf.Area, to be able to run stepwise selection
    sp.H.traits.c<-sp.H.traits[!is.na(sp.H.traits$myc.frac),]
    sp.H.traits.c<-sp.H.traits.c[!is.na(sp.H.traits.c$Leaf.Area),]
    dim(sp.H.traits.c)
    # 35 8
    
    # Adjust species list for abundance data to have matching datasets
    H.abund.c<-H.abund[which(rownames(H.abund)%in%rownames(sp.H.traits.c)),]
    dim(H.abund.c)
    #35 4
    
    vif(H1)
    sp.H.traits.c$Max.Root.Loca<-NULL
    
    
    
    
    
  # 3.1 lm - abundance Ratio - doesn't meet assumptions & only significant bc of outliers ####
  #========================================================================#

  # A) Check trait-trait interactions ####
    ----
  H1<-lm(H.abund.c$abund.ratio~.,data=sp.H.traits.c)
    
  
  # B) Check non-linear relationships ####
    ----


  # C) Fit best minimal model ####
    ----
    
  H1<-update(H1,~.+I(myc.frac^2))
  H0<-lm(H.abund.c$abund.ratio~1,data=sp.H.traits.c)
  
  step.lm.H<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace = F)
  
  summary(step.lm.H)
  # lm(formula = H.abund.c$abund.ratio ~ I(myc.frac^2) + myc.frac + 
  # Min.Root.Loca, data = sp.H.traits.c)
  # Coefficients:
  #                Estimate Std. Error t value Pr(>|t|)   
  # (Intercept)    -0.2491     0.8470  -0.294  0.77066   
  # I(myc.frac^2)   0.7189     0.2318   3.102  0.00408 **
  # myc.frac       -0.7269     0.3657  -1.988  0.05571 . 
  # Min.Root.Loca   0.4164     0.3025   1.377  0.17849 
  # Residual standard error: 1.69 on 31 degrees of freedom
  # Multiple R-squared:  0.4697,	Adjusted R-squared:  0.4183
  
  anova(step.lm.H,update(step.lm.H,~.-Min.Root.Loca))
    # Res.Df    RSS Df Sum of Sq      F Pr(>F)
    # 1     31 88.568                           
    # 2     32 93.983 -1   -5.4144 1.8951 0.1785
  # Min.Root.Loca not necessary
  step.lm.H<-update(step.lm.H,~.-Min.Root.Loca)
  
  AIC(step.lm.H)
  # 141.89
  
  # D) Diagnostic plots ####
  ----
  
  plot(step.lm.H)
  # Variance in residuals decreases with mean - homoscedastic, but an x-y relationship remain unexplained by model
  # QQ-plot -  error not normal - s-shaped - bad
  # Leverage - CASC & LYAN have a lot of influence on model - try model without that datapoint.
  
  # E) Model without outliers ####
  ----
  
  
  H.abund.c.no.outliers<-H.abund.c[!rownames(H.abund.c)%in%c('CASC','LYAN'),]
  dim(H.abund.c.no.outliers) #33 4
  sp.H.traits.c.no.outliers<-sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('CASC','LYAN'),]
  dim(sp.H.traits.c.no.outliers) #33 8
  
  summary(lm(
    H.abund.c.no.outliers$abund.ratio~I(myc.frac^2)+myc.frac,data=sp.H.traits.c.no.outliers
                          ))
  # Coefficients:
  #               Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)    1.15518    0.23141   4.992 2.38e-05 ***
  # I(myc.frac^2) -0.03709    0.21574  -0.172    0.865    
  # myc.frac      -0.09088    0.19465  -0.467    0.644    
  # ---
  # 
  # Residual standard error: 0.9436 on 30 degrees of freedom
  # Multiple R-squared:  0.007838,	Adjusted R-squared:  -0.05831 
  # F-statistic: 0.1185 on 2 and 30 DF,  p-value: 0.8887
  
  #-----------------------------------------------------
  # model was significant only because of CASC and LYAN
  
  
  # 3.2 lm - log abundance ratio - meets assumptions & only significnat bc of outliers ####
  #=================================================================#

  # A) Check trait-trait interactions ####
  ----
  
  # B) Check higher-order effects ####
  ----
  ml<-lm(H.abund.c$log.abund.ratio~.-Max.Root.Loca+I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
           I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
         data=sp.H.traits.c)
  summary(ml)
  # myc.frac^2 is significant.
  
  # C) Fit best model ####
  ----
  H0<-lm(H.abund.c$log.abund.ratio~1,data=sp.H.traits.c)
  summary(step.ml<-step(H0,scope=list(lower=H0,upper=ml),direction='both',trace=F))
  
  anova(step.ml,update(step.ml,~.-Ht.veg)) # difference n.s.
  step.ml<-update(step.ml,~.-Ht.veg)
  
  anova(step.ml,update(step.ml,~.-LDMC)) # marginally sifnificant, so remove
  step.ml<-update(step.ml,~.-LDMC)
  
  summary(step.ml)
  #               Estimate Std. Error t value Pr(>|t|)   
  # (Intercept)    -0.4659     0.2134  -2.183  0.03625 * 
  # I(myc.frac^2)   0.3728     0.1192   3.127  0.00367 **
  # ---
  # 
  # Residual standard error: 1.031 on 33 degrees of freedom
  # Multiple R-squared:  0.2286,	Adjusted R-squared:  0.2052 
  # F-statistic:  9.78 on 1 and 33 DF,  p-value: 0.003672
  
  # D) diagnostic plots ####
  ----
  
  plot(step.ml)
  # Variance homoscedastic
  # Residuals - normal !
  # Leverage - CASC, LYAN & OXMO pulling regression
  
 # E) fit model without outliers ####
  ----
    
  summary(lm(
    H.abund.c[rownames(H.abund.c)!='OXMO','log.abund.ratio']~
      LDMC+I(myc.frac^2), data=sp.H.traits.c[rownames(sp.H.traits.c)!='OXMO',]
  ))
  
  #-----------------------
  # LDMC was only significant because of OXMO - 
  
  # 3.3 lm - PowerTransform response variable - meets assumptions & only significant bc of outliers ####
  #=================================================================================
  
  # A) test trait-trait interactions ####
  ----
  
  m<-lm(H.abund.c$pt.abund.ratio~(Min.Root.Loca+Ht.veg+Min.Root.Loca+Lamina.thck+LMA+
                                    LDMC+Leaf.Area+myc.frac)^2,
        data=sp.H.traits.c)
  summary(m)
  
  #nothing significant
  
  # B) Test higher-level effects ####
  ----
  
  m<-lm(H.abund.c$pt.abund.ratio~.-Max.Root.Loca+I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
          I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
        data=sp.H.traits.c)
  summary(m)
  # myc.frac^2 is singificant. 
  
  # C) fit best minimal model ####
  ----
    
  H0<-lm(H.abund.c$pt.abund.ratio~1,data=sp.H.traits.c)
  summary(step.m<-step(H0,scope=list(lower=H0,upper=m),direction='both',trace=F))
  
  # Coefficients:
  #               Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)    0.95222    0.02904  32.785  < 2e-16 ***
  # I(myc.frac^2)  0.05369    0.01611   3.332  0.00219 ** 
  # LDMC           0.04264    0.02603   1.638  0.11117    
  # ---
  # 
  # Residual standard error: 0.1356 on 32 degrees of freedom
  # Multiple R-squared:  0.3498,	Adjusted R-squared:  0.3091 
  # F-statistic: 8.607 on 2 and 32 DF,  p-value: 0.001021
  
  anova(step.m,update(step.m,~.-LDMC)) # models not significantly different.
  step.m<-update(step.m,~.-LDMC)
  
  # D) Diagnostic plots ####
  ----
    
  plot(step.m)
  # Variance homoscedastic
  # Residuals - close to normal
  # Leverage - CASC, LYAN & OXMO pulling regression
  
  # E) fit model without the outliers ####
  ----
    
  mno<-lm(
    H.abund.c[!rownames(H.abund.c)%in%c('CASC','LYAN','OXMO'),'pt.abund.ratio']~
      I(myc.frac^2), data=sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('CASC','LYAN','OXMO'),]
  )
  
  summary(mno)
  # Coefficients:
  #               Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)    0.98322    0.02614  37.616   <2e-16 ***
  # I(myc.frac^2)  0.01313    0.02428   0.541    0.593
  
  # Residual standard error: 0.1049 on 30 degrees of freedom
  # Multiple R-squared:  0.009657,	Adjusted R-squared:  -0.02335 
  # F-statistic: 0.2925 on 1 and 30 DF,  p-value: 0.5926
  
  # model only fits the outliers !!! :-(
  
  # F) plot regression ####
  ----
    
  plot(H.abund.c$pt.abund.ratio~myc.frac,data=sp.H.traits.c,type='n')
  text(sp.H.traits.c$myc.frac,H.abund.c$pt.abund.ratio,label=rownames(sp.H.traits.c),cex=0.7)
  
  # add curve
  points(x=sp.H.traits.c$myc.frac, y=fitted(step.m), col='red', pch=20)
  lines(x=sort(sp.H.traits.c$myc.frac),y=fitted(step.m)[sort(sp.H.traits.c$myc.frac)],col='red',type='b' ) 
  # doesn't work for some reason
  
  # G) plot plot regression without outliers ####
  ----
    
  plot(H.abund.c[!rownames(H.abund.c)%in%c('CASC','LYAN','OXMO'),'pt.abund.ratio']~
         myc.frac,data=sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('CASC','LYAN','OXMO'),],
       type='n')
  text(x=sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('CASC','LYAN','OXMO'),'myc.frac'],
       y= H.abund.c[!rownames(H.abund.c)%in%c('CASC','LYAN','OXMO'),'pt.abund.ratio'],
       label=rownames(sp.H.traits.c)[!rownames(H.abund.c)%in%c('CASC','LYAN','OXMO')],cex=0.7)
  points(x=sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('CASC','LYAN','OXMO'),'myc.frac'],
         y=fitted(mno),
         col='red', pch=20)

  
  # 3.4 - Ratio of log abundance - does not meet assumptions & only significant bc of outliers ####
  #================================================
  
  # A) Test trait-trait interactions ####
  ----
   
    
  m1<-lm(H.abund.c$pt.abund.ratio~.-max.root.Loca^2+I(myc.frac^2),data=sp.H.traits.c)
  m0<-lm(H.abund.c$pt.abund.ratio~1,data=sp.H.traits.c)
  
  step.m<-step(m0,scope=list(lower=m0,upper=m1),direction='both',trace=F)
  summary(step.m)
  # lm(formula = H.abund.c$pt.abund.ratio ~ I(myc.frac^2) + LDMC, 
  # data = sp.H.traits.c)
  
  # Coefficients:
  #          Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)    0.95222    0.02904  32.785  < 2e-16 ***
  # I(myc.frac^2)  0.05369    0.01611   3.332  0.00219 ** 
  # LDMC           0.04264    0.02603   1.638  0.11117
  
  # Residual standard error: 0.1356 on 32 degrees of freedom
  # Multiple R-squared:  0.3498,	Adjusted R-squared:  0.3091 
  # F-statistic: 8.607 on 2 and 32 DF,  p-value: 0.001021
  
  # interactions not retained  
  
  # B) test higher-order effects ####
  ----
  mrl<-lm(H.abund.c$ratio.log.abund~.-Max.Root.Loca+I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
            I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
          data=sp.H.traits.c)
  summary(mrl) #LMA^2 significant
  
  # C) Fit best model ####
  ----
  
  H0<-lm(H.abund.c$ratio.log.abund~1,data=sp.H.traits.c)
  summary(step.mrl<-step(H0,scope=list(lower=H0,upper=mrl),direction='both',trace=F))
  
  anova(step.mrl,update(step.mrl,~.-Ht.veg))
  #   Res.Df    RSS Df Sum of Sq      F  Pr(>F)  
  # 1     31 303.74                              
  # 2     32 342.50 -1   -38.769 3.9568 0.05557 . 
  # marginally significant. remove
  
  # step.mrl<-update(step.mrl,~.-Ht.veg)
  # summary(step.mrl)
  # Coefficients:
  #             Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)   2.1363     0.6873   3.108 0.003933 ** 
  # I(LMA^2)     -2.0807     0.5151  -4.040 0.000313 ***
  # LMA           1.9552     0.6621   2.953 0.005850 ** 
  # ---
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # Residual standard error: 3.272 on 32 degrees of freedom
  # Multiple R-squared:  0.4219,	Adjusted R-squared:  0.3858 
  # F-statistic: 11.68 on 2 and 32 DF,  p-value: 0.0001556
  
  # D) diagnostic plots ####
  ----
    
  plot(step.mrl)
  # Variance- homoscedastic, negative relationship between residuals and fitted values - non-linear relationship?
  # Residuals - not normal
  # Leverage - LYOB & IMCA pulling regression
  
  # E) fit model without outliers ####
  ----
    
  mno2<-lm(
    H.abund.c[!rownames(H.abund.c)%in%c('LYOB','IMCA','OXMO'),'ratio.log.abund']~
      I(LMA^2)+LMA, data=sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('LYOB','IMCA','OXMO'),]
  )
  
  summary(mno2) #AdjR2=0.008, p-val=0.33
  
  # -------------------------------
  # model only works because of outliers

  
  # 3.5 glm.nb - doesn't meet assumption ####
  #=========================================#
  
  # A) Check Trait-trait interactions #####
  ----
    
  # B) Check higher-order effects ####
  ----
    
  # C) Fit best model ####
  ----
  H1<-glm.nb(H.abund.c$abund.ratio~.+I(myc.frac^2)-Max.Root.Loca,data=sp.H.traits.c)
  H0<-glm.nb(H.abund.c$abund.ratio~1,data=sp.H.traits.c)
  
  step.glm.H<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=F)
  summary(step.glm.H)
  
  # glm.nb(formula = H.abund.c$abund.ratio ~ I(myc.frac^2), data = sp.H.traits.c,
  # init.theta = 5.892508028, link = log)
  
  # Coefficients:
  # Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)   -0.05089    0.20302  -0.251    0.802    
  # I(myc.frac^2)  0.32586    0.07283   4.474 7.67e-06 ***
  
  # Null deviance: 50.976  on 34  degrees of freedom
  # Residual deviance: 32.344  on 33  degrees of freedom
  # AIC: 108.7
  
  # Null deviance = total variance & model deviance = variance explained by model. hence, 
  # D2 = (model$null.deviance - model$deviance) / model$null.deviance
  (step.glm.H$null.deviance-step.glm.H$deviance)/step.glm.H$deviance
  # 0.57
  
  # D) Diagnostic plot ####
  ----
  
  plot(step.glm.H)
  # Variance in residuals is homoscedastic
  # QQ-plot - distribution still is still banana-shaped :-(
  # Leverage - CASC & LYAN have a lot of influence on model - try model without that datapoint.
  
  # E) fit model without outliers
  ----
    
  summary(glm.nb(
    H.abund.c.no.outliers$abund.ratio~I(myc.frac^2),data=sp.H.traits.c.no.outliers
    ))
  
  #Coefficients:
  #Estimate Std. Error z value Pr(>|z|)
  #(Intercept)    0.12861    0.22803   0.564    0.573
  #I(myc.frac^2) -0.02696    0.21878  -0.123    0.902
  
  #    Null deviance: 27.319  on 32  degrees of freedom
  #Residual deviance: 27.304  on 31  degrees of freedom
  #AIC: 88.115
  
  # model isn't significant without high-leverage datapoints
  
          # Testing whether we do need an overdispersion parameter (need for negative binomial distribution)
          # by comparing the negative binomial model with with a poisson regression, which does not assume overdispersion
          
          m3 <- glm(H.abund.c$abund.ratio ~ I(myc.frac^2), family = "poisson", data = sp.H.traits.c)
          pchisq(2 * (logLik(step.glm.H) - logLik(m3)), df = 1, lower.tail = FALSE)
          # Doesn't work bc poisson distribution doesn't work with my data
  
  # 3.6 - glm - lognormal link function? #####
  

  
