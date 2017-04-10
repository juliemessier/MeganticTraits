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

#<<LOAD FILES>>

load(file=paste0(wrk.dir,'list.trait.names.herbivory.layer.Rdata'))
load(file=paste0(wrk.dir,"Herbaceous.Layer.Traits.Standardized.RData"))
load(file=paste0(wrk.dir,'species-level.traits.Herbaceous.layer.Rdata'))

load(file=paste0(wrk.dir,'list.trait.names.canopy.layer.Rdata'))
load(file=paste0(wrk.dir,"Canopy.Layer.Traits.Standardized.RData"))
load(file=paste0(wrk.dir,'species-level.traits.Canopy.layer.Rdata'))

load(file=paste0(wrk.dir,'Species.abundances.Inf.removed.Rdata'))

#============================================================================================

# 0-  Make trait and abundance dataframes correspond based on shared species ####

  # 0A - Herbivory layer
    H.abund<-abund[which(rownames(abund)%in%rownames(sp.H.traits)),c('presence.change','abund.change')]
    H.abund<-H.abund[order(rownames(H.abund)),]
    View(H.abund)
    dim(H.abund) #45 2
    
    dim(sp.H.traits)
    # 53 8
    sp.H.traits<-sp.H.traits[which(rownames(sp.H.traits)%in%rownames(H.abund)),]
    dim(sp.H.traits) 
    #45 8
    all(rownames(H.abund)==rownames(sp.H.traits)) # TRUE

  # 0B - Canopy layer
    C.abund<-abund[which(rownames(abund)%in%rownames(sp.C.traits)),c('presence.change','abund.change')]
    C.abund<-C.abund[order(rownames(C.abund)),]
    View(C.abund)
    dim(C.abund) 
    #21 2

    dim(sp.C.traits)
    #22 4
    dim (C.abund)
    #21 2
    sp.C.traits<-sp.C.traits[which(rownames(sp.C.traits)%in%rownames(C.abund)),]
    all(rownames(sp.C.traits)==rownames(C.abund)) 
    # TRUE

# 1 - Data Exploration ####

  # 1A - Herbivory layer
    # 1A.1 - Abundance
  
      # Any pairwise interactions?
      plot.new()
      pairs(c(sp.H.traits,H.abund[,c('abund.change','presence.change')]),panel=panel.smooth)
      #yes, LMA-LDMC ; Root.Loca-myc.frac ; 
      
      #Any important interactions on response variable?
    
      model<-tree(H.abund$abund.change~.,data=sp.H.traits)
      plot(model)
      text(model)
      # Myc.frac is most important variable, 
      # for those with high myc.frac, LDMC matters,
        
    # 1A.2 - Presence
        
      #Any important interactions on response variable?
      
      model<-tree(H.abund$presence.change~.,data=sp.H.traits)
      plot(model)
      text(model)
      # LMA is most important variable, 
      # for those with low LMA, LDMC matters
      # for those with low LMA and high LDMC, Ht.Veg matters
      # Non-lineariry in LMA bc it appears twice
      
  # 1B.1 - Canopy layer
      
      # 1B.1 - Abundance
      
      # Pairwise interactions?
      plot.new()
      pairs(c(sp.C.traits,C.abund),panel=panel.smooth)
      #yes, LMA-LDMC ; Lamina.thickness-LMA ; 
      
      #Any important interactions on response variable?
      
      model<-tree(C.abund$abund.change~.,data=sp.C.traits)
      plot(model)
      text(model)
      # Lamina.thickness is most important variable, 
      # for those with low Lamina.thickness, Leaf.Area matters,
      
      # 1B.2 - Presence
      
      #Any important interactions on response variable?
      
      model<-tree(C.abund$presence.change~.,data=sp.C.traits)
      plot(model)
      text(model)
      # Leaf Area is the only important variable

# 2 - Abund~Traits, Individual variables - look at diagnostic plots ####
     
  # 1) Residuals vs Fitted - detects non-linear relationships between x & y variables
  # 2) QQ plot - shows whether residuals are normally distributed (will follow straight line)
  # 3) Scale-Location - shows whether variance (residuals) increase with incrasing mean (tests
      # for homoscedasticity)
  # 4) Residuals vs Leverage - shows whether individual datapoints have a lot of 'weight' on 
      # the regression
               
  # 2A - Abundance
    names(sp.H.traits)
    #[1] "Ht.veg"        "Min.Root.Loca" "Max.Root.Loca" "Lamina.thck"   "LMA"           "LDMC"          "Leaf.Area"    
    #[8] "myc.frac"
    
    H.ht<-lm(H.abund$abund.change~Ht.veg,data=sp.H.traits)
    summary(H.ht) #n.s.
    plot(H.ht) 
    # Diagnostic plots - residuals are right skewed
    
    H.min.root<-lm(H.abund$abund.change~Min.Root.Loca,data=sp.H.traits)
    summary(H.min.root)#n.s.
    plot(H.min.root)
    # Diagnostic plots - residuals are right skewed; mini heteroscedasticity; EPHE large weight
    
    
    H.max.root<-lm(H.abund$abund.change~Max.Root.Loca,data=sp.H.traits)
    summary(H.max.root)#n.s.
    plot(H.max.root)
    # Diagnostic plots - residuals are right skewed; mini heteroscedasticity;EPHE large weight
    
    H.lamina<-lm(H.abund$abund.change~Lamina.thck,data=sp.H.traits)
    summary(H.lamina)#n.s.
    plot(H.lamina)
    # Diagnostic plots - residuals are right skewed;mini heteroscedasticity;
    
    H.LMA<-lm(H.abund$abund.change~LMA,data=sp.H.traits)
    summary(H.LMA) #n.s.
    plot(H.LMA)
    # Diagnostic plots - residuals are right skewed;
    
    H.LDMC<-lm(H.abund$abund.change~LDMC,data=sp.H.traits)
    summary(H.LDMC)#n.s.
    # Diagnostic plots normal
    
    H.area<-lm(H.abund$abund.change~Leaf.Area,data=sp.H.traits)
    summary(H.area)#n.s.
    plot(H.area)
    # Diagnostic plots - residuals are right skewed;LYAN large weight
    
    H.myc<-lm(H.abund$abund.change~myc.frac,data=sp.H.traits)
    summary(H.myc) # AdjR2 = 0.23, p-val=0.001
    plot(H.myc)
    # Diagnostic plot - non-linear relationship; heteroscedasticity ; CASC large weight
    
  #2B 

# 3 - log- and power- transform response, linear models ####

powerTransform(H.abund$abund.change+0.01) # 0.1519625
plot(density((H.abund$abund.change+0.01)^0.1519625))
plot(density(log(H.abund$abund.change+0.01)))

# Abund~ Individual traits & Diagnostic plots
H.abund$log.abund.change<-log(H.abund$abund.change+0.01)
  # LMA, LDMC & myc.frac marginally significantly correlated with abundance with log transform, residuals 
  # more normally distributed, but variance heteroscedastic and OXMO has large weight

H.abund$pt.abund.change<-(H.abund$abund.change+0.01)^0.1519625
  # Residuals right-skewed, variance homoscedastic

# 4 - Stepwise selection on H layer ####

  # 4.1 lm - full model, no interactions ####

  # Remove species with NAs in myc.frac or Leaf.Area, to be able to run stepwise selection
  sp.H.traits.c<-sp.H.traits[!is.na(sp.H.traits$myc.frac),]
  sp.H.traits.c<-sp.H.traits.c[!is.na(sp.H.traits.c$Leaf.Area),]
  dim(sp.H.traits.c)
  # 35 8
  
  # Adjust species list for abundance data to have matching datasets
  H.abund.c<-H.abund[which(rownames(H.abund)%in%rownames(sp.H.traits.c)),]
  dim(H.abund.c)
  #35 4
  
  #Set Max model
  H1<-lm(H.abund.c$abund.change~.,data=sp.H.traits.c)
  vif(H1)

  # Min Root Location and Max Root Location highly correlated. 
  # LDMC-LmA correlated.
  H1<-update(H1,~.-Max.Root.Loca)
  H0<-lm(H.abund.c$abund.change~1,data=sp.H.traits.c)
  
  step.lm.H<-step(H0,scope=list(lower=H0,upper=H1),direction='both')
  
  summary(step.lm.H)
  # lm(formula = H.abund.c$abund.change ~ myc.frac + Min.Root.Loca, 
  # data = sp.H.traits.c)
  # Residual standard error: 1.904 on 32 degrees of freedom
  # Multiple R-squared:  0.305,	Adjusted R-squared:  0.2616 
  # F-statistic: 7.022 on 2 and 32 DF,  p-value: 0.002961
  AIC(step.lm.H)
  # 149.2829
  
  plot(step.lm.H)
  # Variance in residuals slightly decreases with mean - heteroscedastic.
  # QQ-plot -  error not normal - right skewed - bad
  # Leverage - CASC has a lot of influence on model - try model without that datapoint. 
  
  summary(lm(
    H.abund.c[rownames(H.abund.c)!='CASC','abund.change']~
    myc.frac+Min.Root.Loca, data=sp.H.traits.c[rownames(sp.H.traits.c)!='CASC',]
                          ))
  #               Estimate Std. Error t value Pr(>|t|)  
  # (Intercept)    1.05313    0.59751   1.763   0.0878 .
  # myc.frac      -0.57460    0.24390  -2.356   0.0250 *
  # Min.Root.Loca  0.09564    0.22023   0.434   0.6671 
  
  # Multiple R-squared:  0.1677,	Adjusted R-squared:  0.114 
  # F-statistic: 3.123 on 2 and 31 DF,  p-value: 0.05812
  
  # Min.Root.Loca was significant only because of CASC


  # 4.2 glm.nb - full model, no interactions ####
  
  H1<-glm.nb(H.abund.c$abund.change~.-Max.Root.Loca,data=sp.H.traits.c)
  H0<-glm.nb(H.abund.c$abund.change~1,data=sp.H.traits.c)
  
  step.glm.H<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=T)
  summary(step.glm.H)
  
  # glm.nb(formula = H.abund.c$abund.change ~ myc.frac + Min.Root.Loca, 
  # data = sp.H.traits.c, init.theta = 10.39193325, link = log)
  
  # Coefficients:
  #               Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)    -0.4679     0.4834  -0.968   0.3331    
  # myc.frac       -0.6465     0.1478  -4.375 1.21e-05 ***
  # Min.Root.Loca   0.2721     0.1625   1.674   0.0941 .
  
  #Null deviance: 56.766  on 34  degrees of freedom
  #residual deviance: 34.821  on 32  degrees of freedom
  #AIC: 110.26
  
  plot(step.glm.H)
  # Variance in residuals slightly decrease with mean - heteroscedastic?
  # QQ-plot -  error wavy - what does it mean?
  # Leverage - CASC has a lot of influence on model - try model without that datapoint.
  
  summary(glm.nb(
    H.abund.c[rownames(H.abund.c)!='CASC','abund.change']~
      myc.frac+Min.Root.Loca, data=sp.H.traits.c[rownames(sp.H.traits.c)!='CASC',]
  ))
  
  #               Estimate Std. Error z value Pr(>|z|)  
  #(Intercept)   -0.008488   0.475910  -0.018   0.9858  
  #myc.frac      -0.402754   0.167047  -2.411   0.0159 *
  #Min.Root.Loca  0.077373   0.171825   0.450   0.6525
  
  #Null deviance: 38.086  on 33  degrees of freedom
  #Residual deviance: 31.668  on 31  degrees of freedom
  #AIC: 98.185
  
  # Min.Root.Loca was significant only because of CASC
  
  # Testing whether we do need an overdispersion parameter (need for negative binomial distribution)
  # by comparing with a poisson regression, which does not assume overdispersion
  

  # 4.3 lm - powerTransform response variable. full model, no interactions ####
  
  H1<-lm(H.abund.c$pt.abund.change~.-Max.Root.Loca,data=sp.H.traits.c)
  H0<-lm(H.abund.c$pt.abund.change~1,data=sp.H.traits.c)
  
  step.lm.pt.H<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=T)
  summary(step.lm.pt.H)
  # lm(formula = H.abund.c$pt.abund.change ~ myc.frac + LDMC + Ht.veg, 
  #    data = sp.H.traits.c)
  
  # Coefficients:
  # Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)  1.01530    0.02548  39.843   <2e-16 ***
  # myc.frac    -0.04922    0.02478  -1.986   0.0559 .  
  # LDMC         0.05856    0.02814   2.081   0.0457 *  
  # Ht.veg       0.06017    0.04205   1.431   0.1625 
   
  # Residual standard error: 0.1465 on 31 degrees of freedom
  # Multiple R-squared:  0.2645,	Adjusted R-squared:  0.1933 
  # F-statistic: 3.716 on 3 and 31 DF,  p-value: 0.02158
  
  plot(step.lm.pt.H)
  # Variance in residuals homoscedastic - good
  # QQ-plot - more linear
  # Leverage - no datapoint with too much leverage
  
  AIC (step.lm.pt.H)
  #-29.35463
  
  step.lm.pt.H.2<-update(step.lm.pt.H,~.-Ht.veg)
  anova(step.lm.pt.H.2,step.lm.pt.H)
    #   Res.Df     RSS Df Sum of Sq      F Pr(>F)
    # 1     32 0.70963                           
    # 2     31 0.66568  1  0.043955 2.0469 0.1625
  
    #not significantly different, so keep the more parsimoneous model. 
  RsquareAdj(step.lm.pt.H.2) # 0.17
  
  plot(step.lm.pt.H.2)
  
  # 4.4 lm - log response variable. full model, no interactions ####
  H1<-lm(H.abund.c$log.abund.change~.-Max.Root.Loca,data=sp.H.traits.c)
  H0<-lm(H.abund.c$log.abund.change~1,data=sp.H.traits.c)
  
  step.lm.log.H<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=T)
  summary(step.lm.log.H)
  # lm(formula = H.abund.c$log.abund.change ~ LDMC + Ht.veg + myc.frac, 
  # data = sp.H.traits.c)
  # Coefficients:
  #             Estimate Std. Error t value Pr(>|t|)  
  # (Intercept)  0.02905    0.18228   0.159   0.8744  
  # LDMC         0.45425    0.20125   2.257   0.0312 *
  # Ht.veg       0.50950    0.30081   1.694   0.1003  
  # myc.frac    -0.26019    0.17725  -1.468   0.1522
  
  # Residual standard error: 1.048 on 31 degrees of freedom
  # Multiple R-squared:  0.2511,	Adjusted R-squared:  0.1786 
  # F-statistic: 3.465 on 3 and 31 DF,  p-value: 0.02796
  
  plot(step.lm.log.H)
  
  # Variance in residuals homoscedastic - good
  # QQ-plot - more linear - good
  # Leverage - OXMO has a lot of influence on model - try model without that datapoint.
  
  summary(lm(
    H.abund.c[rownames(H.abund.c)!='OXMO','log.abund.change']~
      LDMC+Ht.veg+myc.frac, data=sp.H.traits.c[rownames(sp.H.traits.c)!='OXMO',]
  ))
  # LDMC & vegetative height were only significant because of OXMO - 

 
  # 4.5 - exlore non-linearities with GAM & trait-trait interactions with tree ####
  
  library(mgcv)
  par(mfrow=c(2,2))
  m<-gam(H.abund.c$pt.abund.change~s(sp.H.traits.c$myc.frac)+s(sp.H.traits.c$Min.Root.Loca)+
        s(sp.H.traits.c$LDMC),data=sp.H.traits.c)
  plot(m)
  #myc.frac appears non-linear
  #Min.Root.Loca appears non-linear
  m<-gam(H.abund.c$pt.abund.change~s(sp.H.traits.c$Ht.veg)+s(sp.H.traits.c$Lamina.thck)+
  s(sp.H.traits.c$LMA),data=sp.H.traits.c)
  plot(m)
  #LMA appears correlated
  m<-gam(H.abund.c$pt.abund.change~s(sp.H.traits.c$Leaf.Area),data=sp.H.traits.c)
  plot(m)
  
  t<-tree(lm(H.abund.c$pt.abund.change~.-Max.Root.Loca,data=sp.H.traits.c))
  plot(t)
  text(t)
  # myc.frac. and LDMC appear non-linear
  
  # 4.6 model with interactions ####
  m<-lm(H.abund.c$pt.abund.change~(Min.Root.Loca+Ht.veg+Min.Root.Loca+Lamina.thck+LMA+
                                      LDMC+Leaf.Area+myc.frac)^2,
        data=sp.H.traits.c)
  summary(m)
  #no significant traits, nor trait-triat interactions
  
  # curvature model
  m<-lm(H.abund.c$pt.abund.change~I(Min.Root.Loca^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
          I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
        data=sp.H.traits.c)
  summary(m)
  # myc.frac^2 is singificant. 
  
  m1<-lm(H.abund.c$pt.abund.change~myc.frac*LDMC+I(myc.frac^2)+I(LDMC^2),data=sp.H.traits.c)
  m0<-lm(H.abund.c$pt.abund.change~1,data=sp.H.traits.c)
  
  step.m<-step(m0,scope=list(lower=m0,upper=m1),direction='both',trace=T)
  summary(step.m)
  # lm(formula = H.abund.c$pt.abund.change ~ I(myc.frac^2) + LDMC, 
  # data = sp.H.traits.c)
  
  # Coefficients:
  #          Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)    0.95222    0.02904  32.785  < 2e-16 ***
  # I(myc.frac^2)  0.05369    0.01611   3.332  0.00219 ** 
  # LDMC           0.04264    0.02603   1.638  0.11117
  
  # Residual standard error: 0.1356 on 32 degrees of freedom
  # Multiple R-squared:  0.3498,	Adjusted R-squared:  0.3091 
  # F-statistic: 8.607 on 2 and 32 DF,  p-value: 0.001021
  m2<-lm(H.abund.c$pt.abund.change ~ I(myc.frac^2),data = sp.H.traits.c)
  
  anova(step.m,m2)
  #   Res.Df     RSS Df Sum of Sq      F Pr(>F)
  # 1     32 0.58850                           
  # 2     33 0.63785 -1 -0.049355 2.6837 0.1112
  
  # Adding LDMC does not significantly improve the model, so drop it. 
  
  plot(m2)
  # Variance in residuals homoscedastic - good
  # QQ-plot - more linear, with slight tails
  # Leverage - CASC borderline of excessive leverage.
  
  RsquareAdj(m2)
  # 0.27
