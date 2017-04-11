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
#================================================================================

  # 0A - Herbivory layer
    H.abund<-abund[which(rownames(abund)%in%rownames(sp.H.traits)),c('presence.change','abund.change','ratio.labund')]
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
#================================================================================

  # 1A - Herbivory layer
    # 1A.1 - Abundance
  
      # Any pairwise interactions?
      plot.new()
      pairs(c(sp.H.traits,H.abund[,c('abund.change','ratio.labund')]),panel=panel.smooth)
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
      
      # Any important interactions on response variable?
      
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

# 2 - Abund~Traits, Individual variables ####
      #================================================================================
     
  # Diagnostic plots
      
  # 1) Residuals vs Fitted - detects non-linear relationships between x & y variables + heteroscedasticity (wedge)
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
    # Diagnostic plots - residuals corelated with fitted (non-linearity); resids right skewed;EPHE large weight
    
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
    
    ##=========================
    # Residuals are non-normally distributed & non-linearities with predictor variables. 
    
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


    # 4.0 test higher-level effects and trait-trait interactions - strong effect, but only driven by outliers####
    
    # Test Trait-Triat interactions

    m<-lm(H.abund.c$pt.abund.change~(Min.Root.Loca+Ht.veg+Min.Root.Loca+Lamina.thck+LMA+
                                       LDMC+Leaf.Area+myc.frac)^2,
          data=sp.H.traits.c)
    summary(m)
    
    #no significant traits effects, nor trait-trait interactions
    
    # Test higher-level effects 
    
    m<-lm(H.abund.c$pt.abund.change~.-Max.Root.Loca+I(Ht.veg^2)+I(Min.Root.Loca^2)+I(Lamina.thck^2)+
            I(LMA^2)+I(LDMC^2)+I(Leaf.Area^2)+I(myc.frac^2),
          data=sp.H.traits.c)
    summary(m)
    # myc.frac^2 is singificant. 
    summary(m<-update(m,~.-I(LMA^2)))
    summary(m<-update(m,~.-Leaf.Area))
    summary(m<-update(m,~.-Ht.veg))
    summary(m<-update(m,~.-LDMC))
    summary(m<-update(m,~.-I(Ht.veg^2)))
    summary(m<-update(m,~.-I(Lamina.thck^2)))
    summary(m<-update(m,~.-I(LDMC^2)))
    summary(m<-update(m,~.-I(Min.Root.Loca^2)))
    summary(m<-update(m,~.-Lamina.thck))
    summary(m<-update(m,~.-Min.Root.Loca))
    summary(m<-update(m,~.-myc.frac))
    summary(m<-update(m,~.-I(Leaf.Area^2)))
    summary(m<-update(m,~.-LMA))
     
      # Call:
      #   lm(formula = H.abund.c$pt.abund.change ~ I(myc.frac^2), data = sp.H.traits.c)
     
      # Coefficients:
      #                Estimate Std. Error t value Pr(>|t|)    
      # (Intercept)    0.94001    0.02878  32.663  < 2e-16 ***
      # I(myc.frac^2)  0.05977    0.01607   3.718 0.000743 ***
      #   ---
            # 
      # Residual standard error: 0.139 on 33 degrees of freedom
      # Multiple R-squared:  0.2952,	Adjusted R-squared:  0.2739 
      # F-statistic: 13.82 on 1 and 33 DF,  p-value: 0.0007432
    
    plot(m)
    
    # plot regression
    plot(H.abund.c$pt.abund.change~myc.frac,data=sp.H.traits.c,type='n')
    text(sp.H.traits.c$myc.frac,H.abund.c$pt.abund.change,label=rownames(sp.H.traits.c),cex=0.7)
    
    # add curve
    points(x=sp.H.traits.c$myc.frac, y=fitted(m), col='red', pch=20)
    lines(x=sort(sp.H.traits.c$myc.frac),y=fitted(m)[sort(sp.H.traits.c$myc.frac)],col='red',type='b' ) 
    # doesn't work for some reason
    
    # plot without outliers
    mno<-lm(
      H.abund.c[!rownames(H.abund.c)%in%c('CASC','LYAN','OXMO'),'pt.abund.change']~
        I(myc.frac^2), data=sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('CASC','LYAN','OXMO'),]
    )
    plot(H.abund.c[!rownames(H.abund.c)%in%c('CASC','LYAN','OXMO'),'pt.abund.change']~
           myc.frac,data=sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('CASC','LYAN','OXMO'),],
                                       type='n')
    text(x=sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('CASC','LYAN','OXMO'),'myc.frac'],
         y= H.abund.c[!rownames(H.abund.c)%in%c('CASC','LYAN','OXMO'),'pt.abund.change'],
         label=rownames(sp.H.traits.c)[!rownames(H.abund.c)%in%c('CASC','LYAN','OXMO')],cex=0.7)
    points(x=sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('CASC','LYAN','OXMO'),'myc.frac'],
           y=fitted(mno),
           col='red', pch=20)
    
    RsquareAdj(mno) #-0.02
    summary(mno)
    # Coefficients:
    #               Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)    0.98322    0.02614  37.616   <2e-16 ***
    # I(myc.frac^2)  0.01313    0.02428   0.541    0.593
    
    # Residual standard error: 0.1049 on 30 degrees of freedom
    # Multiple R-squared:  0.009657,	Adjusted R-squared:  -0.02335 
    # F-statistic: 0.2925 on 1 and 30 DF,  p-value: 0.5926
    
    # model only fits the outliers !!!
    
    # Test Trait-Trait interactions
    
    m1<-lm(H.abund.c$pt.abund.change~myc.frac*LDMC+I(myc.frac^2),data=sp.H.traits.c)
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
    
  # 4.1 lm - doesn't meet assumptions & only significant bc of outliers ####

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
  H1<-update(H1,~.-Max.Root.Loca+I(myc.frac^2))
  H0<-lm(H.abund.c$abund.change~1,data=sp.H.traits.c)
  
  step.lm.H<-step(H0,scope=list(lower=H0,upper=H1),direction='both')
  
  summary(step.lm.H)
  # lm(formula = H.abund.c$abund.change ~ I(myc.frac^2) + myc.frac + 
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
  
  plot(step.lm.H)
  # Variance in residuals decreases with mean - homoscedastic, but an x-y relationship remain unexplained by model
  # QQ-plot -  error not normal - s-shaped - bad
  # Leverage - CASC & LYAN have a lot of influence on model - try model without that datapoint. 
  
  H.abund.c.no.outliers<-H.abund.c[!rownames(H.abund.c)%in%c('CASC','LYAN'),]
  dim(H.abund.c.no.outliers) #33 4
  sp.H.traits.c.no.outliers<-sp.H.traits.c[!rownames(sp.H.traits.c)%in%c('CASC','LYAN'),]
  dim(sp.H.traits.c.no.outliers) #33 8
  
  summary(lm(
    H.abund.c.no.outliers$abund.change~I(myc.frac^2)+myc.frac,data=sp.H.traits.c.no.outliers
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
  
  # model was significant only because of CASC and LYAN
  
  # 4.2 lm - log-transform response variable - meets assumptions ####
  
  H1<-lm(H.abund.c$log.abund.change~.+I(myc.frac^2)-Max.Root.Loca,data=sp.H.traits.c)
  H0<-lm(H.abund.c$log.abund.change~1,data=sp.H.traits.c)
  
  step.lm.log.H<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=F)
  summary(step.lm.log.H)
  # lm(formula = H.abund.c$log.abund.change ~ LDMC + Ht.veg + myc.frac, 
  # data = sp.H.traits.c)
  # Coefficients:
  #               Estimate Std. Error t value Pr(>|t|)  
  # (Intercept)    -0.3076     0.2124  -1.448   0.1577  
  # I(myc.frac^2)   0.3131     0.1162   2.696   0.0113 *
  # LDMC            0.3855     0.1898   2.032   0.0508 .
  # Ht.veg          0.4558     0.2805   1.625   0.1143  
  # ---
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  # Residual standard error: 0.9757 on 31 degrees of freedom
  # Multiple R-squared:  0.3511,	Adjusted R-squared:  0.2883 
  # F-statistic: 5.592 on 3 and 31 DF,  p-value: 0.003481
  
  
  anova(step.lm.log.H,update(step.lm.log.H,~.-Ht.veg))
  #   Res.Df    RSS Df Sum of Sq      F Pr(>F)
  # 1     31 29.511                           
  # 2     32 32.024 -1    -2.513 2.6398 0.1143
  
  step.lm.log.H<-update(step.lm.log.H,~.-Ht.veg)
  
  plot(step.lm.log.H)
  
  # Variance in residuals homoscedastic - good
  # QQ-plot - more linear - good
  # Leverage - OXMO has a lot of influence on model - try model without that datapoint.
  
  summary(lm(
    H.abund.c[rownames(H.abund.c)!='OXMO','log.abund.change']~
      LDMC+I(myc.frac^2), data=sp.H.traits.c[rownames(sp.H.traits.c)!='OXMO',]
  ))
  # LDMC was only significant because of OXMO - 
  
  # 4.3 lm - powerTransform response variable - meets assumptions ####
  
  H1<-lm(H.abund.c$pt.abund.change~.+I(myc.frac^2)-Max.Root.Loca,data=sp.H.traits.c)
  H0<-lm(H.abund.c$pt.abund.change~1,data=sp.H.traits.c)
  
  step.lm.pt.H<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=T)
  summary(step.lm.pt.H)
  # lm(formula = H.abund.c$pt.abund.change ~ myc.frac + LDMC + Ht.veg, 
  #    data = sp.H.traits.c)
  
  # Coefficients:
  # Estimate Std. Error t value Pr(>|t|)    
  #              Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)    0.95222    0.02904  32.785  < 2e-16 ***
  #  I(myc.frac^2)  0.05369    0.01611   3.332  0.00219 ** 
  #  LDMC           0.04264    0.02603   1.638  0.11117
  
  # Residual standard error: 0.1356 on 32 degrees of freedom
  # Multiple R-squared:  0.3498,	Adjusted R-squared:  0.30 
  # F-statistic: 8.607 on 2 and 32 DF,  p-value: 0.001021
  
  anova(step.lm.pt.H,update(step.lm.pt.H,~.-LDMC))
  # Res.Df     RSS Df Sum of Sq      F Pr(>F)
  # 1     32 0.58850                           
  # 2     33 0.63785 -1 -0.049355 2.6837 0.1112
  
  # Inclusion of LDMC not warranted. Remove
  step.lm.pt.H<-update(step.lm.pt.H,~.-LDMC)
  
  plot(step.lm.pt.H)
  # Variance in residuals homoscedastic - good
  # QQ-plot - quite linear !
  # Leverage - CASC is borderline  with too much leverage
  
  AIC (step.lm.pt.H)
  #-34.85
  
  # 4.4 glm.nb - doesn't meet assumption ####
  
  H1<-glm.nb(H.abund.c$abund.change~.+I(myc.frac^2)-Max.Root.Loca,data=sp.H.traits.c)
  H0<-glm.nb(H.abund.c$abund.change~1,data=sp.H.traits.c)
  
  step.glm.H<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=F)
  summary(step.glm.H)
  
  # glm.nb(formula = H.abund.c$abund.change ~ I(myc.frac^2), data = sp.H.traits.c,
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
  
  plot(step.glm.H)
  # Variance in residuals is homoscedastic
  # QQ-plot - distribution still is still banana-shaped :-(
  # Leverage - CASC & LYAN have a lot of influence on model - try model without that datapoint.
  
  summary(glm.nb(
    H.abund.c.no.outliers$abund.change~I(myc.frac^2),data=sp.H.traits.c.no.outliers
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
  
  m3 <- glm(H.abund.c$abund.change ~ I(myc.frac^2), family = "poisson", data = sp.H.traits.c)
  pchisq(2 * (logLik(step.glm.H) - logLik(m3)), df = 1, lower.tail = FALSE)
  # Doesn't work bc poisson distribution doesn't work with my data
  
  # 4.5 - glm - lognormal link function?
  
  # 5 - exlore non-linearities with GAM & trait-trait interactions with tree ####
  
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
  

  
  # Observe long-term patterns: 
  plot(abund[abund$Layer=='H','avg.abundance.2012']~abund[abund$Layer=='H','avg.abundance.1970'])
  text(abund[abund$Layer=='H','avg.abundance.1970'],abund[abund$Layer=='H','avg.abundance.2012'],
       labels=rownames(abund)[which(abund$Layer=='H')])
  
  #on log-scale, much better
  plot(log(abund[abund$Layer=='H','avg.abundance.2012'])~log(abund[abund$Layer=='H','avg.abundance.1970']))
  text(log(abund[abund$Layer=='H','avg.abundance.1970']),log(abund[abund$Layer=='H','avg.abundance.2012']),
       labels=rownames(abund)[which(abund$Layer=='H')])
  abline(0,1)
  
  # for Canopy layer: 
  plot(abund[abund$Layer=='C','avg.abundance.2012']~abund[abund$Layer=='C','avg.abundance.1970'])
  text(abund[abund$Layer=='C','avg.abundance.1970'],abund[abund$Layer=='C','avg.abundance.2012'],
       labels=rownames(abund)[which(abund$Layer=='C')])
  
  plot(log(abund[abund$Layer=='C','avg.abundance.2012'])~log(abund[abund$Layer=='C','avg.abundance.1970']))
  text(log(abund[abund$Layer=='C','avg.abundance.1970']),log(abund[abund$Layer=='C','avg.abundance.2012']),
       labels=rownames(abund)[which(abund$Layer=='C')])
  abline(0,1)
