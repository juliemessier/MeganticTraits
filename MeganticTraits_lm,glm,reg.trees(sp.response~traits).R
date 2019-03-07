# <<TABLE OF CONTENTS>> ####
# Make trait and abundance dataframes correspond based on shared species
#(A) Abundance vs 8 Traits
#  A1 - Explore Trait Interactions
#  A2 - Explore non-linearity 
#  A3 - Model selection 
#     A3.1 - lm 
#     A3.2 - glm
#  A4 Rerun Tree & GLM models without high-leverage points
#(B) Occurence vs 8 Traits
#  B1 - Explore Trait Interactions
#  B2 - Explore non-linearity
#  B3 - Model selection 
#     B3.2 - glm
#(C) Elevation vs 8 Traits
#  C1 - Explore Trait Interactions
#  C2 - Explore non-linearity
#  C3 - Model selection
#     C3.1 - lm full

#====================================================================================#


#
#<<WORKSPACES>>
wrk.dir<-("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Traits/Workspaces/") # Workspaces
data.dir<-(("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Trait/Data/")) # data
res.dir<-("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Traits/Results/")  # Results
grp.dir<-("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Trait/Graphs/")   # Graphs
fct.dir<-("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Trait/Functions/") # Functions

#<<LIBRARIES>>
library(car) # for powerTransform
library(vegan) # for decostand (standardizing data)
library(tree) # for tree()
library(MASS) # for glm.nb()  
library(mgcv) # for gam()

#<<LOAD FILES>>

load(file=paste0(wrk.dir,'Species-Level.Traits.Understory-Layer.with.MycFrac.and.SeedSize.and.PCs.RData')) 
   # Object = U.traits4.sp
   # Species-level traits
   # Understory layer only
   # Includes mycorhizal fraction, seed size, log-transformed traits and principal components of PCA performed on traits
   dim(U.traits4.sp)
   # dim = 53 30

load(file=paste0(wrk.dir,'All.Response.variables.Understory.Layer.Inf.removed.RData')) 
   # object = U.Ys.c
   # Species temporal responses
   # Understory layer only
   # Includes abundance ratio, occurence ratio, elevation shift, ratio of log(abundances) and ratio of log(occurences)
   dim(U.Ys.c)
   # dim = 47 14

#===============================================================================#
# Make trait and abundance dataframes correspond based on shared species ####
#===============================================================================#

    # subset abundance data to include only the species with trait data 
    sp.response<-U.Ys.c[rownames(U.Ys.c)%in%rownames(U.traits4.sp),]
    dim(sp.response)
      # 47 14
    
    # subset trait data to include only the species with abundance data
    Xs<-U.traits4.sp[rownames(U.traits4.sp)%in%rownames(U.Ys.c),]
    dim(Xs)
      # 47 30
    
    all(rownames(sp.response)==rownames(Xs)) 
      # TRUE
    
    # Create long and short lists of traits 
    Trait.Names.8t<-c("Log.Ht.veg","Max.Root.Loca","Log.Lamina.thck","LDMC","Log.Leaf.Area",  
                        "Leaf.Mass.Frac","SRL","Myc.Frac")
      # 8 traits
    length(which(complete.cases(Xs[,Trait.Names.8t])))
      # 36 species (instead of 35 previously)
    
    Trait.Names.5t<-c("Log.Ht.veg","Log.Lamina.thck","LDMC","Log.Leaf.Area","Leaf.Mass.Frac")
      # 5 traits
    length(which(complete.cases(Xs[,Trait.Names.5t])))
      # 46 species (instead of 44 previously)
    
#============================#
#(A) Abundance vs 8 Traits####
#============================#
    
    # A1 - Explore Trait Interactions #### 
    #==============================#
   
    # A1.1 - Tree Model
    par(mfrow=c(1,1))
    tree.model<-tree(sp.response$abundance.ratio~.,data=Xs[,Trait.Names.8t])
    plot(tree.model); text(tree.model);title('Regression Tree \n Abundance.ratio vs 8 traits')
    # Myc.frac is most important variable, 
    # for those with high myc.frac, SRL matters
    
    tree.model 
    # node), split, n, deviance, yval
    #      * denotes terminal node
    # 
    # 1) root 36 178.1000 1.5740  
    #   2) Myc.Frac < 0.58 5 104.0000 4.3820 *
    #   3) Myc.Frac > 0.58 31  28.3700 1.1220  
    #     6) SRL < 10.7122 12  13.1100 1.5160  
    #      12) Max.Root.Loca < 2.05 7   0.9664 1.0740 *
    #      13) Max.Root.Loca > 2.05 5   8.8620 2.1360 *
    #     7) SRL > 10.7122 19  12.2100 0.8723 *
    
    #  i.e. null deviance = 178.1
 
    summary(tree.model)
   # Regression tree:
   # tree(formula = sp.response$abundance.ratio ~ ., data = Xs[, Trait.Names.8t])
   # Variables actually used in tree construction:
   # [1] "Myc.Frac"      "SRL"           "Max.Root.Loca"
   # Number of terminal nodes:  4 
   # Residual mean deviance:  3.938 = 126 / 32 
   # Distribution of residuals:
   #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # -3.4060 -0.6045 -0.1723  0.0000  0.2187  8.2850 
    
    # R2 = 1 - SSres/SStotal
    # Pseudo-R2
    # 1-(deviance(model)/Total(Null) Deviance)
    
    1-(deviance(tree.model)/178.1) 
      # 0.2924518
    plot(cv.tree(tree.model)) 
    # c.v.  splits the data into 'training' set for model fitting and a 
    # validation set to evaluate goodness of fit 
    # look for how many splits produces the minimum deviance - run multiple times
    
    # lowest deviations at 2-4 branches
    
          # See if regression tree retains Myc.Frac-log(Myc.Frac), like lm and glms do.
          Xs$Log.Myc.Frac<-log(Xs$Myc.Frac)
          tree.model2<-tree(sp.response$abundance.ratio~.,data=Xs[,c(Trait.Names.8t,'Log.Myc.Frac')])
          
          plot(tree.model2)
          text(tree.model2)
          title('Regression Tree2 \n Abundance.ratio vs 8 traits & log(Myc.Frac) \n n=36' )
          
          tree.model2
          # node), split, n, deviance, yval
          #   * denotes terminal node
          #
          # 1) root 36 178.1000 1.5740  
          #   2) Myc.Frac < 0.58 5 104.0000 4.3820 *
          #   3) Myc.Frac > 0.58 31  28.3700 1.1220  
          #     6) SRL < 10.7122 12  13.1100 1.5160  
          #      12) Max.Root.Loca < 2.05 7   0.9664 1.0740 *
          #      13) Max.Root.Loca > 2.05 5   8.8620 2.1360 *
          #     7) SRL > 10.7122 19  12.2100 0.8723 *
          
          summary(tree.model2)
         # Regression tree:
         # tree(formula = sp.response$abundance.ratio ~ ., data = Xs[, c(Trait.Names.8t, 
         #     "Log.Myc.Frac")])
         # Variables actually used in tree construction:
         # [1] "Myc.Frac"      "SRL"           "Max.Root.Loca"
         # Number of terminal nodes:  4 
         # Residual mean deviance:  3.938 = 126 / 32 
         # Distribution of residuals:
         #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
         # -3.4060 -0.6045 -0.1723  0.0000  0.2187  8.2850 
          
               # No, only myc.frac retained
          
          1-(deviance(tree.model2)/178.1) 
            # 0.29
          
          # Check that the model did consider including log(myc.frac)
          tree.model2$terms[[3]]
             # Log.Ht.veg + Max.Root.Loca + Log.Lamina.thck + LDMC + Log.Leaf.Area + 
             # Leaf.Mass.Frac + SRL + Myc.Frac + Log.Myc.Frac
          
       #=====================================================================================#
       # Summary of regression tree
       #
       # Regression Tree retains Myc.Frac, Max.Root.Loca, SRL and their interactions - R2 = 0.29
       # Regression Tree does not retain log(Myc.Frac) like glm does. This is because it does 
       # need to fit the exact shape of the regression. Instead, it pools species into bins
       #=====================================================================================#
    
    # A1.2 - Test trait-trait interactions 
    ----
       
    # need to check in subgroups because sample size too small to test all interactions
    # simultaneously
    summary(lm(sp.response$abundance.ratio~(Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC)^2,
               data=merge(sp.response,Xs,by="row.names",all=T)))
 
      # No significant interactions
    
    summary(lm(sp.response$abundance.ratio~(Log.Leaf.Area+Leaf.Mass.Frac+SRL+Myc.Frac)^2,
               data=merge(sp.response,Xs,by="row.names",all=T)))
    
    # SRL:Myc.Frac significant @p=0.01
    # Log.Leaf.Area:SRL significant @p=0.01
    
    summary(lm(sp.response$abundance.ratio~(Log.Ht.veg+Max.Root.Loca+Log.Leaf.Area+Leaf.Mass.Frac)^2,
               data=merge(sp.response,Xs,by="row.names",all=T)))
    
    # Log.Ht.veg:Max.Root.Loca marginally significant @p=0.0584
    # Max.Root.Loca:Leaf.Mass.Frac marginally significant @p=0.0977
    
    summary(lm(sp.response$abundance.ratio~(Log.Ht.veg+Max.Root.Loca+SRL+Myc.Frac)^2,
               data=merge(sp.response,Xs,by="row.names",all=T)))
    # SRL:Myc.Frac significant @p=0.001
    # Log.Ht.veg:SRL marginally significant @p=0.077
    # Max.Root.Loca:SRL marginally significant @p=0.059

    
    summary(lm(sp.response$abundance.ratio~(Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac)^2,
               data=merge(sp.response,Xs,by="row.names",all=T)))
    
    # No significant interactions
    
    summary(lm(sp.response$abundance.ratio~(Log.Lamina.thck+LDMC+SRL+Myc.Frac)^2,
               data=merge(sp.response,Xs,by="row.names",all=T)))
    
    # SRL:Myc.Frac significant @p=0.0004
    # LDMC:SRL significant @p=0.02
    
    #===============================================================#
    # SUMMARY of trait interactions #
    # significant trait-trait interactions to include in full model: 
    #    SRL:Myc.Frac  
    #    Log.Leaf.Area:SRL
    #    Log.Ht.veg:Max.Root.Loca
    #    Max.Root.Loca:Leaf.Mass.Frac
    #    Log.Ht.veg:SRL
    #    Max.Root.Loca:SRL 
    #    LDMC:SRL
    #===============================================================#
    
    # A2 - Explore non-linearity ####
    #=============================#
    
    # A2.1 GAMs
    ---
    par(mfrow=c(2,2),mar=c(4,4,3,3))
    plot(gam(sp.response$abundance.ratio~s(Log.Ht.veg)+
               +s(SRL)+s(Myc.Frac),data=Xs))
      #myc.frac is non-linear
    
    plot(gam(sp.response$abundance.ratio~s(Max.Root.Loca)+s(Log.Lamina.thck)+
             s(LDMC),data=Xs))
      # No relationships
    
    plot(gam(sp.response$abundance.ratio~s(Log.Leaf.Area)+s(Leaf.Mass.Frac),data=Xs))
      # No relationships
    
    
    # A2.2 Test polynomial relationships 
    ----
    summary(lm(sp.response$abundance.ratio~I(Log.Ht.veg^2)+I(Max.Root.Loca^2)+I(Log.Lamina.thck^2)+
                  I(LDMC^2)+I(Log.Leaf.Area^2)+I(Leaf.Mass.Frac^2)+I(SRL^2)+I(Myc.Frac^2),
               data=merge(sp.response,Xs,by="row.names",all=T)))
      # myc.frac^2 significant @ p=0.01
    
    summary(lm(sp.response$abundance.ratio~log(Myc.Frac),
               data=merge(sp.response,Xs,by="row.names",all=T)))
      # log(myc.frac) significant @ p=0.0004
    
   
# A3 - Model selection ####
#=========================#

    dat.8t<-merge(sp.response[,c("LatinName","abundance.ratio","occurence.ratio",'ElevDif')],
                    Xs[which(complete.cases(Xs[,Trait.Names.8t])),Trait.Names.8t],
                    by="row.names",all.y=T)
    # Cleanup
    dim(dat.8t)
    # 36 13
    rownames(dat.8t)<-dat.8t$Row.names
    dat.8t$Row.names<-NULL
    dim(dat.8t)
    # 36 12
    
    dat.5t<-merge(sp.response[,c("LatinName","abundance.ratio","occurence.ratio",'ElevDif')],
                     Xs[which(complete.cases(Xs[,Trait.Names.5t])),Trait.Names.5t],
                     by="row.names",all.y=T)
    # Cleanup
    dim(dat.5t)
    # 46 10
    rownames(dat.5t)<-dat.5t$Row.names
    dat.5t$Row.names<-NULL
    dim(dat.5t)
    # 46 9
    
   # A3.1 - lm ####
   #===================================#
   
   # A3.1.1 lm full #### 
   #=========================================#
   # Includes all 1st order terms 
   # + significant trait interactions - four 2nd order, one 3rd order
   # + significant non-linear variables.
    
  H1<-lm(abundance.ratio~
            Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac+SRL+Myc.Frac+
            log(Myc.Frac)+
            SRL:Myc.Frac + Log.Leaf.Area:SRL + Max.Root.Loca:SRL + Log.Ht.veg:SRL+
            Myc.Frac:SRL:Log.Ht.veg,
         data=dat.8t)
    # use 3-way interactions with Myc.Frac
  H0<-lm(abundance.ratio~1,
         data=dat.8t)
  
  best.lm.abund.8t<-step(H1,scope=list(lower=H0,upper=H1),direction='backward',trace = T)
  
  summary(best.lm.abund.8t)
   # Call:
   # lm(formula = abundance.ratio ~ Log.Ht.veg + SRL + Myc.Frac + 
   #     log(Myc.Frac) + SRL:Myc.Frac + Log.Ht.veg:SRL + Log.Ht.veg:SRL:Myc.Frac, 
   #     data = dat.8t)
   # 
   # Residuals:
   #     Min      1Q  Median      3Q     Max 
   # -1.8856 -0.6668 -0.0396  0.3563  3.4124 
   # 
   # Coefficients:
   #                         Estimate Std. Error t value Pr(>|t|)    
   # (Intercept)             -26.0082     8.3353  -3.120 0.004163 ** 
   # Log.Ht.veg               -1.5515     0.6471  -2.398 0.023407 *  
   # SRL                      -1.6238     0.3844  -4.225 0.000230 ***
   # Myc.Frac                 33.3629     8.7854   3.798 0.000721 ***
   # log(Myc.Frac)           -23.2078     5.8643  -3.957 0.000471 ***
   # SRL:Myc.Frac              1.7504     0.4750   3.685 0.000971 ***
   # Log.Ht.veg:SRL            0.6201     0.1115   5.560 6.02e-06 ***
   # Log.Ht.veg:SRL:Myc.Frac  -0.6856     0.1379  -4.972 2.99e-05 ***
   # ---
   # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
   # 
   # Residual standard error: 1.121 on 28 degrees of freedom
   # Multiple R-squared:  0.8026,	Adjusted R-squared:  0.7532 
   # F-statistic: 16.26 on 7 and 28 DF,  p-value: 2.43e-08
  
  AIC(best.lm.abund.8t) # 119.319
  
  plot(best.lm.abund.8t) # Abnormal. points clumped to left w strong negative slope 
                        # datapoints 30,8 and 11 are poorly fitted (high residuals vs fitted) -- osmorhiza claytonii, circea alpina & Cornus canadensis
                        # datapoints 30 and 7 have high leverage (>1) -- Carex scabrata and osmorhiza claytonii 
  
  vif(best.lm.abund.8t)
   #               Log.Ht.veg                     SRL                Myc.Frac 
   #                3.953579              456.477998               48.486129 
   #           log(Myc.Frac)            SRL:Myc.Frac          Log.Ht.veg:SRL 
   #               50.038599              434.131104              351.792030 
   # Log.Ht.veg:SRL:Myc.Frac 
   #              318.856627
      # It is normal to have high VIFs for interaction terms. Don't worry about it. 
  
  # Plot residuals vs each variable (in the model)
  plot((resid(best.lm.abund.8t))~dat.8t$Log.Ht.veg)
  plot((resid(best.lm.abund.8t))~dat.8t$SRL)
  plot((resid(best.lm.abund.8t))~dat.8t$Myc.Frac)
  plot((resid(best.lm.abund.8t))~log(dat.8t$Myc.Frac))
  # Plot residuals vs each variable (NOT in the model)
  plot((resid(best.lm.abund.8t))~log(dat.8t$Max.Root.Loca))
  plot((resid(best.lm.abund.8t))~log(dat.8t$Log.Lamina.thck))
  plot((resid(best.lm.abund.8t))~log(dat.8t$LDMC))
  plot((resid(best.lm.abund.8t))~log(dat.8t$Log.Leaf.Area))
  plot((resid(best.lm.abund.8t))~log(dat.8t$Leaf.Mass.Frac))
  
 
  # Why does the model retain myc.frac - log (myc.frac)?
  # Look at shape of that relationship on a graph
  x<-seq(0,1,0.05)
  y<-x-log(x)
  plot(x,y,main='y~x-log(x)\n x=[0-1]')
  
  y<--log(x)
  plot(x,y,main='y~x-log(x)\n x=[0-1]')
  
  # Checkout what abundance ratio vs Myc.Frc look like graphically
   model<-lm(abundance.ratio~Myc.Frac+log(Myc.Frac),data=dat.8t)
   plot(dat.8t$Myc.Frac,dat.8t$abundance.ratio,type='n',main='Abund = -33 + 35 Myc.Frac - 28 log(Myc.Frac)')
   text(dat.8t$Myc.Frac,dat.8t$abundance.ratio,dat.8t$LatinName,cex=0.7)
   lines(sort(dat.8t$Myc.Frac), predict(model)[order(dat.8t$Myc.Frac)], 
        type='b', 
        col='red')
  
  model<-lm(abundance.ratio~log(Myc.Frac),data=dat.8t)
  plot(dat.8t$Myc.Frac,dat.8t$abundance.ratio,type='n',main='Abund = -0.1 - 5.4 log(Myc.Frac)')
  text(dat.8t$Myc.Frac,dat.8t$abundance.ratio,dat.8t$LatinName,cex=0.7)
  lines(sort(dat.8t$Myc.Frac), predict(model)[order(dat.8t$Myc.Frac)], 
        type='b', 
        col='red')

  
  
        # A3.1.2 lm full, with Myc.Frac only (no log(Myc.Frac)) #### 
        #=========================================#
        # Not both untransformed and log transformed Myc.Frac.
        
        H1<-lm(abundance.ratio~
                  Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac+SRL+Myc.Frac+
                  SRL:Myc.Frac + Log.Leaf.Area:SRL + Max.Root.Loca:SRL + Log.Ht.veg:SRL+
                  Myc.Frac:SRL:Log.Ht.veg,
               data=dat.8t)
        
        H0<-lm(abundance.ratio~1,
               data=dat.8t)
        
        best.lm.abund.8t<-step(H1,scope=list(lower=H0,upper=H1),direction='backward',trace = T)
        
        summary(best.lm.abund.8t)
            # Call:
            # lm(formula = abundance.ratio ~ Log.Ht.veg + Max.Root.Loca + Log.Lamina.thck + 
            #     SRL + Myc.Frac + SRL:Myc.Frac + Max.Root.Loca:SRL + Log.Ht.veg:SRL + 
            #     Log.Ht.veg:SRL:Myc.Frac, data = dat.8t)
            # 
            # Residuals:
            #     Min      1Q  Median      3Q     Max 
            # -2.0461 -0.7133 -0.1584  0.5898  2.9374 
            # 
            # Coefficients:
            #                         Estimate Std. Error t value Pr(>|t|)   
            # (Intercept)             10.21304    5.23861   1.950  0.06209 . 
            # Log.Ht.veg              -1.44851    0.81502  -1.777  0.08723 . 
            # Max.Root.Loca           -0.64328    0.57003  -1.129  0.26941   
            # Log.Lamina.thck         -0.99739    0.79312  -1.258  0.21973   
            # SRL                     -1.22919    0.56066  -2.192  0.03750 * 
            # Myc.Frac                 1.91902    3.50042   0.548  0.58821   
            # SRL:Myc.Frac             1.13993    0.72163   1.580  0.12627   
            # Max.Root.Loca:SRL        0.03644    0.02888   1.261  0.21835   
            # Log.Ht.veg:SRL           0.52788    0.16270   3.244  0.00323 **
            # Log.Ht.veg:SRL:Myc.Frac -0.55837    0.20261  -2.756  0.01055 * 
            # ---
            # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
            # 
            # Residual standard error: 1.391 on 26 degrees of freedom
            # Multiple R-squared:  0.7174,	Adjusted R-squared:  0.6196 
            # F-statistic: 7.334 on 9 and 26 DF,  p-value: 2.996e-05
        
        AIC(best.lm.abund.8t) 
         # 136.2296
        
        plot(best.lm.abund.8t) #### abnormal residual vs fitted. very strong slope. 8, 23 and 11 are outliers
                              ## non-normal residuals
                              ### 7 and 30 have large leverage/cooks distance
        
      
        # Plot residuals vs each variable (in the model)
        plot((resid(best.lm.abund.8t))~dat.8t$Log.Ht.veg)
        plot((resid(best.lm.abund.8t))~dat.8t$SRL)
        plot((resid(best.lm.abund.8t))~dat.8t$Myc.Frac)
      
        # Plot residuals vs each variable (NOT in the model)
        plot((resid(best.lm.abund.8t))~log(dat.8t$Max.Root.Loca))
        plot((resid(best.lm.abund.8t))~log(dat.8t$Log.Lamina.thck))
        plot((resid(best.lm.abund.8t))~log(dat.8t$LDMC))
        plot((resid(best.lm.abund.8t))~log(dat.8t$Log.Leaf.Area))
        plot((resid(best.lm.abund.8t))~log(dat.8t$Leaf.Mass.Frac))
  
        # A3.1.3 lm full, with log(Myc.Frac) only (no Myc.Frac) #### 
        #=============================================#
        H1<-lm(abundance.ratio~
                  Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac+SRL+log(Myc.Frac)+
                  SRL:log(Myc.Frac) + Log.Leaf.Area:SRL + Max.Root.Loca:SRL + Log.Ht.veg:SRL+
                  log(Myc.Frac):SRL:Log.Ht.veg,
               data=dat.8t)
        
        H0<-lm(abundance.ratio~1,
               data=dat.8t)
        
        best.lm.abund.8t<-step(H1,scope=list(lower=H0,upper=H1),direction='backward',trace = T)
        
        summary(best.lm.abund.8t)
            # lm(formula = abundance.ratio ~ Log.Ht.veg + SRL + log(Myc.Frac) + 
            #     SRL:log(Myc.Frac) + Log.Ht.veg:SRL + Log.Ht.veg:SRL:log(Myc.Frac), 
            #     data = dat.8t)
            # 
            # Residuals:
            #      Min       1Q   Median       3Q      Max 
            # -2.28408 -0.72965 -0.07581  0.31226  2.96131 
            # 
            # Coefficients:
            #                              Estimate Std. Error t value Pr(>|t|)    
            # (Intercept)                   4.69892    2.39361   1.963 0.059288 .  
            # Log.Ht.veg                   -1.10072    0.71700  -1.535 0.135581    
            # SRL                           0.03011    0.15507   0.194 0.847411    
            # log(Myc.Frac)                -0.44366    1.90943  -0.232 0.817895    
            # SRL:log(Myc.Frac)             0.95670    0.37723   2.536 0.016853 *  
            # Log.Ht.veg:SRL               -0.03588    0.04905  -0.732 0.470325    
            # Log.Ht.veg:SRL:log(Myc.Frac) -0.41110    0.10545  -3.899 0.000527 ***
            # ---
            # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
            # 
            # Residual standard error: 1.258 on 29 degrees of freedom
            # Multiple R-squared:  0.7421,	Adjusted R-squared:  0.6888 
            # F-statistic: 13.91 on 6 and 29 DF,  p-value: 2.148e-07
        
        AIC(best.lm.abund.8t) 
         # 126.9312
        
        plot(best.lm.abund.8t)
        
        # Plot residuals vs each variable (in the model)
        plot((resid(best.lm.abund.8t))~dat.8t$Log.Ht.veg)
        plot((resid(best.lm.abund.8t))~dat.8t$SRL)
        plot((resid(best.lm.abund.8t))~dat.8t$Myc.Frac)
        
        # Plot residuals vs each variable (NOT in the model)
        plot((resid(best.lm.abund.8t))~log(dat.8t$Max.Root.Loca))
        plot((resid(best.lm.abund.8t))~log(dat.8t$Log.Lamina.thck))
        plot((resid(best.lm.abund.8t))~log(dat.8t$LDMC))
        plot((resid(best.lm.abund.8t))~log(dat.8t$Log.Leaf.Area))
        plot((resid(best.lm.abund.8t))~log(dat.8t$Leaf.Mass.Frac))
        
        # A3.1.4 lm, 1st-order effects only ####
        #===========================================#
      
        H1<-lm(abundance.ratio~
                  Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac+SRL+Myc.Frac,
                  data=dat.8t)
        
        H0<-lm(abundance.ratio~1,
               data=dat.8t)
        
        best.lm.abund.8t<-step(H1,scope=list(lower=H0,upper=H1),direction='backward',trace = T)
        
        summary(best.lm.abund.8t)
            # Call:
            # lm(formula = abundance.ratio ~ Max.Root.Loca + Myc.Frac, data = dat.8t)
            # 
            # Residuals:
            #     Min      1Q  Median      3Q     Max 
            # -1.9131 -1.0621 -0.4691  0.5486  8.1567 
            # 
            # Coefficients:
            #               Estimate Std. Error t value Pr(>|t|)    
            # (Intercept)     7.4148     1.6992   4.364 0.000119 ***
            # Max.Root.Loca   0.5056     0.3535   1.430 0.162053    
            # Myc.Frac       -8.8992     2.4352  -3.654 0.000887 ***
            # ---
            # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
            # 
            # Residual standard error: 1.96 on 33 degrees of freedom
            # Multiple R-squared:  0.2883,	Adjusted R-squared:  0.2452 
            # F-statistic: 6.685 on 2 and 33 DF,  p-value: 0.003651
        
        AIC(best.lm.abund.8t)
         # 155.4782
        drop1(best.lm.abund.8t)
        best.lm.abund.8t<-update(best.lm.abund.8t,~.-Max.Root.Loca)
        AIC(best.lm.abund.8t)
         # 155.6432
        drop1(best.lm.abund.8t)
        
        summary(best.lm.abund.8t)
            # Call:
            # lm(formula = abundance.ratio ~ Myc.Frac, data = dat.8t)
            # 
            # Residuals:
            #     Min      1Q  Median      3Q     Max 
            # -2.3888 -1.2434 -0.2599  0.3705  8.7824 
            # 
            # Coefficients:
            #             Estimate Std. Error t value Pr(>|t|)    
            # (Intercept)    7.152      1.715   4.170 0.000198 ***
            # Myc.Frac      -7.426      2.240  -3.315 0.002188 ** 
            # ---
            # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
            # 
            # Residual standard error: 1.99 on 34 degrees of freedom
            # Multiple R-squared:  0.2442,	Adjusted R-squared:  0.222 
            # F-statistic: 10.99 on 1 and 34 DF,  p-value: 0.002188
        
        AIC(best.lm.abund.8t) 
         # 155.6432
        
        plot(best.lm.abund.8t) # negative slope in fitted vs residual (non-independence of data?)
                              
        # data not normal, with sp 7,8 and 11 with high residuals
                              # point 7 has high leverage
        Xs[7,2] 
         # Carex scabrata
        
        # Plot residuals vs each variable (in the model)
        plot((resid(best.lm.abund.8t))~dat.8t$Myc.Frac)
        
        # Plot residuals vs each variable (NOT in the model)
        plot((resid(best.lm.abund.8t))~dat.8t$Log.Ht.veg)
        plot((resid(best.lm.abund.8t))~dat.8t$SRL)
        plot((resid(best.lm.abund.8t))~log(dat.8t$Max.Root.Loca))
        plot((resid(best.lm.abund.8t))~log(dat.8t$Log.Lamina.thck))
        plot((resid(best.lm.abund.8t))~log(dat.8t$LDMC))
        plot((resid(best.lm.abund.8t))~log(dat.8t$Log.Leaf.Area))
        plot((resid(best.lm.abund.8t))~log(dat.8t$Leaf.Mass.Frac))
        
   #===============================================================#
   # SUMMARY OF LMs 
   # assumptions badly violated, so we need to fit a glm
   # by far, bestmodel has both myc.frac and log.myc.frac (adjR2=0.75, aic=119)
   # 2nd best model has myc.frac (no Log(Myc.Frac)) with all other variables (adjR2=0.69, aic=126)
   #===============================================================#
        

   # A3.2 - glm ####
  
   # A3.2.1 glm full with myc.frac & log(myc.frac) ####
   #=========================================#
   # Includes all 1st order terms 
   # + significant trait interactions - four 2nd order, one 3rd order
   # + significant non-linear variables.
   
        H1<-glm(abundance.ratio~
              Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac+SRL+Myc.Frac+
                 + log(Myc.Frac)+
                 SRL:Myc.Frac + Log.Leaf.Area:SRL + Max.Root.Loca:SRL + Log.Ht.veg:SRL+
                 Myc.Frac:SRL:Log.Ht.veg,
              data=dat.8t,
              family=Gamma(link='log'),
              maxit=1000)
        
        H0<-glm(abundance.ratio~1,
           data=dat.8t,
           family=Gamma(link='log'),
           maxit=1000)
        
        best.glm.abund.8t<-step(H1,scope=list(lower=H0,upper=H1),direction='backward',trace = T)
        
        summary(best.glm.abund.8t)
        AIC(best.glm.abund.8t)
        # 102.99
        
        # can I simplify any further?
        drop1(best.glm.abund.8t) # drop Max.Root.Loca decreases AIC by 1.9
        best.glm.abund.8t.minus1<-update(best.glm.abund.8t,~.-Max.Root.Loca)
        drop1(best.glm.abund.8t.minus1) # can't simplify any further
        
        summary(best.glm.abund.8t.minus1)
         # Call:
         # glm(formula = abundance.ratio ~ Log.Ht.veg + SRL + Myc.Frac + 
         #     log(Myc.Frac) + SRL:Myc.Frac + Log.Ht.veg:SRL + Log.Ht.veg:SRL:Myc.Frac, 
         #     family = Gamma(link = "log"), data = dat.8t, maxit = 1000)
         # 
         # Deviance Residuals: 
         #     Min       1Q   Median       3Q      Max  
         # -1.8907  -0.4921  -0.1665   0.1199   1.7460  
         # 
         # Coefficients:
         #                         Estimate Std. Error t value Pr(>|t|)   
         # (Intercept)             -9.22969    6.54073  -1.411  0.16923   
         # Log.Ht.veg              -0.91118    0.50775  -1.795  0.08353 . 
         # SRL                     -0.86938    0.30163  -2.882  0.00750 **
         # Myc.Frac                12.59711    6.89394   1.827  0.07834 . 
         # log(Myc.Frac)           -9.89562    4.60171  -2.150  0.04030 * 
         # SRL:Myc.Frac             0.90387    0.37275   2.425  0.02202 * 
         # Log.Ht.veg:SRL           0.26178    0.08752   2.991  0.00574 **
         # Log.Ht.veg:SRL:Myc.Frac -0.27393    0.10820  -2.532  0.01724 * 
         # ---
         # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
         # 
         # (Dispersion parameter for Gamma family taken to be 0.7733397)
         # 
         #     Null deviance: 37.596  on 35  degrees of freedom
         # Residual deviance: 17.588  on 28  degrees of freedom
         # AIC: 92.371
         # 
         # Number of Fisher Scoring iterations: 8
        
         # Explained Deviance
         1-(best.glm.abund.8t.minus1$deviance/best.glm.abund.8t.minus1$null)
            # 0.5321823
         
         AIC(best.glm.abund.8t.minus1)
            # 92.37

         #library(Mumin)
         library(MuMIn)
         AICc(best.glm.abund.8t.minus1)
            # 99.29
         
         # Validate model
         par(mfrow=c(2,2));plot(best.glm.abund.8t.minus1)
            # Looks awesome
         
         # Plot residuals vs each variable (in the model)
         plot((resid(best.glm.abund.8t.minus1))~dat.8t$Log.Ht.veg)
         plot((resid(best.glm.abund.8t.minus1))~dat.8t$SRL)
         plot((resid(best.glm.abund.8t.minus1))~dat.8t$Myc.Frac)

         # Plot residuals vs each variable (NOT in the model)
         plot((resid(best.glm.abund.8t.minus1))~dat.8t$Max.Root.Loca)
         plot((resid(best.glm.abund.8t.minus1))~dat.8t$Leaf.Mass.Frac)
         plot((resid(best.glm.abund.8t.minus1))~dat.8t$Log.Lamina.thck)
         plot((resid(best.glm.abund.8t.minus1))~dat.8t$LDMC)
         plot((resid(best.glm.abund.8t.minus1))~dat.8t$Log.Leaf.Area)
         
         # Try with MuMIn{} dredge()
         # rewrite H1 with na.action='fail' - required by dredge ()
         H1<-glm(abundance.ratio~
                    Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac+SRL+Myc.Frac+
                    + log(Myc.Frac)+
                    SRL:Myc.Frac + Log.Leaf.Area:SRL + Max.Root.Loca:SRL + Log.Ht.veg:SRL+
                    Myc.Frac:SRL:Log.Ht.veg,
                 data=dat.8t,
                 family=Gamma(link='log'),
                 maxit=1000,na.action='na.fail')
         
         d.abund<-dredge(H1,beta='sd',rank=AIC,trace=F)
         head(d.abund)
         # best model includes:
         # log(Myc.Frac) + Myc.Frac +
         # SRL+ Log.Ht.veg + Log.Lmn.thc +
         # Log.Ht.veg:SRL + Myc.Frc:SRL
         # Log.Ht.veg:Myc.Frc:SRL
         
         d.abund.aicc<-dredge(H1,beta='sd',rank=AICc,trace=F)
         head(d.abund.aicc)
         # best model includes:
         # Myc.Frac + log(Myc.Frac)
         
         
         
   # A3.2.2 glm full, with Myc.Frac only #### 
   #=========================================#
   # Not both Myc.Frac and log(Myc.Frac).
  
    H1<-glm(abundance.ratio~
              Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac+SRL+Myc.Frac+
              SRL:Myc.Frac + Log.Leaf.Area:SRL + Max.Root.Loca:SRL + Log.Ht.veg:SRL+
              Myc.Frac:SRL:Log.Ht.veg,
            data=dat.8t,
            family=Gamma(link='log'),
            maxit=100,na.action=na.omit)

    # use 3-way interactions with Myc.Frac
    H0<-glm(abundance.ratio~1,
           data=dat.8t,
           family=Gamma(link='log'),
           maxit=100)
    
    best.glm.abund.8t<-step(H1,scope=list(lower=H0,upper=H1),direction='backward',trace = T)
    
    # Warning: glm.fit: algorithm did not converge - fixed - add parameter maxit=100 in H0 and H1 models
    # why did I include a dispersion parameter in previous model ?
      # family=Gamma(link='log')),dispersion=1
    
    # Explained Deviance
     1-(best.glm.abund.8t$deviance/best.glm.abund.8t$null)
    # 0.4825224
     
     AIC(best.glm.abund.8t) 
     # 96.30
     library (MuMIn); AICc(best.glm.abund.8t)
     # 103.23
     
     summary(best.glm.abund.8t)
     
     # Call:
     #    glm(formula = abundance.ratio ~ Log.Ht.veg + Leaf.Mass.Frac + 
     #           SRL + Myc.Frac + SRL:Myc.Frac + Log.Ht.veg:SRL + Log.Ht.veg:SRL:Myc.Frac, 
     #        family = Gamma(link = "log"), data = dat.8t, maxit = 100)
     # 
     # Deviance Residuals: 
     #    Min       1Q   Median       3Q      Max  
     # -1.6157  -0.5183  -0.2746   0.1531   1.7109  
     # 
     # Coefficients:
     #                         Estimate Std. Error t value Pr(>|t|)   
     # (Intercept)              3.74528    1.85804   2.016  0.05352 . 
     # Log.Ht.veg              -0.50642    0.50807  -0.997  0.32742   
     # Leaf.Mass.Frac           1.52255    0.72489   2.100  0.04483 * 
     # SRL                     -0.77606    0.29277  -2.651  0.01307 * 
     # Myc.Frac                -3.50276    1.97328  -1.775  0.08676 . 
     # SRL:Myc.Frac             0.81487    0.36596   2.227  0.03419 * 
     # Log.Ht.veg:SRL           0.23476    0.08478   2.769  0.00986 **
     # Log.Ht.veg:SRL:Myc.Frac -0.24574    0.10548  -2.330  0.02726 * 
     #    ---
     #    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
     # 
     # (Dispersion parameter for Gamma family taken to be 0.7247453)
     # 
     # Null deviance: 37.596  on 35  degrees of freedom
     # Residual deviance: 19.455  on 28  degrees of freedom
     # AIC: 96.307
     # 
     # Number of Fisher Scoring iterations: 12
    
    # can I simplify any further?
    summary(best.glm.abund.8t)
    drop1(best.glm.abund.8t) 
      # droping LMF only increases AIC by 1.2, but the term is significant
    
          best.glm.abund.8t.minus1<-update(best.glm.abund.8t,~.-Leaf.Mass.Frac)
          summary(best.glm.abund.8t.minus1)# most variables are now NS
          drop1(best.glm.abund.8t.minus1)# increases AIC to above 2.0 above best model.  
          
          
          AIC(best.glm.abund.8t.minus1) 
            #98.83
          AICc(best.glm.abund.8t.minus1)
            # 104.16
    


    # validate model
    par(mfrow=c(2,2))
    plot(best.glm.abund.8t)
    # YESS !!! All good :-)
    
    # Plot residuals vs each variable (in the model)
    plot((resid(best.glm.abund.8t))~dat.8t$Log.Ht.veg)
    plot((resid(best.glm.abund.8t))~log(dat.8t$Leaf.Mass.Frac))
    plot((resid(best.glm.abund.8t))~dat.8t$Myc.Frac)
    plot((resid(best.glm.abund.8t))~dat.8t$SRL)
    # Plot residuals vs each variable (NOT in the model)
    plot((resid(best.glm.abund.8t))~log(dat.8t$Max.Root.Loca))
    plot((resid(best.glm.abund.8t))~log(dat.8t$Log.Lamina.thck))
    plot((resid(best.glm.abund.8t))~log(dat.8t$LDMC))
    plot((resid(best.glm.abund.8t))~log(dat.8t$Log.Leaf.Area))
    
    # A3.2.3 glm full, with log(Myc.Frac) only #### 
    #=========================================#
    # Not both Myc.Frac and log(Myc.Frac).
    
    H1<-glm(abundance.ratio~
               Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac+SRL+log(Myc.Frac)+
               SRL:log(Myc.Frac) + Log.Leaf.Area:SRL + Max.Root.Loca:SRL + Log.Ht.veg:SRL+
               log(Myc.Frac):SRL:Log.Ht.veg,
            data=dat.8t,
            family=Gamma(link='log'),
            maxit=100,na.action=na.omit)
    
    # use 3-way interactions with Myc.Frac
    H0<-glm(abundance.ratio~1,
            data=dat.8t,
            family=Gamma(link='log'),
            maxit=100)
    
    best.glm.abund.8t<-step(H1,scope=list(lower=H0,upper=H1),direction='backward',trace = T)
    
    # Warning: glm.fit: algorithm did not converge - fixed - add parameter maxit=100 in H0 and H1 models
    # why did I include a dispersion parameter in previous model ?
    # family=Gamma(link='log')),dispersion=1
    
    # Explained Deviance
    1-(best.glm.abund.8t$deviance/best.glm.abund.8t$null)
    # 0.4990414
    
    AIC(best.glm.abund.8t) 
    # 95.03
    library (MuMIn); AICc(best.glm.abund.8t)
    # 101.9612
    
    # can I simplify any further?
    drop1(best.glm.abund.8t) 
      # increase AIC by only 1.1 if I drop Leaf.Mass.Frac, 
      # BUT all terms quickly become NS if I do so. 
    
    summary(best.glm.abund.8t)
       # Call:
       #    glm(formula = abundance.ratio ~ Log.Ht.veg + Leaf.Mass.Frac + 
       #           SRL + log(Myc.Frac) + SRL:log(Myc.Frac) + Log.Ht.veg:SRL + 
       #           Log.Ht.veg:SRL:log(Myc.Frac), family = Gamma(link = "log"), 
       #        data = dat.8t, na.action = na.omit, maxit = 100)
       # 
       # Deviance Residuals: 
       #     Min       1Q   Median       3Q      Max  
       # -1.6865  -0.5066  -0.2450   0.2125   1.7284  
       # 
       # Coefficients:
       # Estimate Std. Error t value Pr(>|t|)  
       # (Intercept)                   0.32640    1.86403   0.175   0.8623  
       # Log.Ht.veg                   -0.49768    0.49257  -1.010   0.3210  
       # Leaf.Mass.Frac                1.44981    0.70421   2.059   0.0489 *
       # SRL                           0.01359    0.10886   0.125   0.9015  
       # log(Myc.Frac)                -2.62064    1.38320  -1.895   0.0685 .
       # SRL:log(Myc.Frac)             0.56513    0.25952   2.178   0.0380 *
       # Log.Ht.veg:SRL               -0.00293    0.03373  -0.087   0.9314  
       # Log.Ht.veg:SRL:log(Myc.Frac) -0.16657    0.07167  -2.324   0.0276 *
       #    ---
       # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
       # 
       # (Dispersion parameter for Gamma family taken to be 0.705269)
       # 
       # Null deviance: 37.596  on 35  degrees of freedom
       # Residual deviance: 18.834  on 28  degrees of freedom
       # AIC: 95.038
       # 
       # Number of Fisher Scoring iterations: 11
    

    
       # SUMMARY OF GLMs ===============================================# 
       # Model w both Myc.Frac and log(Myc.Frac) is best  
       # Pseudo-R2 = 53%, AIC=92
       # Model with either Myc.Frac or log(Myc.Frac) perform equally well
       # and select for same variables
       # Pseudo-R2 = 0.48 and 0.50; AIC= 96 and 95 (with Myc.Frac and log(Myc.Frac, respectively))
       # all 3 models retain same variables
       #================================================================#
    

   # A4 Rerun Tree & GLM models without high-leverage points ####
   #============================================================#
   # (CIAL, LYAN,OSCL) 
  ----
  dim(dat.8t) #36 12
  dat.8t.no.outliers<-dat.8t[!rownames(H.abund.c)%in%c('CIAL','LYAN','OSCL'),]
  dim(dat.8t.no.outliers) #33 12
  
  best.glm.abund.8t.no.outlier<-glm(formula = abundance.ratio ~ Log.Ht.veg + Leaf.Mass.Frac + SRL + Myc.Frac +
         SRL:Myc.Frac + Log.Ht.veg:SRL + Log.Ht.veg:SRL:Myc.Frac, 
      family = Gamma(link = "log"), data = dat.8t.no.outliers, na.action = "na.fail", 
      maxit = 100)
  
  summary(best.glm.abund.8t.no.outlier)
     #Call:
     # glm(formula = abundance.ratio ~ Log.Ht.veg + Leaf.Mass.Frac + 
     #     SRL + Myc.Frac + SRL:Myc.Frac + Log.Ht.veg:SRL + Log.Ht.veg:SRL:Myc.Frac, 
     #     family = Gamma(link = "log"), data = dat.8t.no.outliers, 
     #     na.action = "na.fail", maxit = 100)
     # 
     # Deviance Residuals: 
     #    Min        1Q    Median        3Q       Max  
     # -1.55224  -0.25626  -0.05457   0.14405   1.45426  
     # 
     # Coefficients:
     #                         Estimate Std. Error t value Pr(>|t|)    
     # (Intercept)              2.77280    1.56222   1.775 0.088091 .  
     # Log.Ht.veg              -0.67373    0.41854  -1.610 0.120017    
     # Leaf.Mass.Frac           0.81224    0.58860   1.380 0.179817    
     # SRL                     -0.87213    0.24219  -3.601 0.001369 ** 
     # Myc.Frac                -0.81428    1.63627  -0.498 0.623086    
     # SRL:Myc.Frac             0.89672    0.31627   2.835 0.008936 ** 
     # Log.Ht.veg:SRL           0.28732    0.07225   3.977 0.000526 ***
     # Log.Ht.veg:SRL:Myc.Frac -0.31171    0.09874  -3.157 0.004128 ** 
     # ---
     # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
     # 
     # (Dispersion parameter for Gamma family taken to be 0.4317756)
     # 
     # Null deviance: 33.034  on 32  degrees of freedom
     # Residual deviance: 11.574  on 25  degrees of freedom
     # AIC: 67.519
  
  1-(best.glm.abund.8t.no.outlier$deviance/best.glm.abund.8t.no.outlier$null)
   # 0.6496394
  
  
  tree.model<-tree(dat.8t.no.outliers$abundance.ratio~.,data=dat.8t.no.outliers[,Trait.Names.8t])
  par(mfrow=c(1,1));plot(tree.model);  text(tree.model)
  title('Regression Tree, no outliers \n Abundance.ratio vs 8 traits')
  tree.model
  
   # 1) root 33 152.5000 1.3860  
   #     2) Myc.Frac < 0.595 5 102.9000 3.6480 *
   #     3) Myc.Frac > 0.595 28  19.4500 0.9824  
   #        6) LDMC < 0.198315 18   2.2330 0.6716 *
   #        7) LDMC > 0.198315 10  12.3500 1.5420  
   #           14) Max.Root.Loca < 2.38333 5   8.3010 2.1120 *
   #           15) Max.Root.Loca > 2.38333 5   0.7923 0.9712 *
  
  1-(deviance(tree.model)/152.5) 
  # 0.2511
  plot(cv.tree(tree.model))
  
  # c.v.  splits the data into 'training' set for model fitting and a 
  # validation set to evaluate goodness of fit 
  # look for how many splits produces the minimum deviance - run multiple times
  
  # very unstable. most frequent results are lowest dev at 1-2 or 2-3 splits. 
  
  # SUMMARY of removing outliers ============================# 
  # Model still valid without outliers, but 
  # Different model is selected without the 3 outlier species
  #==========================================================#
  
#==============================#
# (B) Occurence vs 8 Traits ####
#==============================#
  
     # B1 - Explore Trait Interactions #### 
     #==============================#
     
     # B1.1 - Tree Model
     par(mfrow=c(1,1))
     tree.model<-tree(dat.8t$occurence.ratio~.,data=dat.8t[,Trait.Names.8t])
     plot(tree.model); text(tree.model); title('Regression Tree \n Occurence.ratio vs 8 traits')
     # Myc.Frac is best predictor,
     # for those with high Myc.Frac, SRL is best predictor
  
      tree.model
         # 1) root 36 30.7900 1.5530  
         #     2) Myc.Frac < 0.595 6 16.1900 2.4540 *
         #     3) Myc.Frac > 0.595 30  8.7560 1.3730  
         #        6) SRL < 27.7246 25  7.5060 1.4570  
         #           12) Max.Root.Loca < 1.71364 13  5.6240 1.6400  
         #               24) Log.Lamina.thck < 4.67958 5  3.2850 2.1500 *
         #               25) Log.Lamina.thck > 4.67958 8  0.2284 1.3220 *
         #           13) Max.Root.Loca > 1.71364 12  0.9730 1.2580 *
         #        7) SRL > 27.7246 5  0.1851 0.9514 *
      
      1-(deviance(tree.model)/30.79)
      # 0.3223963
      plot(cv.tree(tree.model))
      # c.v.  splits the data into 'training' set for model fitting and a 
      # validation set to evaluate goodness of fit 
      # look for how many splits produces the minimum deviance - run multiple times
      
      # unstable. selects 1-2 or 3-4
      # could stop at 2 branches?
      
         #=====================================================================================#
         # Summary of regression tree
         # Regression Tree retains Myc.Frac, SRL, Max.Root.Loca, Log.Lamian.thck and their interactions - R2 = 0.29
         #=====================================================================================#
         
      # B1.2 - Test trait-trait interactions 
      ----
         
      # need to check in subgroups because sample size too small to test all interactions
      # simultaneously
      summary(lm(dat.8t$occurence.ratio~(Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC)^2,
                 data=dat.8t))
   
      # Max.Root.Loca:LDMC, @p= 0.0461
      # Log.Ht.veg:LDMC,    @p= 0.0467
      
      summary(lm(dat.8t$occurence.ratio~(Log.Leaf.Area+Leaf.Mass.Frac+SRL+Myc.Frac)^2,
                 data=dat.8t))
      
      # Log.Leaf.Area:Myc.Frac     @p=0.00015
      
      summary(lm(dat.8t$occurence.ratio~(Log.Ht.veg+Max.Root.Loca+Log.Leaf.Area+Leaf.Mass.Frac)^2,
                 data=dat.8t))
      
      # Max.Root.Loca:Log.Leaf.Area @p=0.0004
      # Log.Leaf.Area:Leaf.Mass.Frac @p=0.012075
      # Log.Ht.veg:Log.Leaf.Area     @p=0.031010
      
      summary(lm(dat.8t$occurence.ratio~(Log.Ht.veg+Max.Root.Loca+SRL+Myc.Frac)^2,
                 data=dat.8t))
      
      # No significant terms
      
            #=====================================================================================#
            # Summary of trait-trait interactions. 
            # None of the interactions in regression tree are selected by lm. Potential candidates are:
            #  - Log.Leaf.Area:Myc.Frac
            #  - Log.Leaf.Area:Max.Root.Loca
            #  - Log.Leaf.Area:Leaf.Mass.Frac
            #  - Log.Leaf.Area:Log.Ht.veg
            #  - LDMC:Log.Ht.veg
            #  - LDMC:Max.Root.Loca
            #=====================================================================================#
      
      # B2 - Explore non-linearity ####
      #=============================#
      
      # B2.1 GAMs
      ---
      par(mfrow=c(2,2),mar=c(4,4,3,3))
      plot(gam(occurence.ratio~s(Log.Ht.veg)+
                  +s(SRL)+s(Myc.Frac),data=dat.8t))
         # myc.frac is non-linear
      
      plot(gam(occurence.ratio~s(Max.Root.Loca)+s(Log.Lamina.thck)+
                  s(LDMC),data=dat.8t))
         # Linear relationships
      
      plot(gam(occurence.ratio~s(Log.Leaf.Area)+s(Leaf.Mass.Frac),data=dat.8t))
         # Non-linear for Log.Leaf.Area
         # can't transform it further - already a log - gam w log link to transform occurence
         # equivalent to logging occurence.ratio
      
      # B2.2 Test polynomial relationships 
      ----
      summary(lm(occurence.ratio~I(Log.Ht.veg^2)+I(Max.Root.Loca^2)+I(Log.Lamina.thck^2)+
                       I(LDMC^2)+I(Log.Leaf.Area^2)+I(Leaf.Mass.Frac^2)+I(SRL^2)+I(Myc.Frac^2),
                    data=dat.8t))
         # No significant relationships
      
      # Myc.Frac
      summary(lm(occurence.ratio~log(Myc.Frac),
                 data=dat.8t))
         # log(myc.frac) significant @ p=0.0018 ; AdjR2=23%
      
      summary(lm(occurence.ratio~Myc.Frac+log(Myc.Frac),
                 data=dat.8t))
         # log(myc.frac) significant @p=4.91e-07 Myc.Frac significant @p=3.10e-06; AdjR2=59%
      
      # Log.Leaf.Area
      summary(lm(occurence.ratio~Log.Leaf.Area,
                 data=dat.8t))
         # significant @p=0.0087
      
         #=====================================================================================#
         # Summary of non-linearity. 
         # Myc.Frac is non-linear and fits best as Myc.Frac + log(Myc.Frac)
         # Use a gamma with log link to account for non-linear y ~ Log.Leaf.Area
         #=====================================================================================#
         
      
# B3 - Model selection ####
#=========================#
      
      # No lms - response variable lower-bound @ zero - violates normal distribution assumptions 
      
      # B3.2 - glm ####
      
      # B3.2.1 glm full with myc.frac & log(myc.frac) ####
      #=========================================#
      # Includes all 1st order terms 
      # + significant non-linear variables
      # + significant trait interactions - seven 2nd order, one 3rd order
      # .
      
      H1<-glm(occurence.ratio~
                 Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac+SRL+Myc.Frac+
                 log(Myc.Frac)+
                 Myc.Frac:SRL + SRL:Max.Root.Loca + Myc.Frac:SRL:Max.Root.Loca+
                 Log.Leaf.Area:Max.Root.Loca + Log.Leaf.Area:Max.Root.Loca + Log.Leaf.Area:Log.Ht.veg +
                 LDMC:Log.Ht.veg + LDMC:Max.Root.Loca,
              data=dat.8t,
              family=Gamma(link='log'),
              maxit=1000)

      
      H0<-glm(occurence.ratio~1,
              data=dat.8t,
              family=Gamma(link='log'),
              maxit=1000)
      
      best.glm.occur.8t<-step(H1,scope=list(lower=H0,upper=H1),direction='backward',trace = T)
      AIC(best.glm.occur.8t)
         #56.04
      
      # Can I simplify any further?
      drop1(best.glm.occur.8t) # removing Log.Ht.veg:LDMC drops AIC by 0.3
      best.glm.occur.8t.minus1<-update(best.glm.occur.8t,~.-Log.Ht.veg:LDMC)
      AIC(best.glm.occur.8t.minus1) 
         # 56.56
      drop1(best.glm.occur.8t.minus1) # removing LDMC drops AIC by 1.9
      best.glm.occur.8t.minus2<-update(best.glm.occur.8t.minus1,~.-LDMC)
      AIC(best.glm.occur.8t.minus2)
         # 54.74
      drop1(best.glm.occur.8t.minus2) # removing Log.Ht.veg decreases AIC by 1.4
      best.glm.occur.8t.minus3<-update(best.glm.occur.8t.minus2,~.-Log.Ht.veg)
      AIC(best.glm.occur.8t.minus3)
         # 53.60
      drop1(best.glm.occur.8t.minus3) # No
      
      # Lowest Global AIC = minus3 
            
      1-(deviance(best.glm.occur.8t.minus3)/best.glm.occur.8t.minus3$null.deviance)
         # 0.6153814
      AIC( best.glm.occur.8t.minus3)
         # 53.60
      AICc(best.glm.occur.8t.minus3)
         # 60.52
      
      summary(best.glm.occur.8t.minus3)
      
      # Call:
      # glm(formula = occurence.ratio ~ Max.Root.Loca + SRL + Myc.Frac + 
      #        log(Myc.Frac) + SRL:Myc.Frac + Max.Root.Loca:SRL + Max.Root.Loca:SRL:Myc.Frac, 
      #        family = Gamma(link = "log"), data = dat.8t, maxit = 1000)
      # 
      # Deviance Residuals: 
      #      Min        1Q    Median        3Q       Max  
      # -0.46667  -0.28745  -0.00207   0.18468   0.80930  
      # 
      # Coefficients:
      # Estimate Std. Error t value Pr(>|t|)    
      # (Intercept)                -11.187967   2.618448  -4.273 0.000201 ***
      # Max.Root.Loca                0.002531   0.118907   0.021 0.983166    
      # SRL                         -0.128116   0.054978  -2.330 0.027221 *  
      # Myc.Frac                    11.853757   2.689533   4.407 0.000140 ***
      # log(Myc.Frac)               -8.969936   1.911340  -4.693 6.42e-05 ***
      # SRL:Myc.Frac                 0.167918   0.076954   2.182 0.037656 *  
      # Max.Root.Loca:SRL            0.052807   0.022517   2.345 0.026334 *  
      # Max.Root.Loca:SRL:Myc.Frac  -0.072163   0.029321  -2.461 0.020276 *  
      # ---
      # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # 
      # (Dispersion parameter for Gamma family taken to be 0.1089945)
      # 
      # Null deviance: 7.4065  on 35  degrees of freedom
      # Residual deviance: 2.8487  on 28  degrees of freedom
      # AIC: 53.599
      # 
      # Number of Fisher Scoring iterations: 7
      
      # Validate model
      par(mfrow=c(2,2));plot(best.glm.occur.8t.minus5,which = c(1,3:5))
      # LYAN is an outlier with high leverage
      
      # Plot residuals vs each variable (in the model)
      plot((resid(best.glm.occur.8t.minus5))~dat.8t$Max.Root.Loca)
      plot((resid(best.glm.occur.8t.minus5))~dat.8t$Myc.Frac)
      plot((resid(best.glm.occur.8t.minus5))~dat.8t$SRL)
      
      # Plot residuals vs each variable (NOT in the model)
      plot((resid(best.glm.occur.8t.minus5))~dat.8t$Log.Ht.veg)
      plot((resid(best.glm.occur.8t.minus5))~dat.8t$Leaf.Mass.Frac)
      plot((resid(best.glm.occur.8t.minus5))~dat.8t$Log.Lamina.thck)
      plot((resid(best.glm.occur.8t.minus5))~dat.8t$LDMC)
      plot((resid(best.glm.occur.8t.minus5))~dat.8t$Log.Leaf.Area)
      
      # Try with MuMIn {} dredge ()
      # rewrite H1 with na.action='fail' - required by dredge ()
      H1<-glm(occurence.ratio~
                 Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac+SRL+Myc.Frac+
                 + log(Myc.Frac)+
                 Myc.Frac:SRL+SRL:Max.Root.Loca+Myc.Frac:SRL:Max.Root.Loca+
                 Log.Leaf.Area:Max.Root.Loca + Log.Leaf.Area:Max.Root.Loca + Log.Leaf.Area:Log.Ht.veg+
                 LDMC:Log.Ht.veg + LDMC:Max.Root.Loca,
              data=dat.8t,
              family=Gamma(link='log'),
              maxit=1000,na.action='na.fail')
      
      d.occur<-dredge(H1,beta='sd',rank=AIC,trace=F)
      head(d.occur)
         # best model includes:
         # Myc.Frac + log(Myc.Frac) 
         # + Max.Root.Loca + SRL 
         # + Max.Rot.Loc:SRL
         # + Myc.Frc:SRL
         # + Max.Rot.Loc:Myc.Frc:SRL
      
         # 2nd equivalently good model (same AIC) includes:
         # Myc.Frac + log(Myc.Frac) 
      
      d.occur.aicc<-dredge(H1,beta='sd',rank=AICc,trace=F)
      head(d.occur.aicc)
         # # best model includes Myc.Frac + log(Myc.Frac)
      
      # B3.2.2 glm full with myc.frac only ####
      #=========================================#
      # Not Myc.Frac and log(Myc.Frac)
      
      H1<-glm(occurence.ratio~
                 Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac+SRL+Myc.Frac+
                 Myc.Frac:SRL+SRL:Max.Root.Loca+Myc.Frac:SRL:Max.Root.Loca+
                 Log.Leaf.Area:Max.Root.Loca + Log.Leaf.Area:Max.Root.Loca + Log.Leaf.Area:Log.Ht.veg+
                 LDMC:Log.Ht.veg + LDMC:Max.Root.Loca,
              data=dat.8t,
              family=Gamma(link='log'),
              maxit=1000)
      
      H0<-glm(occurence.ratio~1,
              data=dat.8t,
              family=Gamma(link='log'),
              maxit=1000)
      
      best.glm.occur.8t<-step(H1,scope=list(lower=H0,upper=H1),direction='backward',trace = T)
      
      
      1-(deviance(best.glm.occur.8t)/best.glm.occur.8t$null.deviance)
         #  0.5600242
      AIC(best.glm.occur.8t)
         # 66.51
      AICc(best.glm.occur.8t)
         # 83.05
      
      summary(best.glm.occur.8t)
      # Coefficients:
      #                             Estimate Std. Error t value Pr(>|t|)   
      # (Intercept)                  2.49958    1.15085   2.172  0.03997 * 
      # Log.Ht.veg                  -0.21786    0.30940  -0.704  0.48813   
      # Max.Root.Loca               -0.33169    0.14997  -2.212  0.03675 * 
      # LDMC                        -5.48208    3.32229  -1.650  0.11195   
      # Log.Leaf.Area               -0.30023    0.10130  -2.964  0.00677 **
      # SRL                         -0.11447    0.06629  -1.727  0.09707 . 
      # Myc.Frac                    -0.78618    0.96748  -0.813  0.42443   
      # SRL:Myc.Frac                 0.12866    0.09031   1.425  0.16711   
      # Max.Root.Loca:SRL            0.05619    0.02818   1.994  0.05759 . 
      # Max.Root.Loca:Log.Leaf.Area  0.10827    0.04294   2.521  0.01874 * 
      # Log.Ht.veg:LDMC              1.99315    1.28439   1.552  0.13379   
      # Max.Root.Loca:SRL:Myc.Frac  -0.06593    0.03601  -1.831  0.07958 . 
      # ---
      # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # 
      # (Dispersion parameter for Gamma family taken to be 0.1417102)
      # 
      # Null deviance: 7.4065  on 35  degrees of freedom
      # Residual deviance: 3.2587  on 24  degrees of freedom
      # AIC: 66.508
      # 
      # Number of Fisher Scoring iterations: 7
      
      # Can I simplify any further?
      drop1(best.glm.occur.8t) # Dropping Log.Ht.veg:LDMC increases AIC by 0.2 
      best.glm.occur.8t.minus1<-update(best.glm.occur.8t,~.-Log.Ht.veg:LDMC)
      AIC(best.glm.occur.8t.minus1)
         # 67.85
      drop1(best.glm.occur.8t.minus1) # dropping LDMC decreases AIC by 1.7
      best.glm.occur.8t.minus2<-update(best.glm.occur.8t.minus1,~.-LDMC)
      AIC(best.glm.occur.8t.minus2)
         # 66.39
      drop1(best.glm.occur.8t.minus2) # dropping Log.Ht.veg increases AIC by 0.1
      best.glm.occur.8t.minus3<-update(best.glm.occur.8t.minus2,~.-Log.Ht.veg)
      AIC(best.glm.occur.8t.minus3)
         # 67.58
      drop1(best.glm.occur.8t.minus3) # dropping Max.Root.Loca:SRL:Myc.Frac increases AIC by 1.2
      best.glm.occur.8t.minus4<-update(best.glm.occur.8t.minus3,~.-Max.Root.Loca:SRL:Myc.Frac)
      AIC(best.glm.occur.8t.minus4)
         # 69.91
      drop1(best.glm.occur.8t.minus4)# dropping Max.Root.Loca:SRL decreases AIC by 0.7
      best.glm.occur.8t.minus5<-update(best.glm.occur.8t.minus4,~.-Max.Root.Loca:SRL)
      AIC(best.glm.occur.8t.minus5)   
         # 68.62 
      drop1(best.glm.occur.8t.minus5) # dropping SRL:Myc.Frac decreases AIC by 1.3
      best.glm.occur.8t.minus6<-update(best.glm.occur.8t.minus5,~.-SRL:Myc.Frac)
      AIC(best.glm.occur.8t.minus6)
         # 67.58
      drop1(best.glm.occur.8t.minus6) # dropping Myc.Frac decreases AIC by 1.2
      best.glm.occur.8t.minus7<-update(best.glm.occur.8t.minus6,~.-Myc.Frac)
      AIC(best.glm.occur.8t.minus7)
         # 66.70
      drop1(best.glm.occur.8t.minus7) # dropping SRL decreases AIC by 0.7
      best.glm.occur.8t.minus8<-update(best.glm.occur.8t.minus7,~.-SRL)
      AIC(best.glm.occur.8t.minus8)
         # 66.50
      drop1(best.glm.occur.8t.minus8) # none
      
         # two best model with best global minimum: orginal best and minus8
         summary(best.glm.occur.8t.minus8)
         # Coefficients:
         # Estimate Std. Error t value Pr(>|t|)    
         # (Intercept)                  0.86384    0.17226   5.015  1.9e-05 ***
         # Max.Root.Loca               -0.17782    0.08657  -2.054   0.0482 *  
         # Log.Leaf.Area               -0.19494    0.07746  -2.517   0.0171 *  
         # Max.Root.Loca:Log.Leaf.Area  0.06627    0.03829   1.731   0.0931 .  
         # ---
         # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
         # 
         # (Dispersion parameter for Gamma family taken to be 0.1636046)
         # 
         # Null deviance: 7.4065  on 35  degrees of freedom
         # Residual deviance: 5.0398  on 32  degrees of freedom
         # AIC: 66.502 
         # 
         # Number of Fisher Scoring iterations: 8
         
      1-(deviance(best.glm.occur.8t.minus8)/best.glm.occur.8t.minus8$null.deviance)
         # 0.3195456
      AIC( best.glm.occur.8t.minus8)
         # 66.50
      AICc(best.glm.occur.8t.minus8)
         # 68.50
      
      # Validate model
      par(mfrow=c(2,2));plot(best.glm.occur.8t.minus8,which = c(1,3:5))
         # LYAN is an outlier with too much leverage 
      
      # Plot residuals vs each variable (in the model)
      plot((resid(best.glm.occur.8t.minus8))~dat.8t$Myc.Frac)
      plot((resid(best.glm.occur.8t.minus8))~dat.8t$Log.Leaf.Area)
      
      # Plot residuals vs each variable (NOT in the model)
      plot((resid(best.glm.occur.8t.minus8))~dat.8t$Log.Ht.veg)
      plot((resid(best.glm.occur.8t.minus8))~dat.8t$SRL)
      plot((resid(best.glm.occur.8t.minus8))~dat.8t$Max.Root.Loca)
      plot((resid(best.glm.occur.8t.minus8))~dat.8t$Leaf.Mass.Frac)
      plot((resid(best.glm.occur.8t.minus8))~dat.8t$Log.Lamina.thck)
      plot((resid(best.glm.occur.8t.minus8))~dat.8t$LDMC)
   
      # B3.2.3 glm full with log(myc.frac) only ####
      #=========================================#
      # Not Myc.Frac and log(Myc.Frac)
      
      H1<-glm(occurence.ratio~
                 Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac+SRL+log(Myc.Frac)+
                 log(Myc.Frac):SRL+SRL:Max.Root.Loca+log(Myc.Frac):SRL:Max.Root.Loca+
                 Log.Leaf.Area:Max.Root.Loca + Log.Leaf.Area:Max.Root.Loca + Log.Leaf.Area:Log.Ht.veg+
                 LDMC:Log.Ht.veg + LDMC:Max.Root.Loca,
              data=dat.8t,
              family=Gamma(link='log'),
              maxit=1000)
      
      H0<-glm(occurence.ratio~1,
              data=dat.8t,
              family=Gamma(link='log'),
              maxit=1000)
      
      best.glm.occur.8t<-step(H1,scope=list(lower=H0,upper=H1),direction='backward',trace = T)
      
      1-(deviance(best.glm.occur.8t)/best.glm.occur.8t$null.deviance)
      # 0.5561831
      AIC( best.glm.occur.8t)
      # 66.83
      AICc(best.glm.occur.8t)
      # 83.37112
      
      summary(best.glm.occur.8t) # best global minimum AIC
      # Coefficients:
      # Estimate Std. Error t value Pr(>|t|)  
      # (Intercept)                      1.739941   0.863932   2.014   0.0554 .
      # Log.Ht.veg                      -0.243406   0.314005  -0.775   0.4458  
      # Max.Root.Loca                   -0.314525   0.157302  -1.999   0.0570 .
      # LDMC                            -5.782599   3.347172  -1.728   0.0969 .
      # Log.Leaf.Area                   -0.295976   0.110026  -2.690   0.0128 *
      # SRL                              0.009599   0.025719   0.373   0.7122  
      # log(Myc.Frac)                   -0.719863   0.716116  -1.005   0.3248  
      # SRL:log(Myc.Frac)                0.090923   0.068784   1.322   0.1987  
      # Max.Root.Loca:SRL               -0.007180   0.010708  -0.671   0.5089  
      # Max.Root.Loca:Log.Leaf.Area      0.106629   0.045614   2.338   0.0281 *
      # Log.Ht.veg:LDMC                  2.098881   1.293137   1.623   0.1176  
      # Max.Root.Loca:SRL:log(Myc.Frac) -0.044189   0.026776  -1.650   0.1119  
      # ---
      # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # 
      # (Dispersion parameter for Gamma family taken to be 0.1444615)
      # 
      # Null deviance: 7.4065  on 35  degrees of freedom
      # Residual deviance: 3.2871  on 24  degrees of freedom
      # AIC: 66.826

      # Can I simplify any further?
      drop1(best.glm.occur.8t) # dropping Log.Ht.veg:LDMC increases AIC by 0.4
      best.glm.occur.8t.minus1<-update(best.glm.occur.8t,~.-Log.Ht.veg:LDMC)
      AIC(best.glm.occur.8t.minus1)
         #68.53
      drop1(best.glm.occur.8t.minus1)# dropping LDMC decreases AIC by 1.8
      best.glm.occur.8t.minus2<-update(best.glm.occur.8t.minus1,~.-LDMC)
      AIC(best.glm.occur.8t.minus2)
         # 67.12
      drop1(best.glm.occur.8t.minus2) # dropping Log.Ht.veg increases AIC by 0.02
      best.glm.occur.8t.minus3<-update(best.glm.occur.8t.minus2,~.-Log.Ht.veg)
      AIC(best.glm.occur.8t.minus3)
         # 68.26
      drop1(best.glm.occur.8t.minus3) # dropping Max.Root.Loca:SRL:log(Myc.Frac) decreases AIC by 0.03
      best.glm.occur.8t.minus4<-update(best.glm.occur.8t.minus3,~.-Max.Root.Loca:SRL:log(Myc.Frac))
      AIC(best.glm.occur.8t.minus4)
         # 69.12
      drop1(best.glm.occur.8t.minus4) # dropping Max.Root.Loca:SRL decreases AIC by 1.7
      best.glm.occur.8t.minus5<-update(best.glm.occur.8t.minus4,~.-Max.Root.Loca:SRL)
      AIC(best.glm.occur.8t.minus5)
         # 67.60
      drop1(best.glm.occur.8t.minus5) # dropping SRL:log(Myc.Frac) decreases AIC by 0.8
      best.glm.occur.8t.minus6<-update(best.glm.occur.8t.minus5,~.-SRL:log(Myc.Frac))
      AIC(best.glm.occur.8t.minus6)
         # 66.54
      drop1(best.glm.occur.8t.minus6) # dropping SRL decreases AIC by 0.9
      best.glm.occur.8t.minus7<-update(best.glm.occur.8t.minus6,~.-SRL )
      AIC(best.glm.occur.8t.minus7)
         # 66.14
      drop1(best.glm.occur.8t.minus7)# dropping Max.Root.Loca:Log.Leaf.Area decreases AIC by 0.5
      best.glm.occur.8t.minus8<-update(best.glm.occur.8t.minus7,~.-Max.Root.Loca:Log.Leaf.Area)
      AIC(best.glm.occur.8t.minus8)
         # 66.17
      drop1(best.glm.occur.8t.minus8) # dropping Max.Root.Loca decreases AIC by 1.7
      best.glm.occur.8t.minus9<-update(best.glm.occur.8t.minus8,~.-Max.Root.Loca)
      AIC(best.glm.occur.8t.minus9)
         #64.54
      drop1(best.glm.occur.8t.minus9) # dropping Log.Leaf.Area increases AIC by only 0.2
      best.glm.occur.8t.minus10<-update(best.glm.occur.8t.minus9,~.-Log.Leaf.Area)
      AIC(best.glm.occur.8t.minus10)
         # 65.54
      drop1(best.glm.occur.8t.minus10) # No
   
      
      # lowest global AIC is minus 9
      summary(best.glm.occur.8t.minus9)
      # Coefficients:
      # Estimate Std. Error t value Pr(>|t|)  
      # (Intercept)    0.33767    0.16912   1.997   0.0542 .
      # Log.Leaf.Area -0.06205    0.04252  -1.459   0.1539  
      # log(Myc.Frac) -0.63539    0.33566  -1.893   0.0672 .
      # ---
      # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # 
      # (Dispersion parameter for Gamma family taken to be 0.1901384)
      # 
      # Null deviance: 7.4065  on 35  degrees of freedom
      # Residual deviance: 5.0455  on 33  degrees of freedom
      # AIC: 64.544
      # 
      # Number of Fisher Scoring iterations: 8
      
      1-(deviance(best.glm.occur.8t.minus9)/best.glm.occur.8t.minus9$null.deviance)
      # 0.3187757
      
      AIC(best.glm.occur.8t.minus9)
         # 64.54
      AICc(best.glm.occur.8t.minus9)
         # 65.83

      # Validate model
      par(mfrow=c(2,2));plot(best.glm.occur.8t.minus9,which = c(1,3:5))
      # LYAN is an outlier with too much leverage 
      
      # Plot residuals vs each variable (in the model)
      plot((resid(best.glm.occur.8t.minus9))~dat.8t$Myc.Frac)
      plot((resid(best.glm.occur.8t.minus9))~dat.8t$Log.Leaf.Area)
      
      # Plot residuals vs each variable (NOT in the model)
      plot((resid(best.glm.occur.8t.minus9))~dat.8t$Log.Ht.veg)
      plot((resid(best.glm.occur.8t.minus9))~dat.8t$Max.Root.Loca)
      plot((resid(best.glm.occur.8t.minus9))~dat.8t$SRL)
      plot((resid(best.glm.occur.8t.minus9))~dat.8t$Leaf.Mass.Frac)
      plot((resid(best.glm.occur.8t.minus9))~dat.8t$Log.Lamina.thck)
      plot((resid(best.glm.occur.8t.minus9))~dat.8t$LDMC)

      
         #===========SUMMARY OF GLMS============================================#
         # Best model includes both Myc.Frac and log(Myc.Frac) and other variables
         # only glm with both myc.frac + log(Myc.Frac) select same variables as regression tree 
         # Pseudo-R2= 61%, AIC= 53
         # Models with Myc.Frac only and log(Myc.frac) retain different variables
         # Myc.Frac: Pseudo-R2=32%, AIC=66 ; log(Myc.Frac)Pseudo-R2=0.31, AIC=64
         # Regression tree retains Myc.Frac + SRL, pseudo-R2 = 32%
         #======================================================================#

# B4 - Rerun glm and tree models without outliers (LYAN)
      
      
#==============================#
# (C) Elevation vs 8 Traits ####
#==============================#
      
   # C1 - Explore Trait Interactions ####
   #==============================#
      # C1.1 - Tree Model
      par(mfrow=c(1,1))
      tree.model<-tree(dat.8t$ElevDif~.,data=dat.8t[,Trait.Names.8t])
      plot(tree.model); text(tree.model); title('Regression Tree \n Elevation Shift vs 8 traits')
      # Leaf.Mass.Frac is best predictor,
      # for those with high Leaf.Mass.Frac, Max.Root.Loca is best predictor
      
      tree.model
         # 1) root 32 92210  36.94  
         #     2) Leaf.Mass.Frac < 0.302389 5  7741 -25.60 *
         #     3) Leaf.Mass.Frac > 0.302389 27 61290  48.52  
         #        6) Max.Root.Loca < 2.15 17 30940  69.41  
         #           12) Log.Ht.veg < 2.83148 9 13630  89.78 *
         #           13) Log.Ht.veg > 2.83148 8  9376  46.50 *
         #        7) Max.Root.Loca > 2.15 10 10310  13.00  
         #           14) Max.Root.Loca < 2.60256 5  1556  -6.00 *
         #           15) Max.Root.Loca > 2.60256 5  5148  32.00 *
      
      # Pseudo-R2
      1-(deviance(tree.model)/92210)
      # 0.59832
      plot(cv.tree(tree.model))
      
      # c.v.  splits the data into 'training' set for model fitting and a 
      # validation set to evaluate goodness of fit 
      # look for how many splits produces the minimum deviance - run multiple times
      
      # unstable. out of 30 runs, select4-5 splits 15 times
      # 3-4 splits 7 times
      # 1-2 splits 6 times
      # 2-3 splits 1 time
      
      # stop at 3 or 4 branches?
      
         #=====================================================================================#
         # Summary of regression tree
         # Regression Tree retains:
         #  Leaf.Mass.Frac, 
         #  Max.Root.Loca, 
         #  Log.Ht.veg
         #  and their interactions - R2 = 0.59
         #=====================================================================================#
      
      # C1.2 - Test trait-trait interactions 
      ----
         
      # need to check in subgroups because sample size too small to test all interactions
      # simultaneously
      summary(lm(dat.8t$ElevDif~(Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC)^2,
                 data=dat.8t))
      
      # No significant terms
      
      summary(lm(dat.8t$ElevDif~(Log.Leaf.Area+Leaf.Mass.Frac+SRL+Myc.Frac)^2,
                 data=dat.8t))
      
      # Log.Leaf.Area:Myc.Frac     @p=0.0039
      
      summary(lm(dat.8t$ElevDif~(Log.Ht.veg+Max.Root.Loca+Log.Leaf.Area+Leaf.Mass.Frac)^2,
                 data=dat.8t))
      
      # Max.Root.Loca:Log.Ht.veg  @p=0.0895
 
      summary(lm(dat.8t$ElevDif~(Log.Ht.veg+Max.Root.Loca+SRL+Myc.Frac)^2,
                 data=dat.8t))
      
      # Max.Root.Loca:Log.Ht.veg @p=0.0704
      
         #=====Summary of trait-trait interactions =======#
         # include in model: 
         #        Log.Leaf.Area:Myc.Frac        significant
         #        Max.Root.Loca:Log.Ht.veg      marginally significant
         #====================================#
      
      #=====================================================================================#
      # Summary of trait-trait interactions. 
      # Main interactions in regression tree does not appear in glm. Potential candidates are:
      #     Leaf.Mass.Frac:Max.Root.Loca, 
      #     Max.Root.Loca:Log.Ht.veg, 
      #     Leaf.Mass.Frac:Max.Root.Loca:Log.Ht.veg
      #     Log.Leaf.Area:Myc.Frac        
      #     Max.Root.Loca:Log.Ht.veg      
      # ====================================================================================#
  
   # C2 - Explore non-linearity ####
   #==============================#
      
      # C2.1 GAMs
      ---
      par(mfrow=c(2,2),mar=c(4,4,3,3))
      plot(gam(ElevDif~s(Log.Ht.veg)+s(SRL)+s(Myc.Frac),data=dat.8t))
         # SRL is non-linear
         # Log.Ht.veg is wonky
      
      plot(gam(ElevDif~s(Max.Root.Loca)+s(Log.Lamina.thck)+
                  s(LDMC),data=dat.8t))
         # Log.Lamina.thck is non-linear
      
      plot(gam(ElevDif~s(Log.Leaf.Area)+s(Leaf.Mass.Frac),data=dat.8t))
      
      # C2.2 Test polynomial relationships 
      ----
      summary(lm(ElevDif~I(Log.Ht.veg^2)+I(Max.Root.Loca^2)+I(Log.Lamina.thck^2)+
                    I(LDMC^2)+I(Log.Leaf.Area^2)+I(Leaf.Mass.Frac^2)+I(SRL^2)+I(Myc.Frac^2),
                 data=dat.8t))
      # I(Leaf.Mass.Frac^2) significant @p=0.0003
      # I(Max.Root.Loca^2)  significant @p=0.02
      # I(LDMC^2)           significant @p=0.02
      # I(SRL^2)            significant @p=0.02
      # I(Log.Leaf.Area^2)  significant @p=0.04
      
   # C3 - Model selection ####
   #==============================#
      
      # C3.1 - lm full ####
      #=========================================#
      # Start with lm bc values range from negative to positive inf. 
      
      # Includes all 1st order terms 
      # + 5 significant non-linear variables
      # + 5 significant trait interactions
      #====================================#
      
   H1<-lm(ElevDif~
      Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac+SRL+Myc.Frac+
         I(Leaf.Mass.Frac^2) + I(Max.Root.Loca^2) + I(SRL^2) + I(LDMC^2) + I(Log.Leaf.Area^2)+
         Log.Leaf.Area:Myc.Frac + Max.Root.Loca:Log.Ht.veg  +
         Leaf.Mass.Frac:Max.Root.Loca + Max.Root.Loca:Log.Ht.veg + Leaf.Mass.Frac:Max.Root.Loca:Log.Ht.veg,
      data=dat.8t)


   H0<-lm(ElevDif~1,
           data=dat.8t)

   best.lm.elev.8t<-step(H1,scope=list(lower=H0,upper=H1),direction='backward',trace = T)
   summary(best.lm.elev.8t)
      #AdjR2 = 0.65
   AIC(best.lm.elev.8t)
      # 325.00
   
   # Can the model be simplified?
   drop1(best.lm.elev.8t) # dropping I(SRL^2) only increases AIC by 0.6
   best.lm.elev.8t.minus1<-update(best.lm.elev.8t,~.-I(SRL^2) )
   AIC(best.lm.elev.8t.minus1) 
      # 325.61
   drop1(best.lm.elev.8t.minus1) # dropping Log.Lamina.thck increases AIC by 1.4
   best.lm.elev.8t.minus2<-update(best.lm.elev.8t.minus1,~.-Log.Lamina.thck )
   AIC(best.lm.elev.8t.minus2)
      # 327.00
   drop1(best.lm.elev.8t.minus2) # dropping SRL increases AIc by 0.6
   best.lm.elev.8t.minus3<-update(best.lm.elev.8t.minus2,~.-SRL)
   AIC(best.lm.elev.8t.minus3)
      # 327.53
   drop1(best.lm.elev.8t.minus3) # dropping Log.Ht.veg:Max.Root.Loca:Leaf.Mass.Frac decreases AIC by 0.1
   best.lm.elev.8t.minus4<-update(best.lm.elev.8t.minus3,~.-Log.Ht.veg:Max.Root.Loca:Leaf.Mass.Frac)
   AIC(best.lm.elev.8t.minus4)
      # 327.41
   drop1(best.lm.elev.8t.minus4) # dropping Max.Root.Loca:Leaf.Mass.Frac decreases AIC by 1.4 
   best.lm.elev.8t.minus5<-update(best.lm.elev.8t.minus4,~.-Max.Root.Loca:Leaf.Mass.Frac)
   AIC(best.lm.elev.8t.minus5)
      # 325.68
   
   drop1(best.lm.elev.8t.minus5) # dropping I(Log.Leaf.Area^2) increases AIC by 0.4
   best.lm.elev.8t.minus6<-update(best.lm.elev.8t.minus5,~.-I(Log.Leaf.Area^2))
   AIC(best.lm.elev.8t.minus6)
      # 326.0845 
   
      # all parameters of .minus6 are statistically significant
      # globally lowest AIC = best.lm.elev.8t 
   
   # Best model (globally lowest AIC)
   summary(best.lm.elev.8t)
      # Coefficients:
      #                                          Estimate Std. Error t value Pr(>|t|)    
      # (Intercept)                              514.14674  207.56027   2.477 0.023394 *  
      # Log.Ht.veg                                49.71068   28.88259   1.721 0.102370    
      # Max.Root.Loca                            234.42043  113.02893   2.074 0.052703 .  
      # Log.Lamina.thck                          -43.09295   27.10131  -1.590 0.129229    
      # Log.Leaf.Area                           -210.54889   48.24919  -4.364 0.000374 ***
      # Leaf.Mass.Frac                           217.82081   65.82081   3.309 0.003900 ** 
      # SRL                                       -5.41214    3.18663  -1.698 0.106653    
      # Myc.Frac                                -547.70525  165.30436  -3.313 0.003866 ** 
      # I(SRL^2)                                   0.08798    0.07122   1.235 0.232625    
      # I(Log.Leaf.Area^2)                         8.13330    4.39313   1.851 0.080595 .  
      # Log.Leaf.Area:Myc.Frac                   239.21141   63.77892   3.751 0.001464 ** 
      # Log.Ht.veg:Max.Root.Loca                 -80.57210   37.00948  -2.177 0.043027 *  
      # Max.Root.Loca:Leaf.Mass.Frac            -277.76071  181.92097  -1.527 0.144186    
      # Log.Ht.veg:Max.Root.Loca:Leaf.Mass.Frac   78.31334   54.71178   1.431 0.169456    
      # ---
      #    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      # 
      # Residual standard error: 32.4 on 18 degrees of freedom
      # (4 observations deleted due to missingness)
      # Multiple R-squared:  0.795,	Adjusted R-squared:  0.647 
      # F-statistic: 5.371 on 13 and 18 DF,  p-value: 0.0006709
   
   # Validation plots
   plot(best.lm.elev.8t)
      # # All looks good, except that CLBO has a large leverage (>0.5, < 1)
   
   # minimal model - smallest with 2 AIC of best global model
   summary(best.lm.elev.8t.minus6)
      # Coefficients:
      #                          Estimate Std. Error t value Pr(>|t|)    
      # (Intercept)                237.10     109.28   2.170 0.040152 *  
      # Log.Ht.veg                  61.94      21.01   2.949 0.007009 ** 
      # Max.Root.Loca               74.04      34.32   2.157 0.041205 *  
      # Log.Leaf.Area             -155.61      40.77  -3.817 0.000836 ***
      # Leaf.Mass.Frac             164.92      32.95   5.005  4.1e-05 ***
      # Myc.Frac                  -567.17     143.27  -3.959 0.000585 ***
      # Log.Leaf.Area:Myc.Frac     207.15      51.01   4.061 0.000451 ***
      # Log.Ht.veg:Max.Root.Loca   -32.98      11.59  -2.845 0.008944 ** 
      # ---
      # Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      # 
      # Residual standard error: 34.42 on 24 degrees of freedom
      # (4 observations deleted due to missingness)
      # Multiple R-squared:  0.6916,	Adjusted R-squared:  0.6016 
      # F-statistic: 7.688 on 7 and 24 DF,  p-value: 6.754e-05

   #library(Mumin)
   library(MuMIn)
   AICc(best.lm.elev.8t)
   # 355.01
   AICc(best.lm.elev.8t.minus6)
   # 334.26
   
   # Validate model
   par(mfrow=c(2,2)); plot(best.lm.elev.8t.minus6)
   # Looks awesome
   
   AQUI
   # Plot residuals vs each variable (in the model)
   plot((resid(best.lm.elev.8t.minus6))~dat.8t[!is.na(dat.8t$ElevDif),'Log.Ht.veg'])
   plot((resid(best.lm.elev.8t.minus6))~dat.8t[!is.na(dat.8t$ElevDif),'Max.Root.Loca'])
   plot((resid(best.lm.elev.8t.minus6))~dat.8t[!is.na(dat.8t$ElevDif),'Log.Leaf.Area'])
   plot((resid(best.lm.elev.8t.minus6))~dat.8t[!is.na(dat.8t$ElevDif),'Leaf.Mass.Frac'])
   plot((resid(best.lm.elev.8t.minus6))~dat.8t[!is.na(dat.8t$ElevDif),'Myc.Frac'])
 
   
   # Plot residuals vs each variable (NOT in the model)
   plot((resid(best.lm.elev.8t.minus6))~dat.8t[!is.na(dat.8t$ElevDif),'Log.Lamina.thck'])
   plot((resid(best.lm.elev.8t.minus6))~dat.8t[!is.na(dat.8t$ElevDif),'LDMC'])
   plot((resid(best.lm.elev.8t.minus6))~dat.8t[!is.na(dat.8t$ElevDif),'SRL'])

   # Try with MuMIn{} dredge()
   # rewrite H1 with na.action='fail' - required by dredge (), add dat8[!is.na(ElevDif),] to remove rows with NA
   dat.8t.nona<-dat.8t[!is.na(dat.8t$ElevDif),]
   H1<-lm(ElevDif~
             Log.Ht.veg+Max.Root.Loca+Log.Lamina.thck+LDMC+Log.Leaf.Area+Leaf.Mass.Frac+SRL+Myc.Frac+
             I(Leaf.Mass.Frac^2) + I(Max.Root.Loca^2) + I(SRL^2) + I(LDMC^2) + I(Log.Leaf.Area^2)+
             Log.Leaf.Area:Myc.Frac + Max.Root.Loca:Log.Ht.veg  +
             Leaf.Mass.Frac:Max.Root.Loca + Max.Root.Loca:Log.Ht.veg + Leaf.Mass.Frac:Max.Root.Loca:Log.Ht.veg,
          data=dat.8t.nona,
          na.action='na.fail')
   
   d.abund<-dredge(H1,beta='sd',rank=AIC,trace=F)
   head(d.abund)
   # 40 models within 2 AIC points... 
   
   # best model includes (AIC = 325):
   # Lef.Mss.Frc 
   # Log.Ht.veg 
   # Log.Lmn.thc 
   # Log.Lef.Are 
   # Log.Lef.Are^2 
   # Max.Rot.Loc
   # Myc.Frc    
   # SRL  
   # SRL^2 
   # Log.Ht.veg:Max.Rot.Loc
   # Log.Lef.Are:Myc.Frc
   
   d.abund.aicc<-dredge(H1,beta='sd',rank=AICc,trace=F)
   head(d.abund.aicc)
   # 7 models within 3 AICc points of each other 

   # best model includes (AIC = 333.8):
   # LDMC
   # Lef.Mss.Frc
   # Log.Lef.Are^2
   # Max.Rot.Loc^2
   # SRL
   # SRL^2
   
#   ====================================

     
      
  # Diagnostic plots
  # 1) Residuals vs Fitted - detects residual non-linear relationships between x & y variables 
  #    and heteroscedasticity (wedge)
  # 2) QQ plot - shows whether residuals are normally distributed (will follow straight line)
  # 3) Scale-Location - shows whether variance (residuals) increase with incrasing mean (also tests
  # for homoscedasticity). Another version of plot 1
  # 4) Residuals vs Leverage - shows whether individual datapoints have a lot of 'weight' on 
  # the regression parameters
  
  
  
  # fitting Gamma with dispersion=1 is supposed to be equivalent to fitting lognormal. 
  # see: http://stats.stackexchange.com/questions/21447/how-to-specify-a-lognormal-distribution-in-the-glm-family-argument-in-r
  # see: https://stats.stackexchange.com/questions/67547/when-to-use-gamma-glms
  # see: http://seananderson.ca/2014/04/08/gamma-glms.html
  

  # Null deviance = total variance & model deviance = variance explained by model. hence, 

  