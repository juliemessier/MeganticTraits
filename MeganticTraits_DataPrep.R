# <<TABLE OF CONTENTS>>
# A - Setup trait data
  # A1 - Herbivory layer
    # A1.1 - Remove date effect on traits
    # A1.2 - Turn individual-level table into species-level table.
    # A1.3 - Add species mean traits (Seed size and mychorizae)
    # A1.4 - Data Exploration
    # A1.5 - Transform trait data as needed
    

  # A2 - Canopy layer
    # A2.1 - Transform trait data for normality
    # A2.2 - Remove date effect on traits
    # A2.3 - Turn individual-level table into species-level table.
    # A2.4 - Add species mean traits (Seed size and mychorizae)
# B- Setup abundance data
# B2- Create abundance dataframe & new variables
# B3 - Problem with infinity values ####

#<<WORKSPACES>>
wrk.dir<-("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Traits/Workspaces/") # Workspaces
data.dir<-(("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Traits/Data/")) # data
res.dir<-("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Traits/Results/")  # Results
grp.dir<-("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Traits/Graphs/")   # Graphs
fct.dir<-("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Traits/Functions/") # Functions

#<<LIBRARIES>>
library(car) # for powerTransform
library(vegan) # for decostand (standardizing data)
library(lattice) # For fancy multipanel graphs
source(paste0(wrk.dir,"HighstatLibV10.R")) # to make fancy graphs

# ==================================================================================#

# A - Setup trait data ####
# ========================#

  # A1 - HERBACEOUS LAYER ####
  #=======================#

    Megtraits<-read.csv(paste0(data.dir,'MegTraits_20181106.csv'))
    dim(Megtraits) #640 47
    str(Megtraits)
    Megtraits$Min.Root.Loca<-as.ordered(Megtraits$Min.Root.Loca) # remains a factor, but is ordered (ordinal)
    Megtraits$Max.Root.Loca<-as.ordered(Megtraits$Max.Root.Loca) # remains a factor, but is ordered (ordinal)
    
      # Simplify names
      names(Megtraits)[names(Megtraits)=="Lamina.thck.æm."] <- "Lamina.thck"
      names(Megtraits)[names(Megtraits)=='Vein.thck.æm.']<-'Vein.thck'
      names(Megtraits)[names(Megtraits)=='LMA.g.cm2.']<-'LMA'
      names(Megtraits)[names(Megtraits)=='LDMC.g.g.']<-'LDMC'
      names(Megtraits)[names(Megtraits)=='Leaf.Area.cm2.']<-'Leaf.Area'
      names(Megtraits)[names(Megtraits)=='Fine.Root.Diam.mm.']<-'F.Root.Diam'
      names(Megtraits)[names(Megtraits)=='SRL.cm.mg.']<-'SRL'
      
      # create herbaceous trait dataframe with traits I will work with
      
      H.traits<-Megtraits[Megtraits$Layer=='H',c("Plant.ID","Species","Layer","Date.recolte","Ht.veg","Min.Root.Loca",
                                                 "Max.Root.Loca","Lamina.thck","LMA","LDMC","Leaf.Area","Leaf.Mass.Frac",
                                                 "Supp.Mass.Frac","Rep.Mass.Frac","Stor.Mass.Frac","F.Root.Diam","SRL")]
      
      dim(H.traits)
      # 459 17
      H.traits<-droplevels(H.traits)
      names(H.traits)
      # [1] "Plant.ID"       "Species"        "Layer"          "Date.recolte"   "Ht.veg"         "Min.Root.Loca"  "Max.Root.Loca" 
      # [8] "Lamina.thck"    "LMA"            "LDMC"           "Leaf.Area"      "Leaf.Mass.Frac" "Supp.Mass.Frac" "Rep.Mass.Frac" 
      # [15] "Stor.Mass.Frac" "F.Root.Diam"    "SRL"
      
      str(H.traits)
      
        # 'data.frame':	459 obs. of  17 variables:
        # $ Plant.ID      : Factor w/ 459 levels "ARNU1","ARNU2",..: 1 2 3 4 5 6 7 8 9 10 ...
        # $ Species       : Factor w/ 51 levels "ARNU","ARTR",..: 1 1 1 1 1 1 1 1 1 2 ...
        # $ Layer         : Factor w/ 1 level "H": 1 1 1 1 1 1 1 1 1 1 ...
        # $ Date.recolte  : Factor w/ 49 levels "01-06-2016","01-08-2016",..: 47 49 49 49 1 1 1 1 24 44 ...
        # $ Ht.veg        : num  32 24 34.5 42 28 47.5 26 14 56 15.5 ...
        # $ Min.Root.Loca : Ord.factor w/ 6 levels "0"<"1"<"2"<"3"<..: 1 1 1 1 1 3 1 2 2 5 ...
        # $ Max.Root.Loca : Ord.factor w/ 6 levels "0"<"1"<"2"<"3"<..: 1 2 1 1 1 3 1 2 2 5 ...
        # $ Lamina.thck   : num  71.3 73.8 70.5 67.8 73.8 ...
        # $ LMA           : num  0.0012 0.0013 0.0011 0.0014 0.0015 0.0017 0.0011 0.0012 0.0015 0.0014 ...
        # $ LDMC          : num  0.16 0.182 0.153 0.175 0.193 0.203 0.181 0.18 0.194 0.102 ...
        # $ Leaf.Area     : num  18 14 21.2 37 18.3 ...
        # $ Leaf.Mass.Frac: num  0.21 0.3 0.09 0.3 0.5 0.41 0.59 0.22 0.29 0.42 ...
        # $ Supp.Mass.Frac: num  0.17 0.2 0.16 0.27 0.28 0.31 0.36 0.17 0.22 0.11 ...
        # $ Rep.Mass.Frac : num  0 0 0 0 0 0 0 0 0.11 0 ...
        # $ Stor.Mass.Frac: num  0.62 0.49 0.75 0.43 0.22 0.28 0.04 0.61 0.38 0.47 ...
        # $ F.Root.Diam   : num  NA 0.567 NA 0.74 0.548 ...
        # $ SRL           : num  NA 2.37 NA 1.42 2.38 ...      
      
      # Reponse Variables
      names(H.abund)
      head(H.abund)
      
      # A1.1 - Remove date effect on traits  ####
      # =======================================#
      
      names(H.traits)
      # [1] "Plant.ID"       "Species"        "Layer"          "Date.recolte"   "Ht.veg"        
      # [6] "Min.Root.Loca"  "Max.Root.Loca"  "Lamina.thck"    "LMA"            "LDMC"          
      # [11] "Leaf.Area"      "Leaf.Mass.Frac" "Supp.Mass.Frac" "Rep.Mass.Frac"  "Stor.Mass.Frac"
      # [16] "F.Root.Diam"    "SRL"
      
      Trait.Names <- c("Ht.veg","Min.Root.Loca","Max.Root.Loca","Lamina.thck","LMA","LDMC","Leaf.Area",
                       "Leaf.Mass.Frac","Supp.Mass.Frac","Rep.Mass.Frac","Stor.Mass.Frac","F.Root.Diam",
                       "SRL")
      
      Trait.Names         # list of trait names
      # [1] "Ht.veg"         "Min.Root.Loca"  "Max.Root.Loca"  "Lamina.thck"    "LMA"           
      # [6] "LDMC"           "Leaf.Area"      "Leaf.Mass.Frac" "Supp.Mass.Frac" "Rep.Mass.Frac" 
      # [11] "Stor.Mass.Frac" "F.Root.Diam"    "SRL"
      length(Trait.Names) #13
      
      save(Trait.Names,file=paste0(wrk.dir,'list.trait.names.herbivory.layer.Rdata'))
      
      # Change calendar dates into julian dates
      H.traits$Date.recolte<- as.numeric(format(as.Date(H.traits$Date.recolte, format = "%d-%m-%Y"),"%j"))


      # loop testing for date effects and replacing with regression residuals
      
      for (t in Trait.Names[-c(2:3)]){  # not for min. and max. root location, because they are factors
        x<-lm(H.traits[[t]]~H.traits$Date.recolte,na.action=na.exclude)  # # regress the trait against the julian date
        
        if (summary(x)$adj.r.squared > 0.02 & summary(x)[["coefficients"]][2,4] < 0.05) # if regression R2 > 0.02 AND it is statistically significant (with P-value<0.05), go through this next loop
        
          # show me which traits
            print(c(names(H.traits[t]),
                  paste( "R2= ",round(summary(x)$adj.r.squared,digits=3)),
                  paste('Sign= ',sign(summary(x)[["coefficients"]][2,1])))) 
        
          # replace with residuals
        # for (i in 1:nrow(H.traits[t])){                                                   # if each row is numeric (not an NA)
        #   if (is.numeric(H.traits[t][i,]))
        #   {H.traits[t][i,]<-resid(x)[i]}                                              # then replace the value with the residual
        # }                                                                           # else, do nothing to the cells with value of NA
      }
      
      # [1] "Ht.veg"     "R2=  0.163" "Sign=  1"  
      # [1] "Lamina.thck" "R2=  0.274"  "Sign=  -1"  
      # [1] "LDMC"       "R2=  0.143" "Sign=  1"  
      # [1] "Leaf.Area" "R2=  0.03" "Sign=  1" 
      # [1] "Leaf.Mass.Frac" "R2=  0.097"     "Sign=  1"      
      # [1] "Supp.Mass.Frac" "R2=  0.089"     "Sign=  1"      
      # [1] "Stor.Mass.Frac" "R2=  0.215"     "Sign=  -1" 
      
      # Okay, but look visually to see if these relationships are (a) linear and (b) due to species-date association? 
      
      # (1) Ht.veg - # NO DATE EFFECT
      ---
      plot(Ht.veg~Date.recolte, data=H.traits) 
      # funnel pattern, but is this due to species effect ?
      
      scatterplot(Ht.veg~Date.recolte|Species, data=H.traits, legend=list(cex=0.8),
                  xlab="Sampling Date",
                  ylab="Vegetative Height",
                  smooth=F,
                  col=c(1:51))
      
        # Some taller species sampled later
        # Ht.veg doesn't appear to be correlated with sampling date within species.
      
      m0<-lm(Ht.veg~Species, data=H.traits)
      summary(m0)
      
        # Residual standard error: 13.1 on 406 degrees of freedom
        # (2 observations deleted due to missingness)
        # Multiple R-squared:  0.8664,	Adjusted R-squared:  0.8499 
        # F-statistic: 52.65 on 50 and 406 DF,  p-value: < 2.2e-16
      
      # Height varies by species
      
      # Is the Sampling Date parameter significant, once we account for Species?
      
      m1<-lm(Ht.veg~Species+Date.recolte, data=H.traits)
      summary(m1)
      #Date.recolte is not significant in this model ! 
      
      anova(m0,m1)
      # m1 not significantly better. 
      
      # NO DATE EFFECT
      
      
      # (2) Lamina.thck - NO DATE EFFECT
      ---
      plot(Lamina.thck~Date.recolte, data=H.traits) 
      # negative exp. pattern, but is this due to species effect ?
      
      scatterplot(Lamina.thck~Date.recolte|Species, data=H.traits, legend=list(cex=0.8),
                  xlab="Sampling Date",
                  ylab="Lamina thickness",
                  smooth=F,
                  col=c(1:51))
      
      # Some thicker species sampled earlier
      # Lamina.thck doesn't appear to be correlated with sampling date within species.
      
      summary(m0<-lm(Lamina.thck~Species, data=H.traits))
        # Multiple R-squared:  0.9328,	Adjusted R-squared:  0.9244 
        # F-statistic: 111.1 on 50 and 400 DF,  p-value: < 2.2e-16
      
      # Lamina.thck varies by species
      
      # Is the Sampling Date parameter significant, once we account for Species?
      
      summary(m1<-lm(Lamina.thck~Species+Date.recolte, data=H.traits))
        # Multiple R-squared:  0.9328,	Adjusted R-squared:  0.9242 
        # F-statistic: 108.7 on 51 and 399 DF,  p-value: < 2.2e-16
      
      # Sampling date parameter not significant in this model. 
      anova (m0,m1)
      # Not significant
      
      # NO DATE EFFECT
        
      # (3) LDMC - significant sampling effect, but effect size is 0.3%  
      ---
        
      plot(LDMC~Date.recolte, data=H.traits) 
      # positive pattern, but is this due to species effect ?
      
      scatterplot(LDMC~Date.recolte|Species, data=H.traits, legend=list(cex=0.8),
                  xlab="Sampling Date",
                  ylab="LDMC",
                  smooth=F,
                  col=c(1:51))
      
      summary(m0<-lm(LDMC~Species, data=H.traits))
      # Multiple R-squared:  0.9259,	Adjusted R-squared:  0.9168 
      # F-statistic:   102 on 50 and 408 DF,  p-value: < 2.2e-16
      
      # LDMC varies by species
      
      # Is the Sampling Date parameter significant, once we account for Species?
      
      summary(m1<-lm(LDMC~Species+Date.recolte, data=H.traits))
      # Date.recolte  0.0005974  0.0001391   4.296 2.18e-05 ***
      # ---
      # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # 
      # Residual standard error: 0.02767 on 407 degrees of freedom
      # Multiple R-squared:  0.9291,	Adjusted R-squared:  0.9203 
      # F-statistic: 104.6 on 51 and 407 DF,  p-value: < 2.2e-16
      
      # Sampling date parameter IS significant in this model. 
      anova (m0,m1)
      as.numeric(RsquareAdj(m1)[2])-as.numeric(RsquareAdj(m0)[2])
      # [1] 0.003411269
      
      # Yes, better with Date, BUT, R2 explained by sampling date is 0.0035
        
      # (4) Leaf.Area
      ---
      
      plot(Leaf.Area~Date.recolte, data=H.traits) 
      # no visible pattern
      
      scatterplot(Leaf.Area~Date.recolte|Species, data=H.traits, legend=list(cex=0.8),
                  xlab="Sampling Date",
                  ylab="Leaf.Area",
                  smooth=F,
                  col=c(1:51))
      
      summary(m0<-lm(Leaf.Area~Species, data=H.traits))
      # Multiple R-squared:  0.7494,	Adjusted R-squared:  0.7166 
      # F-statistic: 22.88 on 49 and 375 DF,  p-value: < 2.2e-1
      
      # Leaf.Area varies by species
      
      # Is the Sampling Date parameter significant, once we account for Species?
      
      summary(m1<-lm(Leaf.Area~Species+Date.recolte, data=H.traits))
      # Multiple R-squared:  0.7496,	Adjusted R-squared:  0.7162 
      # F-statistic:  22.4 on 50 and 374 DF,  p-value: < 2.2e-16
      
      # Sampling date parameter not significant in this model. 
      anova (m0,m1)
      
      # Sampling date does not explain any variance after species effect taken into account.
        
      # (5) Leaf.Mass.Frac - significant sampling effect, but effect size is 0.1%
      ---
      
      plot(Leaf.Mass.Frac~Date.recolte, data=H.traits) 
      # positive trend
      
      scatterplot(Leaf.Mass.Frac~Date.recolte|Species, data=H.traits, legend=list(cex=0.8),
                  xlab="Sampling Date",
                  ylab="Leaf.Mass.Fraction",
                  smooth=F,
                  col=c(1:51))
      
      summary(m0<-lm(Leaf.Mass.Frac~Species, data=H.traits))
      # Multiple R-squared:  0.8619,	Adjusted R-squared:  0.8449 
      # F-statistic:  50.7 on 50 and 406 DF,  p-value: < 2.2e-16
      
      
      # Leaf.Area varies by species
      
      # Is the Sampling Date parameter significant, once we account for Species?
      
      summary(m1<-lm(Leaf.Mass.Frac~Species+Date.recolte, data=H.traits))
      # Date.recolte -0.0009841  0.0004599  -2.140 0.032954 *  
      # ---
      # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # 
      # Residual standard error: 0.09108 on 405 degrees of freedom
      # (2 observations deleted due to missingness)
      # Multiple R-squared:  0.8635,	Adjusted R-squared:  0.8463 
      # F-statistic: 50.23 on 51 and 405 DF,  p-value: < 2.2e-16
      
      # Sampling date significant, but small effect size.
      
      # Sampling date parameter not significant in this model. 
      anova (m0,m1)
      as.numeric(RsquareAdj(m1)[2])-as.numeric(RsquareAdj(m0)[2])
      # [1] 0.001355086
          
      # (6) Supp.Mass.Frac - significant sampling effect, but effect size is 0.2%
      ---
      plot(Supp.Mass.Frac~Date.recolte, data=H.traits) 
      # positive trend
      
      scatterplot(Supp.Mass.Frac~Date.recolte|Species, data=H.traits, legend=list(cex=0.8),
                  xlab="Sampling Date",
                  ylab="Support Mass Fraction",
                  smooth=F,
                  col=c(1:51))
      
      summary(m0<-lm(Supp.Mass.Frac~Species, data=H.traits))
        # Multiple R-squared:  0.7578,	Adjusted R-squared:  0.7279 
        # F-statistic: 25.35 on 48 and 389 DF,  p-value: < 2.2e-16
      
      
      # Leaf.Area varies by species
      
      # Is the Sampling Date parameter significant, once we account for Species?
      
      summary(m1<-lm(Supp.Mass.Frac~Species+Date.recolte, data=H.traits))
        # Date.recolte -6.759e-04  3.188e-04  -2.120 0.034614 *  
        # ---
        # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
        # 
        # Residual standard error: 0.06263 on 388 degrees of freedom
        # (21 observations deleted due to missingness)
        # Multiple R-squared:  0.7606,	Adjusted R-squared:  0.7303 
        # F-statistic: 25.15 on 49 and 388 DF,  p-value: < 2.2e-16
        
      # Sampling date significant, but small effect size.
      
      # Sampling date parameter not significant in this model. 
      anova (m0,m1)
      as.numeric(RsquareAdj(m1)[2])-as.numeric(RsquareAdj(m0)[2])
      # [1] 0.00242344
      
      # (7) Stor.Mass.Frac - significant sampling effect, but effect size is 0.1%
      ---
      
      plot(Stor.Mass.Frac~Date.recolte, data=H.traits) 
      # negative trend
      
      scatterplot(Stor.Mass.Frac~Date.recolte|Species, data=H.traits, legend=list(cex=0.8),
                  xlab="Sampling Date",
                  ylab="Storage Mass Fraction",
                  smooth=F,
                  col=c(1:51))
      
      summary(m0<-lm(Stor.Mass.Frac~Species, data=H.traits))
        # Multiple R-squared:  0.9042,	Adjusted R-squared:  0.8927 
        # F-statistic: 78.65 on 42 and 350 DF,  p-value: < 2.2e-16
      
      # Leaf.Area varies by species
      
      # Is the Sampling Date parameter significant, once we account for Species?
      
      summary(m1<-lm(Stor.Mass.Frac~Species+Date.recolte, data=H.traits))
        # Date.recolte  0.0010184  0.0004688   2.173 0.030484 *  
        # ---
        # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
        # 
        # Residual standard error: 0.08477 on 349 degrees of freedom
        # (66 observations deleted due to missingness)
        # Multiple R-squared:  0.9055,	Adjusted R-squared:  0.8938 
        # F-statistic: 77.75 on 43 and 349 DF,  p-value: < 2.2e-16
        
      # Sampling date significant, but small effect size.
      
      # Sampling date parameter not significant in this model. 
      anova (m0,m1)
      as.numeric(RsquareAdj(m1)[2])-as.numeric(RsquareAdj(m0)[2])  
      # [1] 0.001128434
        
      save(H.traits,file=paste0(wrk.dir,"Herbaceous.Layer.Traits.Unstandardized.RData"))
      
      load(paste0(wrk.dir,"Herbaceous.Layer.Traits.Unstandardized.RData"))
      
      # Standardize all variables
      H.traits.stand<-decostand(H.traits[,c(5,8:17)],method='standardize', margin=2,na.rm=TRUE)
      save(H.traits.stand,file=paste0(wrk.dir,"Herbaceous.Layer.Traits.Standardized.RData"))
      
      # A1.2 - Turn individual-level table into species-level table  ####
      #=================================================================#    
      #apply(H.traits[,5:11], 2, function(x) tapply(x, H.traits$Species, mean,na.rm=T))
      
      H.traits$Min.Root.Loca<-as.numeric(H.traits$Min.Root.Loca) # can't take mean or median of ordinal var. 
      H.traits$Max.Root.Loca<-as.numeric(H.traits$Max.Root.Loca) # can't take mean or median of ordinal var. 
      
      H.traits.sp<-as.data.frame(aggregate(H.traits[,5:17],by=list(H.traits$Species),mean,na.rm=T))
      head(H.traits.sp)
      #H.traits.sp[,2:8]<-round(sp.H.traits[2:8],digits=3)
      H.traits.sp[is.na(H.traits.sp)] <-NA # replace NaN with Na
      rownames(H.traits.sp)<-H.traits.sp$Group.1 
      colnames(H.traits.sp)[1]<-'Species'
      #H.traits.sp<-H.traits.sp[,-1]
      H.traits.sp<-H.traits.sp[order(rownames(H.traits.sp)),]
      dim(H.traits.sp) # 51  14
      
      save(H.traits.sp,file=paste0(wrk.dir,"Species-Level.Herbaceous.Layer.Traits.RData"))
      
      # A1.3 - Add species mean traits (Seed size and mychorizae) ####
      #=================================================================#
      
      load(paste0(wrk.dir,'Species-Level.Herbaceous.Layer.Traits.RData')) #H.traits.sp
      dim(H.traits.sp)
        #51 14
      
      myc<-read.csv(paste0(data.dir,'myc.csv'))
      head(myc)
      dim(myc) #43 2
      rownames(myc)<-myc$species
      myc<-myc[-1]
      #myc<-decostand(myc,method='standardize',margin=2)
      dim(myc) #43 1
      names(myc)<-"Myc.Frac"
      # rows TRGR and TNCO don't exist. Replace TRGR with TRER. Delete TNCO bc we already have a value 
      # for TICO
      row.names(myc)[which(row.names(myc)=='TRGR')]<-'TRER'
      which(row.names(myc)=='TNCO')
        #39
      myc<-myc[-c(39),,drop = FALSE]
      dim(myc)
        # 42 1
      
      H.traits2.sp<-merge(H.traits.sp,myc, by="row.names",all=T)
      head(H.traits2.sp)
      H.traits2.sp$Row.names<-NULL
      rownames(H.traits2.sp)<-H.traits2.sp$Species
      
      dim(H.traits2.sp) # 51 15
      
      save(H.traits2.sp,file=paste0(wrk.dir,'Species-Level.Herbaceous.Layer.Traits.withMycFrac.Rdata'))
      
      # Add in Seed size some other time
      
      # seed.file<-read.delim('C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Traits/Data/Seed size/2950_05042017080522/Seed.size.2950.txt')
      # dim(seed.file)
      # View (seed.file)
      
      # A1.4 - Data Exploration (Following Highlands Stats course) ####
      #=================================================================#
      
        #A1.4.1 - Outliers on Y and X 
        
            #Y1 - Abundance - 1 outliers (CASC), variance homogeneous
            ---
            par(mfrow = c(1, 2))  
            boxplot(H.abund$abund.ratio,
                    main = "Abundance Ratio")
            dotchart(H.abund$abund.ratio, #i.e. cleveland dot chart
                     labels=rownames(H.abund),
                     xlab='range of data', ylab='Order of the data', main='abundance ratios') # 1 outlier

            #Y2 - Elevation
            ---
            # ... to do
            
            #X1 - Ht.veg - 2 outlier species - COCO and VIAL are tall (3 & 6 sd)
            ---
            par(mfrow=c(1,2))
            boxplot(H.traits2.sp$Ht.veg,
                    main = "vegetative height")
            dotchart(H.traits2.sp$Ht.veg, #i.e. cleveland dot chart
                     xlab='range of data', ylab='Order of the data',main='Ht.veg',
                     labels=H.traits2.sp$Species)
            
            # X2 - Min.Root.Loca - no outliers - only 6 categories. 7 species with NAs
            ---
            boxplot(H.traits2.sp$Min.Root.Loca,
                    main = "Minimum Root Location")
            
            H.traits2.sp[,c('Species','Min.Root.Loca')] 
              # All VIAL, SAPU, RUID, LOCA, COCO & COAL are 'NA'
            
            dotchart(as.numeric(H.traits2.sp$Min.Root.Loca), #i.e. cleveland dot chart
                     xlab='range of data',
                     ylab='Order of the data',
                     main='Min.root.Loca',
                     labels=H.traits2.sp$Species)
            
            plot(as.numeric(H.traits2.sp$Max.Root.Loca)~as.numeric(H.traits2.sp$Min.Root.Loca))
            # pretty well correlated, but not fully redundant.
            
            # X3 - Max.Root.Loca - no outliers - only 6 categories. 7 species with NAs
            ---
            boxplot(H.traits2.sp$Max.Root.Loca,
                    main = "Maximum Root Location")
            
            H.traits2.sp[,c('Species','Max.Root.Loca')]  
              # All VIAL, SAPU, RUID, LOCA, COCO & COAL are 'NA'
            
            # X4 - Lamina.thck - CLBO, GAPR and ERAM are over 2 SD thicker than other species
            ---
            boxplot(H.traits2.sp$Lamina.thck,
                    main = "Lamina.thck") # a few outliers
            
            dotchart(H.traits2.sp$Lamina.thck, #i.e. cleveland dot chart
                     xlab='range of data',
                     ylab='Order of the data',
                     main='Lamina.thickness',
                     labels=H.traits2.sp$Species)  # two species are very thick
            
            dotchart(log(H.traits2.sp$Lamina.thck), #i.e. cleveland dot chart
                     xlab='range of data',
                     ylab='Order of the data',
                     main='Lamina.thickness',
                     labels=H.traits2.sp$Species)
            
            dotchart(decostand(H.traits2.sp$Lamina.thck,method='standardize',margin=2),
                        xlab='range of data',
                        ylab='Order of the data',
                        main='Lamina.thickness',
                        labels=H.traits2.sp$Species)
            
            # X5 - LMA - GAPR and LYOB have LMA over 2 SD higher than other species
            ---
            boxplot(H.traits2.sp$LMA,
                    main = "LMA") # 3 outliers
            
            dotchart(H.traits2.sp$LMA, #i.e. cleveland dot chart
                     xlab='range of data',
                     ylab='Order of the data',
                     main='LMA',
                     labels=H.traits2.sp$Species) # two outlier species
            
            H.traits2.sp[(H.traits2.sp$Species=='SAPU'|H.traits2.sp$Species=='COAL'),'LMA']
            
            # X6 - LDMC - okay
            boxplot(H.traits2.sp$LDMC,
                    main = "LDMC")
            
            dotchart(H.traits2.sp$LDMC, #i.e. cleveland dot chart
                     xlab='range of data',
                     ylab='Order of the data',
                     main='LDMC') # fine
            
            # X7 - Leaf.Area # VEVI and SAPU leaf area over 2 SD larger than other sp. 
            ---
              
            boxplot(H.traits2.sp$Leaf.Area,
                    main = "Leaf.Area") # 5 outliers
            
            dotchart(H.traits2.sp$Leaf.Area, #i.e. cleveland dot chart
                     xlab='range of data',
                     ylab='Order of the data',
                     labels=H.traits2.sp$Species,
                     main='Leaf.Area') # one species outlier
            
            # X8 - Leaf.Mass.Frac - okay
            ---
            boxplot(H.traits2.sp$Leaf.Mass.Frac,
                    main = "Leaf.Mass.Frac") # fine
            
            dotchart(H.traits2.sp$Leaf.Mass.Frac, #i.e. cleveland dot chart
                     xlab='range of data',
                     ylab='Order of the data',
                     labels=H.traits2.sp$Species,
                     main='Leaf.Mass.Frac') # fine
            
            # X9 - Supp.Mass.Frac - Okay
            ---
            boxplot(H.traits2.sp$Supp.Mass.Frac,
                    main = "Supp.Mass.Frac") # fine
            
            dotchart(H.traits2.sp$Supp.Mass.Frac, #i.e. cleveland dot chart
                     xlab='range of data',
                     ylab='Order of the data',
                     labels=H.traits2.sp$Species,
                     main='Supp.Mass.Frac')
            
            # X10 - Rep.Mass.Frac # skewed distribution (lots of very low values), CYAC and CASC over 2SD larger Rep.Mass.Frac than other sp.
            ---
            boxplot(H.traits2.sp$Rep.Mass.Frac,
                    main = "Rep.Mass.Frac") # 0 inflated / Skewed
            
            dotchart(H.traits2.sp$Rep.Mass.Frac, #i.e. cleveland dot chart
                     xlab='range of data',
                     ylab='Order of the data',
                     labels=H.traits2.sp$Species,
                     main='Rep.Mass.Frac') # 0 inflated
            
            # X11 - Stor.Mass.Frac - 0 inflated at species level (turn into presence/absence?)/ Skewed
            ---
            boxplot(H.traits2.sp$Stor.Mass.Frac,
                    main = "Rep.Mass.Frac") # 0 inflated / Skewed
            
            dotchart(H.traits2.sp$Stor.Mass.Frac, #i.e. cleveland dot chart
                     xlab='range of data',
                     ylab='Order of the data',
                     labels=H.traits2.sp$Species,
                     main='Stor.Mass.Frac') # 0 inflated - no outliers
            
            # X12 - F.Root.Diam - CYAC & EPHE have fine roots 6 SD thicker than other spp.
            ---
            boxplot(H.traits2.sp$F.Root.Diam,
                    main = "F.Root.Diam")
            
            dotchart(H.traits2.sp$F.Root.Diam, #i.e. cleveland dot chart
                     xlab='range of data',
                     ylab='Order of the data',
                     labels=H.traits2.sp$Species,
                     main="Fine Root thickness") # 2 sp with large diameter
            
            # X13 - SRL - okay
            ---
            boxplot(H.traits2.sp$SRL,
                    main = "SRL") # 2 outliers
            
            dotchart(H.traits2.sp$SRL, #i.e. cleveland dot chart
                     xlab='range of data',
                     ylab='Order of the data',
                     labels=H.traits2.sp$Species,
                     main="SRL") # No outliers
            
            
           # X14 Myc.Frac - Okay
            ---
            boxplot(H.traits2.sp$Myc.Frac,
                      main = "Myccorhizal Fraction") # 2 outliers
            
            dotchart(H.traits2.sp$Myc.Frac, #i.e. cleveland dot chart
                     xlab='range of data',
                     ylab='Order of the data',
                     labels=H.traits2.sp$Species,
                     main="Myccorhizal Fraction") # No outliers
           
           # TRY HighlandStats code
           
           # Source Highland library v.10 for Highland's own wrapper function for Cleveland dotplots
           source("C:/Users/Julie/Desktop/Postdoc/Workshops/Highland Stats_GLM, GAMS/HighstatLibV10.R")   #<---
           
           Mydotplot(H.traits2.sp[ ,Trait.Names])
          
           
           # Summary
           # log transform leaf area, vegetative height and fine root thickness for sure 
           # because of outliers
           
           # Transformed Leaf thickness and LMA because log forms work better latter on in the models. 
           
           
           # Pay attention to Storage mass Fraction: 0-inflated. Transform to presence/absence?
           # Pay attention to Reproductive mass Fraction: skewed distribution. May need to be transf.
           
           
           # Can't log transform standardized values because values smaller than 1.
           H.traits2.sp$Log.Ht.veg<-log(H.traits2.sp$Ht.veg)
           H.traits2.sp$Log.Leaf.Area<-log(H.traits2.sp$Leaf.Area)
           H.traits2.sp$Log.F.Root.Diam<-log(H.traits2.sp$F.Root.Diam)
           H.traits2.sp$Log.LMA<-log(H.traits2.sp$LMA)
           H.traits2.sp$Log.Lamina.thck<-log(H.traits2.sp$Lamina.thck)
           
          # Update list of traits to use
           Trait.Names<-c("Log.Ht.veg","Min.Root.Loca","Max.Root.Loca","Log.Lamina.thck","Log.LMA",
                          "LDMC","Log.Leaf.Area","Leaf.Mass.Frac","Supp.Mass.Frac","Rep.Mass.Frac",
                          "Stor.Mass.Frac","Log.F.Root.Diam","SRL",'Myc.Frac')
           
        #A1.0.2 - Homogeneity (homoscedasticity) of Y - residuals vs fitted for each var.
         
            # Test after model is built, as validation step
      
        #A1.0.3 - Normality of Y - (least important - use link function instead). test for normality of residuals
           
            # Test after model is built, as validation step
      
        #A1.0.4 - Zero trouble (Y) ?
           
            sort(H.abund$abund.ratio)
            # only one 0 value - Okay! 
            
        #A1.0.5 - Collinearity  (X)
        
            
        # Try PairPlots        
        pairs(H.traits2.sp[,Trait.Names],
              lower.panel = panel.cor) # Requires to source HighstatlibV10.R
        
          # Summary
          # Problem with Min and Max root location: R2 = 0.97 - Pick one
          # Problem with LMA and LDMC: R2=0.8 - Pick one
          # Problem with Leaf Mass Fraction and Storage Mass Fraction: R2=-0.9 - Pick one
          # Problem with Fine Root Diameter and SRL: R2=-0.6 - Pick one
        
          # Maybe problem with Storage Mass fraction and Support Mass Fraction: R2=-0.5 - Pick one
          # Maybe problem with Veg.Height and Leaf Area: R2 = 0.4
          # Maybe problem with Support Mass Fraction and Veg.height: R2 =0.4
          # Maybe problem with LMA and Leaf thickness: R2=0.4
          # Maybe problem with LDMC and Leaf Mass Fraction: R2=0.4
          # Maybe problem with LDMC and Storage Mass Fraction: R2=-0.4
          # Maybe problem with Fine Root Diameter and Leaf Area: R2=0.4
          # Maybe problem with Fine Root Diameter and Rep. Mass Fraction: R2=0.4
        
          # Worry about y1-y2 covariance between 0.5-0.8 if relationshps between X and each of Y1
          # and y2 are weak. 
        
         summary(lm(abund.ratio~Log.F.Root.Diam,
                    data=merge(H.abund,H.traits2.sp,by="row.names",all=T)
                    ))
         # NS
         summary(lm(abund.ratio~SRL,
                    data=merge(H.abund,H.traits2.sp,by="row.names",all=T)
         ))
         # NS   
         
         # Both X-Y correlations are week, so remove one of Fine Root Diameter or SRL
         
         summary(lm(abund.ratio~Supp.Mass.Frac,
                    data=merge(H.abund,H.traits2.sp,by="row.names",all=T)
         ))
         # NS
         
         summary(lm(abund.ratio~Stor.Mass.Frac,
                    data=merge(H.abund,H.traits2.sp,by="row.names",all=T)
         ))
        # NS
         
        # Both X-Y correlations are week, so remove one of Storage or Support Mass Fraction
         
        # Now look at the VIFs
         
        # Do it using dataset without myccorhizal fraction because it this trait contains only 42 species
        # instead of 51, which decreases our sample size a lot. 
        
        # sequentially eliminate the variables with the highest VIFs  
        vif(lm(abund.ratio~Log.Ht.veg+Min.Root.Loca+Max.Root.Loca+Log.Lamina.thck+Log.LMA+
                 LDMC+Log.Leaf.Area+Leaf.Mass.Frac+Supp.Mass.Frac+Rep.Mass.Frac+
                 Stor.Mass.Frac+Log.F.Root.Diam+SRL,
               na.action=na.omit,
            data=merge(H.abund,H.traits2.sp,by="row.names",all=T)
            ))
        
        #     Log.Ht.veg   Min.Root.Loca   Max.Root.Loca Log.Lamina.thck         Log.LMA            LDMC 
        #       2.088027       64.270612       64.708128        4.427030       12.094997       10.354246 
        #  Log.Leaf.Area  Leaf.Mass.Frac  Supp.Mass.Frac   Rep.Mass.Frac  Stor.Mass.Frac Log.F.Root.Diam 
        #       3.439911      157.738467       33.286408       12.404597      179.907708        3.310469 
        #            SRL 
        #       4.605975 
        # 
        # Storage Mass Fraction is the worst
        
        vif(lm(abund.ratio~Log.Ht.veg+Min.Root.Loca+Max.Root.Loca+Log.Lamina.thck+Log.LMA+
                 LDMC+Log.Leaf.Area+Leaf.Mass.Frac+Supp.Mass.Frac+Rep.Mass.Frac+
                 Log.F.Root.Diam+SRL,
               na.action=na.omit,
               data=merge(H.abund,H.traits2.sp,by="row.names",all=T)
        ))
        
        #     Log.Ht.veg   Min.Root.Loca   Max.Root.Loca Log.Lamina.thck         Log.LMA            LDMC 
        #       1.877799       63.093456       63.744531        4.523619       10.608937        9.930118 
        #  Log.Leaf.Area  Leaf.Mass.Frac  Supp.Mass.Frac   Rep.Mass.Frac Log.F.Root.Diam             SRL 
        #       2.175915        3.078716        1.670315        1.683049        3.154330        4.057710
        
        # Max.Root.Loca is the worst
        
        vif(lm(abund.ratio~Log.Ht.veg+Min.Root.Loca+Log.Lamina.thck+Log.LMA+
                 LDMC+Log.Leaf.Area+Leaf.Mass.Frac+Supp.Mass.Frac+Rep.Mass.Frac+
                 Log.F.Root.Diam+SRL,
               na.action=na.omit,
               data=merge(H.abund,H.traits2.sp,by="row.names",all=T)
        ))
        
        #     Log.Ht.veg   Min.Root.Loca Log.Lamina.thck         Log.LMA            LDMC   Log.Leaf.Area 
        #       1.872089        1.850274        3.864454        9.107514        8.212370        2.126254 
        # Leaf.Mass.Frac  Supp.Mass.Frac   Rep.Mass.Frac Log.F.Root.Diam             SRL 
        #       2.808127        1.670286        1.668985        3.122450        3.837607 
         
        # TRY without LMA
        vif(lm(abund.ratio~Log.Ht.veg+Min.Root.Loca+Log.Lamina.thck+
                 LDMC+Log.Leaf.Area+Leaf.Mass.Frac+Supp.Mass.Frac+Rep.Mass.Frac+
                 Log.F.Root.Diam+SRL,
               na.action=na.omit,
               data=merge(H.abund,H.traits2.sp,by="row.names",all=T)
        ))
        
        #     Log.Ht.veg   Min.Root.Loca Log.Lamina.thck            LDMC   Log.Leaf.Area  Leaf.Mass.Frac 
        #       1.724775        1.792609        1.905555        2.490835        2.046712        2.507318 
        # Supp.Mass.Frac   Rep.Mass.Frac Log.F.Root.Diam             SRL 
        #       1.661494        1.579522        3.100699        3.782878 
        
        # Then try without Fine root Diameter, because I don't want to loose SRL
        
        vif(lm(abund.ratio~Log.Ht.veg+Min.Root.Loca+Log.Lamina.thck+
                 LDMC+Log.Leaf.Area+Leaf.Mass.Frac+Supp.Mass.Frac+Rep.Mass.Frac+
                 SRL,
               na.action=na.omit,
               data=merge(H.abund,H.traits2.sp,by="row.names",all=T)
        ))
        
        #     Log.Ht.veg   Min.Root.Loca Log.Lamina.thck            LDMC   Log.Leaf.Area  Leaf.Mass.Frac 
        #       1.658229        1.637638        1.882275        1.988683        2.024708        2.205449 
        # Supp.Mass.Frac   Rep.Mass.Frac             SRL 
        #       1.347331        1.517007        1.761267
        
        # Good to go! All VIFs below 3
        
        
        
        # A1.0.6 Relationship between Y and X. Linear?
      
        # A1.0.7 Interactions?
      
      
      
      # A1.5- Transform trait data Transform trait data as needed ####
      #===========================================================#
      
      # Test best transform for each trait ####
      
      # Following Kerkhoff & Enquist (2009), all size- and growth-related traits are natural-log transformed, because 
      # they follow multiplicative processes: Ht.veg, Lamina.thck, LMA, LDMC, Leaf Area).

      # 1- Ht.veg
        plot(density(na.omit(H.traits$Ht.veg)))
        qqnorm(na.omit(H.traits$Ht.veg)) # super banana
        shapiro.test(na.omit(H.traits$Ht.veg)) # p=2.2e-16
                                                                  # Shapiro test null is that data is normal
                                                                  # i.e. non-significance means data is normal
        
        shapiro.test(log(H.traits$Ht.veg))                      # 0.0004355
        powerTransform(H.traits$Ht.veg) # -0.1247701
        shapiro.test(H.traits$Ht.veg^-0.1247701) #p-value = 0.03754
        shapiro.test(1/H.traits$Ht.veg^-0.1247701)# p = 2.026e-08
        shapiro.test(H.traits$Ht.veg^0.123255)#2.311e-08
        # If I transform with a negative exponent, I can't use 1/values to regain original order
        shapiro.test(-1*H.traits$Ht.veg^-0.1247701)             # 0.03754
        plot(density(log(na.omit(H.traits$Ht.veg))))
        plot(density(log(na.omit(H.traits$Ht.veg^-0.123255))))
        
            t<-na.omit(H.traits[,c("Species","Ht.veg")])
            t<-droplevels(t)
            medi<-tapply(t$Ht.veg,t$Species,FUN= median) # calculate means of species
            boxplot(t$Ht.veg~factor(t$Species,levels=levels(t$Species)[order(medi)]),
                    las=2, main='ht.veg')
            #all looks good
        
        # Best Transform - size-based traits use log & no large difference with PowerTransform exponent
        H.traits$Ht.veg<-log(H.traits$Ht.veg)
      
      # 2- Min.Root.Loca
        plot(density(na.omit(H.traits$Min.Root.Loca))) 
        # Totally bimodal. No transformation will fix that. 
        # Now is an ordered factor
      
      # 3- Max.Root.Loca
        plot(density(na.omit(H.traits$Max.Root.Loca))) 
        # Totally bimodal. No transformation will fix that. 
        # Now is an ordered factor
        
            # Why bimodal?
            t1<-na.omit(H.traits[,c("Species","Min.Root.Loca")])
            t1<-droplevels(t1)
            medi<-tapply(t1$Min.Root.Loca,t1$Species,FUN= median) # calculate means of species
            boxplot(t1$Min.Root.Loca~factor(t1$Species,levels=levels(t1$Species)[order(medi)]),
                    las=2)
      
      # 4- Lamina.thck
        plot(density(na.omit(H.traits$Lamina.thck)))
        shapiro.test(H.traits$Lamina.thck) # 2.2e-16
        shapiro.test(log(na.omit(H.traits$Lamina.thck)))       # 5.668e-10
        powerTransform(na.omit(H.traits$Lamina.thck))# -0.6403293
        shapiro.test(H.traits$Lamina.thck^-0.6403293)          # 9.266e-10
        plot(density(na.omit(log(H.traits$Lamina.thck))))
        plot(density(na.omit(H.traits$Lamina.thck^-0.6403293))) 
        
        # Best Transform - log bc size-related & no large difference with PowerTransform exponent
        H.traits$Lamina.thck<-log(H.traits$Lamina.thck)
      
      # 5 - LMA
        plot(density(na.omit(H.traits$LMA)))
        qqnorm(na.omit(H.traits$LMA))
        shapiro.test(H.traits$LMA)         # 2.2e-16
        shapiro.test(log(na.omit(H.traits$LMA)))            # 3.644e-09
        powerTransform(na.omit(H.traits$LMA))# -0.3581721 
        shapiro.test(H.traits$LMA^-0.3581721)              # 0.0006073
        
        plot(density(na.omit(log(H.traits$LMA))))
        plot(density(na.omit(H.traits$LMA^-0.3581721)))
        
            t2<-na.omit(H.traits[,c("Species","LMA")])
            t2<-droplevels(t2)
            md<-tapply(t2$LMA,t2$Species,FUN= median) # calculate means of species
            boxplot(t2$LMA~factor(t2$Species,levels=levels(t2$Species)[order(md)]),
                    las=2, main="LMA")
        
        # Best Transform - PowerTransform is WAY better, but log transform because it is a multiplicative trait. 
        H.traits$LMA<-log(H.traits$LMA)    
        
      # 6 - LDMC
        plot(density(na.omit(H.traits$LDMC))) # right skewed
        qqnorm(H.traits$LDMC) # S shaped
        shapiro.test(H.traits$LDMC) #  1.026e-15
        shapiro.test(log(na.omit(H.traits$LDMC)))          #1.529e-05
        powerTransform(H.traits$LDMC)# -0.07986793 
        shapiro.test(H.traits$LDMC^-0.07986793 )             #2.582e-05
        plot(density(na.omit(H.traits$LDMC^-0.07986793 )))
        plot(density(na.omit(log(H.traits$LDMC))))
        qqnorm(log(H.traits$LDMC)) 
        qqnorm(H.traits$LDMC^-0.07986793 ) 
        
            # species boxplots
            t3<-na.omit(H.traits[,c("Species","LDMC")])
            t3<-droplevels(t3)
            md<-tapply(t3$LDMC,t3$Species,FUN= median) # calculate means of species
            boxplot(t3$LDMC~factor(t3$Species,levels=levels(t3$Species)[order(md)]),
                    las=2, main='LDMC')
            
        # Best Transform - PowerTransform is way better, but it iss a ratio of weights therefore use log transform
          H.traits$LDMC<-log(H.traits$LDMC)
        
      # 7-Leaf.Area
          plot(density(na.omit(H.traits$Leaf.Area)))
          qqnorm(H.traits$Leaf.Area) # Super Banana due to a few outliers
          shapiro.test(H.traits$Leaf.Area) # 2.2e-16
          shapiro.test(log(na.omit(H.traits$Leaf.Area)))         #3.741e-07
          powerTransform(H.traits$Leaf.Area)# 0.0757846
          shapiro.test(H.traits$Leaf.Area^0.0757846 )            # 0.0002567
          plot(density(na.omit(H.traits$Leaf.Area^0.0757846 )))
          plot(density(na.omit(log(H.traits$Leaf.Area))))
          qqnorm(log(H.traits$Leaf.Area)) # 3 outliers
          qqnorm(H.traits$Leaf.Area^0.0757846 ) 
          
            t4<-na.omit(H.traits[,c("Species","Leaf.Area")])
            t4<-droplevels(t4)
            md<-tapply(t4$Leaf.Area,t4$Species,FUN= median) # calculate means of species
            boxplot(t4$Leaf.Area~factor(t4$Species,levels=levels(t4$Species)[order(md)]),
                    las=2,main="Leaf.Area")
      
          # Best Transform - size trait therefore log
          H.traits$Leaf.Area<-log(H.traits$Leaf.Area)
        
      # 8-Leaf.Mass.Frac
          plot(density(na.omit(H.traits$Leaf.Mass.Frac)))
          qqnorm(H.traits$Leaf.Mass.Frac) 
          shapiro.test(H.traits$Leaf.Mass.Frac) # 2.293e-06
          shapiro.test(log(na.omit(H.traits$Leaf.Mass.Frac)))         #2.2e-16 - worse with log transform
          powerTransform(H.traits$Leaf.Mass.Frac)# 1.024067 # exponent 1! i.e. don't transform it....
          shapiro.test(H.traits$Leaf.Mass.Frac^1.024067)            # 1.212e-06
          plot(density(na.omit(H.traits$Leaf.Mass.Frac^1.024067)))
      
          
          t5<-na.omit(H.traits[,c("Species","Leaf.Mass.Frac")])
          t5<-droplevels(t5)
          md<-tapply(t5$Leaf.Mass.Frac,t5$Species,FUN= median) # calculate means of species
          boxplot(t5$Leaf.Mass.Frac~factor(t5$Species,levels=levels(t5$Species)[order(md)]),
                  las=2,main="Leaf.Mass.Frac")
        
          # do NOT transform
      
      # 9-Supp.Mass.Frac
          plot(density(na.omit(H.traits$Supp.Mass.Frac)))
          qqnorm(H.traits$Supp.Mass.Frac) 
          shapiro.test(H.traits$Supp.Mass.Frac) #  1.775e-06
          shapiro.test(log(na.omit(H.traits$Supp.Mass.Frac)))         #3.108e-14 - worse with log transform
          powerTransform(H.traits$Supp.Mass.Frac)# 0.6124825
          shapiro.test(H.traits$Supp.Mass.Frac^0.6124825)            # 0.009796
          plot(density(na.omit(H.traits$Supp.Mass.Frac^0.6124825)))
          
          
          t6<-na.omit(H.traits[,c("Species","Supp.Mass.Frac")])
          t6<-droplevels(t6)
          md<-tapply(t6$Supp.Mass.Frac,t6$Species,FUN= median) # calculate means of species
          boxplot(t6$Supp.Mass.Frac~factor(t6$Species,levels=levels(t6$Species)[order(md)]),
                  las=2,main="Supp.Mass.Frac")
        
          # Best Transform
          H.traits$Supp.Mass.Frac<-H.traits$Supp.Mass.Frac^0.6124825
      
      # 10-Rep.Mass.Frac
        # H.traits$Rep.Mass.Frac<-Megtraits[Megtraits$Layer=='H','Rep.Mass.Frac'] 
        summary(H.traits$Rep.Mass.Frac[which(H.traits$Rep.Mass.Frac!=0)]) # smallest non-zero Rep.Mass.Frac is 0.01000
        
        H.traits$Rep.Mass.Frac<-H.traits$Rep.Mass.Frac+0.001            # Add 0.001 to all my data to get rid of 0s
        plot(density(na.omit(H.traits$Rep.Mass.Frac)))
        qqnorm(H.traits$Rep.Mass.Frac) 
        shapiro.test(H.traits$Rep.Mass.Frac) # 2.2e-16
        shapiro.test(log(na.omit(H.traits$Rep.Mass.Frac)))         # 2.2e-16
        powerTransform(H.traits$Rep.Mass.Frac)# -0.5502378
        shapiro.test(H.traits$Rep.Mass.Frac^-0.564259)            # 2.2e-16
        plot(density(na.omit(H.traits$Rep.Mass.Frac^-0.564259))) # two normal curves
        plot(density(na.omit(log(H.traits$Rep.Mass.Frac))))
        
        
        t7<-na.omit(H.traits[,c("Species","Rep.Mass.Frac")])
        t7<-droplevels(t7)
        md<-tapply(t7$Rep.Mass.Frac,t7$Species,FUN= median) # calculate means of species
        boxplot(t7$Rep.Mass.Frac~factor(t7$Species,levels=levels(t7$Species)[order(md)]),
                las=2,main="Rep.Mass.Frac")
        
        # Best Transform 
        H.traits$Rep.Mass.Frac<-(H.traits$Rep.Mass.Frac+0.001)^-0.564259
      
      
      # 11-Stor.Mass.Frac
        
        # H.traits$Stor.Mass.Frac<-Megtraits[Megtraits$Layer=='H','Stor.Mass.Frac']
        summary(H.traits$Stor.Mass.Frac[which(H.traits$Stor.Mass.Frac!=0)]) #  # smallest non-zero Rep.Mass.Frac is 0.04000
        
        H.traits$Stor.Mass.Frac<-H.traits$Stor.Mass.Frac+0.001            # Add 0.001 to all my data to get rid of 0s
        plot(density(na.omit(H.traits$Stor.Mass.Frac)))
        qqnorm(H.traits$Stor.Mass.Frac) 
        shapiro.test(H.traits$Stor.Mass.Frac) # 2.2e-16
        shapiro.test(log(na.omit(H.traits$Stor.Mass.Frac)))         # 2.2e-16
        powerTransform(H.traits$Stor.Mass.Frac)# -0.07815023
        shapiro.test(H.traits$Stor.Mass.Frac^-0.07569593)            # 2.2e-16
        plot(density(na.omit(H.traits$Stor.Mass.Frac^-0.07569593))) # two normal curves
        plot(density(na.omit(log(H.traits$Stor.Mass.Frac)))) # two normal curves
        
        
        t8<-na.omit(H.traits[,c("Species","Stor.Mass.Frac")])
        t8<-droplevels(t8)
        md<-tapply(t8$Stor.Mass.Frac,t8$Species,FUN= median) # calculate means of species
        boxplot(t8$Stor.Mass.Frac~factor(t8$Species,levels=levels(t8$Species)[order(md)]),
                las=2,main="Stor.Mass.Frac")
        
        # Best Transform 
        H.traits$Stor.Mass.Frac<-(H.traits$Stor.Mass.Frac+0.001)^-0.07569593
        
      # Summary of trait transforms ####
      #================================#
        
        H.traits$Ht.veg<-log(H.traits$Ht.veg)
        # Min.Root.Loca - no transforms
        # Max.Root.Loca - no transforms
        H.traits$Lamina.thck<-log(H.traits$Lamina.thck)
        H.traits$LMA<-log(H.traits$LMA)
        H.traits$LDMC<-log(H.traits$LDMC)
        H.traits$Leaf.Area<-log(H.traits$Leaf.Area)
        #Leaf.Mass.Frac - no transforms
        H.traits$Supp.Mass.Frac<-H.traits$Supp.Mass.Frac^0.6124825
        H.traits$Rep.Mass.Frac<-(H.traits$Rep.Mass.Frac+0.001)^-0.564259
        H.traits$Stor.Mass.Frac<-(H.traits$Stor.Mass.Frac+0.001)^-0.07569593
        
  

        
  # A2 - Canopy layer ####
  #=======================#
        
    # A2.1- Transform trait data for normality ####
    #=============================================#
        
        # create trait dataframe with traits I will work with
        
          C.traits<-Megtraits[Megtraits$Layer=='C',c("Plant.ID","Species","Layer","Date.recolte",
          "Lamina.thck","LMA","LDMC","Leaf.Area")]
      
      dim(C.traits)
      # 181   8       
      
      # Following Kerkhoff & Enquist (2009), all size- and growth-related traits are natural-log transformed,  
      # because they follow multiplicative processes: Lamina.thck, LMA, LDMC, Leaf Area.
      
      # Test best transform for each trait ####
      
      # 4 - Lamina thickness
        plot(density(na.omit(C.traits$Lamina.thck))) # Bimodal
        shapiro.test(C.traits$Lamina.thck) # 9.971e-16
        shapiro.test(log(na.omit(C.traits$Lamina.thck)))       # 6.935e-09
        powerTransform(na.omit(C.traits$Lamina.thck))# -0.3236387
        shapiro.test(C.traits$Lamina.thck^-0.3236387)   
        
        # Best Transform - log bc size-related & no large difference with PowerTransform exponent
        C.traits$Lamina.thck<-log(C.traits$Lamina.thck)
      
      # 5 - LMA
        plot(density(na.omit(C.traits$LMA))) #bimodal
        qqnorm(na.omit(C.traits$LMA))
        shapiro.test(C.traits$LMA)         # 8.594e-11
        shapiro.test(log(na.omit(C.traits$LMA)))            # 0.00413
        powerTransform(na.omit(C.traits$LMA))# 0.08602909 
        shapiro.test(C.traits$LMA^0.08602909)               # 0.004137
        
        plot(density(na.omit(log(C.traits$LMA))))
        plot(density(na.omit(C.traits$LMA^0.08602909)))
        
        # Best Transform - log transform because it is a multiplicative trait, and just as good as power Transform
        C.traits$LMA<-log(C.traits$LMA)
      
      # 6 - LDMC
        plot(density(na.omit(C.traits$LDMC))) # weirdly shaped
        shapiro.test(C.traits$LDMC) #   9.555e-05
        shapiro.test(log(na.omit(C.traits$LDMC)))          #2.94e-09
        powerTransform(C.traits$LDMC)#1.276166 
        shapiro.test(C.traits$LDMC^1.276166)             # 0.000209
        plot(density(na.omit(C.traits$LDMC^1.276166)))
        plot(density(na.omit(log(C.traits$LDMC))))
        
        # Best Transform - log transform bc it is a multiplicative trait, but powerTransform is way better
        C.traits$LDMC<-log(C.traits$LDMC)
        
      # 7 - Leaf Area
        plot(density(na.omit(C.traits$Leaf.Area)))
        shapiro.test(C.traits$Leaf.Area) # 5.992e-16
        shapiro.test(log(na.omit(C.traits$Leaf.Area)))         #0.002092
        powerTransform(C.traits$Leaf.Area)#  -0.2679859
        shapiro.test(C.traits$Leaf.Area^ -0.2679859 )            #0.2017
        plot(density(na.omit(C.traits$Leaf.Area^-0.2679859 )))
        plot(density(na.omit(log(C.traits$Leaf.Area))))
        
        # Best Transform - size trait therefore log
        C.traits$Leaf.Area<-log(C.traits$Leaf.Area)
        
      # Summary of trait transforms ####
        C.traits$Lamina.thck<-log(C.traits$Lamina.thck)
        C.traits$LMA<-log(C.traits$LMA)
        C.traits$LDMC<-log(C.traits$LDMC)
        C.traits$Leaf.Area<-log(C.traits$Leaf.Area)
        
    # A2.2- Remove date effect on traits ####
    #=======================================#
        
      names(C.traits)
      C.trait.names<-c("Lamina.thck","LMA","LDMC","Leaf.Area")
      length(C.trait.names) #4
      
      save(C.trait.names,file=paste0(wrk.dir,'list.trait.names.canopy.layer.Rdata'))
      
      # Change caledar dates into julian dates
      C.traits$Date.recolte<- as.numeric(format(as.Date(C.traits$Date.recolte, format = "%m/%d/%Y"),"%j"))
      
      # loop testing for date effects and replacing with regression residuals
      for (t in C.trait.names){ 
        x<-lm(C.traits[[t]]~C.traits$Date.recolte,na.action=na.exclude)  # # regress the trait against the julian date
        
        if (summary(x)$adj.r.squared > 0.02 & summary(x)[["coefficients"]][2,4] < 0.05) # if regression R2 > 0.02 AND it is statistically significant (with P-value<0.05), go through this next loop
          print(c(names(C.traits[t]),
                  paste( "R2= ",round(summary(x)$adj.r.squared,digits=3)),
                  paste('Sign= ',sign(summary(x)[["coefficients"]][2,1])))) 
        
        for (i in 1:nrow(C.traits[t])){                                                   # if each row is numeric (not an NA)
          if (is.numeric(C.traits[t][i,]))
          {C.traits[t][i,]<-resid(x)[i]}                                              # then replace the value with the residual
        }                                                                           # else, do nothing to the cells with value of NA
      }
      
      # [1] "Lamina.thck" "R2=  0.159"  "Sign=  -1"  
      # [1] "LMA"         "R2=  0.072"  "Sign=  -1" 
      # [1] "Leaf.Area"   "R2=  0.052"  "Sign=  -1" Why is Leaf area decreasing with season?   
      
      # Standardize all variables
      C.traits[5:8]<-decostand(C.traits[5:8],method='standardize',margin=2)
      
      save(C.traits,file=paste0(wrk.dir,"Canopy.Layer.Traits.Standardized.RData"))
      
    # A2.3 - Turn individual-level table into species-level table ####
    # ===============================================================#  
      
      sp.C.traits<-as.data.frame(aggregate(C.traits[,5:8],by=list(C.traits$Species),mean,na.rm=T))
      sp.C.traits[,2:5]<-round(sp.C.traits[2:5],digits=3)
      sp.C.traits[is.na(sp.C.traits)] <-NA
      rownames(sp.C.traits)<-sp.C.traits$Group.1
      sp.C.traits<-sp.C.traits[,-1]
      sp.C.traits<-sp.C.traits[order(rownames(sp.C.traits)),]
      dim(sp.C.traits) # 24  4
      
      save(sp.C.traits,file=paste0(wrk.dir,'species-level.traits.Canopy.layer.Rdata'))
      
    # A2.4 - Add species mean traits (Seed size) ####
      
        
# B- Setup abundance data ####
# ===========================#

      # B1 - Look at temporal patterns ####
      # ====================================#
      
      plot(abund[abund$Layer=='H','avg.abundance.2012']~abund[abund$Layer=='H','avg.abundance.1970'],type='n',
           xlab='av abundance 1970', ylab='avg abundance 2012',main='herbaceous')
      text(abund[abund$Layer=='H','avg.abundance.1970'],abund[abund$Layer=='H','avg.abundance.2012'],
           labels=rownames(abund)[which(abund$Layer=='H')])
      abline(0,1)
      
      #on log-scale, much better
      plot(log(abund[abund$Layer=='H','avg.abundance.2012'])~log(abund[abund$Layer=='H','avg.abundance.1970']),type='n',
           xlab='log(av abundance 1970)', ylab='log(avg abundance 2012)',main='log herbaceous')
      text(log(abund[abund$Layer=='H','avg.abundance.1970']),log(abund[abund$Layer=='H','avg.abundance.2012']),
           labels=rownames(abund)[which(abund$Layer=='H')])
      abline(0,1)
      
      # for Canopy layer: 
      plot(abund[abund$Layer=='C','avg.abundance.2012']~abund[abund$Layer=='C','avg.abundance.1970'],type='n',
           xlab='av abundance 1970', ylab='avg abundance 2012',main='canopy')
      text(abund[abund$Layer=='C','avg.abundance.1970'],abund[abund$Layer=='C','avg.abundance.2012'],
           labels=rownames(abund)[which(abund$Layer=='C')])
      abline(0,1)
      
      plot(log(abund[abund$Layer=='C','avg.abundance.2012'])~log(abund[abund$Layer=='C','avg.abundance.1970']),type='n',
           xlab='log(av abundance 1970)', ylab='log(avg abundance 2012)',main='log canopy')
      text(log(abund[abund$Layer=='C','avg.abundance.1970']),log(abund[abund$Layer=='C','avg.abundance.2012']),
           labels=rownames(abund)[which(abund$Layer=='C')])
      abline(0,1)
        
    # B2- Create abundance dataframe & new variables ####
    #=====================================================#
        
    abund<-read.csv(paste0(data.dir,'sp.abund.csv'))
    rownames(abund)<-abund$Species
    abund$Species<-NULL
    
      # Abundance/Presence Ratios
      
      abund[,'presence.ratio']<-round(abund$pct.plot.present.2012/abund$pct.plot.present.1970,digits=1)
      abund[,'abund.ratio']<-round(abund$avg.abundance.2012/abund$avg.abundance.1970,digits=1)
      abund$abund.ratio
      rownames(abund[abund$abund.ratio=='Inf',]) 
      # [1] "CAAL" "CACR" "CAPE" "CATH" "COMA" "DEPU" "GAPR" "GATE" "HYAM" "JUTE" "LICO" "LYCL" "LYUN" "MOUN" "PAQU" "POGR"
      # [17] "RARE" "TAOF
      # loosing 18 data points to "Inf"
      plot(density(abund$abund.ratio))
      plot(density(log(abund$abund.ratio+0.01)))
      
      # Log of Abundance/Presence Ratios
      
      abund[,'log.presence.ratio']<-log(abund$presence.ratio+0.01)
      abund[,'log.abund.ratio']<-log(abund$abund.ratio+0.01)
      abund$log.abund.ratio
      rownames(abund[abund$log.abund.ratio=='Inf',])
      # [1] "CAAL" "CACR" "CAPE" "CATH" "COMA" "DEPU" "GAPR" "GATE" "HYAM" "JUTE" "LICO" "LYCL" "LYUN" "MOUN" "PAQU" "POGR"
      # [17] "RARE" "TAOF"
      #loosing 18 data points to Inf
      plot(density(abund$log.abund.ratio))
      rownames(abund[abund$log.abund.ratio<=-3,])
      # [1] "ALTR"  "DRGO"  "EPAN"  "ERST"  "EUPMA" "FRNI"  "GORE"  "LIBO"  "OXMO"  "TSCA" 
      # sencond density peak has 10 species. 
      
      # Ratios of log Abundance/Presence
      
      abund[,'ratio.log.presence']<-round(log(abund$pct.plot.present.2012)/log(abund$pct.plot.present.1970),digits=3)
      abund[,'ratio.log.abund']<-round(log(abund$avg.abundance.2012)/log(abund$avg.abundance.1970),digits=3)
      abund$ratio.log.abund
      plot(density(abund$ratio.log.abund))
      rownames((abund)[abund$ratio.log.abund<=-10,]) # IMCA, POTR
      rownames((abund)[abund$ratio.log.abund>=20,]) # "ACPE"  "ALTR"  "DRGO"  "EPAN"  "ERST"  "EUPMA" "FRNI"  "GORE"  "LIBO"  "TSCA"
      
      # powerTransform Abundance/Presence ####
      
      a<-abund$abund.ratio
      a[which(a=='Inf')]<-NA
      powerTransform(a+0.01) # 0.181293
      abund[,'PT.abund.ratio']<-(a+0.01)^0.181293
      b<-abund$presence.ratio
      b[which(b==Inf)]<-NA
      powerTransform(b+0.01)# -.3390714
      abund[,'PT.presence.ratio']<-(b+0.01)^0.3390714
    
    str(abund)
      # $ Layer                : Factor w/ 2 levels "C","H": 1 1 1 1 2 2 2 2 2 2 ...
      # $ pct.plot.present.1970: num  87.5 52.1 18.8 54.2 89.6 14.6 2.1 58.3 4.2 12.5 ...
      # $ pct.plot.present.2012: num  85.4 70.8 29.2 62.5 79.2 4.2 0 75 6.3 27.1 ...
      # $ avg.abundance.1970   : num  27.2 1.03 5.12 25.97 6.73 ...
      # $ avg.abundance.2012   : num  20.27 4.06 3.82 24.3 4.85 ...
      # $ presence.ratio       : num  1 1.4 1.6 1.2 0.9 0.3 0 1.3 1.5 2.2 ...
      # $ abund.ratio          : num  0.7 3.9 0.7 0.9 0.7 0.3 0 1.1 1.5 2.1 ...
      # $ log.presence.ratio   : num  0.00995 0.34359 0.47623 0.19062 -0.09431 ...
      # $ log.abund.ratio      : num  -0.3425 1.3635 -0.3425 -0.0943 -0.3425 ...
      # $ ratio.log.abund      : num  0.911 45.921 0.821 0.98 0.829 ...
      # $ ratio.log.presence   : num  0.995 1.078 1.15 1.036 0.973 ...
      # $ PT.abund.ratio       : num  0.94 1.28 0.94 0.983 0.94 ...
      # $ PT.presence.ratio    : num  1.003 1.124 1.175 1.067 0.969 ...
    
    save(abund,file=paste0(wrk.dir,'Species.abundances.full.data.Rdata'))
  
    # link function of glms (negative binomial) made to deal with lack of normality of response variable. 

    # B3 - Problem with infinity values ####
    #========================================#
    
      # ABUNDANCE
    
      # Which species were absent in 1970?
      rownames(abund)[which(abund$avg.abundance.1970==0)]  
      # [1] "CAAL" "CACR" "CAPE" "CATH" "COMA" "DEPU" "GAPR" "GATE" "HYAM" "JUTE" "LICO" "LYCL" "LYUN"  "MOUN" "PAQU" "POGR" 
      # [17]"RARE" "TAOF"
      
      # What is the minimum value in dataset?
      abund$avg.abundance.1970[order(abund$avg.abundance.1970)][15:25]
      #0.00 0.00 0.00 0.00 0.01 0.01 0.01 0.01 0.01 0.01 0.01
      
      # replace 0s with minimum value
      abund$non0.avg.abundance.1970<-abund$avg.abundance.1970
      abund[c("CAAL","CACR","CAPE","CATH","COMA","DEPU","GAPR","GATE",
              "HYAM","JUTE","LICO","LYCL","LYUN"),'non0.avg.abundance.1970']<-0.01 
      
      #Re-calculate abundance
      non0.abund.change<-round(abund$avg.abundance.2012/abund$non0.avg.abundance.1970,digits=2)
      plot(density(non0.abund.change)) 
      # now ranges from 0-200, with a single value at 200. Bad idea to replace 1970's 0s.
      
      abund$non0.avg.abundance.1970<-NULL
      
      #save new abundance dataset without INF values
      abund.c<-abund[!rownames(abund)%in% c("CAAL","CACR","CAPE","CATH","COMA","DEPU","GAPR","GATE",
                                     "HYAM","JUTE","LICO","LYCL","LYUN","MOUN","PAQU","POGR",
                                     "RARE","TAOF"),]
    
    dim(abund.c) # 107 13
    
    save(abund.c,file=paste0(wrk.dir,'Species.abundances.Inf.removed.Rdata')) 
    
    
    
    ###############################################
    #Scrap Pile ####
    
    library(nlme)
    date.effect<-lme(Ht.veg~1,random=~1|Species/Date.recolte,data=H.traits,na.action=na.exclude)
    summary(date.effect)
    
    # to check if the date effect is significant, compare with a model with only species effect
    species.effect<-lme(Ht.veg~1,random=~1|Species,data=H.traits,na.action=na.exclude)
    anova(species.effect,date.effect)
    #             Model df      AIC      BIC    logLik   Test  L.Ratio p-value
    # species.effect     1  3 3858.244 3870.612 -1926.122                        
    # date.effect        2  4 3856.969 3873.458 -1924.484 1 vs 2 3.275978  0.0703
    
    # date effect within species (interaction between species and date) is marginally significant
    
    
    
    # Look at effect of date independent of effect of species. Species and date are correlated. 
    
    m<-lm(Date.recolte~Species, data=H.traits)
    summary(m)
    
    lm(formula = Date.recolte ~ Species, data = H.traits)
    
    # Residuals:
    #   Min      1Q  Median      3Q     Max 
    # -38.077  -2.286   0.000   1.600  38.500 
    # 
    # Coefficients:
    #             Estimate Std. Error t value Pr(>|t|)    
    # (Intercept) 153.7778     3.2838  46.829  < 2e-16 ***
    # SpeciesARTR  28.4444     4.6440   6.125 2.14e-09 ***
    # SpeciesATFI  47.2222     4.9647   9.512  < 2e-16 ***
    # SpeciesCAIN  16.6222     4.5264   3.672 0.000272 ***
    # SpeciesCASC  27.2222     4.5264   6.014 4.02e-09 ***
    # SpeciesCIAL  31.2222     4.7869   6.522 2.05e-10 ***
    # SpeciesCLBO -10.1778     4.5264  -2.249 0.025076 *  
    # SpeciesCLCA -27.5051     4.4279  -6.212 1.29e-09 ***
    # SpeciesCOAL  26.3333     4.6440   5.670 2.70e-08 ***
    # SpeciesCOCA   7.4222     4.5264   1.640 0.101827    
    # SpeciesCOCO  75.3651     4.9647  15.180  < 2e-16 ***
    # SpeciesCOTR -18.6778     4.5264  -4.126 4.47e-05 ***
    # SpeciesCYAC  -5.2778     7.7012  -0.685 0.493534    
    # SpeciesDEAC  60.9365     4.9647  12.274  < 2e-16 ***
    # SpeciesDEDI  28.2222     4.5264   6.235 1.13e-09 ***
    # SpeciesDEPU  36.7222     4.7869   7.671 1.26e-13 ***
    # SpeciesDRIN  36.5556     4.6440   7.872 3.19e-14 ***
    # SpeciesEPHE  39.8472     4.7869   8.324 1.29e-15 ***
    # SpeciesERAM -29.3778     4.5264  -6.490 2.49e-10 ***
    # SpeciesGAPR -27.7778     4.5264  -6.137 2.00e-09 ***
    # SpeciesGATE  22.0000     4.6440   4.737 3.00e-06 ***
    # SpeciesGATR  35.7222     5.1922   6.880 2.26e-11 ***
    # SpeciesHULU  69.6508     4.9647  14.029  < 2e-16 ***
    # SpeciesIMCA  18.4222     4.5264   4.070 5.65e-05 ***
    # SpeciesLOCA  26.9495     4.4279   6.086 2.67e-09 ***
    # SpeciesLYAN  67.5079     4.9647  13.598  < 2e-16 ***
    # SpeciesLYOB  63.9222     4.5264  14.122  < 2e-16 ***
    # SpeciesMACA -18.2778     4.5264  -4.038 6.44e-05 ***
    # SpeciesMEVI  -7.5556     4.6440  -1.627 0.104521    
    # SpeciesMIRE  41.7222     4.7869   8.716  < 2e-16 ***
    # SpeciesOCAC  43.2991     4.2719  10.136  < 2e-16 ***
    # SpeciesOSCI  68.9365     4.9647  13.885  < 2e-16 ***
    # SpeciesOSCL  16.9495     4.4279   3.828 0.000150 ***
    # SpeciesOXMO  -2.2778     4.5264  -0.503 0.615083    
    # SpeciesPHCO  35.8222     4.5264   7.914 2.37e-14 ***
    # SpeciesPOPU  -2.0778     4.5264  -0.459 0.646455    
    # SpeciesPRAL  23.7222     4.5264   5.241 2.57e-07 ***
    # SpeciesPYSE   6.3131     4.4279   1.426 0.154701    
    # SpeciesRUID  66.6508     4.9647  13.425  < 2e-16 ***
    # SpeciesRUPU  48.7937     4.9647   9.828  < 2e-16 ***
    # SpeciesSAPU  43.6508     4.9647   8.792  < 2e-16 ***
    # SpeciesSMRA  -1.2323     4.4279  -0.278 0.780916    
    # SpeciesSTAM  16.1222     4.5264   3.562 0.000412 ***
    # SpeciesSTLA  -0.9596     4.4279  -0.217 0.828538    
    # SpeciesTHNO  35.2222     4.9647   7.095 5.77e-12 ***
    # SpeciesTICO  21.3333     4.6440   4.594 5.81e-06 ***
    # SpeciesTRBO  -7.0778     4.5264  -1.564 0.118674    
    # SpeciesTRER -14.9778     4.5264  -3.309 0.001020 ** 
    # SpeciesTRUN -15.1778     4.5264  -3.353 0.000874 ***
    # SpeciesVEVI  21.5222     4.5264   4.755 2.76e-06 ***
    # SpeciesVIAL  63.7937     4.9647  12.850  < 2e-16 ***
    # ---
    # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 9.851 on 408 degrees of freedom
    # Multiple R-squared:  0.8989,	Adjusted R-squared:  0.8865 
    # F-statistic: 72.56 on 50 and 408 DF,  p-value: < 2.2e-16
    
    # Get intraspecific variation in date by extracting residuals of Date~Species
    length(resid(m)) #459
    resid.Date<-resid(m) 
    
    
    
    # does sampling date have an effect, once we remove the species effect?
    library(lme4)
    df<-H.traits[,c('Ht.veg','Species','Date.recolte')]
    df[,c(1,3)]<- scale(df[,c(1,3)])
    m0<-lmer(Ht.veg~(1|Species), data=df)
    m1<-lmer(Ht.veg~(Date.recolte|Species), data=df)
    anova(m0,m1)
    
    #Data: df
    # Models:
    #   m0: Ht.veg ~ (1 | Species)
    # m1: Ht.veg ~ (Date.recolte | Species)
    # Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
    # m0  3 1301.9 1314.3 -647.95   1295.9                        
    # m1  5 1305.9 1326.5 -647.95   1295.9     0      2          1
    
    # No effect of date, once species effect is accounted for. 
    
    # Test with residuals
    summary(lm(H.traits$Ht.veg~resid.Date))
    # Not significant
    
    
    # Looking for outliers in individual-level dataset
    ######################################################
    
    #X1 - Ht.veg - 2 outlier species - COCO and VIAL are tall
    boxplot(H.traits$Ht.veg,
            main = "vegetative height")
    dotchart(H.traits$Ht.veg, #i.e. cleveland dot chart
             xlab='range of data',
             ylab='Order of the data')
    
    dotchart(H.traits$Ht.veg, #i.e. cleveland dot chart - conditional on species
             groups=H.traits$Species,
             xlab='range of data',
             ylab='Order of the data')
    
    # X2 - Min.Root.Loca - no outliers - only 6 categories. 7 species with NAs
    boxplot(H.traits$Min.Root.Loca,
            main = "Minimum Root Location")
    
    H.traits[,c('Species','Min.Root.Loca')] # All VIAL, SAPU, RUID, LOCA, COCO & COAL are 'NA'
    
    dotchart(as.numeric(H.traits$Min.Root.Loca), #i.e. cleveland dot chart
             xlab='range of data',
             ylab='Order of the data')
    
    dotchart(as.numeric(H.traits$Min.Root.Loca), #i.e. cleveland dot chart - conditional on species
             groups=H.traits$Species,
             xlab='range of data',
             ylab='Order of the data')
    
    plot(as.numeric(H.traits$Max.Root.Loca)~as.numeric(H.traits$Min.Root.Loca))
    
    
    # X3 - Max.Root.Loca - no outliers - only 6 categories. 7 species with NAs
    boxplot(H.traits$Max.Root.Loca,
            main = "Maximum Root Location")
    
    H.traits[,c('Species','Max.Root.Loca')]  # All VIAL, SAPU, RUID, LOCA, COCO & COAL are 'NA'
    
    # X4 - Lamina.thck - CLCA, CLBO and GAPR, ERAM are 3x thicker than other species
    boxplot(H.traits$Lamina.thck,
            main = "Lamina.thck") # many outliers
    
    dotchart(H.traits$Lamina.thck, #i.e. cleveland dot chart
             xlab='range of data',
             ylab='Order of the data')  # two species are very thick
    
    dotchart(H.traits$Lamina.thck, #i.e. cleveland dot chart - conditional on species
             groups=H.traits$Species,
             xlab='range of data',
             ylab='Order of the data')
    
    # X5 - LMA - LYOB and GAPR have 3x higher LMA than other species
    boxplot(H.traits$LMA,
            main = "LmA") # lots of outliers
    
    dotchart(H.traits$LMA, #i.e. cleveland dot chart
             xlab='range of data',
             ylab='Order of the data',
             main='LMA') # two outlier species
    
    dotchart(H.traits$LMA, #i.e. cleveland dot chart - conditional on species
             groups=H.traits$Species,
             xlab='range of data',
             ylab='Order of the data') # LYOB, GAPR have high LMA
    
    H.traits[(H.traits$Species=='SAPU'|H.traits$Species=='COAL'),c('Plant.ID','LMA')]
    
    # X6 - LDMC - okay
    boxplot(H.traits$LDMC,
            main = "LDMC")
    
    dotchart(H.traits$LDMC, #i.e. cleveland dot chart
             xlab='range of data',
             ylab='Order of the data') # fine
    
    # X7 - Leaf.Area # SAPU has 20x larger leaf area than other sp. 
    
    boxplot(H.traits$Leaf.Area,
            main = "Leaf.Area") # outliers
    
    dotchart(H.traits$Leaf.Area, #i.e. cleveland dot chart
             xlab='range of data',
             ylab='Order of the data') # one species outlier
    
    dotchart(H.traits$Leaf.Area, #i.e. cleveland dot chart - conditional on species
             groups=H.traits$Species,
             xlab='range of data',
             ylab='Order of the data') # SAPU is outlier
    
    # X8 - Leaf.Mass.Frac - okay
    boxplot(H.traits$Leaf.Mass.Frac,
            main = "Supp.Mass.Frac") # fine
    
    dotchart(H.traits$Leaf.Mass.Frac, #i.e. cleveland dot chart
             xlab='range of data',
             ylab='Order of the data') # fine
    
    # X9 - Supp.Mass.Frac - Okay
    
    boxplot(H.traits$Supp.Mass.Frac,
            main = "Supp.Mass.Frac") # fine
    
    dotchart(H.traits$Supp.Mass.Frac, #i.e. cleveland dot chart
             xlab='range of data',
             ylab='Order of the data')
    
    # X10 - Rep.Mass.Frac # 0 inflated at individual level (not clustered by species)
    
    boxplot(H.traits$Rep.Mass.Frac,
            main = "Rep.Mass.Frac") # 0 inflated
    
    dotchart(H.traits$Rep.Mass.Frac, #i.e. cleveland dot chart
             xlab='range of data',
             ylab='Order of the data') # 0 inflated
    
    dotchart(H.traits$Rep.Mass.Frac, #i.e. cleveland dot chart - conditional on species
             groups=H.traits$Species,
             xlab='range of data',
             ylab='Order of the data')
    
    # X11 - Stor.Mass.Frac - 0 inflated at species level (turn into presence/absence?)
    boxplot(H.traits$Stor.Mass.Frac,
            main = "Rep.Mass.Frac") # 0 inflated
    
    dotchart(H.traits$Stor.Mass.Frac, #i.e. cleveland dot chart
             xlab='range of data',
             ylab='Order of the data') # 0 inflated
    
    dotchart(H.traits$Stor.Mass.Frac, #i.e. cleveland dot chart - conditional on species
             groups=H.traits$Species,
             xlab='range of data',
             ylab='Order of the data') # some species at 0, others not. 
    
    # X12 - F.Root.Diam - CYAC & EPHE have 4x thicker fine roots
    boxplot(H.traits$F.Root.Diam,
            main = "F.Root.Diam")
    
    dotchart(H.traits$F.Root.Diam, #i.e. cleveland dot chart
             xlab='range of data',
             ylab='Order of the data',
             main="fine Root thickness") # 2 sp with large diameter
    
    dotchart(H.traits$F.Root.Diam, #i.e. cleveland dot chart - conditional on species
             groups=H.traits$Species,
             xlab='range of data',
             ylab='Order of the data',
             main='SRL') # CYAC & EPHE are outliers 
    
    # X13 - SRL - okay
    
    boxplot(H.traits$SRL,
            main = "SRL") # many outliers
    
    dotchart(H.traits$SRL, #i.e. cleveland dot chart
             xlab='range of data',
             ylab='Order of the data') # 3 species with outliers
    
    dotchart(H.traits$SRL, #i.e. cleveland dot chart - conditional on species
             groups=H.traits$Species,
             xlab='range of data',
             ylab='Order of the data') 
    
    H.traits[(H.traits$Species=='STAM'|H.traits$Species=='STLA'),c('Plant.ID','SRL')]
    
    
    