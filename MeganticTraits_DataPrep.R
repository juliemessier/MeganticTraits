# <<TABLE OF CONTENTS>>
# A - Setup trait data
  # A1 - Herbivory layer
    # A1.1 - Transform trait data for normality
    # A1.2 - Remove date effect on traits
    # A1.3 - Turn individual-level table into species-level table.
    # A1.4 - Add species mean traits (Seed size and mychorizae)
  # A2 - Canopy layer
    # A2.1 - Transform trait data for normality
    # A2.2 - Remove date effect on traits
    # A2.3 - Turn individual-level table into species-level table.
    # A2.4 - Add species mean traits (Seed size and mychorizae)
# B- Setup abundance data
# B2- Create abundance dataframe & new variables
# B3 - Problem with infinity values ####



#<<WORKSPACES>>
wrk.dir<-("C:/Users/Julie/Desktop/Postdoc/Megantic Traits/Workspaces/") # Workspaces
data.dir<-(("C:/Users/Julie/Desktop/Postdoc/Megantic Traits/Data/")) # data
res.dir<-("C:/Users/Julie/Desktop/Postdoc/Megantic Traits/Results/")  # Results
grp.dir<-("C:/Users/Julie/Desktop/Postdoc/Megantic Traits/Graphs/")   # Graphs
fct.dir<-("C:/Users/Julie/Desktop/Postdoc/Megantic Traits/Functions/") # Functions

#<<LIBRARIES>>
library(car) # for powerTransform
library(vegan) # for decostand (standardizing data)

# ==================================================================================#

# A - Setup trait data ####
# ========================#

  # A1 - HERBIVORY LAYER ####
  #=======================#

    # A1.1- Transform trait data for normality ####

    Megtraits<-read.csv(paste0(data.dir,'MegTraits_20170406.csv'))
    dim(Megtraits) #640 46
    str(Megtraits)
    Megtraits$Min.Root.Loca<-as.ordered(Megtraits$Min.Root.Loca) # remains a factor, but is ordered (ordinal)
    Megtraits$Max.Root.Loca<-as.ordered(Megtraits$Max.Root.Loca) # remains a factor, but is ordered (ordinal)
    
      # Simplify names
      names(Megtraits)[names(Megtraits)=="Lamina.thck.æm."] <- "Lamina.thck"
      names(Megtraits)[names(Megtraits)=='Vein.thck.æm.']<-'Vein.thck'
      names(Megtraits)[names(Megtraits)=='LMA.g.cm2.']<-'LMA'
      names(Megtraits)[names(Megtraits)=='LDMC.g.g.']<-'LDMC'
      names(Megtraits)[names(Megtraits)=='Leaf.Area.cm2.']<-'Leaf.Area'
      
      # create trait dataframe with traits I will work with
      
      H.traits<-Megtraits[Megtraits$Layer=='H',c("Plant.ID","Species","Layer","Date.recolte","Ht.veg","Min.Root.Loca","Max.Root.Loca",
                           "Lamina.thck","LMA","LDMC","Leaf.Area","Leaf.Mass.Frac","Supp.Mass.Frac",        
                           "Rep.Mass.Frac","Stor.Mass.Frac")]
      
      dim(H.traits)
      # 459 15
      H.traits<-droplevels(H.traits)
      str(H.traits)
      
        # 'data.frame':	459 obs. of  15 variables:
        #   $ Plant.ID      : Factor w/ 459 levels "ARNU1","ARNU2",..: 155 64 70 159 160 68 158 69 157 67 ...
        # $ Species       : Factor w/ 51 levels "ARNU","ARTR",..: 19 8 8 19 19 8 19 8 19 8 ...
        # $ Layer         : Factor w/ 1 level "H": 1 1 1 1 1 1 1 1 1 1 ...
        # $ Date.recolte  : Factor w/ 49 levels "5/11/2016","5/12/2016",..: 7 7 11 11 11 11 11 11 11 11 ...
        # $ Ht.veg        : num  NA 5 4.5 5.5 6 4.5 8.5 4.5 4.5 7 ...
        # $ Min.Root.Loca : Ord.factor w/ 6 levels "0"<"1"<"2"<"3"<..: 6 1 2 2 5 2 5 1 2 4 ...
        # $ Max.Root.Loca : Ord.factor w/ 6 levels "0"<"1"<"2"<"3"<..: 6 1 2 2 5 2 5 1 2 5 ...
        # $ Lamina.thck   : num  315 305 325 402 432 ...
        # $ LMA           : num  0.0029 0.0029 0.0017 0.0031 0.0042 0.0024 0.0031 0.0023 0.0035 0.0023 ...
        # $ LDMC          : num  0.1 0.098 0.051 0.097 0.113 0.087 0.111 0.125 0.103 0.106 ...
        # $ Leaf.Area     : num  8.67 1.34 1.36 13.07 23.2 ...
        # $ Leaf.Mass.Frac: num  0.48 0.06 0.05 0.69 0.38 0.1 0.54 0.08 0.6 0.06 ...
        # $ Supp.Mass.Frac: num  0.14 0.05 0.02 0.1 0.07 0.05 0.09 0.06 0.04 0.06 ...
        # $ Rep.Mass.Frac : num  0 0.08 0.04 0 0.11 0.07 0 0.09 0 0.12 ...
        # $ Stor.Mass.Frac: num  0.37 0.8 0.89 0.21 0.44 0.78 0.37 0.76 0.36 0.76 ...
      
      
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
        
    # A1.2- Remove date effect on traits  ####
    # =======================================#
    
        names(H.traits)
        #  [1] "Plant.ID"       "Species"        "Layer"          "Date.recolte"   "Ht.veg"        
        #  [6] "Min.Root.Loca"  "Max.Root.Loca"  "Lamina.thck"    "LMA"            "LDMC"          
        # [11] "Leaf.Area"      "Leaf.Mass.Frac" "Supp.Mass.Frac" "Rep.Mass.Frac"  "Stor.Mass.Frac"
    
        # list of trait names
        H.trait.names<-c("Ht.veg","Min.Root.Loca","Max.Root.Loca","Lamina.thck","LMA","LDMC","Leaf.Area",
                       "Leaf.Mass.Frac","Supp.Mass.Frac","Rep.Mass.Frac","Stor.Mass.Frac")
        length(H.trait.names) #11
        
        save(H.trait.names,file=paste0(wrk.dir,'list.trait.names.herbivory.layer.Rdata'))
        
        # Change calendar dates into julian dates
        H.traits$Date.recolte<- as.numeric(format(as.Date(H.traits$Date.recolte, format = "%m/%d/%Y"),"%j"))
        
        # loop testing for date effects and replacing with regression residuals
        for (t in H.trait.names[-c(2:3)]){  # not for min. and max. root location, because they are factors
          x<-lm(H.traits[[t]]~H.traits$Date.recolte,na.action=na.exclude)  # # regress the trait against the julian date
          
          if (summary(x)$adj.r.squared > 0.02 & summary(x)[["coefficients"]][2,4] < 0.05) # if regression R2 > 0.02 AND it is statistically significant (with P-value<0.05), go through this next loop
            print(c(names(H.traits[t]),
                    paste( "R2= ",round(summary(x)$adj.r.squared,digits=3)),
                    paste('Sign= ',sign(summary(x)[["coefficients"]][2,1])))) 
          
          for (i in 1:nrow(H.traits[t])){                                                   # if each row is numeric (not an NA)
            if (is.numeric(H.traits[t][i,]))
            {H.traits[t][i,]<-resid(x)[i]}                                              # then replace the value with the residual
          }                                                                           # else, do nothing to the cells with value of NA
        }
        
        # [1] "Ht.veg"                 "R2=  0.203" "Sign=  1"              
        # [1] "Lamina.thck"            "R2=  0.273" "Sign=  -1"             
        # [1] "LDMC"                   "R2=  0.179" "Sign=  1"              
        # [1] "Leaf.Mass.Frac"         "R2=  0.097" "Sign=  1"              
        # [1] "Supp.Mass.Frac"         "R2=  0.092" "Sign=  1" 
        # [1] "Stor.Mass.Frac"         "R2=  0.178" "Sign=  1"
        
        str(H.traits)
          # 'data.frame':	459 obs. of  15 variables:
          # $ Plant.ID      : Factor w/ 640 levels "ABBA1","ABBA10",..: 226 135 141 230 231 139 229 140 228 138 ...
          # $ Species       : Factor w/ 75 levels "ABBA","ACPE",..: 27 16 16 27 27 16 27 16 27 16 ...
          # $ Layer         : w/ 1 level "H": 1 1 1 1 1 1 1 1 1 1 ...
          # $ Date.recolte  : num  123 123 124 124 124 124 124 124 124 124 ...
          # $ Ht.veg        : num  NA -0.68 -0.799 -0.598 -0.511 ...
          # $ Min.Root.Loca : Factor w/ 6 levels "0","1","2","3",..: 6 1 2 2 5 2 5 1 2 4 ...
          # $ Max.Root.Loca : Factor w/ 6 levels "0","1","2","3",..: 6 1 2 2 5 2 5 1 2 5 ...
          # $ Lamina.thck   : num  0.613 0.582 0.652 0.864 0.936 ...
          # $ LMA           : num  0.409 0.409 -0.127 0.473 0.777 ...
          # $ LDMC          : num  -0.1908 -0.211 -0.8713 -0.2284 -0.0757 ...
          # $ Leaf.Area     : num  0.111 -1.756 -1.747 0.516 1.09 ...
          # $ Leaf.Mass.Frac: num  0.0392 -0.3808 -0.3933 0.2467 -0.0633 ...
          # $ Supp.Mass.Frac: num  -0.0361 -0.1764 -0.2463 -0.0933 -0.1412 ...
          # $ Rep.Mass.Frac : num  15.2 -29.9 -28 15.2 -30.6 ...
          # $ Stor.Mass.Frac: num  -0.1009 -0.1619 -0.1748 -0.0586 -0.1196 ...
        
        save(H.traits,file=paste0(wrk.dir,"Herbaceous.Layer.Traits.Unstandardized.RData"))
        
        # Standardize all variables
        H.traits[8:15]<-decostand(H.traits[8:15],method='standardize', margin=2)
        save(H.traits,file=paste0(wrk.dir,"Herbaceous.Layer.Traits.Standardized.RData"))
     
    # A1.3 - Turn individual-level table into species-level table  ####
    #=================================================================#    
        #apply(H.traits[,5:11], 2, function(x) tapply(x, H.traits$Species, mean,na.rm=T))

        H.traits$Min.Root.Loca<-as.numeric(H.traits$Min.Root.Loca) # can't take mean or median of ordinal var. 
        H.traits$Max.Root.Loca<-as.numeric(H.traits$Max.Root.Loca) # can't take mean or median of ordinal var. 
        
        sp.H.traits<-as.data.frame(aggregate(H.traits[,5:11],by=list(H.traits$Species),mean,na.rm=T))
        sp.H.traits[,2:8]<-round(sp.H.traits[2:8],digits=3)
        sp.H.traits[is.na(sp.H.traits)] <-NA
        rownames(sp.H.traits)<-sp.H.traits$Group.1
        sp.H.traits<-sp.H.traits[,-1]
        sp.H.traits<-sp.H.traits[order(rownames(sp.H.traits)),]
        dim(sp.H.traits) # 51  7
        
    # A1.4 - Add species mean traits (Seed size and mychorizae) ####
    #=================================================================#
        
        # Deleted "FRAM" row by hand - We did not collect roots on trees... what is this row?
        myc<-read.csv(paste0(data.dir,'myc.csv'))
        head(myc)
        dim(myc) #43 2
        rownames(myc)<-myc$species
        myc<-myc[-1]
        plot(density(myc$myc.frac)) # left skewed
        myc<-decostand(myc,method='standardize',margin=2)
        dim(myc) #43 1
    
        sp.H.traits<-merge(sp.H.traits,myc, by="row.names",all=T)
        head(sp.H.traits)
        rownames(sp.H.traits)<-sp.H.traits$Row.names
        sp.H.traits$Row.names<-NULL
        dim(sp.H.traits) # 53 8
        
        save(sp.H.traits,file=paste0(wrk.dir,'species-level.traits.Herbaceous.layer.Rdata'))
        
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
    