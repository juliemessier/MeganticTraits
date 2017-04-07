# <<TABLE OF CONTENTS>>
# A - Setup trait data
  # A1 - Herbivory layer
    # A1.1 - Transform trait data for normality
    # A1.2 - Remove date effect on traits
    # A1.3 - Turn individual-level table into species-level table.
    # A1.4 - Add species mean traits (Seed size and mychorizae)
  # A2 - Canopy layer

# B- Setup abundance data

#<<WORKSPACES>>
wrk.dir<-("C:/Users/Julie/Desktop/Postdoc/Mégantic Traits/Workspaces/") # Workspaces
data.dir<-(("C:/Users/Julie/Desktop/Postdoc/Mégantic Traits/Data/")) # data
res.dir<-("C:/Users/Julie/Desktop/Postdoc/Mégantic Traits//Results/")  # Results
grp.dir<-("C:/Users/Julie/Desktop/Postdoc/Mégantic Traits//Graphs/")   # Graphs
fct.dir<-("C:/Users/Julie/Desktop/Postdoc/Mégantic Traits//Functions/") # Functions

#<<LIBRARIES>>
library(car) # for 
library(vegan) # for decostand (standardizing data)

# ==================================================================================

# A - Setup trait data ####

  # A1 - Herbivory layer
    # A1.1- Transform trait data for normality ####

    Megtraits<-read.csv(paste0(data.dir,'MegTraits_20170406.csv'))
    dim(Megtraits) #640 46
    str(Megtraits)
    
      # Simplify names
      names(Megtraits)[names(Megtraits)=="Lamina.thck.æm."] <- "Lamina.thck"
      names(Megtraits)[names(Megtraits)=='Vein.thck.æm.']<-'Vein.thck'
      names(Megtraits)[names(Megtraits)=='LMA.g.cm2.']<-'LMA'
      names(Megtraits)[names(Megtraits)=='LDMC.g.g.']<-'LDMC'
      names(Megtraits)[names(Megtraits)=='Leaf.Area.cm2.']<-'Leaf.Area'
      
      # create dataframes with traits I will work with
      
      H.traits<-Megtraits[Megtraits$Layer=='H',c("Plant.ID","Species","Layer","Date.recolte","Ht.veg","Min.Root.Loca","Max.Root.Loca",
                           "Lamina.thck","LMA","LDMC","Leaf.Area","Leaf.Mass.Frac","Supp.Mass.Frac",        
                           "Rep.Mass.Frac","Stor.Mass.Frac")]
      dim(H.traits)
      # 459 15
      
      C.traits<-Megtraits[Megtraits$Layer=='C',c("Plant.ID","Species","Layer","Date.recolte",
                          "Lamina.thck","LMA","LDMC","Leaf.Area")]
      dim(C.traits)
      # 181   8
      
      # Following Kerkhoff & Enquist (2009), all size- and growth-related traits are natural-log transformed, because 
      # they follow multiplicative processes:
      #  Ht.veg, Lamina.thck, LMA, LDMC, Leaf Area).
    
      # Test best transform for each trait ####
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
      
      # 3- Max.Root.Loca
      plot(density(na.omit(H.traits$Max.Root.Loca))) 
      # Totally bimodal. No transformation will fix that. 
      
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
      
      plot(density(na.omit(C.traits$Lamina.thck))) # Bimodal
      shapiro.test(C.traits$Lamina.thck) # 9.971e-16
      shapiro.test(log(na.omit(C.traits$Lamina.thck)))       # 6.935e-09
      powerTransform(na.omit(C.traits$Lamina.thck))# -0.3236387
      shapiro.test(C.traits$Lamina.thck^-0.3236387)   
      
      # Best Transform - log bc size-related & no large difference with PowerTransform exponent
      C.traits$Lamina.thck<-log(C.traits$Lamina.thck)
      
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
          
      plot(density(na.omit(C.traits$LDMC))) # weirdly shaped
      shapiro.test(C.traits$LDMC) #   9.555e-05
      shapiro.test(log(na.omit(C.traits$LDMC)))          #2.94e-09
      powerTransform(C.traits$LDMC)#1.276166 
      shapiro.test(C.traits$LDMC^1.276166)             # 0.000209
      plot(density(na.omit(C.traits$LDMC^1.276166)))
      plot(density(na.omit(log(C.traits$LDMC))))
      
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
        
        plot(density(na.omit(C.traits$Leaf.Area)))
        shapiro.test(C.traits$Leaf.Area) # 5.992e-16
        shapiro.test(log(na.omit(C.traits$Leaf.Area)))         #0.002092
        powerTransform(C.traits$Leaf.Area)#  -0.2679859
        shapiro.test(C.traits$Leaf.Area^ -0.2679859 )            #0.2017
        plot(density(na.omit(C.traits$Leaf.Area^-0.2679859 )))
        plot(density(na.omit(log(C.traits$Leaf.Area))))
         
        # Best Transform - size trait therefore log
        C.traits$Leaf.Area<-log(C.traits$Leaf.Area)
        
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
        
        H.traits$Ht.veg<-log(H.traits$Ht.veg)
        # Min.Root.Loca - no transforms
        # Max.Root.Loca - no transforms
        H.traits$Lamina.thck<-log(H.traits$Lamina.thck)
        C.traits$Lamina.thck<-log(C.traits$Lamina.thck)
        H.traits$LMA<-log(H.traits$LMA)
        C.traits$LMA<-log(C.traits$LMA)
        H.traits$LDMC<-log(H.traits$LDMC)
        H.traits$Leaf.Area<-log(H.traits$Leaf.Area)
        C.traits$Leaf.Area<-log(C.traits$Leaf.Area)
        #Leaf.Mass.Frac - no transforms
        H.traits$Supp.Mass.Frac<-H.traits$Supp.Mass.Frac^0.6124825
        H.traits$Rep.Mass.Frac<-(H.traits$Rep.Mass.Frac+0.001)^-0.564259
        H.traits$Stor.Mass.Frac<-(H.traits$Stor.Mass.Frac+0.001)^-0.07569593
        
    # A1.2- Remove date effect on traits  ####
    
        names(H.traits)
        #  [1] "Plant.ID"       "Species"        "Layer"          "Date.recolte"   "Ht.veg"        
        #  [6] "Min.Root.Loca"  "Max.Root.Loca"  "Lamina.thck"    "LMA"            "LDMC"          
        # [11] "Leaf.Area"      "Leaf.Mass.Frac" "Supp.Mass.Frac" "Rep.Mass.Frac"  "Stor.Mass.Frac"
    
        # list of trait names
        H.trait.names<-c("Ht.veg","Min.Root.Loca","Max.Root.Loca","Lamina.thck","LMA","LDMC","Leaf.Area",
                       "Leaf.Mass.Frac","Supp.Mass.Frac","Rep.Mass.Frac","Stor.Mass.Frac")
        length(H.trait.names) #11
        
        C.trait.names<-c("Lamina.thck","LMA","LDMC","Leaf.Area")
        length(C.trait.names) #4
        
        # Change caledar dates into julian dates
        H.traits$Date.recolte<- as.numeric(format(as.Date(H.traits$Date.recolte, format = "%m/%d/%Y"),"%j"))
        C.traits$Date.recolte<- as.numeric(format(as.Date(C.traits$Date.recolte, format = "%m/%d/%Y"),"%j"))
    
        # loop testing for date effects and replacing with regression residuals
        for (t in H.trait.names){ 
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
        
        # [1] "Ht.veg"                 "R2=  0.202951260763896" "Sign=  1"              
        # [1] "Min.Root.Loca"          "R2=  0.0435787342801673" "Sign=  1"               
        # [1] "Max.Root.Loca"          "R2=  0.0338570003191077" "Sign=  1"               
        # [1] "Lamina.thck"            "R2=  0.272643266243074" "Sign=  -1"             
        # [1] "LDMC"                   "R2=  0.179448937492485" "Sign=  1"              
        # [1] "Leaf.Area"              "R2=  0.0303995275437033" "Sign=  1"               
        # [1] "Leaf.Mass.Frac"         "R2=  0.096874667091267" "Sign=  1"              
        # [1] "Supp.Mass.Frac"         "R2=  0.0916507749900383" "Sign=  1" 
        # [1] "Stor.Mass.Frac"         "R2=  0.178232187256645" "Sign=  1"
        
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
        # [1] "LMA"        "R2=  0.072" "Sign=  -1" 
        # [1] "LDMC"       "R2=  0.043" "Sign=  -1" 
        # [1] "Leaf.Area"  "R2=  0.052" "Sign=  -1"
        
        str(H.traits)
          # 'data.frame':	459 obs. of  15 variables:
          # $ Plant.ID      : Factor w/ 640 levels "ABBA1","ABBA10",..: 226 135 141 230 231 139 229 140 228 138 ...
          # $ Species       : Factor w/ 75 levels "ABBA","ACPE",..: 27 16 16 27 27 16 27 16 27 16 ...
          # $ Layer         : Factor w/ 2 levels "C","H": 2 2 2 2 2 2 2 2 2 2 ...
          # $ Date.recolte  : Factor w/ 60 levels "5/11/2016","5/12/2016",..: 7 7 11 11 11 11 11 11 11 11 ...
          # $ Ht.veg        : num  NA -0.681 -0.696 -0.666 -0.651 ...
          # $ Min.Root.Loca : num  2.026 -0.971 -0.372 -0.372 1.427 ...
          # $ Max.Root.Loca : num  1.975 -1.025 -0.425 -0.425 1.375 ...
          # $ Lamina.thck   : num  2.65 2.52 2.79 3.88 4.3 ...
          # $ LMA           : num  0.218 0.218 -0.463 0.331 0.955 ...
          # $ LDMC          : num  -1.123 -1.164 -2.499 -1.185 -0.874 ...
          # $ Leaf.Area     : num  -0.3847 -0.5896 -0.589 -0.2618 0.0213 ...
          # $ Leaf.Mass.Frac: num  -0.372 -2.179 -2.222 0.532 -0.802 ...
          # $ Supp.Mass.Frac: num  -0.85 -1.6 -1.85 -1.18 -1.43 ...
          # $ Rep.Mass.Frac : num  -0.446 0.682 0.118 -0.446 1.105 ...
          # $ Stor.Mass.Frac: num  0.6656 2.3183 2.6643 0.0506 0.9346 ...
        H.traits$Min.Root.Loca<-as.factor(H.traits$Min.Root.Loca)
        H.traits$Max.Root.Loca<-as.factor(H.traits$Max.Root.Loca)
     
    # A1.3 - Turn individual-level table into species-level table.  ####
        #apply(H.traits[,5:11], 2, function(x) tapply(x, H.traits$Species, mean,na.rm=T))
        
        #first standardize all variables
        H.traits[8:15]<-decostand(H.traits[8:15],method='standardize', margin=2)
        C.traits[5:8]<-decostand(C.traits[5:8],method='standardize',margin=2)
        
        sp.H.traits<-as.data.frame(aggregate(H.traits[,5:11],by=list(H.traits$Species),mean,na.rm=T))
        sp.H.traits[,2:8]<-round(sp.H.traits[2:8],digits=3)
        sp.H.traits[is.na(sp.H.traits)] <-NA
        rownames(sp.H.traits)<-sp.H.traits$Group.1
        sp.H.traits<-sp.H.traits[,-1]
        sp.H.traits<-sp.H.traits[order(rownames(sp.H.traits)),]
        dim(sp.H.traits) # 54  8
        
        sp.C.traits<-as.data.frame(aggregate(C.traits[,5:8],by=list(C.traits$Species),mean,na.rm=T))
        sp.C.traits[,2:5]<-round(sp.C.traits[2:5],digits=3)
        sp.C.traits[is.na(sp.C.traits)] <-NA
        rownames(sp.C.traits)<-sp.C.traits$Group.1
        sp.C.traits<-sp.C.traits[,-1]
        sp.C.traits<-sp.C.traits[order(rownames(sp.C.traits)),]
        dim(sp.C.traits) # 24  4
        
    # A1.4 - Add species mean traits (Seed size and mychorizae) ####
        # Deleted "FRAM" row by hand - We did not collect roots on trees... what is this row?
        myc<-read.csv(paste0(data.dir,'myc.csv'))
        head(myc)
        dim(myc) #43 2
        rownames(myc)<-myc$species
        myc<-myc[-1]
        myc<-decostand(myc,method='standardize',margin=2)
        dim(myc) #43 1
    
        sp.H.traits<-merge(sp.H.traits,myc, by="row.names",all=T)
        head(sp.H.traits)
        rownames(sp.H.traits)<-sp.H.traits$Row.names
        sp.H.traits$Row.names<-NULL
        dim(sp.H.traits) # 53 8
        
        
    
# B- Setup abundance data ####

abund<-read.csv(paste0(data.dir,'sp.abund.csv'))
rownames(abund)<-abund$Species
abund$Species<-NULL
abund[,'presence.change']<-round(abund$pct.plot.present.2012/abund$pct.plot.present.1970,digits=1)
abund[,'abund.change']<-round(abund$avg.abundance.2012/abund$avg.abundance.1970,digits=1)
str(abund)
dim(abund) # 125 7

#Look at Frequency distribution of response variables
# transform reponse variable bc model
plot(density(abund$presence.change)) # lots of small values
plot(density(log(abund$presence.change))) # looks way better
plot(density(abund$abund.change)) # lots of small values
plot(density(log(abund$presence.change)))# looks way better

  # problem with infinity values
  rownames(abund)[which(abund$avg.abundance.1970==0)] # replace infinity values
  # [1] "CAAL" "CACR" "CAPE" "CATH" "COMA" "DEPU" "GAPR" "GATE" "HYAM" "JUTE" "LICO" "LYCL" "LYUN" 
  # [17] "MOUN" "PAQU" "POGR" "RARE" "TAOF"
  abund$avg.abundance.1970[order(abund$avg.abundance.1970)][15:25]# What is the minimum value in dataset?
  #0.00 0.00 0.00 0.00 0.01 0.01 0.01 0.01 0.01 0.01 0.01
  abund[c("CAAL","CACR","CAPE","CATH","COMA","DEPU","GAPR","GATE",
          "HYAM","JUTE","LICO","LYCL","LYUN"),'avg.abundance.1970']<-0.01 # replace 0s with minimum value
  
  rownames(abund)[which(abund$pct.plot.present.1970==0)] 
  # [1] "CAAL" "CACR" "CAPE" "CATH" "COMA" "DEPU" "GAPR" "GATE" "HYAM" "JUTE" "LICO" "LYCL" "LYUN" 
  # [17] "MOUN" "PAQU" "POGR" "RARE" "TAOF"
  abund$pct.plot.present.1970[order(abund$pct.plot.present.1970)][15:25]
  # 0.0 0.0 0.0 0.0 2.1 2.1 2.1 2.1 2.1 2.1 2.1 
  abund[c("CAAL","CACR","CAPE","CATH","COMA","DEPU","GAPR","GATE",
          "HYAM","JUTE","LICO","LYCL","LYUN"),'pct.plot.present.1970']<-2.1
  
  #Re-calculate presence and abundance changes
  abund[,'presence.change']<-round(abund$pct.plot.present.2012/abund$pct.plot.present.1970,digits=1)
  abund[,'abund.change']<-round(abund$avg.abundance.2012/abund$avg.abundance.1970,digits=1)
  
  #Look at Frequency distribution of response variables
  plot(density(abund$presence.change)) # lots of small values
  plot(density(abund$abund.change)) # lots of small values
  
  #DEPU has values that are way too extreme - absent in 1970 but present in 2012. Delete. 
  abund<-abund[rownames(abund)!='DEPU',]
  abund<-abund[rownames(abund)!="GAPR",]
  abund<-abund[rownames(abund)!="GATE",]
  abund<-abund[rownames(abund)!="POGR",]
  
  dim(abund) # 122 7
  
# A6 - Make trait and abundance dataframes correspond based on shared species

# Herbivory layer
H.abund<-abund[which(rownames(abund)%in%rownames(sp.H.traits)),c('presence.change','abund.change')]
H.abund<-H.abund[order(rownames(H.abund)),]
H.abund
dim(H.abund) #45 2

# Canopy layer
C.abund<-abund[which(rownames(abund)%in%rownames(sp.C.traits)),c('presence.change','abund.change')]
C.abund<-C.abund[order(rownames(C.abund)),]
C.abund
dim(C.abund) #21 2

length(rownames(sp.H.traits))
# 53
length(rownames(H.abund))
# 45
sp.H.traits<-sp.H.traits[which(rownames(sp.H.traits)%in%rownames(H.abund)),]
dim(sp.H.traits) #45 8
all(rownames(H.abund)==rownames(sp.H.traits)) # TRUE

length(rownames(sp.C.traits))
#22
length(rownames(C.abund))
#21
sp.C.traits<-sp.C.traits[which(rownames(sp.C.traits)%in%rownames(C.abund)),]
all(rownames(sp.C.traits)==rownames(C.abund)) # TRUE

#Any pairwise interactions?
pairs(c(sp.H.traits,H.abund),panel=panel.smooth)
  #yes, the two rooting depths, and rooting depth-myc.frac. and LMA-LDMC

#Any important interactions on response variable?

library(tree)
model<-tree(H.abund$abund.change~.,data=sp.H.traits)
plot(model)
text(model)
  # Myc.frac is most important variable, 
  # for those with high myc.frac, LDMC matters,
  # for those with high LDMC, root.location matters.

