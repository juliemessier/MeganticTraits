# <<TABLE OF CONTENTS>> ####
# A- Setup abundance data
  # A1- Create abundance dataframe & new variables 
  # A2 - Problem with infinity values
  # A3 - Look at temporal patterns
  
# B - Setup trait data    
  # B1 - HERBACEOUS LAYER
    # B1.1 - Remove date effect on traits
    # B1.2 - Turn individual-level table into species-level table
    # B1.3 - Add species mean traits (Seed size and mychorizae)
  # B2 - CANOPY LAYER 

# C - Data Exploration (Following Highlands Stats course) 
  # C1 - Outliers on Y and X
  # C2 - Homogeneity (homoscedasticity) of Y 
  # C3 - Normality of Y
  # C4 - Zero trouble (Y) ?
  # C5 - Collinearity  (Xs)
    # C5.1 Deal with colinearity using VIFs
    # C5.2 Deal with colinearity using PCs 
  # C6 Relationship between Y and X linear? 
  # C7 Interactions?


#<<WORKSPACES>>
wrk.dir<-("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Traits/Workspaces/") # Workspaces
data.dir<-(("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Traits/Data/")) # data
res.dir<-("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Traits/Results/")  # Results
grp.dir<-("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Traits/Graphs/")   # Graphs
fct.dir<-("C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Traits/Functions/") # Functions

#<<LIBRARIES>>
library(car) # for powerTransform
library(vegan) # for decostand(), rda(), etc. 
library(lattice) # For fancy multipanel graphs
source(paste0(wrk.dir,"HighstatLibV10.R")) # to make fancy graphs

# ==================================================================================#

# A- Setup response variables ####
# ===========================#

# A1 - Setup Abundance & Occurence Data ####
#-----------------------------#

# A1.1- Create abundance dataframe & new variables ####
#=====================================================#

abund<-read.csv(paste0(data.dir,'sp.abund.csv'))
str(abund)

rownames(abund)<-abund$Species
colnames(abund)[1]<-'SpCode'

# Abundance and Occurence ratios

abund[,'occurence.ratio']<-round(abund$pct.plot.present.2012/abund$pct.plot.present.1970,digits=1)
abund[,'abundance.ratio']<-round(abund$avg.abundance.2012/abund$avg.abundance.1970,digits=1)
abund$abundance.ratio


# Ratios of log Abundance & Occurence

# abund[,'ratio.log.coverage']<-round(log(abund$pct.plot.present.2012)/log(abund$pct.plot.present.1970),digits=3)
# abund[,'ratio.log.abund']<-round(log(abund$avg.abundance.2012)/log(abund$avg.abundance.1970),digits=3)
# abund$ratio.log.abund
# plot(density(abund$ratio.log.abund))
# # Outliers?
# rownames((abund)[abund$ratio.log.abund<=-10,])# IMCA, POTR
# rownames((abund)[abund$ratio.log.abund>=20,]) # "ACPE"  "ALTR"  "DRGO"  "EPAN"  "ERST"  "EUPMA" "FRNI"  "GORE"  "LIBO"  "TSCA"

# If this is an issue - link function of glms (negative binomial) made to deal with lack of normality of response variable. 

str(abund)
# data.frame':	125 obs. of  8 variables:
# $ Species              : Factor w/ 125 levels "ABBA","ACPE",..: 1 2 3 4 5 6 7 8 9 10 ...
# $ Layer                : Factor w/ 2 levels "C","H": 1 1 1 1 2 2 2 2 2 2 ...
# $ pct.plot.present.1970: num  87.5 52.1 18.8 54.2 89.6 14.6 2.1 58.3 4.2 12.5 ...
# $ pct.plot.present.2012: num  85.4 70.8 29.2 62.5 79.2 4.2 0 75 6.3 27.1 ...
# $ avg.abundance.1970   : num  27.2 1.03 5.12 25.97 6.73 ...
# $ avg.abundance.2012   : num  20.27 4.06 3.82 24.3 4.85 ...
# $ occurence.ratio      : num  1 1.4 1.6 1.2 0.9 0.3 0 1.3 1.5 2.2 ...
# $ abund.ratio          : num  0.7 3.9 0.7 0.9 0.7 0.3 0 1.1 1.5 2.1 ...



save(abund,file=paste0(wrk.dir,'Species.abundances.full.data.Rdata'))


# A1.2 - Remove infinity values ####
#===============================#

# Species absent in 1974 that are present in 2012 give a ratio of infinity
# Which species have ratios of infinity?

rownames(abund[abund$abund.ratio=='Inf',]) 
# [1] "CAAL" "CACR" "CAPE" "CATH" "COMA" "DEPU" "GAPR" "GATE" "HYAM" "JUTE" "LICO" "LYCL" "LYUN" "MOUN" "PAQU" "POGR"
# [17] "RARE" "TAOF
# loosing 18 data points to "Inf"

dim(abund[abund$occurence.ratio=='Inf',])
#18 
# we loose 18 species to infinity

# are the infinity values for abundance and occurence the same? 
all(rownames(abund[abund$occurence.ratio=='Inf',])==rownames(abund[abund$abund.ratio=='Inf',])) 
# TRUE


# Try replacing inf values with a number, 
# where value for 1970 is half of the smallest observed

# What is the minimum value in dataset?
abund$avg.abundance.1970[order(abund$avg.abundance.1970)][15:25]
#0.00 0.00 0.00 0.00 0.01 0.01 0.01 0.01 0.01 0.01 0.01

# try replacing 0s with half of minimum value
abund$non0.avg.abundance.1970<-abund$avg.abundance.1970
abund[c("CAAL","CACR","CAPE","CATH","COMA","DEPU","GAPR","GATE",
        "HYAM","JUTE","LICO","LYCL","LYUN"),'non0.avg.abundance.1970']<-0.01/2 

# Re-calculate abundance
non0.abund.change<-round(abund$avg.abundance.2012/abund$non0.avg.abundance.1970,digits=2)
plot(density(non0.abund.change)) 
# now ranges from 0-400, with a single value at 400. 
# Therefore, it is a bad idea to replace 1970's 0s with a tiny number.

# Delete this column
abund$non0.avg.abundance.1970<-NULL

#save new abundance dataset without INF values
abund.c<-abund[abund$abundance.ratio!='Inf',]

dim(abund.c) # 107 9

# look at distribution of values
plot(density(abund$abund.ratio)) # right skewed
plot(density(log(abund$abund.ratio+0.01))) # bimodal

plot(density(abund$occurence.ratio)) # right skewed
plot(density(log(abund$occurence.ratio))) # normal-ish. 

save(abund.c,file=paste0(wrk.dir,'Species.abundances.Inf.removed.Rdata')) 


# A1.3 - Look at temporal patterns ####
# ====================================#

par(mfrow=c(1,2))

plot(abund.c[abund.c$Layer=='H','avg.abundance.2012']~abund.c[abund.c$Layer=='H','avg.abundance.1970'],type='n',
     xlab='mean abundance 1970', ylab='mean abundance 2012',main='Abundance, herbaceous')
text(abund.c[abund.c$Layer=='H','avg.abundance.1970'],abund.c[abund.c$Layer=='H','avg.abundance.2012'],
     labels=rownames(abund.c)[which(abund.c$Layer=='H')],cex=0.8)
abline(0,1)

# DRIN and OXMO are "outliers" : abundance ratios much larger/smaller than other species 

#on log-scale, much better
plot(log(abund.c[abund.c$Layer=='H','avg.abundance.2012'])~log(abund.c[abund.c$Layer=='H','avg.abundance.1970']),type='n',
     xlab='log(av abundance 1970)', ylab='log(avg abundance 2012)',main='Abundance, log herbaceous')
text(log(abund.c[abund.c$Layer=='H','avg.abundance.1970']),log(abund.c[abund.c$Layer=='H','avg.abundance.2012']),
     labels=rownames(abund.c)[which(abund.c$Layer=='H')],cex=0.8)
abline(0,1)

# for coverage : 

plot(abund.c[abund.c$Layer=='H','pct.plot.present.2012']~abund.c[abund.c$Layer=='H','pct.plot.present.1970'],type='n',
     xlab='pct.plot.present.1970', ylab='pct.plot.present.2012',main='coverage herbaceous')
text(abund.c[abund.c$Layer=='H','pct.plot.present.1970'],abund.c[abund.c$Layer=='H','pct.plot.present.2012'],
     labels=rownames(abund.c)[which(abund.c$Layer=='H')],cex=0.8)
abline(0,1)

# DRIN and OXMO are "outliers" : abundance ratios much larger/smaller than other species 

#on log-scale, much better
plot(log(abund.c[abund.c$Layer=='H','pct.plot.present.2012'])~log(abund.c[abund.c$Layer=='H','pct.plot.present.1970']),type='n',
     xlab='log(pct.plot.present.1970)', ylab='log(pct.plot.present.2012)',main='coverage, log herbaceous')
text(log(abund.c[abund.c$Layer=='H','pct.plot.present.1970']),log(abund.c[abund.c$Layer=='H','pct.plot.present.2012']),
     labels=rownames(abund.c)[which(abund.c$Layer=='H')],cex=0.8)
abline(0,1)

# A2 - Setup elevation change data ####
#=====================================#

  # A2.1 - Create elevation shift data ####
  #=====================================================#
  elev<-read.csv(paste0(data.dir,'Elevation_shifts.csv'))
  str(elev)
  # 'data.frame':	58 obs. of  8 variables:
  # $ CodeSP    : Factor w/ 58 levels "ABIBAL","ACEPEN",..: 6 7 8 9 10 13 14 15 16 18 ...
  # $ latin.name: Factor w/ 58 levels "abies balsamea",..: 6 7 31 8 19 11 12 13 14 16 ...
  # $ my.code   : Factor w/ 58 levels "ABBA","ACPE",..: 6 7 31 8 18 11 12 13 17 15 ...
  # $ PotStrata : Factor w/ 2 levels "Canopy","Understory": 2 2 2 2 2 2 2 2 2 2 ...
  # $ SpElev1970: int  613 632 699 616 605 587 590 706 650 803 ...
  # $ SpElev2012: int  692 611 680 623 584 655 655 870 774 897 ...
  # $ ElevMean  : num  652 622 690 620 594 ...
  # $ ElevDif   : int  79 -21 -19 7 -21 68 65 164 124 94 ...
  dim(elev)
  # 58 8
  
  # Remove unnecessary columns
  elev$CodeSP<-NULL
  elev$ElevMean<-NULL
  rownames(elev)<-elev$my.code
  
  # A2.2 - Remove infinity values #### 
  # ====================================#
  # Cleaning already done by José Savage - no 0s or Inf in dataset
  
  # A2.3 - Look at temporal patterns ####
  # ====================================#
  
  par(mfrow=c(1,2))
  
  plot(elev[elev$PotStrata=='Understory','SpElev2012']~elev[elev$PotStrata=='Understory','SpElev1970'],type='n',
       xlab='mean elevation 1970', ylab='mean elevation 2012',main='Elevation Shift,\ Understorey')
  text(elev[elev$PotStrata=='Understory','SpElev1970'],elev[elev$PotStrata=='Understory','SpElev2012'],
       labels=rownames(elev)[which(elev$PotStrata=='Understory')],cex=0.8)
  abline(0,1)
  
  # No outliers ! 
  
  save(elev,file=paste0(wrk.dir,'Species.elevation.shifts.Rdata'))

# A3 - Combine Abundance, Occurence and Elevation shifts into one dataset ####
#============================================================================#
  
  colnames(abund)
  # [1] "SpCode"               "Layer"                 "pct.plot.present.1970" "pct.plot.present.2012"
  # [5] "avg.abundance.1970"    "avg.abundance.2012"    "occurence.ratio"       "abundance.ratio"  
  colnames(elev)
  # [1] "latin.name" "my.code"    "PotStrata"  "SpElev1970" "SpElev2012" "ElevDif" 

  # Match colnames
  colnames(elev)[1:3]<-c('LatinName','SpCode','Layer')
  
  # Match layer names
  levels(abund$Layer)[levels(abund$Layer)=="H"] <- "Understory"
  levels(abund$Layer)[levels(abund$Layer)=="C"] <- "Canopy"
  
  # How many species in 'abund' dataset absent from 'elev' dataset?
  length(abund$SpCode[which(!abund$SpCode%in%elev$SpCode)])
  # 68 species! 
  # How many Understory species in abund dataset absent from elev dataset?
  dim(abund[!(abund$SpCode%in%elev$SpCode) & abund$Layer=='Understory',])
  # 58 species
  
  # How many species in "elev" dataset absent from 'abund' dataset?
  length(elev$SpCode[which(!elev$SpCode%in%abund$SpCode)])
  # 1 
  elev$SpCode[which(!elev$SpCode%in%abund$SpCode)]
  # BEPO
  
  Ys<-merge(abund,elev,by="SpCode",all=T)
  dim(Ys)
  # 126 13
  dim(abund)
  # 125 13 - one new column from elev that wasn't in abund
  
  # reorder columns
  colnames(Ys)
  Ys<-Ys[,c(1,9,2,10,5,6,3,4,11,12,8,7,13)]
  View(Ys)
  # I checked there are no contradictions between the two layer columns
  
  # Replace NAs in Layer.x
  Ys$Layer.x[126]<-Ys$Layer.y[126]

  # Remove duplicate column
  Ys$Layer.y<-NULL
  colnames(Ys)[3]<-'Layer'
  
  save(Ys,file=paste0(wrk.dir,'All.Response.variables.RData')) #Ys
  
  Ys.c<-Ys[Ys$abundance.ratio!='Inf',]
  dim(Ys.c)
  # 108.12
  
  save(Ys.c,file=paste0(wrk.dir,'All.Response.variables.Inf.removed.RData')) 
    
# B - Setup trait data ####
# ========================#

  # B1 - HERBACEOUS LAYER ####
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
      save(H.abund,file=paste0(wrk.dir,'Herbaceous.layer.species.abundance.change.for.sp.with.traits.RData'))
      
      # B1.1 - Remove date effect on traits  ####
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
        
      # (4) Leaf.Area - No Date effect
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
      #NS
      
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
      
      # Supp.Mass.Frac  varies by species
      
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
      
      # B1.2 - Turn individual-level table into species-level table  ####
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
      
      # B1.3 - Add species mean traits (Seed size and mychorizae) ####
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
      
      # save(H.traits2.sp,file=paste0(wrk.dir,'Species-Level.Herbaceous.Layer.Traits.withMycFrac.Rdata'))
      
      # Add in Seed Size 
      
      load(paste0(wrk.dir,'Seed.Size.TRY.and.TOPIC.RData'))
      dim(Seed.Size)
      # 60 4
      str(Seed.Size)
      # 'data.frame':	60 obs. of  4 variables:
      #   $ SpeciesName      : Factor w/ 60 levels "Abies balsamea",..: 1 2 3 4 5 6 7 8 9 10 ...
      # $ Seed.Weight.TOPIC: num  7.59 23.23 69.44 4.19 1.02 ...
      # $ Seed.Weight.TRY  : num  7.37 20.36 NA 6.37 1.85 ...
      # $ Combined.DataSets: num  7.59 23.23 69.44 4.19 1.02 ...
      
      H.traits2.sp<-merge(H.traits2.sp,Seed.Size[,-c(1:3),drop=FALSE],
                  by="row.names",all.x=T)
      head(H.traits2.sp)
      H.traits2.sp$Row.names<-NULL
      rownames(H.traits2.sp)<-H.traits2.sp$Species
      dim(H.traits2.sp)
      #51 16
      
      colnames(H.traits2.sp)[colnames(H.traits2.sp)=='Combined.DataSets']<-'Seed.Size'
      
      
      save(H.traits2.sp,file=paste0(wrk.dir,'Species-Level.Herbaceous.Layer.Traits.withMycFrac.Rdata'))
      
    # B2 - CANOPY LAYER ####
    #=======================#
      
      # TO DO 
      
# C - Data Exploration (Following Highlands Stats course) ####
#=================================================================#

  #C1 - Outliers on Y and X  ----
  
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
        
      
      #Y3 - coverage
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
               main='Stor.Mass.Frac',
               cex=0.6) # 0 inflated - no outliers
      
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
     
  #C2 - Homogeneity (homoscedasticity) of Y ----
   
      # Test after model is built, as validation step
      # residuals vs fitted for each var.

  #C3 - Normality of Y ----
     
     # Test after model is built, as validation step
     # (least important assumption - use link function instead of transforming Y).
     # test for normality of residuals

  #C4 - Zero trouble (Y) ? ----
     
      sort(H.abund$abund.ratio)
      # only one 0 value - Okay! 
      
  #C5 - Collinearity  (Xs) ----
  
      
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
   
  # C5.1 Deal with colinearity using VIFs ####
  #-----------------------------------------#
   
  # Do VIF tests using dataset without myccorhizal fraction because it this trait contains only 42 species
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
   
  #1 TRY without LMA 
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
  
  #2 Then try without Fine root Diameter, because I don't want to loose SRL
  
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
  
  # C5.2 Deal with colinearity using PCs ####
  #-----------------------------------------#
  
  # Is it better to collapse the 13 traits into PCs?
  MyPCATraits<-c('Log.Ht.veg','Min.Root.Loca','Max.Root.Loca','Log.Lamina.thck','Log.LMA',
                 'LDMC','Log.Leaf.Area','Leaf.Mass.Frac','Supp.Mass.Frac','Rep.Mass.Frac',
                 'Stor.Mass.Frac','Log.F.Root.Diam','SRL')
  
  MyPCAData<-H.traits2.sp[,MyPCATraits][complete.cases(H.traits2.sp[,MyPCATraits]),]
  
  library(vegan)
  
  my.pca<-rda(MyPCAData,scale=T)
  biplot(my.pca,
         display = c("sites", 
                     "species"),
         type = c("text",
                  "text"))
  
  biplot(my.pca,
         choices=c(3,4),
         display = c("sites", 
                     "species"),
         type = c("text",
                  "text"))
  
  # Plot with vegan functions to test for significance of axes and variables
  
  # "scaling = 1" means that we are looking at correlation biplot
  #  angles between lines reflect their correlations
  
  ordi12<-ordiplot(my.pca, scaling=1, type="t",
                   main = "PCA 1-2")
  ordiequilibriumcircle(my.pca, ordi12) 
  
  ordi34<-ordiplot(my.pca, scaling=1, type="t",choices=c(3,4),
                   main = "PCA 3-4")
  ordiequilibriumcircle(my.pca, ordi34) 
  
  ordi56<-ordiplot(my.pca, scaling=1, type="t",choices=c(5,6),
                   main = "PCA 3-4")
  ordiequilibriumcircle(my.pca, ordi56) 
  
  summary(my.pca)
  # Call:
  # rda(X = MyPCAData, scale = T) 
  # 
  # Partitioning of correlations:
  #               Inertia Proportion
  # Total              13          1
  # Unconstrained      13          1
  # 
  # Eigenvalues, and their contribution to the correlations 
  # 
  # Importance of components:
  #                          PC1    PC2    PC3    PC4     PC5     PC6     PC7     PC8     PC9    PC10    PC11
  # Eigenvalue            3.3980 2.6571 1.9697 1.4526 1.29771 0.84580 0.54751 0.40279 0.23704 0.13444 0.04764
  # Proportion Explained  0.2614 0.2044 0.1515 0.1117 0.09982 0.06506 0.04212 0.03098 0.01823 0.01034 0.00366
  # Cumulative Proportion 0.2614 0.4658 0.6173 0.7290 0.82886 0.89392 0.93604 0.96702 0.98525 0.99560 0.99926
  # 
  #                           PC12     PC13
  # Eigenvalue            0.007333 0.002292
  # Proportion Explained  0.000560 0.000180
  # Cumulative Proportion 0.999820 1.000000
  # 
  # Scaling 2 for species and site scores
  # * Species are scaled proportional to eigenvalues
  # * Sites are unscaled: weighted dispersion equal on all dimensions
  # * General scaling constant of scores:  4.745172 
  # 
  # 
  # Species scores
  # 
  #                      PC1      PC2     PC3       PC4      PC5       PC6
  # Log.Ht.veg      -0.02985 -0.21835 -0.9195 -0.272806  0.49282  0.559430
  # Min.Root.Loca   -0.26507 -1.06552 -0.3410  0.533190 -0.27254 -0.043970
  # Max.Root.Loca   -0.27311 -1.10149 -0.3119  0.484425 -0.24021 -0.082546
  # Log.Lamina.thck  0.66651 -0.38192  0.8019  0.316441  0.18111 -0.423685
  # Log.LMA         -0.69319 -0.56407  0.7795  0.031087  0.43419  0.009899
  # LDMC            -1.05592 -0.10757  0.4676 -0.209673  0.32976  0.231548
  # Log.Leaf.Area    0.83523 -0.32294 -0.5881 -0.002711  0.10727 -0.190455
  # Leaf.Mass.Frac  -1.11506 -0.14794 -0.2758 -0.249071 -0.38860 -0.240577
  # Supp.Mass.Frac  -0.33950  0.58991 -0.4020  0.164217  0.76681 -0.609162
  # Rep.Mass.Frac   -0.15623  0.07729  0.0272  0.954380  0.61519  0.434712
  # Stor.Mass.Frac   1.12833 -0.13501  0.4361 -0.037680 -0.08079  0.397475
  # Log.F.Root.Diam  0.58770 -0.67438 -0.2335 -0.335773  0.57536 -0.302895
  # SRL              0.02943  0.88160 -0.2334  0.814536 -0.27616 -0.078143
  ---
  
  # which PCs are statistically significant?
    
  library(BiodiversityR) # for PCAsignificance() and ordiequilibriumcircle()
  PCAsignificance(my.pca)
  # First 5 PCs are significant
  
  # save those values to trait data frame
  PCA.scores<-my.pca$CA$u[,1:5]
  H.traits3.sp<-merge(H.traits2.sp,PCA.scores,by='row.names',all=T)
  head(H.traits3.sp)
  H.traits3.sp$Row.names<-NULL
  dim(H.traits3.sp)
  #51 21
  
  
  save(H.traits3.sp,file=paste0(wrk.dir,'Species-Level.Herbaceous.Layer.Traits.withMycFrac.PCs.RData'))
  
  # C6 Relationship between Y and X linear? ----
  
  Myxyplot(merge(H.abund,H.traits2.sp,by="row.names",all=T),
           Trait.Names,'abund.ratio',MyYlab="Abundance Ratio")
  
  # Everything linear - good to go. 
  # Two "outliers" - species with high abundance ratios. 
  
  # Try with ratio of log abundances as response variable instead.Remove IMCA (row 19 bc is problematic)
  Myxyplot(merge(H.abund,H.traits2.sp,by="row.names",all=T)[-24,],
          Trait.Names,'ratio.log.abund',MyYlab="Log-Abundance Ratio")

  # C7 Interactions? ----

  coplot(abund.ratio ~ Rep.Mass.Frac |SRL,
         data = merge(H.abund,H.traits2.sp,by="row.names",all=T), 
         xlab = "Rep.Mass.Frac",
         ylab = "Abundance Ratio")
  
  MyTraits<-c("Log.Ht.veg","Log.Leaf.Area",'Supp.Mass.Frac',"Rep.Mass.Frac","Log.Lamina.thck","LDMC",
              "Leaf.Mass.Frac","Min.Root.Loca","SRL",'Myc.Frac')
  
  
  for (i in 1:(length(MyTraits)-1)){
    for (j in (i+1):length(MyTraits)){
      
      x1<-MyTraits[i]
      x2<-MyTraits[j]
      
      pdf(file=paste0(data.dir,'Data Exploration/X-interaction plots/Interactions btw ',x1,' and ',x2,'.pdf'))
      
      data <- merge(H.abund,H.traits2.sp,by="row.names",all=T)
      
      coplot(data$abund.ratio ~ data[,x1]|data[,x2],
           xlab = c(x1,x2),
           ylab = "Abundance Ratio",
           number=4,overlap=0.3)
      dev.off()
      
       }
  }
  
  # Summary
  # 45 interactions
  # All graphs printed in C:/Users/Julie/Desktop/Postdoc/PROJECT - Megantic Traits/Data/Data Exploration/X-interaction plots/
  summary(lm(abund.ratio~Rep.Mass.Frac*Log.Lamina.thck,data=merge(H.abund,H.traits2.sp,by="row.names",all=T)))
  # Rep.Mass.Fraction x Log.Lamina.thck significant (for sp with thick lamina, abundance increases 
  # with increasing reproductive mass fraction)
  summary(lm(abund.ratio~Rep.Mass.Frac*Min.Root.Loca,data=merge(H.abund,H.traits2.sp,by="row.names",all=T)))
  # Rep.Mass.Fraction x Min.Root.Loca significant (for sp with shallow minimum root location, abundance increases 
  # with increasing reproductive mass fraction)
  summary(lm(abund.ratio~Supp.Mass.Frac*Min.Root.Loca,data=merge(H.abund,H.traits2.sp,by="row.names",all=T)))
  # Supp.Mass.Fraction x Min.Root.Loca significant
  
  # A1.0.8 Independance?
    # No spatial or temporal variables to test against the response variable
    # Look at independence of residuals after model fit. 
  



  
# Canopy layer ? ####
#=======================#
  
# To do
  

    
      

    
    
    
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
             ylab='Order of the data',
             cex=0.6,
             main="Storage Mass Fraction") # some species at 0, others not. 
    
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
    
    
    