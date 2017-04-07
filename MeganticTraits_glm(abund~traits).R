# B1 - GLMs ####

# B1.1 - Try individual variables
names(sp.H.traits)
#[1] "Ht.veg"        "Min.Root.Loca" "Max.Root.Loca" "Lamina.thck"   "LMA"           "LDMC"          "Leaf.Area"    
#[8] "myc.frac"

H.ht<-lm(H.abund$abund.change~Ht.veg,data=sp.H.traits)
summary(H.ht) #n.s.
H.min.root<-glm(H.abund$abund.change~Min.Root.Loca,data=sp.H.traits)
summary(H.min.root)#n.s.
H.max.root<-glm(H.abund$abund.change~Max.Root.Loca,data=sp.H.traits)
summary(H.max.root)#n.s.
H.lamina<-lm(H.abund$abund.change~Lamina.thck,data=sp.H.traits)
summary(H.lamina)#n.s.
H.LMA<-lm(H.abund$abund.change~LMA,data=sp.H.traits)
summary(H.LMA) #n.s.
H.LDMC<-lm(H.abund$abund.change~LDMC,data=sp.H.traits)
summary(H.LDMC)#n.s.
H.area<-lm(H.abund$abund.change~Leaf.Area,data=sp.H.traits)
summary(H.area)#n.s.
H.myc<-lm(H.abund$abund.change~myc.frac,data=sp.H.traits)
summary(H.myc) # AdjR2 = 0.23, p-val=0.001

#B1.2

#Stepwise selection on H layer

#full model, no interactions
H1<-glm(H.abund$abund.change~.,data=sp.H.traits,na.action=na.exclude)
summary(H1,correlation=T) # see ?summary.glm

# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     1.8132     0.4006   4.527 0.000109 ***
# Ht.veg          1.4288     1.0711   1.334 0.193343    
# Min.Root.Loca   0.7317     0.6596   1.109 0.277073    
# Lamina.thck     0.3589     0.6484   0.554 0.584456    
# LMA            -0.3568     0.8177  -0.436 0.666061    
# LDMC            0.3841     0.7161   0.536 0.596086    
# Leaf.Area      -0.6979     0.6448  -1.082 0.288654    
# myc.frac       -1.2558     0.3739  -3.359 0.002344 ** 

RsquareAdj(H1)
#$r.squared
#[1] 0.3826448

#$adj.r.squared
#[1] 0.1926894

vif(H1)

#Min Root Location and Max Root Location highly correlated. 
# LDMC correlated with something, but less-so LMA. 

# Remove species with NAs in myc.frac or Leaf.Area, to be able to run stepwise selection
sp.H.traits.c<-sp.H.traits[!is.na(sp.H.traits$myc.frac),]
sp.H.traits.c<-sp.H.traits.c[!is.na(sp.H.traits.c$Leaf.Area),]
dim(sp.H.traits.c)
# 35 8

# Adjust species list for abundance data to have matching datasets
H.abund.c<-H.abund[which(rownames(H.abund)%in%rownames(sp.H.traits.c)),]
dim(H.abund.c)
#35 2

#Stepwise selection
H1<-glm(H.abund.c$abund.change~.,data=sp.H.traits.c,na.action=na.exclude)
H0<-glm(H.abund.c$abund.change~1,data=sp.H.traits.c,na.action=na.exclude)

step.H<-step(H0,scope=list(lower=H0,upper=H1),direction='both',trace=T)
summary(step.H)

# glm(formula = H.abund.c$abund.change ~ myc.frac + Min.Root.Loca, 
#     data = sp.H.traits.c, na.action = na.exclude)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     1.5841     0.3220   4.919 2.52e-05 ***
#   myc.frac       -1.3174     0.3518  -3.745 0.000713 ***
#   Min.Root.Loca   0.8475     0.5657   1.498 0.143878 

# Null deviance: 167.00  on 34  degrees of freedom
# Residual deviance: 116.06  on 32  degrees of freedom
# AIC: 149.28

plot(step.H)

step.H.2<-glm(H.abund.c$abund.change ~ myc.frac + Min.Root.Loca+myc.frac*LDMC*Min.Root.Loca,
              data = sp.H.traits.c, na.action = na.exclude)

summary(step.H.2) #AIC=149.54; 3-way interaction least significant
step.H.3<-update(step.H.2,~.-myc.frac:Min.Root.Loca:LDMC)
summary(step.H.3) #AIC=149.59; Min.Root.Loca:LDMC interaction least significant
step.H.4<-update(step.H.3,~.-Min.Root.Loca:LDMC)
summary(step.H.4) # AIC= 145.84; LDMC is least significant

#                        Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              1.9254     0.3456   5.572 5.19e-06 ***
# myc.frac                -1.3522     0.3511  -3.851 0.000598 ***
# Min.Root.Loca            1.0585     0.5645   1.875 0.070858 .  
# LDMC                     0.3926     0.3664   1.072 0.292728    
# myc.frac:LDMC           -0.6422     0.3518  -1.826 0.078223 .  
# myc.frac:Min.Root.Loca  -1.3029     0.5772  -2.257 0.031696 *

# Can I remove LDMC? Interactin significant, but not term alone.
# Group non-different categories of Min.Root.Location (see crawley, p.433)
RsquareAdj(step.H.4) # 0.38
plot(step.H.4)

# What's up with CASC?
# Transform response variable ?
# test curvature +|(ht.veg^2)

#full model, interactions
H1x<-glm(H.abund$abund.change~(.-Min.Root.Loca)^2,data=sp.H.traits,na.action=na.exclude)
summary(H1x)

