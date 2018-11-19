# Load TRY Data ####
#----------------------------------#

seed.file<-read.delim(paste0(data.dir,'TRY.Seed.Size.txt'),header=T)
dim(seed.file)
# 81420    27

str(seed.file)

seed.file[1:5,1:5]
View(seed.file)[1:30,]

# Create seed weights data.frame ####
#----------------------------------#

# 1- Remove the many rows which are descriptors, not data. ####
#---
unique(seed.file$TraitID)
# NA  26  18 131
unique(seed.file$TraitName)
# 18 - Plant Height (m)
# 26 - Seed Dry Mass (mg)
# 131 - Seed number per plant

all.seed.wt<-seed.file[which(seed.file$TraitID==26),]
dim(all.seed.wt)
# 671  27


# 2- Many rows have min and max seed values. Delete those. ####
#---
View(all.seed.wt)

unique(all.seed.wt$DataName)
  # [1] Seed dry mass                  Seed mass original value: mean Seed mass min                 
  # [4] Seed mass max                  Seed mass original value: min  Seed mass original value: max    

# Remove 'Seed mass max' values
length(which(all.seed.wt$DataName=='Seed mass max')) 
# 9
dim(all.seed.wt[-which(all.seed.wt$DataName=='Seed mass max'),])
# 662  27
all.seed.wt<-all.seed.wt[-which(all.seed.wt$DataName=='Seed mass max'),]

# Remove 'Seed mass min' values
length(which(all.seed.wt$DataName=='Seed mass min')) 
#9
dim(all.seed.wt[-which(all.seed.wt$DataName=='Seed mass min'),])
# 653  27
all.seed.wt<-all.seed.wt[-which(all.seed.wt$DataName=='Seed mass min'),]


# Remove 'Seed mass original value: min' values
length(which(all.seed.wt$DataName=='Seed mass original value: min')) 
# 9
dim(all.seed.wt[-which(all.seed.wt$DataName=='Seed mass original value: min'),])
# 644  27
all.seed.wt<-all.seed.wt[-which(all.seed.wt$DataName=='Seed mass original value: min'),] 

# Remove 'Seed mass original value: max' values
length(which(all.seed.wt$DataName=='Seed mass original value: max')) 
# 9
dim(all.seed.wt[-which(all.seed.wt$DataName=='Seed mass original value: max'),])
# 635  27
all.seed.wt<-all.seed.wt[-which(all.seed.wt$DataName=='Seed mass original value: max'),] 

all.seed.wt<-droplevels(all.seed.wt)
unique(all.seed.wt$DataName)
View(all.seed.wt)

# 3- Sort dataset by Species name ####
#---
all.seed.wt<-all.seed.wt[order(all.seed.wt$AccSpeciesName),]

unique(all.seed.wt$UnitName)
# [1] mg
# Levels: mg

# 4- Take means across datasets per species####
#---
# First, take mean per species within each dataset - often multiple measurements reported per species
Unique.Sp.Dataset<-as.data.frame(aggregate(StdValue ~ AccSpeciesName + Dataset,
                              data=all.seed.wt,
                              FUN=mean))
colnames(Unique.Sp.Dataset)
dim(Unique.Sp.Dataset)
# 163 3

# Then take mean per species across datasets
Seed.size.TRY<-as.data.frame(aggregate(StdValue ~ AccSpeciesName,
                                   data=Unique.Sp.Dataset,
                                   FUN=mean))
dim(Seed.size.TRY)
# 45 2

#Seed.size.TRY$DataSet<-'TRY'

Seed.size.TRY
save(Seed.size.TRY,file=paste0(wrk.dir,'Seed.Size.TRY.RData'))
# Create plant height data.frame ####
#----------------------------------#

# 1- Remove the many rows which are descriptors, not data. ####
#---
unique(seed.file$TraitID)
# NA  26  18 131
unique(seed.file$TraitName)
# 18 - Plant Height (m)
# 26 - Seed Dry Mass (mg)
# 131 - Seed number per plant

all.Plant.Ht<-seed.file[which(seed.file$TraitID==18),]
dim(all.Plant.Ht)
# 5426   2

all.Plant.Ht<- droplevels(all.Plant.Ht) # cleans out memory of unused factors

# 2- Many rows have min and max seed values. Delete those. ####
#---
View(all.Plant.Ht)

unique(all.Plant.Ht$DataName)
# [1] Plant height vegetative                                 
# [2] Plant height reproductive / generative                  
# [3] Coefficient of variation of plant height                
# [4] Height at 20 Years                                      
# [5] Vegetative plant height                                 
# [6] Plant height (unspecified if vegetative or reproductive)
# [7] Maximum plant height                                    
# [8] Plant height observed                                   
# [9] Maximum height                                          
# [10] Maximum height min                                      
# [11] Maximum height max                                      
# [12] Maximum height extreme   

# Start with 5426 rows
# Remove [2] 'Plant height reproductive / generative' values
length(which(all.Plant.Ht$DataName=='Plant height reproductive / generative')) 
# 55
all.Plant.Ht<-all.Plant.Ht[-which(all.Plant.Ht$DataName=='Plant height reproductive / generative'),]
dim(all.Plant.Ht)
# 5371   27

# Remove [3]'Coefficient of variation of plant height' 
length(which(all.Plant.Ht$DataName=='Coefficient of variation of plant height')) 
# 1
all.Plant.Ht<-all.Plant.Ht[-which(all.Plant.Ht$DataName=='Coefficient of variation of plant height'),]
dim(all.Plant.Ht)
# 5370   27

# Remove 'Height at 20 Years'
length(which(all.Plant.Ht$DataName=='Height at 20 Years')) 
# 30
all.Plant.Ht<-all.Plant.Ht[-which(all.Plant.Ht$DataName=='Height at 20 Years'),]
dim(all.Plant.Ht)
# 5340   27

# Remove 'Plant height (unspecified if vegetative or reproductive)'
length(which(all.Plant.Ht$DataName=='Plant height (unspecified if vegetative or reproductive)')) 
# 77
all.Plant.Ht<-all.Plant.Ht[-which(all.Plant.Ht$DataName=='Plant height (unspecified if vegetative or reproductive)'),]
dim(all.Plant.Ht)
# 5263   27

unique(all.Plant.Ht$DataName)

max.Plant.Ht<-all.Plant.Ht[which(all.Plant.Ht$DataName== 'Maximum plant height'|
                                   all.Plant.Ht$DataName== 'Maximum height'|
                                   all.Plant.Ht$DataName== 'Maximum height max'),]

dim(max.Plant.Ht)
# 149 27
View(max.Plant.Ht)

# 3- Sort dataset by Species name ####
#---
max.Plant.Ht<-max.Plant.Ht[order(max.Plant.Ht$AccSpeciesName),]

unique(max.Plant.Ht$UnitName)
# [1] m
# Levels: g

# 4- Take means across datasets per species####
#---
# First, take mean per species within each dataset - often multiple measurements reported per species
Unique.Sp.Ht.Dataset<-as.data.frame(aggregate(StdValue ~ AccSpeciesName + Dataset,
                                           data=max.Plant.Ht,
                                           FUN=mean))
colnames(Unique.Sp.Ht.Dataset)
dim(Unique.Sp.Ht.Dataset)
# 50 3

# Then take mean per species across datasets
Max.Height.TRY<-as.data.frame(aggregate(StdValue ~ AccSpeciesName,
                                       data=Unique.Sp.Ht.Dataset,
                                       FUN=mean))
dim(Max.Height.TRY)
# 27 2

Max.Height.TRY
# only 7 species are understorey species :-(
save(Max.Height.TRY,file=paste0(data.dir,'Max.Height.TRY.RData'))

# Load TOPIC Data ####
#----------------------------------#
Topic.seed.size<-read.csv(paste0(data.dir,'Seed size/TOPIC/TOPIC_Seed.Size.csv'),
         header=T)

str(Topic.seed.size)

# Convert # of seeds per kg of TOPIC into mg per seed of TRY
Topic.seed.size$SdWeight<-1/Topic.seed.size$TraitValue*1000*1000
#Topic.seed.size$TraitValue<-as.numeric(Topic.seed.size$TraitValue)
head(Topic.seed.size$SdWeight)
# 1 - Sort dataset by Species name ####
#---
Topic.seed.size<-Topic.seed.size[order(Topic.seed.size$OfficialName),]
dim(Topic.seed.size)
# 55 12

# 2- Take means across datasets per species####
#---
# First, take mean per species within each dataset - often multiple measurements reported per species
Unique.Sp.seed<-as.data.frame(aggregate(SdWeight ~ OfficialName + Documenter,
                                           data=Topic.seed.size,
                                           FUN=mean))
colnames(Unique.Sp.seed)
dim(Unique.Sp.seed)
# 51 3

# Then take mean per species across datasets 
Seed.size.TOPIC<-as.data.frame(aggregate(SdWeight ~ OfficialName,
                                       data=Unique.Sp.seed,
                                       FUN=mean))
dim(Seed.size.TOPIC)
# 39 2

View(Seed.size.TOPIC)
save(Seed.size.TOPIC,file=paste0(wrk.dir,'Seed.Size.TOPIC.RData'))

# Merge Seed Datasets ####
#=========================#

# Match column names to allow merging
colnames(Seed.size.TOPIC)
# "OfficialName" "SdWeight"   
colnames(Seed.size.TRY)
# "AccSpeciesName" "StdValue"  

colnames(Seed.size.TOPIC)<-c('SpeciesName','Seed.Weight.TOPIC')
colnames(Seed.size.TRY)<-c('SpeciesName','Seed.Weight.TRY')

# Change name of factor 'Acer saccharum var.saccharum' in TOPIC
levels(Seed.size.TOPIC$SpeciesName)[levels(Seed.size.TOPIC$SpeciesName)=="Acer saccharum var. saccharum"] <- "Acer saccharum"

# merge
Seed.Size<-merge(Seed.size.TOPIC,Seed.size.TRY,by='SpeciesName',all=T)
dim(Seed.Size) # 59  3

Seed.Size$Combined.DataSets<-c(Seed.Size$Seed.Weight.TOPIC[1:39],
                          Seed.Size$Seed.Weight.TRY[40:59])

# Add Species Codes as Row.names
Seed.Size$SpeciesName
sp.codes<-c('ABBA','ACRU','ACSA','ARNU','BEAL','BEPA','CACO',
            'CLBO','COAL','COCA','COCO','FAGR','FRAM','FRNI',
            'GAAP','IMGL','LOMA','LOOB','LOTA','MACA','MEVI',
            'MIRE','POAN','PODE','POGR','POHE','POTE','PRSE',
            'PRVI','RUAL','RUCA','RUFL','RUID','SMTR','THOC',
            'TIAM','TRBO','TSCA','VIAL','ACPE','ACSP',
            'ARTR','ATFI','CAIN','CASC','CIAL','EPHE','GATE',
            'GATR','IMCA','OSVI','PIGL','POBA','PRPU','SAPU',
            'SMRA','STAM','THNO','VEVI')
row.names(Seed.Size)<-sp.codes

save(Seed.Size,file=paste0(wrk.dir,'Seed.Size.TRY.and.TOPIC.RData'))
