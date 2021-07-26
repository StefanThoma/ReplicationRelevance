############Many Labs 5#################
###Replication of Payne, Burkley, & Stokes###
###Charlie Ebersole, shortened to data preparation by Stefan Thoma###
###Begin 4-7-2017###

#setwd("rmd/Analysis/Data/Data_payne")
###Reading Data and initial formatting###
Data1<-read.csv(file="Genova_ML5.csv",header=TRUE,stringsAsFactors=FALSE)
Data1$subject<-as.integer(Data1$subject)
Data1$subject2<-Data1$subject
Data1$Site<-"Genova"
Data1$USvIT<-"IT"
length(unique(Data1$subject2))

Data2<-read.csv(file="Milan_ML5.csv",header=TRUE,stringsAsFactors=FALSE)
Data2$subject<-as.integer(Data2$subject)
Data2$subject2<-Data2$subject+200
Data2$Site<-"Milan"
Data2$USvIT<-"IT"
length(unique(Data2$subject2))

Data3<-read.csv(file="Rome_ML5.csv",header=TRUE,stringsAsFactors=FALSE)
Data3$subject<-as.integer(Data3$subject)
Data3$subject2<-Data3$subject+400
Data3$Site<-"Rome"
Data3$USvIT<-"IT"
length(unique(Data3$subject2))

Data4<-read.csv(file="UCF_ML5.csv",header=TRUE,stringsAsFactors=FALSE)
Data4$subject<-as.integer(Data4$subject)
Data4$subject2<-Data4$subject+600
Data4$Site<-"UCF"
Data4$USvIT<-"US"
length(unique(Data4$subject2))

Data5<-read.csv(file="UVA_ML5.csv",header=TRUE,stringsAsFactors=FALSE)
Data5$subject<-as.integer(Data5$subject)
Data5$subject2<-Data5$subject+800
Data5$Site<-"UVA"
Data5$USvIT<-"US"
#Dropping 34 and 162 because they could read the Chinese characters
Data5<-subset(Data5,subject2!=834)
Data5<-subset(Data5,subject2!=962)
length(unique(Data5$subject2))

Data6<-read.csv(file="Padova_ML5.csv",header=TRUE,stringsAsFactors=FALSE)
Data6$subject<-as.integer(Data6$subject)
Data6$subject2<-Data6$subject+1000
Data6$Site<-"Padova"
Data6$USvIT<-"IT"
length(unique(Data6$subject2))

Data7<-read.csv(file="VCU_ML5.csv",header=TRUE,stringsAsFactors=FALSE)
Data7$subject<-as.integer(Data7$subject)
Data7$subject2<-Data7$subject+1200
Data7$Site<-"VCU"
Data7$USvIT<-"US"
length(unique(Data7$subject2))

Data8<-read.csv(file="PLU_ML5.csv",header=TRUE,stringsAsFactors=FALSE)
Data8$subject<-as.integer(Data8$subject)
Data8$subject2<-Data8$subject+1600
Data8$Site<-"PLU"
Data8$USvIT<-"US"
length(unique(Data8$subject2))

Data<-rbind(Data1,Data2,Data3,Data4,Data5,Data6,Data7,Data8)
Data<-subset(Data,subject2<5000)
head(Data)
str(Data)
length(Data$subject2)
length(unique(Data$subject2))

###Recoding DV###
Data$response[Data$response=="minus2"]<-"1"
head(Data)
Data$response[Data$response=="minus1"]<-"2"
Data$response[Data$response=="plus1"]<-"3"
Data$response[Data$response=="plus2"]<-"4"

###Creating DVs###

WhiteDirect<-subset(Data,trialcode=="whitedirect")
length(WhiteDirect$subject2)
length(unique(WhiteDirect$subject2))
head(WhiteDirect)
WhiteDirect$response<-as.integer(WhiteDirect$response)
WhiteDirect$subject2<-as.integer(WhiteDirect$subject2)

WhiteDirectM<-summaryBy(response~subject2,data=WhiteDirect,FUN=mean,na.rm=TRUE)
colnames(WhiteDirectM)[2]<-"WhiteDirectMean"
WhiteDirectM

BlackDirect<-subset(Data,trialcode=="blackdirect")
length(BlackDirect$subject2)
length(unique(BlackDirect$subject2))
head(BlackDirect)
BlackDirect$response<-as.integer(BlackDirect$response)
BlackDirect$subject<-as.integer(BlackDirect$subject2)

BlackDirectM<-summaryBy(response~subject2,data=BlackDirect,FUN=mean,na.rm=TRUE)
colnames(BlackDirectM)[2]<-"BlackDirectMean"
BlackDirectM

WhiteIndirect<-subset(Data,trialcode=="whiteindirect")
length(WhiteIndirect$subject2)
length(unique(WhiteIndirect$subject2))
head(WhiteIndirect)
WhiteIndirect$response<-as.integer(WhiteIndirect$response)
WhiteIndirect$subject<-as.integer(WhiteIndirect$subject2)

WhiteIndirectM<-summaryBy(response~subject2,data=WhiteIndirect,FUN=mean,na.rm=TRUE)
colnames(WhiteIndirectM)[2]<-"WhiteIndirectMean"
WhiteIndirectM

BlackIndirect<-subset(Data,trialcode=="blackindirect")
length(BlackIndirect$subject2)
length(unique(BlackIndirect$subject2))
head(BlackIndirect)
BlackIndirect$response<-as.integer(BlackIndirect$response)
BlackIndirect$subject<-as.integer(BlackIndirect$subject2)

BlackIndirectM<-summaryBy(response~subject2,data=BlackIndirect,FUN=mean,na.rm=TRUE)
colnames(BlackIndirectM)[2]<-"BlackIndirectMean"
BlackIndirectM

###Merging together###

Clean1<-merge(WhiteDirectM, BlackDirectM,by="subject2")
Clean1
Clean2<-merge(Clean1, WhiteIndirectM,by="subject2")
Clean<-merge(Clean2, BlackIndirectM,by="subject2")
head(Clean)
length(Clean$subject2)

Clean$Indirect<-Clean$WhiteIndirectMean-Clean$BlackIndirectMean
Clean$Direct<-Clean$WhiteDirectMean-Clean$BlackDirectMean

mean(Clean$Indirect)
sd(Clean$Indirect)
mean(Clean$Direct)
sd(Clean$Direct)
t.test(Clean$Indirect,mu=0)
t.test(Clean$Direct,mu=0)
cor(Clean$Indirect, Clean$Direct)

###Ratio Scoring###

Clean$IndirectR<-Clean$WhiteIndirectMean/Clean$BlackIndirectMean
Clean$DirectR<-Clean$WhiteDirectMean/Clean$BlackDirectMean
mean(Clean$IndirectR)
sd(Clean$IndirectR)
hist(Clean$IndirectR)
mean(Clean$DirectR)
sd(Clean$DirectR)
hist(Clean$DirectR)

LowSP<-subset(Clean,subject2 %% 2 != 0)
LowSP$subject2
LowSP$Condition<-"-1"

HighSP<-subset(Clean,subject2 %% 2 == 0)
HighSP$subject2
HighSP$Condition<-"1"

DataFull<-rbind(LowSP,HighSP)
head(DataFull)
length(DataFull$subject2)

###Adding Demographics
str(Data)
Data$trialcode<-as.factor(Data$trialcode)
list(levels(Data$trialcode))
length(unique(Data$subject2))

Age<-subset(Data,trialcode=="age")
Age<-Age[,c("subject2","response")]
colnames(Age)<-c("subject2","age")

Anni<-subset(Data,trialcode=="anni")
Anni<-Anni[,c("subject2","response")]
colnames(Anni)<-c("subject2","age")

Ethnicity<-subset(Data,trialcode=="ethnicity")
Ethnicity <-Ethnicity[,c("subject2","response")]
colnames(Ethnicity)<-c("subject2","ethnicity")

Etnia<-subset(Data,trialcode=="etnia")
Etnia <-Etnia[,c("subject2","response")]
colnames(Etnia)<-c("subject2","ethnicity")

Gender<-subset(Data,trialcode=="gender")
Gender <-Gender[,c("subject2","response")]
colnames(Gender)<-c("subject2","gender")

Sesso<-subset(Data,trialcode=="sesso")
Sesso <-Sesso[,c("subject2","response")]
colnames(Sesso)<-c("subject2","gender")

USdem<-merge(Age,Ethnicity,by="subject2")
USdem<-merge(USdem,Gender,by="subject2")
ITdem<-merge(Anni,Etnia,by="subject2")
ITdem<-merge(ITdem,Sesso,by="subject2")
Dem<-rbind(USdem,ITdem)
head(Dem)
length(Dem$subject2)
length(unique(Dem$subject2))
Dem2<-Dem[!duplicated(Dem[c('subject2')]),]
head(Dem2)
length(Dem2$subject2)

DataFull<-merge(DataFull,Dem2,by="subject2",all.x=TRUE)
head(DataFull)
length(DataFull$subject2)
length(unique(DataFull$subject2))

###Adding site variables

Genova<-subset(DataFull,subject2<200)
list(Genova$subject2)
Genova$Site<-"Genova"
Genova$USvIT<-"1"

Milan<-subset(DataFull,subject2>200 & subject2<400)
list(Milan$subject2)
Milan$Site<-"Milan"
Milan$USvIT<-"1"

Rome<-subset(DataFull,subject2>400 & subject2<600)
list(Rome$subject2)
Rome$Site<-"Rome"
Rome$USvIT<-"1"

UCF<-subset(DataFull,subject2>600 & subject2<800)
list(UCF$subject2)
UCF$Site<-"UCF"
UCF$USvIT<-"-1"

UVA<-subset(DataFull,subject2>800 & subject2<1000)
list(UVA$subject2)
UVA$Site<-"UVA"
UVA$USvIT<-"-1"

Padova<-subset(DataFull,subject2>1000 & subject2<1200)
list(Padova$subject2)
Padova$Site<-"Padova"
Padova$USvIT<-"1"

VCU<-subset(DataFull,subject2>1200 & subject2<1600)
list(VCU$subject2)
VCU$Site<-"VCU"
VCU$USvIT<-"-1"

PLU<-subset(DataFull,subject2>1600)
list(PLU$subject2)
PLU$Site<-"PLU"
PLU$USvIT<-"-1"

Total<-rbind(Genova,Milan,Rome,UCF,UVA,Padova,VCU,PLU)
str(Total)
Total$Site<-as.factor(Total$Site)
Total$USvIT<-as.integer(Total$USvIT)
Total$Condition<-as.integer(Total$Condition)

###Sample Descriptives###
UStotal<-length(PLU$subject2)+length(VCU$subject2)+length(UVA$subject2)+length(UCF$subject2)
UStotal

ITtotal<-length(Genova$subject2)+length(Milan$subject2)+length(Rome$subject2)+length(Padova$subject2)
ITtotal

Total$age<-as.integer(Total$age)
Total$gender<-as.factor(Total$gender)
Total$ethnicity<-as.factor(Total$ethnicity)
Total$ConditionF<-as.factor(Total$Condition)
Total$USvITF<-as.factor(Total$USvIT)


write_csv(Total, "payne.csv")
rm(list=ls())



#US<-subset(Total,USvITF=="-1")
#IT<-subset(Total,USvITF=="1")
#
#mean(US$age,na.rm=TRUE)
#sd(US$age,na.rm=TRUE)
#summary(US$gender)
#343/(343+208+3)
#summary(US$ethnicity)
#65/(65+73+54+41+5+299+15+2)
#73/(65+73+54+41+5+299+15+2)
#54/(65+73+54+41+5+299+15+2)
#41/(65+73+54+41+5+299+15+2)
#5/(65+73+54+41+5+299+15+2)
#299/(65+73+54+41+5+299+15+2)
#15/(65+73+54+41+5+299+15+2)
#2/(65+73+54+41+5+299+15+2)
#
#mean(IT$age,na.rm=TRUE)
#sd(IT$age,na.rm=TRUE)
#summary(IT$gender)
#375/(385+170)
#summary(IT$ethnicity)
#(131+339)/(131+339+5+69+1)
#5/(131+339+5+69+1)
#1/(131+339+5+69+1)
#69/(131+339+5+69+1)
#
############Primary Analyses############
#
####Variable Guide###
##Direct: direct ratings
##Indirect: indirect ratings
##Condition: whether the participant was in the high-pressure (1) or low-pressure condition (-1)
##Site: denotes which data collection site the participant came from, in this simulation there are 10 sites
##USvIT: indicates whether the participant is from the United States (-1) or Italy (1); there are 4 collection sites in each country
#
##if you do not have the lme4 package, run the command: install.packages("lme4")
#require(lme4)
#require(lmerTest)
#require(MuMIn)
#
##Step 1 Model 
#
#Model1<-lmer(Direct~Indirect+Condition+(Condition|Site),data=Total)
#summary(Model1)
#
##Testing primary replication effect
#
#Model2<-lmer(Direct~Indirect*Condition+(Condition|Site),data=Total)
#summary(Model2)
#anova(Model1,Model2)
#r.squaredLR(Model2,Model1)
#
##Testing moderation by sample source
#
##First creating comparison model with new two way interactions
#Model3.0<-lmer(Direct~Indirect*Condition+Indirect*USvIT+Condition*USvIT+(Condition|Site),data=Total)
#summary(Model3.0)
#str(Total$USvIT)
#Model3.1<-lmer(Direct~Indirect*Condition*USvIT+(Condition|Site),data=Total)
#summary(Model3.1)
#anova(Model3.0,Model3.1)
#r.squaredLR(Model3.1,Model3.0)
#
##For Visualization
#Model4<-lmer(Direct~Indirect*ConditionF*USvITF+(Condition|Site),data=Total)
#summary(Model4)
#plot(allEffects(Model4))
#
#####Testing for the primary effect in each sample source####
#
###US##
#
#Model1<-lmer(Direct~Indirect+Condition+(Condition|Site),data=US)
#summary(Model1)
#
#Model2<-lmer(Direct~Indirect*Condition+(Condition|Site),data=US)
#summary(Model2)
#anova(Model1,Model2)
#r.squaredLR(Model2,Model1)
#plot(allEffects(Model2))
#
##For visualization
#
#Model3<-lmer(Direct~Indirect*ConditionF+(Condition|Site),data=US)
#plot(allEffects(Model3))
#
###Italy##
#
#Model1<-lmer(Direct~Indirect+Condition+(Condition|Site),data=IT)
#summary(Model1)
#
#Model2<-lmer(Direct~Indirect*Condition+(Condition|Site),data=IT)
#summary(Model2)
#anova(Model1,Model2)
#plot(allEffects(Model2))
#r.squaredLR(Model2,Model1)
#
##For visualization
#
#Model3<-lmer(Direct~Indirect*ConditionF+(Condition|Site),data=IT)
#plot(allEffects(Model3))
#
####Further Explorations within Country###
#require(compute.es)
#
#t.test(Direct~Condition,data=Total,var.equal=TRUE)
#t.test(Indirect~Condition,data=Total,var.equal=TRUE)
#
#LowP<-subset(US,Condition=="-1")
#HighP<-subset(US,Condition=="1")
#cor.test(LowP$Direct,LowP$Indirect,method="pearson",use="completeobs")
#cor.test(HighP$Direct,HighP$Indirect,method="pearson",use="completeobs")
#
#t.test(Direct~Condition,data=US,var.equal=TRUE)
#tes(0.40047, 286,272)
#t.test(Indirect~Condition,data=US,var.equal=TRUE)
#tes(1.4643, 286,272)
#
#LowP<-subset(IT,Condition=="-1")
#HighP<-subset(IT,Condition=="1")
#cor.test(LowP$Direct,LowP$Indirect,method="pearson",use="completeobs")
#cor.test(HighP$Direct,HighP$Indirect,method="pearson",use="completeobs")
#
#t.test(Direct~Condition,data=IT,var.equal=TRUE)
#tes(2.1582, 273,272)
#t.test(Indirect~Condition,data=IT,var.equal=TRUE)
#tes(0.91919, 273,272)
#
####Breaking apart the Direct and Indirect Scores###
#
#summaryBy(BlackDirectMean~Condition,data=US,FUN=c(mean,sd),na.rm=TRUE)
#summaryBy(BlackIndirectMean~Condition,data=US,FUN=c(mean,sd),na.rm=TRUE)
#summaryBy(WhiteDirectMean~Condition,data=US,FUN=c(mean,sd),na.rm=TRUE)
#summaryBy(WhiteIndirectMean~Condition,data=US,FUN=c(mean,sd),na.rm=TRUE)
#
#t.test(BlackDirectMean ~Condition,data=US,var.equal=TRUE)
#tes(-0.20778, 286,272)
#t.test(BlackIndirectMean ~Condition,data=US,var.equal=TRUE)
#tes(0.81665, 286,272)
#t.test(WhiteDirectMean ~Condition,data=US,var.equal=TRUE)
#tes(0.10473, 286,272)
#t.test(WhiteIndirectMean ~Condition,data=US,var.equal=TRUE)
#tes(1.7878, 286,272)
#
#t.test(BlackDirectMean~Condition,data=IT,var.equal=TRUE)
#tes(-1.334, 273,272)
#t.test(BlackIndirectMean~Condition,data=IT,var.equal=TRUE)
#tes(-0.14681, 273,272)
#t.test(WhiteDirectMean~Condition,data=IT,var.equal=TRUE)
#tes(0.8316, 273,272)
#t.test(WhiteIndirectMean~Condition,data=IT,var.equal=TRUE)
#tes(0.6384, 273,272)
#
##All data
#t.test(BlackDirectMean ~Condition,data=Total,var.equal=TRUE)
#t.test(BlackIndirectMean ~Condition,data=Total,var.equal=TRUE)
#t.test(WhiteDirectMean ~Condition,data=Total,var.equal=TRUE)
#t.test(WhiteIndirectMean ~Condition,data=Total,var.equal=TRUE)
#
####Break down by ethnicity###
#
#list(levels(Total$ethnicity))
#summary(Total$ethnicity)
#USsub<-subset(US,ethnicity== "White, non-Hispanic or Latino" )
#ModelE<-lm(Direct~Indirect*Condition,data=USsub)
#summary(ModelE)
#ModelEV<-lm(Direct~Indirect*ConditionF,data=USsub)
#plot(allEffects(ModelEV))
#
#LowP<-subset(USsub,Condition=="-1")
#HighP<-subset(USsub,Condition=="1")
#
#cor.test(USsub$Direct,USsub$Indirect,method="pearson",use="completeobs")
#cor.test(USsub$DirectR,USsub$IndirectR,method="pearson",use="completeobs")
#
#cor.test(LowP$Direct,LowP$Indirect,method="pearson",use="completeobs")
#cor.test(HighP$Direct,HighP$Indirect,method="pearson",use="completeobs")
#
#cor.test(LowP$DirectR,LowP$IndirectR,method="pearson",use="completeobs")
#cor.test(HighP$DirectR,HighP$IndirectR,method="pearson",use="completeobs")
#
####Simple effects based on high vs. low implicit scores
#
#length(USsub$subject2)
#mean(USsub$Indirect)
#sd(USsub$Indirect)
#
#-0.007466082 - 0.2446654
#-0.007466082 + 0.2446654
#
#Low<-subset(USsub,Indirect < -0.2521315)
#Low$HighLow<- -1
#
#High<-subset(USsub,Indirect > 0.2371993)
#High$HighLow<- 1
#
#SimpleEffectData<-rbind(Low,High)
#SimpleEffectData$HighLowF<-as.character(SimpleEffectData$HighLow)
#SimpleEffectData$HighLowF<-as.factor(SimpleEffectData$HighLowF)
#
#SimpleModel<-lm(Direct~Condition*HighLow,data=SimpleEffectData)
#summary(SimpleModel)
#plot(allEffects(SimpleModel))
#
#require(car)
#Anova(SimpleModel,type="III")
#
##Visualizing
#SimpleModel2<-lm(Direct~ConditionF*HighLowF,data=SimpleEffectData)
#summary(SimpleModel2)
#plot(allEffects(SimpleModel2))
#
#t.test(Direct~Condition,data=High,var.equal=TRUE)
#High$Condition
#tes(-0.93525, 22, 15)
#
#t.test(Direct~Condition,data=Low,var.equal=TRUE)
#Low$Condition
#tes(0.035491, 16, 13)
#
##Examing condition effects among all White participants
#t.test(Indirect~Condition,data=USsub,var.equal=TRUE)
#summary(USsub$ConditionF)
#tes(0.38579, 154, 145)
#
#t.test(Direct~Condition,data=USsub,var.equal=TRUE)
#tes(-0.28265, 154, 145)
#
#t.test(IndirectR~Condition,data=USsub,var.equal=TRUE)
#tes(0.41455, 154, 145)
#
#t.test(DirectR~Condition,data=USsub,var.equal=TRUE)
#tes(-0.60505, 154, 145)
#
#t.test(WhiteDirectMean~Condition,data=USsub,var.equal=TRUE)
#tes(-1.2702, 154, 145)
#
#t.test(BlackDirectMean~Condition,data=USsub,var.equal=TRUE)
#tes(-1.0274, 154, 145)
#
#t.test(WhiteIndirectMean~Condition,data=USsub,var.equal=TRUE)
#tes(-0.92208, 154, 145)
#
#t.test(BlackIndirectMean~Condition,data=USsub,var.equal=TRUE)
#tes(-1.0794, 154, 145)
#
####Breaking down further into component pieces of scores###
#cor.test(LowP$BlackDirectMean,LowP$BlackIndirectMean,method="pearson",use="completeobs")
#cor.test(LowP$WhiteDirectMean,LowP$WhiteIndirectMean,method="pearson",use="completeobs")
#
#cor.test(HighP$BlackDirectMean,HighP$BlackIndirectMean,method="pearson",use="completeobs")
#cor.test(HighP$WhiteDirectMean,HighP$WhiteIndirectMean,method="pearson",use="completeobs")
#
#cor.test(LowP$DirectR,LowP$IndirectR,method="pearson",use="completeobs")
#cor.test(HighP$DirectR,HighP$IndirectR,method="pearson",use="completeobs")
#
#t.test(Direct~Condition,data=USsub,var.equal=TRUE)
#t.test(Indirect~Condition,data=USsub,var.equal=TRUE)
#
#t.test(Direct~ethnicity,data=USsub,var.equal=TRUE)
#tes(-5.2265, 66, 299)
#t.test(Indirect~ethnicity,data=USsub,var.equal=TRUE)
#tes(-3.8747, 66, 299)
#
#t.test(DirectR~ethnicity,data=USsub,var.equal=TRUE)
#tes(-4.2038, 66, 299)
#t.test(IndirectR~ethnicity,data=USsub,var.equal=TRUE)
#tes(-3.5722, 66, 299)
#
#summaryBy(Direct~ethnicity,data=USsub,FUN=c(mean,sd))
#summaryBy(Indirect~ethnicity,data=USsub,FUN=c(mean,sd))
#summaryBy(DirectR~ethnicity,data=USsub,FUN=c(mean,sd))
#summaryBy(IndirectR~ethnicity,data=USsub,FUN=c(mean,sd))
#
##########Ratio Analysis###########
#str(Total)
#
#mean(Total$IndirectR,na.rm=TRUE)
#sd(Total$IndirectR,na.rm=TRUE)
#
#mean(Total$DirectR,na.rm=TRUE)
#sd(Total$DirectR,na.rm=TRUE)
#
#cor(Total$IndirectR, Total$DirectR)
#
#t.test(Total$IndirectR,mu=1)
#t.test(Total$DirectR,mu=1)
#
#
####Primary Analyses###
#
#Model1<-lmer(DirectR~IndirectR+Condition+(Condition|Site),data=Total)
#summary(Model1)
#
##Testing primary replication effect
#
#Model2<-lmer(DirectR~IndirectR*Condition+(Condition|Site),data=Total)
#summary(Model2)
#anova(Model1,Model2)
#r.squaredLR(Model2,Model1)
#
##Testing moderation by sample source
#
##First creating comparison model with new two way interactions
#Model3.0<-lmer(DirectR~IndirectR*Condition+IndirectR*USvIT+Condition*USvIT+(Condition|Site),data=Total)
#summary(Model3.0)
#
#Model3.1<-lmer(DirectR~IndirectR*Condition*USvIT+(Condition|Site),data=Total)
#summary(Model3.1)
#anova(Model3.0,Model3.1)
#r.squaredLR(Model3.1,Model3.0)
#
##For Visualization
#Model4<-lmer(DirectR~IndirectR*ConditionF*USvITF+(Condition|Site),data=Total)
#summary(Model4)
#plot(allEffects(Model4))
#
#####Testing for the primary effect in each sample source####
#
#US<-subset(Total,USvITF=="-1")
#IT<-subset(Total,USvITF=="1")
#
###US##
#
#Model1<-lmer(DirectR~IndirectR+Condition+(Condition|Site),data=US)
#summary(Model1)
#
#Model2<-lmer(DirectR~IndirectR*Condition+(Condition|Site),data=US)
#summary(Model2)
#anova(Model1,Model2)
#plot(allEffects(Model2))
#r.squaredLR(Model2,Model1)
#
##For Visualization
#
#Model3<-lmer(DirectR~IndirectR*ConditionF+(Condition|Site),data=US)
#plot(allEffects(Model3))
#
###IT##
#
#Model1<-lmer(DirectR~IndirectR+Condition+(Condition|Site),data=IT)
#summary(Model1)
#
#Model2<-lmer(DirectR~IndirectR*Condition+(Condition|Site),data=IT)
#summary(Model2)
#anova(Model1,Model2)
#r.squaredLR(Model2,Model1)
#plot(allEffects(Model2))
#
##For Visualization
#
#Model3<-lmer(DirectR~IndirectR*ConditionF+(Condition|Site),data=IT)
#plot(allEffects(Model3))
#
####Breaking down further###
#
#LowP<-subset(US,Condition=="-1")
#HighP<-subset(US,Condition=="1")
#cor.test(LowP$DirectR,LowP$IndirectR,method="pearson",use="completeobs")
#cor.test(HighP$DirectR,HighP$IndirectR,method="pearson",use="completeobs")
#cor.test(US$DirectR,US$IndirectR,method="pearson",use="completeobs")
#
#summary(US$ConditionF)
#t.test(DirectR~Condition,data=US,var.equal=TRUE)
#tes(-0.046602, 286,272)
#t.test(IndirectR~Condition,data=US,var.equal=TRUE)
#tes(1.1741, 286,272)
#
#LowP<-subset(IT,Condition=="-1")
#HighP<-subset(IT,Condition=="1")
#cor.test(LowP$DirectR,LowP$IndirectR,method="pearson",use="completeobs")
#cor.test(HighP$DirectR,HighP$IndirectR,method="pearson",use="completeobs")
#cor.test(IT$DirectR,IT$IndirectR,method="pearson",use="completeobs")
#
#summary(IT$ConditionF)
#t.test(DirectR~Condition,data=IT,var.equal=TRUE)
#tes(2.2057, 273,272)
#t.test(IndirectR~Condition,data=IT,var.equal=TRUE)
#tes(0.40982, 273,272)
#
####Break down by ethnicity###
#
#list(levels(Total$ethnicity))
#
####Merging two levels of Caucasico
#length(Total$subject2)
#summary(Total$ethnicity)
#
#Caucasico<-subset(Total,ethnicity!="Caucasico/Europeo")
#Europeo<-subset(Total,ethnicity=="Caucasico/Europeo")
#Europeo$ethnicity<-"Caucasico"
#Eth<-rbind(Caucasico, Europeo)
#length(Eth$subject2)
#summary(Eth$ethnicity)
#str(Eth)
#USE<-subset(Eth,USvITF=="-1")
#ITE<-subset(Eth,USvITF=="1")
#
#USsub<-subset(USE,ethnicity=="White, non-Hispanic or Latino")
#ModelE<-lm(Direct~Indirect*Condition,data=USsub)
#summary(ModelE)
#ModelEV<-lm(DirectR~IndirectR*ConditionF,data=USsub)
#plot(allEffects(ModelEV))
#
#LowP<-subset(USsub,Condition=="-1")
#HighP<-subset(USsub,Condition=="1")
#
#cor.test(LowP$DirectR,LowP$IndirectR,method="pearson",use="completeobs")
#cor.test(HighP$DirectR,HighP$IndirectR,method="pearson",use="completeobs")
#
#t.test(DirectR~Condition,data=USsub,var.equal=TRUE)
#t.test(IndirectR~Condition,data=USsub,var.equal=TRUE)
#
##Checking basic info within group
#
#mean(USsub$IndirectR,na.rm=TRUE)
#sd(USsub$IndirectR,na.rm=TRUE)
#
#mean(USsub$DirectR,na.rm=TRUE)
#sd(USsub$DirectR,na.rm=TRUE)
#
#cor(USsub$IndirectR, USsub$DirectR)
#
#t.test(USsub$IndirectR,mu=1)
#t.test(USsub$DirectR,mu=1)
#
####Descriptives###
#
#summary(Eth$ethnicity)
#summaryBy(Indirect~ethnicity,data=Eth,FUN=c(mean,sd),na.rm=TRUE)
#summaryBy(Direct~ethnicity,data=Eth,FUN=c(mean,sd),na.rm=TRUE)
#summaryBy(IndirectR~ethnicity,data=Eth,FUN=c(mean,sd),na.rm=TRUE)
#summaryBy(DirectR~ethnicity,data=Eth,FUN=c(mean,sd),na.rm=TRUE)
#
#
#####Site by Site####
#require(lsr)
#list(levels(Total$Site))
##Replace Site == "" with whichever site you'd like to see results from. Names of sites provided by the preceeding command.
#Site<-subset(Total,Site=="Genova")
#SiteModel<-lm(Direct~Indirect*Condition,data=Site)
#summary(SiteModel)
#etaSquared(SiteModel)
#plot(allEffects(SiteModel))
#sqrt(4.436101e-04)
#
##Checking direction of effect
#LowP<-subset(Site,Condition=="-1")
#HighP<-subset(Site,Condition=="1")
#
#cor.test(LowP$Direct,LowP$Indirect,method="pearson",use="completeobs")
#cor.test(HighP$Direct,HighP$Indirect,method="pearson",use="completeobs")
#
##Summarizing data by site
#summary(Total$Site)
#summaryBy(Indirect~Site,data=Total,FUN=c(mean,sd),na.rm=TRUE)
#summaryBy(Direct~Site,data=Total,FUN=c(mean,sd),na.rm=TRUE)
#summaryBy(IndirectR~Site,data=Total,FUN=c(mean,sd),na.rm=TRUE)
#summaryBy(DirectR~Site,data=Total,FUN=c(mean,sd),na.rm=TRUE)
#
#
#summaryBy(Indirect~Site*Condition,data=Total,FUN=c(mean,sd),na.rm=TRUE)
#summaryBy(Direct~Site*Condition,data=Total,FUN=c(mean,sd),na.rm=TRUE)
#summaryBy(IndirectR~Site*Condition,data=Total,FUN=c(mean,sd),na.rm=TRUE)
#summaryBy(DirectR~Site*Condition,data=Total,FUN=c(mean,sd),na.rm=TRUE)
#
#
####Plotting Effects###
#require(metafor)
#
#R<-c(-.021, -.075, .134, .141, -.268, -.363, .063, -.076)
#N<-c(135, 136, 139, 135, 146, 106, 160, 146)
#SiteName<-c("University of Genova", "University of Milan-Bicocca", "University of Padova", "University of Rome Sapienza", "Pacific Lutheran University", "University of Central Florida", "University of Virginia", "Virginia Commonwealth University")
#Country<-c("Italy","Italy","Italy","Italy","United States","United States", "United States", "United States")
#
#GraphData<-as.data.frame(cbind(R,N,SiteName,Country),stringsAsFactors=FALSE)
#GraphData
#str(GraphData)
#GraphData$R<-as.numeric(GraphData$R)
#GraphData$N<-as.numeric(GraphData$N)
#
#MetaData<-escalc(measure="COR",ri= GraphData$R,ni= GraphData$N)
#MetaData
#MetaData<-as.data.frame(cbind(MetaData$yi,MetaData$vi,SiteName,Country),stringsAsFactors=FALSE)
#str(MetaData)
#MetaData$V1<-as.numeric(MetaData$V1)
#MetaData$V2<-as.numeric(MetaData$V2)
#MetaData$Country<-as.factor(MetaData$Country)
#
#Model<-rma(V1,V2,data=MetaData)
#summary(Model)
#forest(Model,slab=MetaData$SiteName)
#
#Model2<-rma(V1,V2,mods=~MetaData$Country,data=MetaData)
#summary(Model2)
#forest.rma(Model2,xlab="Replication Effect Size (partial correlation)",slab=MetaData$SiteName)
#
####Simulated Version###
#
#R<-c(.19, .14, .16, .09, .25, .31, .35, .22)
#N<-c(135, 135, 135, 135, 135, 135, 135, 135)
#SiteName<-c("University of Genova", "University of Milan-Bicocca", "University of Padova", "University of Rome Sapienza", "Pacific Lutheran University", "University of Central Florida", "University of Virginia", "Virginia Commonwealth University")
#Country<-c("Italy","Italy","Italy","Italy","United States","United States", "United States", "United States")
#
#GraphData<-as.data.frame(cbind(R,N,SiteName,Country),stringsAsFactors=FALSE)
#GraphData
#str(GraphData)
#GraphData$R<-as.numeric(GraphData$R)
#GraphData$N<-as.numeric(GraphData$N)
#
#MetaData<-escalc(measure="COR",ri= GraphData$R,ni= GraphData$N)
#MetaData
#MetaData<-as.data.frame(cbind(MetaData$yi,MetaData$vi,SiteName,Country),stringsAsFactors=FALSE)
#str(MetaData)
#MetaData$V1<-as.numeric(MetaData$V1)
#MetaData$V2<-as.numeric(MetaData$V2)
#MetaData$Country<-as.factor(MetaData$Country)
#
#Model<-rma(V1,V2,data=MetaData)
#summary(Model)
#forest(Model,slab=MetaData$SiteName)
#
#Model2<-rma(V1,V2,mods=~MetaData$Country,data=MetaData)
#summary(Model2)
#forest.rma(Model2,xlab="Indirect Ratings x Condition Interaction",slab=MetaData$SiteName)
#
#