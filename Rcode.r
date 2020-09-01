#Daily temperature variation lowers the lethal and sublethal impact of a pesticide pulse due to a higher degradation rate 
#Vienna Delnat, Jonathan Verborgt, Lizanne Janssens & Robby Stoks
#R code tested on 01/09/2020

######Packages######
#install packages
install.packages("lme4")     
install.packages("afex")     
install.packages("car")     
install.packages("effects")     
install.packages("emmeans")     
install.packages("drc")     
install.packages("doBy")     
install.packages("dplyr")     
install.packages("lmerTest")     

#load packages
library(lme4)
library(afex)
library(car)
library(effects)
library(emmeans)
library(drc)
library(doBy)
library(dplyr)
library(lmerTest) 
#lmerTest needed to test effect of random factors with ranova; 
#problem: gives errors when running rbind for contrast tests; 
#restart needed and no reloading of this package when running contrast analysis

######Import datasets######

###Chlorpyrifos

#Chlorpyrifos concentration after 24 hours
dataCPF=read.csv("./Delnat-et-al_ChlorpyrifosConcentration.csv", sep=",", na.strings=c(""))
dataCPF$DTV=as.factor(dataCPF$DTV)
dataCPF$Competition=factor(dataCPF$Competition, levels=c("daphnia","combi","mug"))
dataCPF24=subset(dataCPF, Hour=="24")
str(dataCPF24)

###Culex pipiens data

#Range finder Chlorpyrifos
dataRF=read.csv("./Delnat-et-al_CulexRangefinder_190111.csv", sep=",", na.strings=c(""))
str(dataRF)

#Lethal effect - Survival
dataLethal=read.csv("./Delnat-et-al_CulexTotalSurvivalBinomial.csv", sep=",", na.strings=c(""))
str(dataLethal)

#Sublethal effect - Development time, Pupae mass
dataSublethal=read.csv("./Delnat-et-al_CulexSublethal.csv", sep=",", na.strings=c(""))
str(dataSublethal)

#Pre-pesticide exposure - Appendix H
dataPrePest=read.csv("./Delnat-et-al_CulexPrePesticideExposure.csv", sep=",", na.strings=c(""))
str(dataPrePest)

#Survival after 24h, 48h and 72h
dataSurvival=read.csv("./Delnat-et-al_CulexSurvival-24h-48h-72h.csv", sep=",", na.strings=c(""))
str(dataSurvival)


###Daphnia magna data

#Intrinisc population growth rate calculations using the Euler Lotka equation
dataPopGrowth=read.csv("./Delnat-et-al_DaphniaEulerLotka.csv", sep=",", na.strings=c(""))
dataPopGrowth = subset(dataPopGrowth, Day != "NA") #exclude NA 
dataPopGrowth$Adults = as.numeric(as.character(dataPopGrowth$Adults))
dataPopGrowth$Juveniles = as.numeric(as.character(dataPopGrowth$Juveniles))
dataPopGrowth$Day = as.numeric(as.character(dataPopGrowth$Day))
str(dataPopGrowth)

#Intrinsic population growth rate (lambda), Adult survival, Total sizes of first and second broods per mother, Times till first and second brood release
dataDaphnia=read.csv("./Delnat-et-al_Daphnia.csv", sep=",", na.strings=c(""))
str(dataDaphnia)
dataBrood1=subset(dataDaphnia,Death1 == "no")
dataBrood2=subset(dataDaphnia,Death2 == "no")


######Culex pipiens - Chlorpyrifos range finder######
#drc::drm --> function drm from the package drc
#Concentration needs to be numeric (nominal concentrations are used and not measured concentrations!)
#Mortality in percentage

#Dose response curve based on the survival 72 hours after the pesticide exposure
EC72L4 <- drc::drm(Survival72 ~ Chlorpyrifos, data = dataRF, na.action = na.omit, fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "EC50")))
summary(EC72L4)

#lack-of-fit test
drc::modelFit(EC72L4) 
#F11,109 = 1.45, P = 0.16 --> not significant, OK
#the dose response curve is a good fit to the data 

#LC50 and steepness of dose-response curve
ED(EC72L4, c(10,50,90),interval='delta')
#<10% mortality up to 0.282 µg/L (95% CI [0.275; 0.289]) 
#>90% mortality starting at 0.309 µg/L (95% CI [0.302; 0.315])
#The LC50,72h value of chlorpyrifos was 0.295 µg/L (95% CI [0.292; 0.298]) 

#Figure 2 (exported as svg and dot area manually adapted in Figma to visualise the number of replicate vials)
EC72=plot(EC72L4, type="obs", log="", pch=16,cex=0.8)
EC72=plot(EC72L4, type="confidence", log="", add=TRUE)

#Table D.2 in Appendix D
EC24L4 <- drc::drm(Survival24 ~ Chlorpyrifos, data = dataRF, na.action = na.omit, fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "EC50")))
ED(EC24L4, c(10,50),interval='delta')
EC48L4 <- drc::drm(Survival48 ~ Chlorpyrifos, data = dataRF, na.action = na.omit, fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "EC50")))
ED(EC48L4, c(10,50),interval='delta')
EC72L4 <- drc::drm(Survival72 ~ Chlorpyrifos, data = dataRF, na.action = na.omit, fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "EC50")))
ED(EC72L4, c(10,50),interval='delta')
EC96L4 <- drc::drm(Survival96 ~ Chlorpyrifos, data = dataRF, na.action = na.omit, fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "EC50")))
ED(EC96L4, c(10,50),interval='delta')
EC120L4 <- drc::drm(Survival120 ~ Chlorpyrifos, data = dataRF, na.action = na.omit, fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "EC50")))
ED(EC120L4, c(10,50),interval='delta')
EC144L4 <- drc::drm(Survival144 ~ Chlorpyrifos, data = dataRF, na.action = na.omit, fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "EC50")))
ED(EC144L4, c(10,50),interval='delta')
EC168L4 <- drc::drm(Survival168 ~ Chlorpyrifos, data = dataRF, na.action = na.omit, fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "EC50")))
ED(EC168L4, c(10,50),interval='delta')

#Figure D.1 in Appendix D (exported as svg and lay-out done in Figma)
plotTime=plot(EC24L4, type="confidence", log="", pch=16,cex=0.8)
plotTime=plot(EC48L4, type="confidence", log="", add=TRUE)
plotTime=plot(EC72L4, type="confidence", log="", add=TRUE)
plotTime=plot(EC96L4, type="confidence", log="", add=TRUE)
plotTime=plot(EC120L4, type="confidence", log="", add=TRUE)
plotTime=plot(EC144L4, type="confidence", log="", add=TRUE)
plotTime=plot(EC168L4, type="confidence", log="", add=TRUE)


######Chlorpyrifos concentration after 24 hours######

set_sum_contrasts()
lmCPFint=lm(Concentration ~ DTV*Competition, data=dataCPF24, na.action=na.omit) 
Anova(lmCPFint, type="III") 
#interaction between DTV and interspecific competition was not significant (F2,12 = 0.30, P = 0.74)

#DTV × Competition interaction was removed from the model testing for effects on the chlorpyrifos concentration after 24 h due to limited amount of replicates
set_treatment_contrasts()
lmCPFmain=lm(Concentration ~ DTV+Competition, data=dataCPF24, na.action=na.omit) 
Anova(lmCPFmain, type="II") 

AIC(lmCPFint, lmCPFmain)
#model without the interaction also had a lower AIC score
#lmCPFmain is used as a final model in the article

#emmeans and s.e. used to make Figure 3
emmeans(lmCPFmain, ~ DTV+Competition, type="response")

#assumpties ok?
shapiro.test(resid(lmCPFmain)) #OK                
hist(resid(lmCPFmain))    
qqnorm(resid(lmCPFmain))    
qqline(resid(lmCPFmain))     
leveneTest(Concentration ~ DTV, data = dataCPF24) #OK
leveneTest(Concentration ~ Competition, data = dataCPF24) #OK


######Culex pipiens - Total Survival######

#effect coding due to interaction effects (also use Anova type III instead of type II)
set_sum_contrasts() 

#assumption 
SurvivalBin=glm(TotalSurvival~CPF*DTV*Competition, data=dataLethal, na.action=na.omit, family=quasibinomial(link=logit))
summary(SurvivalBin) 
#Dispersion parameter for quasibinomial family taken to be 1.007197 --> OK

#Generalized linear mixed model with a binomial error structure and the logit link
#individuals as the unit of replication, yet took into account that animals from the same vial were not independent by adding vial as a random factor
SurvivalBin=glmer(TotalSurvival ~ CPF*DTV*Competition + (1|Vial), data=dataLethal, na.action=na.omit, family=binomial(link=logit),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5))) 
Anova(SurvivalBin, type="III") 

#contrasts with false discovery rate (fdr) correction for pairwise posthoc comparisons
interact1<-pairs(emmeans(SurvivalBin, ~DTV|CPF, adjust="none"))
interact2<-pairs(emmeans(SurvivalBin, ~CPF|DTV, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

#effect plots of full model and of significant interactions
plot(effect(mod=SurvivalBin, term="CPF*DTV*Competition"),type = "response")
plot(effect(mod=SurvivalBin, term="CPF*DTV"),type = "response")

#emmeans and s.e. used to make Figure 4A
emmeans(SurvivalBin, ~ DTV*CPF*Competition, type="response")


######Culex pipiens - Survival after 24 h, 48 h and 72 h - Appendix E######

#effect coding due to interaction effects (also use Anova type III instead of type II)
set_sum_contrasts() 


##survival after 24 hours

#general linear model
Survival24=lm(Survival24h~CPF*DTV*Competition, data=dataSurvival) 
Anova(Survival24, type="III")

#effect plots of full model 
plot(effect(mod=Survival24, term="CPF*DTV*Competition"))

#emmeans and s.e. used to make Figure E.1 in Appendix E
emmeans(Survival24, ~DTV*CPF*Competition, adjust="fdr")

#assumptions
shapiro.test(residuals(Survival24)) #not ok
hist(resid(Survival24)) 
#almost no mortality, transforming does not work, do not interpret the Anova
#survival after 24 hours is not analysed in Appendix E
leveneTest(Survival24h~DTV*CPF*Competition, data = dataSurvival) #ok


##survival after 48 hours

#general linear model
Survival48=lm(Survival48h~CPF*DTV*Competition, data=dataSurvival) 
Anova(Survival48, type="III", white.adjust=TRUE)

#emmeans and s.e. used to make Figure E.2 (A) in Appendix E
emmeans(Survival48, ~DTV*CPF*Competition, adjust="fdr")

#assumptions
shapiro.test(residuals(Survival48)) #not OK
hist(resid(Survival48)) 
leveneTest(Survival48h~DTV*CPF*Competition, data = dataSurvival) #not OK
aggregate(Survival48h~DTV*CPF*Competition, data = dataSurvival, var) #not OK

#transformation of the response variable survival using a Box Cox transformation
boxcox=lm(Survival48h ~ CPF*DTV*Competition, data=dataSurvival, na.action=na.omit) 
car::boxCox(boxcox, family="yjPower", plotit=TRUE, lambda = seq(9,10 , length = 10))
bcSurvival48h <- car::yjPower(dataSurvival$Survival48h, 9.6)
dataSurvival$bcSurvival48h=bcSurvival48h
str(dataSurvival)

#general linear model with box cox transformation of survival
Survival48TR=lm(bcSurvival48h~CPF*DTV*Competition, data=dataSurvival) 
Anova(Survival48TR, type="III", white.adjust=TRUE)

#contrasts with false discovery rate (fdr) correction for pairwise posthoc comparisons
interact1<-pairs(emmeans(Survival48TR, ~DTV|CPF, adjust="none"))
interact2<-pairs(emmeans(Survival48TR, ~CPF|DTV, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

#effect plots of full model 
plot(effect(mod=Survival48TR, term="CPF*DTV"))

#assumptions
shapiro.test(residuals(Survival48TR)) #(OK)
hist(resid(Survival48TR)) 
leveneTest(bcSurvival48h~DTV*CPF*Competition, data = dataSurvival) #not OK
aggregate(bcSurvival48h~DTV*CPF*Competition, data = dataSurvival, var) #not OK, use white.adjust in Anova


##Survival after 72 hours

#general linear model
Survival72=lm(Survival72h~CPF*DTV*Competition, data=dataSurvival) 
Anova(Survival72, type="III", white.adjust=TRUE)

#emmeans and s.e. used to make Figure E.2 (B) in Appendix E
emmeans(Survival72, ~DTV*CPF*Competition, adjust="fdr")

#assumptions
shapiro.test(residuals(Survival72)) #not OK
hist(resid(Survival72)) 
leveneTest(Survival72h~DTV*CPF*Competition, data = dataSurvival) #not OK
aggregate(Survival72h~DTV*CPF*Competition, data = dataSurvival, var) #not OK

#transformation of the response variable survival using a Box Cox transformation
boxcox=lm(Survival72h ~ CPF*DTV*Competition, data=dataSurvival, na.action=na.omit) 
car::boxCox(boxcox, family="yjPower", plotit=TRUE, lambda = seq(5,6 , length = 10))
bcSurvival72h <- car::yjPower(dataSurvival$Survival72h, 5.5)
dataSurvival$bcSurvival72h=bcSurvival72h
str(dataSurvival)

#general linear model with box cox transformation of survival
Survival72TR=lm(bcSurvival72h~CPF*DTV*Competition, data=dataSurvival) 
Anova(Survival72TR, type="III", white.adjust=TRUE)

#contrasts with false discovery rate (fdr) correction for pairwise posthoc comparisons
interact1<-pairs(emmeans(Survival72TR, ~DTV|CPF, adjust="none"))
interact2<-pairs(emmeans(Survival72TR, ~CPF|DTV, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

#effect plots of full model 
plot(effect(mod=Survival72TR, term="CPF*DTV"))

#assumptions
shapiro.test(residuals(Survival72TR)) #(OK)
hist(resid(Survival72TR)) 
leveneTest(bcSurvival72h~DTV*CPF*Competition, data = dataSurvival) #not OK
aggregate(bcSurvival72h~DTV*CPF*Competition, data = dataSurvival, var) #not OK, use white.adjust in Anova


######Culex pipiens - Development time######
#Development time from L4 till pupae 

#effect coding due to interaction effects (also use Anova type III instead of type II)
set_sum_contrasts() 

#General linear mixed models with a normal error structure and the identity link
#individuals as the unit of replication, yet took into account that animals from the same vial were not independent by adding vial as a random factor
fitDev=lmer(DevelopmentTimeL4~CPF*DTV*Competition + (1|Vial), data=dataSublethal)
Anova(fitDev, type="III")

#contrasts with false discovery rate (fdr) correction for pairwise posthoc comparisons
interact1<-pairs(emmeans(fitDev, ~DTV|CPF, adjust="none"))
interact2<-pairs(emmeans(fitDev, ~CPF|DTV, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")
interact1<-pairs(emmeans(fitDev, ~Competition|CPF, adjust="none"))
interact2<-pairs(emmeans(fitDev, ~CPF|Competition, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

#effect plots of full model and of significant interactions
plot(effect(mod=fitDev, term="CPF*DTV*Competition"))
plot(effect(mod=fitDev, term="CPF*DTV"))
plot(effect(mod=fitDev, term="CPF*Competition"))

#emmeans and s.e. used to make Figure 4B
emmeans(fitDev, ~DTV*CPF*Competition, adjust="fdr")

#assumption
shapiro.test(residuals(fitDev)) #OK
hist(resid(fitDev)) 
qqnorm(resid((fitDev)))    
qqline(resid((fitDev)))     
leveneTest(DevelopmentTimeL4~DTV*CPF*Competition, data = dataSublethal) #OK


######Culex pipiens - Pupae mass######

#effect coding due to interaction effects (also use Anova type III instead of type II)
set_sum_contrasts() 

#General linear mixed models with a normal error structure and the identity link
#corrected for development time by adding it as a covariate
#individuals as the unit of replication, yet took into account that animals from the same vial were not independent by adding vial as a random factor
fitMass=lmer(PupaeMass~CPF*DTV*Competition + DevelopmentTimeEgg + (1|Vial), data=dataSublethal)
Anova(fitMass, type="III")

#contrasts with false discovery rate (fdr) correction for pairwise posthoc comparisons
interact1<-pairs(emmeans(fitMass, ~DTV|CPF, adjust="none"))
interact2<-pairs(emmeans(fitMass, ~CPF|DTV, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")
interact1<-pairs(emmeans(fitMass, ~DTV|Competition, adjust="none"))
interact2<-pairs(emmeans(fitMass, ~Competition|DTV, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

#effect plots of full model and of significant interactions
plot(effect(mod=fitMass, term="CPF*DTV*Competition"))
plot(effect(mod=fitMass, term="CPF*DTV"))
plot(effect(mod=fitMass, term="DTV*Competition"))

#emmeans and s.e. used to make Figure 4C
emmeans(fitMass, ~CPF*DTV*Competition, adjust="fdr")

#assumption
shapiro.test(residuals(fitMass)) #OK
hist(resid(fitMass)) 
qqnorm(resid(fitMass))    
qqline(resid(fitMass))     
leveneTest(PupaeMass~CPF*DTV*Competition, data = dataSublethal) #not OK
aggregate(PupaeMass~CPF*DTV*Competition, data = dataSublethal, var) #however, no difference of a factor 5 between max (0.075) and min (0.030) variance --> OK


######Culex pipiens - Pre-pesticide exposure - Appendix H######

#effect coding due to interaction effects (also use Anova type III instead of type II)
set_sum_contrasts() 

#general linear model
#Egg clutches were collected on 5 March 2019 and on 15 March 2019 (EggDate) from the same laboratory culture
fitPrePest=lm(Larvae~DTV*EggDate, data=dataPrePest) 
Anova(fitPrePest, type="III")

#effect plots of full model and of significant interactions
plot(effect(mod=fitPrePest, term="DTV*EggDate"))

#emmeans and s.e. used to make Figure H.1 in Appendix H
emmeans(fitPrePest, ~DTV*EggDate, adjust="fdr")

#assumptions
shapiro.test(residuals(fitPrePest)) #OK
hist(resid(fitPrePest)) 
leveneTest(Larvae~DTV*EggDate, data = dataPrePest) #OK


######Daphnia magna - Euler Lotka calculations for intrinsic population growth rate (r) - Appendix B######

# calculate lx and mx are needed to calculate rm
daphnia_growth2 <- dataPopGrowth %>%
  group_by(DTV, Chlorpyrifos, Competition,Vial) %>%
  mutate(lx = Adults/max(Adults)) %>%
  mutate(mx = (Juveniles/Adults)) %>%
  ungroup()
str(daphnia_growth2)

#Remark: If all Adults (=Mothers) were dead, this was noted as '0' adults in the original dataset
#But when 'mx' is calculated this results in 'NaN' values
#These NaN values were removed because it gave a problem when calculating 'r'
daphnia_growth3 <- daphnia_growth2 %>% filter(!is.na(mx))
str(daphnia_growth3)

#make subsets per replicate per treatment to calculate r
daphnia_growth = subset(daphnia_growth3, DTV == "no" & Chlorpyrifos == "no" & Competition == "no") 
daphnia_growth = subset(daphnia_growth3, DTV == "yes" & Chlorpyrifos == "no" & Competition == "no") 
daphnia_growth = subset(daphnia_growth3, DTV == "no" & Chlorpyrifos == "yes" & Competition == "no" & Vial != "8") 
daphnia_growth = subset(daphnia_growth3, DTV == "no" & Chlorpyrifos == "no" & Competition == "yes") 
daphnia_growth = subset(daphnia_growth3, DTV == "yes" & Chlorpyrifos == "yes" & Competition == "no" & Vial != "2" & Vial != "11" & Vial != "13")
daphnia_growth = subset(daphnia_growth3, DTV == "no" & Chlorpyrifos == "yes" & Competition == "yes") 
daphnia_growth = subset(daphnia_growth3, DTV == "yes" & Chlorpyrifos == "no" & Competition == "yes") 
daphnia_growth = subset(daphnia_growth3, DTV == "yes" & Chlorpyrifos == "yes" & Competition == "yes") 
str(daphnia_growth)

#run the following code seperatly for each subset above
rList <- vector("numeric",max(daphnia_growth$Vial))
for(i in unique(daphnia_growth$Vial)){
  daphnia_growthFor = subset(daphnia_growth, Vial == i) 
  
  ##### Calculate rm (intrinsic growth rate) using Euler-Lotka equation: 
  x <- c(daphnia_growthFor$Day) 
  L <- c(daphnia_growthFor$lx) 
  m <- c(daphnia_growthFor$mx) 
  
  r.range<- c(-3, 3) # define lower and upper limit of r (e.g. r.range (0,5))
  eulerlotka <- function(r) sum(L * m * exp(-r*x)) - 1 # introduce the eulerlotka equation for solution (root) 
  eulerlotka
  res <- uniroot(f = eulerlotka, interval = r.range, extendInt = "yes", tol = 1e-8) # define uniroot and set acuracy of estimation (tol)
  rList[i] <- print(res$root)
}

rList
unique(daphnia_growth$Vial)
#collect data in a new excel file


######Daphnia magna - Lambda - Appendix B######
#lambda calculated as lambda = exp(r)
#This has as advantage that it is possible to calculate lambda 
#even when no adult females survive in a certain vial (lambda = 0), 
#while it is not possible to calculate r under that condition (r = -infinity) 

#effect coding due to interaction effects (also use Anova type III instead of type II)
set_sum_contrasts() 

#general linear mixed model including the start date of a vial as random factor
fitPopGrowth=lmer(Lambda~Chlorpyrifos*DTV*Competition + (1|StartDateDaphnia), data=dataDaphnia) 
Anova(fitPopGrowth, type="III", white.adjust = TRUE)
ranova(fitPopGrowth) #load lmerTest

#contrasts with false discovery rate (fdr) correction for pairwise posthoc comparisons
#Remark: make sure lmerTest for ranova is not loaded
interact1<-pairs(emmeans(fitPopGrowth, ~Competition|Chlorpyrifos, adjust="none"))
interact2<-pairs(emmeans(fitPopGrowth, ~Chlorpyrifos|Competition, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

#effect plots of full model and of significant interactions
plot(effect(mod=fitPopGrowth, term="Chlorpyrifos*DTV*Competition"))
plot(effect(mod=fitPopGrowth, term="Chlorpyrifos*Competition"))

#emmeans and s.e. used to make Figure B.2 in Appendix B
emmeans(fitPopGrowth, ~Chlorpyrifos*DTV*Competition, adjust="fdr")

#assumptions
shapiro.test(residuals(fitPopGrowth)) 
hist(resid(fitPopGrowth)) #(OK)
qqnorm(resid(fitPopGrowth))    
qqline(resid(fitPopGrowth))     
leveneTest(Lambda~Chlorpyrifos*DTV*Competition, data = dataDaphnia) #not OK
aggregate(Lambda~Chlorpyrifos*DTV*Competition, data = dataDaphnia, var) #not OK --> use white.adjust=TRUE in Anova


######Daphnia magna - Adult survival - Appendix B######

#effect coding due to interaction effects (also use Anova type III instead of type II)
set_sum_contrasts() 

#general linear mixed model including the start date of a vial as random factor
SurvivalDaphnia=lmer(Survival8DaysPercent ~ Chlorpyrifos*DTV*Competition + (1|StartDateDaphnia), data=dataDaphnia, na.action=na.omit)
Anova(SurvivalDaphnia, type="III",white.adjust=TRUE) 
ranova(SurvivalDaphnia) #load lmerTest

#contrasts with false discovery rate (fdr) correction for pairwise posthoc comparisons
#Remark: make sure lmerTest for ranova is not loaded
interact1<-pairs(emmeans(SurvivalDaphnia, ~Competition|Chlorpyrifos, adjust="none"))
interact2<-pairs(emmeans(SurvivalDaphnia, ~Chlorpyrifos|Competition, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

interact1<-pairs(emmeans(SurvivalDaphnia, ~DTV|Chlorpyrifos, adjust="none"))
interact2<-pairs(emmeans(SurvivalDaphnia, ~Chlorpyrifos|DTV, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

interact1<-pairs(emmeans(SurvivalDaphnia, ~DTV|Competition, adjust="none"))
interact2<-pairs(emmeans(SurvivalDaphnia, ~Competition|DTV, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

#effect plots of full model and of significant interactions
plot(effect(mod=SurvivalDaphnia, term="Chlorpyrifos*DTV*Competition"))
plot(effect(mod=SurvivalDaphnia, term="Chlorpyrifos*Competition"))
plot(effect(mod=SurvivalDaphnia, term="DTV*Competition"))
plot(effect(mod=SurvivalDaphnia, term="Chlorpyrifos*DTV"))

#emmeans and s.e. used to make Figure B.3 in Appendix B
emmeans(SurvivalDaphnia, ~Chlorpyrifos*DTV*Competition, adjust="fdr")

#assumptions
shapiro.test(residuals(SurvivalDaphnia)) #OK
hist(resid(SurvivalDaphnia))
qqnorm(resid(SurvivalDaphnia))    
qqline(resid(SurvivalDaphnia))     
leveneTest(Survival8DaysPercent~Chlorpyrifos*DTV*Competition, data = dataDaphnia) #not OK
aggregate(Survival8DaysPercent~Chlorpyrifos*DTV*Competition, data = dataDaphnia, var) #not OK --> use white.adjust=TRUE in Anova


######Daphnia magna - Total size of first brood per mother - Appendix B######

#effect coding due to interaction effects (also use Anova type III instead of type II)
set_sum_contrasts() 

#general linear mixed model including the start date of a vial as random factor
SizeFirstBrood=lmer(NumberFirstBroodPerMother ~ Chlorpyrifos*DTV*Competition + (1|StartDateDaphnia),data=dataBrood1)
Anova(SizeFirstBrood, type="III")
ranova(SizeFirstBrood) #load lmerTest

#effect plots of full model and of significant interactions
plot(effect(mod=SizeFirstBrood, term="Chlorpyrifos*DTV*Competition"))

#emmeans and s.e. used to make Figure B.4A in Appendix B
emmeans(SizeFirstBrood, ~Chlorpyrifos*DTV*Competition, adjust="fdr")

#assumptions
shapiro.test(residuals(SizeFirstBrood)) #(OK)
hist(resid(SizeFirstBrood))
qqnorm(resid(SizeFirstBrood))    
qqline(resid(SizeFirstBrood))     
leveneTest(NumberFirstBrood~Chlorpyrifos*DTV*Competition, data = dataBrood1) #OK


######Daphnia magna - Total size of second brood per mother - Appendix B######

#effect coding due to interaction effects (also use Anova type III instead of type II)
set_sum_contrasts() 

#general linear mixed model including the start date of a vial as random factor
SizeSecondBrood = lmer(NumberSecondBroodPerMother ~ Chlorpyrifos*DTV*Competition + (1|StartDateDaphnia),data=dataBrood2)
Anova(SizeSecondBrood, type="III")
ranova(SizeSecondBrood) #load lmerTest

#contrasts with false discovery rate (fdr) correction for pairwise posthoc comparisons
#Remark: make sure lmerTest for ranova is not loaded
interact1<-pairs(emmeans(SizeSecondBrood, ~Chlorpyrifos|DTV*Competition, adjust="none"))
interact2<-pairs(emmeans(SizeSecondBrood, ~DTV|Chlorpyrifos*Competition, adjust="none"))
interact3<-pairs(emmeans(SizeSecondBrood, ~Competition|Chlorpyrifos*DTV, adjust="none"))
test(rbind(interact1,interact2,interact3), adjust="fdr")

#effect plots of full model and of significant interactions
plot(effect(mod=SizeSecondBrood, term="Chlorpyrifos*DTV*Competition"))

#emmeans and s.e. used to make Figure B.4B in Appendix B
emmeans(SizeSecondBrood, ~Chlorpyrifos*DTV*Competition, adjust="fdr")

#assumptions
shapiro.test(residuals(SizeSecondBrood)) #OK
hist(resid(SizeSecondBrood))
qqnorm(resid(SizeSecondBrood))    
qqline(resid(SizeSecondBrood))     
leveneTest(NumberSecondBrood~Chlorpyrifos*DTV*Competition, data = dataBrood2) #OK


######Daphnia magna - Time till first brood release - Appendix B######

#effect coding due to interaction effects (also use Anova type III instead of type II)
set_sum_contrasts() 

#general linear mixed model including the start date of a vial as random factor
TimeFirstBrood=lmer(TimeTillFirstBrood ~ Chlorpyrifos*DTV*Competition + (1|StartDateDaphnia),data=dataBrood1)
Anova(TimeFirstBrood, type="III")
ranova(TimeFirstBrood) #load lmerTest

#effect plots of full model and of significant interactions
plot(effect(mod=TimeFirstBrood, term="Chlorpyrifos*DTV*Competition"))

#emmeans and s.e. used to make Figure B.5A in Appendix B
emmeans(TimeFirstBrood, ~Chlorpyrifos*DTV*Competition, adjust="fdr")

#assumptions
shapiro.test(residuals(TimeFirstBrood)) #OK
hist(resid(TimeFirstBrood))
qqnorm(resid(TimeFirstBrood))    
qqline(resid(TimeFirstBrood))     
leveneTest(TimeTillFirstBrood~Chlorpyrifos*DTV*Competition, data = dataBrood1) #OK


######Daphnia magna - Time till second brood release - Appendix B######

#effect coding due to interaction effects (also use Anova type III instead of type II)
set_sum_contrasts() 

#general linear mixed model including the start date of a vial as random factor
TimeSecondBrood=lmer(TimeTillSecondBrood ~ Chlorpyrifos*DTV*Competition + (1|StartDateDaphnia),data=dataBrood2)
Anova(TimeSecondBrood, type="III")
ranova(TimeSecondBrood) #load lmerTest

#effect plots of full model and of significant interactions
plot(effect(mod=TimeSecondBrood, term="Chlorpyrifos*DTV*Competition"))

#emmeans and s.e. used to make Figure B.5A in Appendix B
emmeans(TimeSecondBrood, ~Chlorpyrifos*DTV*Competition, adjust="fdr")

#assumptions
shapiro.test(residuals(TimeSecondBrood)) #OK
hist(resid(TimeSecondBrood))
qqnorm(resid(TimeSecondBrood))    
qqline(resid(TimeSecondBrood))     
leveneTest(TimeTillSecondBrood~Chlorpyrifos*DTV*Competition, data = dataBrood2) #OK


######Save Rdata######
save.image(file="Delnat et al_Impact of Pesticide under DTV_20200901.Rdata")
