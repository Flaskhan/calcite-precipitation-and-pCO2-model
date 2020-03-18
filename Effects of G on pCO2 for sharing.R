# Author: Hares Khan 2020, (hkhan.ch@gmail.com)
# This is the R code used in the study "Khan H., Marcé R., Laas A., Obrador B., "Pelagic calcite precipitation as a major component of the global carbon cycle in lakes and reservoirs"
# The code was used for assessing the effect of calcite precipitation on pCO2. The outputs of this model were then used with the global distribution of lakes per alkalinity level (data from Marcé et al., 2015) for upscaling at a global level the CO2 emitted by calcification.
# The model is based on the carbonate equilibra using the AquaEnv package.

# This code can be used freely and redistributed.
# This code comes with no warranty

### required packages
library(AquaEnv)
library(ggplot2)

### Inputs for Aquaenv that define the initial condition.
S      <- 0.1
t      <- 25
p      <- gauge_p(d=10)  
d      <- 1
pH     <- 8
TA     <- seq(from=0.00005, to=0.00995, by=0.0001) # (mol/kg-soln)
fCO2   <- 0.00041 # atmospheric equilibrium


# Aquaenv uses and returns units as follow:
# DIC and TA expressed in mol/kg-solution
# pCO2 expressed in atm
# We assume 1L water = 1kg


###### Calcification rates (G) are in mol/L/day. This was used to calculate daily effects of G on pCO2 at a local scale (Figure 1B in the study)
### Each equation for calculating G5, G50, and G95, were obtained from the quantile regressions of the relationship between alkalinity and calcification rates, using the 5%, 50% and 95% quantiles.
### The alkalinity and G values from where these equations are obtained are presented in Table 1 in the supplementary methods of this study.
G5<-(1.318182*TA*1000-2.985455)/1000000
G5[which(G5<0)]=0
G50<-(2.277778*TA*1000-3.512222)/1000000
G50[which(G50<0)]=0
G95<-(3.555556*TA*1000-3.843333)/1000000
G95[which(G95<0)]=0

### create dataframe with either the daily calcification rates, or the rates per summer period (choose between the two above) 
data=data.frame(TA,G5,G50,G95)

### Calculate the lower range (5% quantile) of the final values of alkalinity, dissolved inorganic carbon (DIC) and pCO2 using aquaenv function
ae <- aquaenv(S, t, p, TA=data$TA, fCO2=fCO2, SumCO2 = NULL)
DICini = data.frame(matrix(ae$SumCO2))
colnames(DICini) <- c(paste("DICini")) 
TAfin<-data$TA-(2*data$G5) ### 2 equivalences lost per mole of calcite that precipitates
TAfin<-data.frame(TAfin)
colnames(TAfin) <- c(paste("TAfin5")) 
DICfin<- DICini-data$G5
DICfin<-data.frame(DICfin)
colnames(DICfin) <- c(paste("DICfin5")) 
ae2 <- aquaenv(S, t, p, SumCO2=DICfin$DICfin5, TA=TAfin$TAfin5)
ae2 <- convert(ae2, "atm", "uatm", 1000000)
pCO2fin <- data.frame(ae2$fCO2)
colnames(pCO2fin) <- c(paste("pCO2fin5"))
CO2fin <- data.frame(ae2$CO2)
colnames(CO2fin) <- c(paste("CO2fin5"))

data<-cbind(data,TAfin, DICfin, pCO2fin, CO2fin)

### Calculate the mid range (50% quantile) of the final values of alkalinity, dissolved inorganic carbon (DIC) and pCO2 using aquaenv function
ae <- aquaenv(S, t, p, TA=data$TA, fCO2=fCO2, SumCO2 = NULL)
DICini = data.frame(matrix(ae$SumCO2))
colnames(DICini) <- c(paste("DICini")) 
TAfin<-data$TA-(2*data$G50) ### 2 equivalences lost per mole of calcite that precipitates
TAfin<-data.frame(TAfin)
colnames(TAfin) <- c(paste("TAfin50")) 
DICfin<- DICini-data$G50
DICfin<-data.frame(DICfin)
colnames(DICfin) <- c(paste("DICfin50")) 
ae2 <- aquaenv(S, t, p, SumCO2=DICfin$DICfin50, TA=TAfin$TAfin50)
ae2 <- convert(ae2, "atm", "uatm", 1000000)
pCO2fin <- data.frame(ae2$fCO2)
colnames(pCO2fin) <- c(paste("pCO2fin50"))
CO2fin <- data.frame(ae2$CO2)
colnames(CO2fin) <- c(paste("CO2fin50"))

data<-cbind(data,TAfin, DICfin, pCO2fin, CO2fin)

### Calculate the upper range (95% quantile) of the final values of alkalinity, dissolved inorganic carbon (DIC) and pCO2 using aquaenv function
ae <- aquaenv(S, t, p, TA=data$TA, fCO2=fCO2, SumCO2 = NULL)
DICini = data.frame(matrix(ae$SumCO2))
colnames(DICini) <- c(paste("DICini")) 
TAfin<-data$TA-(2*data$G95) ### 2 equivalences lost per mole of calcite that precipitates
TAfin<-data.frame(TAfin)
colnames(TAfin) <- c(paste("TAfin95")) 
DICfin<- DICini-data$G95
DICfin<-data.frame(DICfin)
colnames(DICfin) <- c(paste("DICfin95")) 
ae2 <- aquaenv(S, t, p, SumCO2=DICfin$DICfin95, TA=TAfin$TAfin95)
ae2 <- convert(ae2, "atm", "uatm", 1000000)
pCO2fin <- data.frame(ae2$fCO2)
colnames(pCO2fin) <- c(paste("pCO2fin95"))
CO2fin <- data.frame(ae2$CO2)
colnames(CO2fin) <- c(paste("CO2fin95"))

data<-cbind(data,TAfin, DICfin, pCO2fin, CO2fin)

### convert alkalinity (TA) from mol/l to meq/L
data$TA<-data$TA*1000

### plot the effect of calcification (G) on pCO2. The black line represents the 50% quantile from the regression between alkalinity and G, and the shaded area represents the range between 5% and 95% quantiles
ggplot(data=data, aes(x=TA, y=pCO2fin50)) + xlim(0,6)+geom_line(size=1.2) + 
  xlab("Alkalinity (meq/L)") + ylab("pCO2 (uatm)") + 
  ggtitle(expression(paste("Effects of calcification on pCO"[2])))+
  theme_light()+ theme(text = element_text(size = 13))+
  geom_ribbon(aes(x=TA, ymin=pCO2fin5, ymax=pCO2fin95), alpha=0.3)

### Create a csv output data file.
write.csv(data, "upscaling CO2400_perday.csv")
