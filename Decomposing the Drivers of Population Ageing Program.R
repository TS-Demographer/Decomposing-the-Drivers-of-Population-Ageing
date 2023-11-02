# Decomposing the Drivers of Population Ageing

library(HMDHFDplus)
# Please enter your Human Mortality Database username and password below. No other changes to the program are required.
username.hmd = ""
password.hmd = ""
# Please note that this program will take a few minutes to run.

# This sets up the data frame for the final result
Result<-data.frame(matrix(0,nrow=11, ncol=13))
colnames(Result)=c("Country","OADR_t","OADR_th","OADR change", "Birth component", "Survivorship component", "Migration component","Old-age birth component", "Old-age survivorship component", "Old-age migration component","Working-age birth component", "Working-age survivorship component", "Working-age migration component")

Names<-c("NOR", "CHE", "GBRTENW", "GBR_SCO", "NLD", "SWE", "FIN", "DNK", "FRATNP", "ESP", "ITA")
# This runs the program for the following countries:
# NOR=Norway, CHE=Switzerland, GBRTENW=England & Wales, GBR_SCO=Scotland, NLD=Netherlands, SWE=Sweden, FIN=Finland, DNK=Denmark, FRATNP=France, ESP=Spain, ITA=Italy, 
# AUS=Australia and US=USA are computed separately

# This runs the code 
for (ii in 1:11){
## Computing... ####
  ## Reading and formatting data #### 
# First load data on population
popdata = readHMDweb(CNTRY = Names[ii], item = "Population", username =
                         username.hmd, password = password.hmd)
popdata<-popdata[,c(1,2,4,5,6)]
colnames(popdata)=c("Year", "Age", "Female", "Male", "Total")
  
# Now load births data
birthsdata = readHMDweb(CNTRY = Names[ii], item = "Births", username =
                       username.hmd, password = password.hmd)

# Load deaths data
deathsdata = readHMDweb(CNTRY = Names[ii], item = "Mx_1x1", username =
                          username.hmd, password = password.hmd)
deathsdata=deathsdata[,-6]

# This sets years t and t+h as 2010 and 2020, respectively
hvalue <- 10
timeth <- 2020
timet <- timeth-hvalue

# This restricts the data to the years we are concerned with
popdata2<-popdata[(as.numeric(as.character(popdata$Year))<(timeth+1))&(as.numeric(as.character(popdata$Year))>(timet-101)),]
popdata1<-na.omit(popdata2[as.numeric(as.character(popdata2$Age))<101,])
birthsdata1<-birthsdata[(birthsdata$Year<(timeth))&(birthsdata$Year>(timet-102)),]
deathsdata2<-deathsdata[(deathsdata$Year<(timeth))&(deathsdata$Year>(timet-102)),]
deathsdata1<-na.omit(deathsdata2[as.numeric(as.character(deathsdata2$Age))<101,])

## Population ####
# CHECK: PROBABLY DON'T NEED THESE AS WE'VE ALREADY EXCLUDED 110+
popdata1$Age <- as.character(popdata1$Age)
popdata1$Age[popdata1$Age=="110+"] <- "110"
popdata1$Age <- as.numeric(popdata1$Age)

deathsdata1$Age <- as.character(deathsdata1$Age)
deathsdata1$Age[deathsdata1$Age=="110+"] <- "110"
deathsdata1$Age <- as.numeric(deathsdata1$Age)

# Extrapolating population data at times t and t+h
popdatat=popdata1[popdata1$Year==c(timet),]
popdatath=popdata1[popdata1$Year==c(timeth),]

# Age starts at 0, so this extracts total population counts for ages 65-110+ and 15-64
P65_t <- sum(popdatat$Total[66:101])
P65_th <- sum(popdatath$Total[66:101])
P1564_t <- sum(popdatat$Total[16:65])
P1564_th <- sum(popdatath$Total[16:65])

# This function takes the values of a variable, v, 
# at two times, t and t+h, and returns the derivative over that interval
change<-function(vt,vth,h){
  r<-log(vth/vt)/h
  vth2<-vt*exp(r*h/2)
  vdot=r*vth2
  return(vdot)
}

# This function also takes the values of a variable, v, 
# at two times, t and t+h, but this time returns the midpoint,
# that is, the value of v at time t+(h/2)
midp<-function(vt,vth,h){
  r<-log(vth/vt)/h
  vth2<-vt*exp(r*h/2)
  vth2[is.na(vth2)]<-0
  return(vth2)
}

# Taking growth rates for ages 65-110+ and 15-64 over the interval from t to t+h
midP65<- midp(P65_t,P65_th,hvalue)
midP1564<- midp(P1564_t,P1564_th,hvalue)

# This can be used as a check that the final change in OADR formula is correct
# Taking midpoints for ages 65-110+ and 15-64 over the interval from t to t+h
r_1564=log(P1564_th/P1564_t)/hvalue
r_65=log(P65_th/P65_t)/hvalue
# Taking the OADR at times t and t+h
OADR_t=P65_t/P1564_t
OADR_th=P65_th/P1564_th

# These are the c values at all ages x at times t and t+h
cdatat<-popdatat$Total/sum(popdatat$Total)
cdatath<-popdatath$Total/sum(popdatath$Total)

# This takes the above c values to the midpoint
cmidpointdata<-midp(cdatat,cdatath,hvalue)

# This is the growth rate for all ages x from time t to t+h
rx<-log(popdatath$Total/popdatat$Total)/hvalue
rx[is.infinite(rx)]<-0
rx[is.nan(rx)]<-0

ctotal65<-sum(cmidpointdata[66:101])
ctotal1564<-sum(cmidpointdata[16:65])
cindv65 <- sum(cmidpointdata[66:101]*rx[66:101])
cindv1564 <- sum(cmidpointdata[16:65]*rx[16:65])
OADRdot2 <- (midP65/midP1564)*(((1/ctotal65)*cindv65)-((1/ctotal1564)*cindv1564))

## Births ####
n=length(birthsdata1$Total)
# This takes the log of the ratio of birth rates h years apart all divided by h. This is the change in births.
Rbirths=log(birthsdata1$Total[-(1:hvalue)]/birthsdata1$Total[-((n-(hvalue-1)):n)])/hvalue
# This makes a data frame of the change in birth rates by year
Birthrates<-cbind(birthsdata1$Year[-((n-(hvalue-1)):n)],Rbirths)
colnames(Birthrates) <-c("Year","Rbirths")
# This is a vector of the change in births from youngest to oldest cohort as required
Rbirthsback<-rev(Rbirths)

## Deaths ####
# This is a function that makes a matrix of cohort death rates where rows are age and columns are cohort birth year. The input to the function is the end year
cohortdeathrate<-function(endyear){ 
  endyear=endyear-1
  cohortmxth2 <-matrix(0,101,101)
  dimnames(cohortmxth2) <-list(0:100, (endyear-100):(endyear)) 
  for (t in 1:101){
    cohortmxth <-c()
    for (i in 0:(t-1)){
      new1<-deathsdata1[deathsdata1$Year==c(endyear+1-t+i),]
      new2<-as.numeric(as.character(new1[new1$Age==c(i),5]))
      new2[is.na(new2)]<-0
      cohortmxth<-c(cohortmxth,new2)
    }
    cohortmxth2[1:t,102-t]<-cohortmxth
  }
  return(cohortmxth2)
}
# This uses the above function to get the cohort deaths rates for times t and t+h
cdrtimeth<-cohortdeathrate(timeth)
cdrtimet<-cohortdeathrate(timet)
# This finds the change in survivorship rates (it becomes survivorship as we want improvements (i.e. declines) in death rates to be expressed as positive change)
# The relative derivative of the survival function becomes the difference of the accumulated mortality.
Survrates<-(colSums(cdrtimet)-colSums(cdrtimeth))/hvalue
# As with the births, this is a vector of the change in death rates from youngest to oldest cohort as required
Survratesback<-rev(Survrates)

## Migration (as residual) ####
# This takes the change in migration as the residual of change in population not explained by changes in births or deaths
Migrates <- rx-Rbirthsback - Survratesback
# This creates a table of combined population, birth, deaths and migration growth rates
combinedrates<-cbind(c(0:100), rx, Rbirthsback,Survratesback, Migrates)
colnames(combinedrates) <-c("Year","Rx","Birthrates (back)","Survrates (back)", "Migrationrates")

## Calculating components ####
# This takes the sum of weighted changes in birth, death and migration for the older age group and then for the working age group.
birthcomp65 <- sum(cmidpointdata[66:101]*Rbirthsback[66:101])
survcomp65 <- sum(cmidpointdata[66:101]*Survratesback[66:101])
migcomp65   <- sum(cmidpointdata[66:101]*Migrates[66:101])

birthcomp1564 <- sum(cmidpointdata[16:65]*Rbirthsback[16:65])
survcomp1564 <- sum(cmidpointdata[16:65]*Survratesback[16:65])
migcomp1564   <- sum(cmidpointdata[16:65]*Migrates[16:65])

# This produces the birth, death and migration components of the OADR (split into the two age groups)
oldbirthcomp <-   (midP65/midP1564)*((1/ctotal65)*birthcomp65)
workingbirthcomp <-   (midP65/midP1564)*(-(1/ctotal1564)*birthcomp1564)

oldsurvcomp <-   (midP65/midP1564)*((1/ctotal65)*survcomp65)
workingsurvcomp <-   (midP65/midP1564)*(-(1/ctotal1564)*survcomp1564)

oldmigcomp   <-   (midP65/midP1564)*((1/ctotal65)*migcomp65)
workingmigcomp   <-   (midP65/midP1564)*(-(1/ctotal1564)*migcomp1564)

# This produces the overall birth, death and migration components of the OADR
combinedbirthcomp <-   (midP65/midP1564)*((1/ctotal65)*birthcomp65-(1/ctotal1564)*birthcomp1564)
combinedsurvcomp <-   (midP65/midP1564)*((1/ctotal65)*survcomp65-(1/ctotal1564)*survcomp1564)
combinedmigcomp   <-   (midP65/midP1564)*((1/ctotal65)*migcomp65-(1/ctotal1564)*migcomp1564)

## Useful Figures ####
# This plots the proportion of the population at each age at the midpoint between the cohorts
plot(c(0,100),c(0,0.03),col=0,xlab="Age", main=paste("Population Proportions by Cohort Age in",Names[ii]), ylab="Proportion of population at age x")
lines(c(0:100),cmidpointdata, col="red")

# This is a plot of the growth rates (population, birth, death and migration). They are not weighted by population proportions. 
plot(c(0,100),c(-0.1,0.1),col=0,xlab="Age", ylab="Rate of Change of Cohort Components", main=paste(Names[ii]," (unweighted)"))
lines(c(0:100),combinedrates[,2], col="red")
lines(c(0:100),combinedrates[,3], col="green")
lines(c(0:100),combinedrates[,4], col="blue")
lines(c(0:100),combinedrates[,5], col="orange")
abline(h=0,lty=2)
legend(0,0.105,box.lwd=0,legend=c("Population growth rate","Birth growth rate","Survivorship growth rate","Migration growth rate"),col=c("red","green","blue","orange"),lty=1)

# This is a plot of the growth rates (population, birth, death and migration). This time they are weighted by population proportions. 
plot(c(0,100),c(-0.001,0.001),col=0,xlab="Age", ylab="Rate of Change of Cohort Components", main=paste(Names[ii]," (weighted)"))
lines(c(0:100),combinedrates[,2]*cmidpointdata, col="red")
lines(c(0:100),combinedrates[,3]*cmidpointdata, col="green")
lines(c(0:100),combinedrates[,4]*cmidpointdata, col="blue")
lines(c(0:100),combinedrates[,5]*cmidpointdata, col="orange")
abline(h=0,lty=2)
legend(0,0.00105,box.lwd=0,legend=c("Population growth rate","Birth growth rate","Survivorship growth rate","Migration growth rate"),col=c("red","green","blue","orange"),lty=1)

## Result for the table ####
# This gives OADR at time t, OADR at time t+h, change in OADR and the birth, survivorship and migration components of the change in OADR for Table 2
Result[ii,]<-cbind(Names[ii], OADR_t*100,OADR_th*100,OADRdot2*100,combinedbirthcomp*100,combinedsurvcomp*100,combinedmigcomp*100, oldbirthcomp*100, oldsurvcomp*100,oldmigcomp*100, workingbirthcomp*100, workingsurvcomp*100,workingmigcomp*100)
rm(list=setdiff(ls(), c("Result", "Names", "username.hmd", "password.hmd")))
}
View(Result)

