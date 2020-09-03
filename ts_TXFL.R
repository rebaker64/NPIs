############################
# This script makes Fig. 2c i.e. results for 1 year control in Florida and Texas
############################

source('~/Dropbox/RSV_Global/Scripts/predtsirControl.R', encoding = 'UTF-8')

# Florida
load(file="prelimRSVstatefit.RData")

sub <- outall[outall$state=="Florida",]
pop <- 21480000
births <-221542/52
controlWeekStart = 11
controlWeekLength =  52

# run model
times <- seq(1,55,1/52)
controlStart = 43*52 + controlWeekStart
controlEnd = 43*52 + controlWeekStart + controlWeekLength
pred <- predtsirControl(times = times, births = rep(births, length = length(times)), beta = sub$sea_beta, alpha = 0.97, 
                        S0 =floor(0.8*pop), I0 = floor(0.2*pop), nsim = 10, stochastic = F, 
                        controlStart = controlStart, controlEnd =controlEnd, betachange = 0.8)
subsetI <- pred$I$mean[times > 40]
subsetS <- pred$S$mean[times > 40]

# load florida surveillance data
require("readxl")
dat <- read_excel("RSV_data_US.xlsx")
flo <- dat[dat$state=="Florida",]
flo <- flo[order(flo$year, flo$week),]
dfcorrect <- data.frame(week = rep(seq(1,52), 5), year = rep(seq(2016,2020,1), each = 52))
dffin <- merge(dfcorrect,flo, by=c("week","year"), all.x=T)
dffin <- dffin[order(dffin$year, dffin$week),]
flots <- dffin$percent_specimen_positive
flocases <- (flots/(mean(flots, na.rm=T)))*mean(subsetI[1:52]/pop) # scale % pos to match simulated incidence
timeuse <- seq(2016,2100,1/52)[1:length(subsetI)]

# plot model result
pdf(paste0("FloridaSim.pdf"),width=6,height=4)
par(mar=c(3,3,1,3))
plot(timeuse,subsetI/pop,type="n",col="#F64740",lwd =2,ylim=c(0,0.0035),xlab="",ylab=
       "",bty = "n",xaxs="i",xlim=c(2016,2030))
polygon(c(timeuse[52*4+controlWeekStart],timeuse[52*4+controlWeekStart],timeuse[52*4 + controlWeekStart + controlWeekLength],
          timeuse[52*4 + controlWeekStart + controlWeekLength]), c(-1,1,1,-1), border = NA,col="#D1D3DF")
lines(timeuse,subsetI/pop,col="#F64740",lwd =2,ylim=c(0,0.01),xlab="",ylab=
        "")
points(timeuse[1:length(flocases)],flocases,cex=0.5,col="grey32",pch=18)
abline(v = seq(2016,2100,1),lty=3,col="gray")
par(new = TRUE)
plot(timeuse, subsetS/pop, col="navy",type="l",xlab="",ylab="", 
     axes = F,ylim=c(0.7*min(subsetS/pop),max(subsetS/pop)), lwd = 2, lty = 2)
axis(side = 4)
title(xlab="Year", line = 2)
title(ylab="I/N", line = 2,col.lab="#F64740")
mtext(side = 4, line = 2, 'S/N',col="navy")
text(2016,0.038,"Florida",pos = 4,cex = 1.2)
dev.off()



######## Texas

load(file="prelimRSVstatefitTEXAS.RData")
sub <- outall[outall$state=="Texas",]

pop <- 29000000
births = 378624/52
controlWeekStart = 11
controlWeekLength =  52

# run model
times <- seq(1,55,1/52)
controlStart = 43*52 + controlWeekStart
controlEnd = 43*52 + controlWeekStart + controlWeekLength
pred <- predtsirControl(times = times, births = rep(births, length = length(times)), beta = sub$sea_beta, alpha = 0.97, 
                        S0 =floor(0.8*pop), I0 = floor(0.2*pop), nsim = 10, stochastic = F, 
                        controlStart = controlStart, controlEnd =controlEnd, betachange = 0.8)

subsetI <- pred$I$mean[times > 40]
subsetS <- pred$S$mean[times > 40]

# load RSV surveillance data for Texas
require("readxl")
dat <- read_excel("RSV_data_US.xlsx")
flo <- dat[dat$state=="Texas",]
flo <- flo[order(flo$year, flo$week),]
dfcorrect <- data.frame(week = rep(seq(1,52), 5), year = rep(seq(2016,2020,1), each = 52))
dffin <- merge(dfcorrect,flo, by=c("week","year"), all.x=T)
dffin <- dffin[order(dffin$year, dffin$week),]
flots <- dffin$percent_specimen_positive
flocases <- (flots/(mean(flots, na.rm=T)))*mean(subsetI[1:52]/pop) # scale % pos to match simulated incidence
timeuse <- seq(2016,2100,1/52)[1:length(subsetI)]

# plot Texas
pdf(paste0("TexasSim.pdf"),width=6,height=4)
par(mar=c(3,3,1,3))
plot(timeuse,subsetI/pop,type="n",col="#F64740",lwd =2,ylim=c(0,0.0035),xlab="",ylab=
       "",bty = "n",xaxs="i",xlim=c(2016,2030))
polygon(c(timeuse[52*4+controlWeekStart],timeuse[52*4+controlWeekStart],timeuse[52*4 + controlWeekStart + controlWeekLength],
          timeuse[52*4 + controlWeekStart + controlWeekLength]), c(-1,1,1,-1), border = NA,col="#D1D3DF")
lines(timeuse,subsetI/pop,col="#F64740",lwd =2,ylim=c(0,0.01),xlab="",ylab=
        "")
points(timeuse[1:length(flocases)],flocases,cex=0.5,col="grey32",pch=18)
abline(v = seq(2016,2100,1),lty=3,col="gray")
par(new = TRUE)
plot(timeuse, subsetS/pop, col="navy",type="l",xlab="",ylab="", 
     axes = F,ylim=c(0.7*min(subsetS/pop),max(subsetS/pop)), lwd = 2, lty = 2)
axis(side = 4)
title(xlab="Year", line = 2)
title(ylab="I/N", line = 2,col.lab="#F64740")
mtext(side = 4, line = 2, 'S/N',col="navy")
text(2016,0.16,"Texas",pos = 4,cex = 1.2)
dev.off()






