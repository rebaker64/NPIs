
############################
# This script makes the surface plots in Fig. 2a/b
############################


require("fields")

# source predtsirControl - this function is adapted from the tsiR package
source('predtsirControl.R', encoding = 'UTF-8')


load(file="Data/prelimRSVstatefitFlorida.RData") # load seasonal betas for Florida

sub <- outall[outall$state=="Florida",]
pop <- 21480000 # total population
births <-221542/52 # total weekly births
lengthvec <- seq(4,52,1)
betachangevec <- seq(0.98,0.02,-0.02)
matout <- matrix(NA, ncol = length(lengthvec), nrow = length(betachangevec))
matoutS <- matrix(NA, ncol = length(lengthvec), nrow = length(betachangevec))
matoutTiming <- matrix(NA, ncol = length(lengthvec), nrow = length(betachangevec))
minI <- matrix(NA, ncol = length(lengthvec), nrow = length(betachangevec))

for(l in 1:length(lengthvec)){
  for(b in 1:length(betachangevec)){
    
controlWeekStart = 11
controlWeekLength =  lengthvec[l]
times <- seq(1,55,1/52)
controlStart = 43*52 + controlWeekStart
controlEnd = 43*52 + controlWeekStart + controlWeekLength
pred <- predtsirControl(times = times, births = rep(births, length = length(times)), beta = sub$sea_beta, alpha = 0.97, 
                        S0 =floor(0.8*pop), I0 = floor(0.2*pop), nsim = 10, stochastic = F, 
                        controlStart = controlStart, controlEnd =controlEnd, betachange = betachangevec[b])
subsetI <- pred$I$mean[times > 40]
subsetS <- pred$S$mean[times > 40]
timeuse <- seq(2016,2100,1/52)[1:length(subsetI)]
subsetI <- subsetI[timeuse > 2020 & timeuse < 2025]
subsetS <- subsetS[timeuse > 2020 & timeuse < 2025]
matout[b,l] <- max(subsetI/pop)
matoutTiming[b,l] <- timeuse[timeuse > 2020 & timeuse < 2025][which.max(subsetI/pop)[1]]
matoutS[b,l] <- max(subsetS/pop)
image.plot(lengthvec,(1-betachangevec), t(matout))
minI[b,l] <- 0
if(min(subsetI) < 1){minI[b,l] <- 1}

  }
}


# plot results
library("fields")
require("pals")

#plot I
subpre2020 <-  max(pred$I$mean[times > 40 & times < 42]/pop)
matoutrel <- matout/subpre2020 
matoutrel[matoutrel > 9] <- 9
image.plot(lengthvec,(1-betachangevec)*100,t(matoutrel),ylab="",xlab="",yaxt="n", col =rev(brewer.spectral(100)),
           zlim=c(1,9), 
           axis.args = list(at = 1:9, labels=c(1:8,"> 9")))
contour(lengthvec,(1-betachangevec)*100,t(minI), nlevels = 1, add = T, lty = 2, drawlabels = F)
axis(2)
title(xlab="Length of control (weeks)",line = 2)
title(ylab="% Reduction in transmission",line = 2)


#plot S
subpre2020S <-  max(pred$S$mean[times > 40 & times < 42]/pop)
matoutrel <- matoutS/subpre2020S
image.plot(lengthvec,(1-betachangevec)*100,t(matoutrel),ylab="",xlab="",yaxt="n",
           col =rev(brewer.spectral(100)))
axis(2)
title(xlab="Length of control (weeks)",line = 2)
title(ylab="% Reduction in transmission",line = 2)

#plot time
subpre2020S <-  max(pred$S$mean[times > 40 & times < 42]/pop)
matoutrel <- matoutS/subpre2020S
image.plot(lengthvec,(1-betachangevec)*100,t(matoutTiming),ylab="",xlab="",yaxt="n",
           col =rev(brewer.spectral(100)),zlim=c(2020.827,2025))
axis(2)
title(xlab="Length of control (weeks)",line = 2)
title(ylab="% Reduction in transmission",line = 2)



###### TEXAS


load(file="Data/prelimRSVstatefitTEXAS.RData")

sub <- outall[outall$state=="Texas",] # load seasonal betas for Florida


pop <- 29000000
births = 378624/52


lengthvec <- seq(4,52,1)
betachangevec <- seq(0.98,0.02,-0.02)
matout <- matrix(NA, ncol = length(lengthvec), nrow = length(betachangevec))
matoutS <- matrix(NA, ncol = length(lengthvec), nrow = length(betachangevec))
matoutTiming <- matrix(NA, ncol = length(lengthvec), nrow = length(betachangevec))
minI <- matrix(NA, ncol = length(lengthvec), nrow = length(betachangevec))
for(l in 1:length(lengthvec)){
  for(b in 1:length(betachangevec)){
    
    controlWeekStart = 11
    controlWeekLength =  lengthvec[l]
    times <- seq(1,55,1/52)
    controlStart = 43*52 + controlWeekStart
    controlEnd = 43*52 + controlWeekStart + controlWeekLength
    pred <- predtsirControl(times = times, births = rep(births, length = length(times)), beta = sub$sea_beta, alpha = 0.97, 
                            S0 =floor(0.8*pop), I0 = floor(0.2*pop), nsim = 10, stochastic = F, 
                            controlStart = controlStart, controlEnd =controlEnd, betachange = betachangevec[b])
    subsetI <- pred$I$mean[times > 40]
    subsetS <- pred$S$mean[times > 40]
    timeuse <- seq(2016,2100,1/52)[1:length(subsetI)]
    subsetI <- subsetI[timeuse > 2020 & timeuse < 2025]
    subsetS <- subsetS[timeuse > 2020 & timeuse < 2025]
    matout[b,l] <- max(subsetI/pop)
    matoutTiming[b,l] <- timeuse[timeuse > 2020 & timeuse < 2025][which.max(subsetI/pop)[1]]
    matoutS[b,l] <- max(subsetS/pop)
    image.plot(lengthvec,(1-betachangevec), t(matout))
    minI[b,l] <- 0
    if(min(subsetI) < 1){minI[b,l] <- 1}
    
    
  }
}


# plot Texas

#plot I
subpre2020 <-  max(pred$I$mean[times > 40 & times < 42]/pop)
matoutrel <- matout/subpre2020 
require(pals)
image.plot(lengthvec,(1-betachangevec)*100,t(matoutrel),ylab="",xlab="",yaxt="n",
           col =rev(brewer.spectral(100)), zlim=c(1,9), 
           axis.args = list(at = 1:9, labels=c(1:9)))
contour(lengthvec,(1-betachangevec)*100,t(minI), nlevels = 1, add = T, lty = 2, drawlabels = F)
axis(2)
title(xlab="Length of control",line = 2)
title(ylab="% Reduction in transmission",line = 2)



#plot S
subpre2020S <-  max(pred$S$mean[times > 40 & times < 42]/pop)
matoutrel <- matoutS/subpre2020S
pdf("TXsensS.pdf",height=4.5,width=5)
image.plot(lengthvec,(1-betachangevec)*100,t(matoutrel),ylab="",xlab="",yaxt="n",col =rev(brewer.spectral(100)))
axis(2)
title(xlab="Length of control",line = 2)
title(ylab="% Reduction in transmission",line = 2)

#plot time
subpre2020S <-  max(pred$S$mean[times > 40 & times < 42]/pop)
matoutrel <- matoutS/subpre2020S
image.plot(lengthvec,(1-betachangevec)*100,t(matoutTiming),ylab="",xlab="",yaxt="n",
           col =rev(brewer.spectral(100)), 
           zlim=c(2020.827,2025))
axis(2)
title(xlab="Length of control",line = 2)
title(ylab="% Reduction in transmission",line = 2)

