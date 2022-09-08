no.rem <- function(R,t){
  
  # R: vector containing ordered removal data
  # t: vector containing time points where number of removals to be evaluated
  
  output <- rep(0,times=length(t))
  for (i in 1:length(t)) {
    output[i] <- sum(R<=t[i]) 
  }
  return(output)
}

plot.removal.curves.obs.rep.center.matched <- function(R.obs,R.STAR,shift,infe.period.distr,distance,
                                                       mar=c(4.1, 4.3, 2.1, 2.1)/1.1,mgp=c(3, 1, 0)/1.1,
                                                       las=0,cex.lab=1.2,cex.legend=1.2,width=4.50,height=3.30) {
  
  # R.obs: ordered observed removal times. Class "numeric"
  # R.STAR: matched predicted removal data from the model. Class "matrix"
  # shift: the shift to be applied on R.STAR. Class "numeric"
  # infe.period.distr: the infectious period distribution of the homogeneously mixing SIR. Only used for labelling plots. Class "character".
  # distance: the distance used. Only used for labelling plots. Class "character".
  
  # set some variables #
  
  n <- length(R.obs) # length of all realizations
  S <- nrow(R.STAR) # number of replicated realizations
  R.STAR <- R.STAR + shift # apply the shift on the replicated realizations
  R.bar <- apply(X=R.STAR,MARGIN=2,FUN=mean) # find the "center" (after shift has been applied)
  
  
  
  # plot #
  
  xlim <- c(min(R.obs,R.bar),max(R.obs,R.bar)) # reasonable window to plot
  ylim <- c(0,n)
  t <- R.STAR[1,] # only run -no.rem- function at points of change
  
  par(mar=mar, mgp=mgp, las=las)
  plot(t,no.rem(R=R.STAR[1,],t),col="cyan",lty=5,main="",xlab=expression(bold(t)),ylim=ylim,
       ylab = expression(bold(z[t])),type="s",lwd=0.5,xlim = xlim,cex.lab=cex.lab) # replicated
  
  for (k in 2:S) {
    t <- R.STAR[k,] # only run -no.rem- function at points of change
    lines(t,no.rem(R=R.STAR[k,],t),col="cyan",lty=5,lwd=0.5,type="s") # replicated
  }
  
  t <- R.obs # only run -no.rem- function at points of change
  lines(t,no.rem(R=R.obs,t),lwd=3,type="s") # observed
  t <- R.bar # only run -no.rem- function at points of change
  lines(t,no.rem(R=R.bar,t),col="red",lwd=3,type="s",lty=2) # averaged
  
  labels <- c(expression(z[t]^obs),expression(z[t]^rep),expression(bar(z)[t]^rep))
  legend("topleft",legend=labels,col=c("black","cyan","red"),lty=c(1,5,2),lwd=c(3,0.5,3),
         cex=cex.legend,bty="n")
  
}

distance.L2.matched <- function(R1,R2)  {
  
  #R1: the ordered time points of first removal sample path. Class "numeric". Must have same length as R2
  #R2: the ordered time points of the second removal sample path. Class "numeric". Must have same length as R1
  
  common.set <- sort(union(R1,R2)) # the ordered union of all removal points. Any possible change at any of the z_t occurs at these points
  z1_t <- no.rem(R=R1,t=common.set) # calculate z1_1 at the ordered union of all points
  z2_t <- no.rem(R=R2,t=common.set) # calculate z2_2 at the ordered union of all points
  interv.length <- diff(common.set)
  interv.length <- c(interv.length,0)
  distance <- sqrt(sum(interv.length*(z1_t-z2_t)^2))
  
  return(distance=distance)
  
}

calculate.shift.L2.matched <- function(R.obs,R.STAR) {
  
  # R.obs: ordered observed removal times. Class "numeric"
  # R.STAR: matched predicted removal data from the model. Class "matrix"
  
  # set some variables #
  
  n <- length(R.obs) # length of all realizations
  S <- nrow(R.STAR) # number of replicated realizations
  
  
  # calculate shift #
  
  shift <- rep(NA,times=S) # the shift is a vector; different for all realizations
  
  for (k in 1:S) {
    wrapper.distance <- function(c) {
      distance.L2.matched(R.obs,R.STAR[k,]+c)  
    }
    
    output.optim <- optimize(f=wrapper.distance,interval = c(-max(R.obs),max(R.obs)))
    shift[k] <- output.optim$minimum
  }
  
  return(shift=shift)
  
}

calculate.distances.obs.rep.center.L2.matched <- function(R.obs,R.STAR,shift) {
  
  # R.obs: ordered observed removal times. Class "numeric"
  # R.STAR: matched predicted removal data from the model. Class "matrix"
  # shift: the shift to be applied on R.STAR. Class "numeric"
  
  # set some variables #
  
  n <- length(R.obs) # length of all realizations
  S <- nrow(R.STAR) # number of replicated realizations
  R.STAR <- R.STAR + shift # apply the shift on the replicated realizations
  R.bar <- apply(X=R.STAR,MARGIN=2,FUN=mean) # find the "center" (after shift has been applied)
  
  # distance(z.obs,z.bar) #
  
  distance.z.obs.z.bar <- distance.L2.matched(R.obs,R.bar)
  
  
  # distance(z.obs,z.rep^(s)) and distance(z.bar,z.rep^(s)) #
  
  distance.z.obs.z.rep <- distance.z.bar.z.rep <- rep(NA,times=S)
  
  for (k in 1:S) {
    distance.z.obs.z.rep[k] <- distance.L2.matched(R.obs,R.STAR[k,])
    distance.z.bar.z.rep[k] <- distance.L2.matched(R.bar,R.STAR[k,])
  }
  
  output.list <- list(distance.z.obs.z.bar=distance.z.obs.z.bar,distance.z.obs.z.rep=distance.z.obs.z.rep,
                      distance.z.bar.z.rep=distance.z.bar.z.rep)
  
  return(output.list=output.list)
  
}

calculate.mid.p.value <- function(Y,y) {
  
  # Y: a sample (vector) from the one dimensional discrete distribution Y. Class "numeric"
  # y: the value at which the mid-p-value to be evaluated. Class "numeric"
  
  mid.p.value <- mean(Y<y) + mean(Y==y)*0.5
  
  return(mid.p.value)
  
}

calculate.position.time.distribution.z.obs.wrt.z.rep.matched <- function(R.obs,R.STAR,shift,t) {
  
  # R.obs: ordered observed removal times. Class "numeric"
  # R.STAR: matched predicted removal data from the model. Class "matrix"
  # shift: the shift to be applied on R.STAR. Class "numeric"
  # t: a vector of time points at each of which the position of z.obs wrt z.rep to be calculated. Class "numeric
  
  
  
  # set some variables #
  
  n <- length(R.obs) # length of all realizations
  S <- nrow(R.STAR) # number of replicated realizations
  R.STAR <- R.STAR + shift # apply the shift on the replicated realizations
  
  
  
  # calculate z.rep for each t #
  
  Z.REP <- matrix(nrow = S,ncol = length(t)) # create a matrix that will store number of removals for all replications
  
  for (k in 1:S) {
    Z.REP[k,] <- no.rem(R.STAR[k,],t)
  }
  
  
  
  # calculate z.obs #
  
  z.obs <- no.rem(R=R.obs,t=t)
  
  
  
  # calculate mid-ppp-value of z.obs wrt z.rep for each t #
  
  z.obs.mid.ppp.value <- rep(NA,times=length(z.obs))
  
  for (k in 1:length(t)) {
    z.obs.mid.ppp.value[k] <- calculate.mid.p.value(Y=Z.REP[,k], y=z.obs[k])  
  }
  
  
  
  # output #
  
  output.list <- list(z.obs.mid.ppp.value=z.obs.mid.ppp.value)
  
  return(output.list)
  
}

calculate.proportion.of.time.z.obs.spends.in.given.position.interval <- function(z.obs.position.time.distribution,
                                                                                 prob.lower,prob.upper,interval.type,S) {
  
  
  # z.obs.time.distribution: a vector giving the position of z.obs wrt z.rep for each time point; use output of 
  #                          -calculate.position.time.distribution.z.obs.wrt.z.rep.matched- function. Class "list"
  # prob.lower: the probability (inverse "quantile") specifying the lower limit of the position interval. Class "numeric"
  # prob.upper: the probability (inverse "quantile") specifying the uppper limit of the position interval. Class "numeric"  
  # interval.type: a character specifying the type of position interval; "closed" for [prob.lower,prob.upper], "open" for (prob.lower,prob.upper),
  #                "left.open" for (prob.lower,prob.upper] and "right.open" for [prob.lower,prob.upper). Class "character"

  if (interval.type=="closed") {
    prop.of.time.mid.ppp.value.based  <- mean(z.obs.position.time.distribution$z.obs.mid.ppp.value >= prob.lower & z.obs.position.time.distribution$z.obs.mid.ppp.value <= prob.upper)
  }
  
  if (interval.type=="open") {
    prop.of.time.mid.ppp.value.based  <- mean(z.obs.position.time.distribution$z.obs.mid.ppp.value > prob.lower & z.obs.position.time.distribution$z.obs.mid.ppp.value < prob.upper)
  }
  
  if (interval.type=="left.open") {
    prop.of.time.mid.ppp.value.based  <- mean(z.obs.position.time.distribution$z.obs.mid.ppp.value > prob.lower & z.obs.position.time.distribution$z.obs.mid.ppp.value <= prob.upper)
  }
  
  if (interval.type=="right.open") {
    prop.of.time.mid.ppp.value.based  <- mean(z.obs.position.time.distribution$z.obs.mid.ppp.value >= prob.lower & z.obs.position.time.distribution$z.obs.mid.ppp.value < prob.upper)
  }
  
  
  # output #
  
  output <- prop.of.time.mid.ppp.value.based
  names(output) <- "prop.of.time.mid.ppp.value.based"
  
  return(output)
  
}

distance.method <- function(R.obs, R.STAR,
                            mar=c(4.1, 4.3, 2.1, 2.1)/1.1,mgp=c(3, 1, 0)/1.1,
                            las=0,cex.lab=1.2,cex.legend=1.2,width=4.50,height=3.30) {
  
  # R.obs: ordered observed removal times. Class "numeric"
  # R.STAR: matched predicted removal data from the model. Class "matrix"
  
  
  # calculate the shift #
  
  shift <- calculate.shift.L2.matched(R.obs=R.obs,R.STAR = R.STAR)
  
  
  # removal curve visual assessment plot #
  
  plot.removal.curves.obs.rep.center.matched(R.obs=R.obs,R.STAR = R.STAR,shift = shift,
                                             infe.period.distr = "",distance="",mar=mar,mgp=mgp,
                                             las=las,cex.lab=cex.lab,cex.legend=cex.legend,width=width,height=height)
  
  
  # distance method #
  
  distances.obs.rep.center <- calculate.distances.obs.rep.center.L2.matched(R.obs=R.obs,R.STAR = R.STAR,shift = shift)
  
  xlim <- range(c(distances.obs.rep.center$distance.z.bar.z.rep,distances.obs.rep.center$distance.z.obs.z.bar))
  par(mar=mar, mgp=mgp, las=las)
  hist(distances.obs.rep.center$distance.z.bar.z.rep,xlab=expression(bold(T[d]^rep)), prob=T,breaks="fd",main="",xlim=xlim,cex.lab=cex.lab,font.lab=2)
  abline(v=distances.obs.rep.center$distance.z.obs.z.bar,col="black",lty=2,lwd=2)
  labels <- c(expression(T[d]^obs))
  legend("topright",legend=labels,col="black",lty=2,lwd=2,cex=cex.legend,bty="n")
  fpppv <- mean((distances.obs.rep.center$distance.z.bar.z.rep <= distances.obs.rep.center$distance.z.obs.z.bar))
  mtext(paste("fppp-value=",round(fpppv,2),sep = ""),side=4,cex = cex.legend)
  print(paste("fppp-value=",round(fpppv,2),sep = ""))
}

position.time.method <- function(R.obs,R.STAR,
                                 mar=c(4.1, 4.3, 2.1, 2.1)/1.1,mgp=c(3, 1, 0)/1.1,
                                 las=0,cex.lab=1.2,cex.legend=1.2,width=4.50,height=3.30) {
  
  
  # calculate the shift #
  
  shift <- calculate.shift.L2.matched(R.obs=R.obs,R.STAR = R.STAR)
  
  
  # removal curve visual assessment plot #
  
  plot.removal.curves.obs.rep.center.matched(R.obs=R.obs,R.STAR = R.STAR,shift = shift,
                                             infe.period.distr = "",distance="",mar=mar,mgp=mgp,
                                             las=las,cex.lab=cex.lab,cex.legend=cex.legend,width=width,height=height)
  
  
  # position-time method #
  
  n.obs <- length(R.obs)
  t <- seq(R.obs[1],R.obs[n.obs],length.out = 500) 
  z.obs.position.time.distribution <- calculate.position.time.distribution.z.obs.wrt.z.rep.matched(R=R.obs,R.STAR = R.STAR,shift = shift,t=t)
  par(mar=mar, mgp=mgp, las=las)
  plot(x=t,y=z.obs.position.time.distribution$z.obs.mid.ppp.value,type="l",ylim = c(0,1),xlab="t",ylab="ppp-value(t)",font.lab=2,cex.lab=cex.lab) # mid.ppp.value history plot for observed path
  abline(h=c(0.5,0.95,0.05),col="red",lty=2)
  
  position.prob <- seq(0,1,by=0.1)
  prop.of.time.mid.ppp.value.based <- rep(NA,times=length(diff(position.prob)))
  
  prop.of.time.mid.ppp.value.based[1] <- calculate.proportion.of.time.z.obs.spends.in.given.position.interval(prob.lower = 0,prob.upper = 0.1,z.obs.position.time.distribution = z.obs.position.time.distribution,
                                                                                                              interval.type = "closed",S=S)[1]
  for (k in 2:length(prop.of.time.mid.ppp.value.based) ) {
    prop.of.time.mid.ppp.value.based[k] <- calculate.proportion.of.time.z.obs.spends.in.given.position.interval(prob.lower = position.prob[k],prob.upper = position.prob[k+1],
                                                                                                                z.obs.position.time.distribution = z.obs.position.time.distribution,
                                                                                                                interval.type = "left.open",S=S)[1]
  }
  
  names(prop.of.time.mid.ppp.value.based) <- c("[0,0.1]","(0.1,0.2]","(0.2,0.3]","(0.3,0.4]","(0.4,0.5]","(0.5,0.6]","(0.6,0.7]","(0.7,0.8]","(0.8,0.9]","(0.9,1]")
  print(prop.of.time.mid.ppp.value.based)
  sqrtmse <- sqrt(mean((z.obs.position.time.distribution$z.obs.mid.ppp.value-0.5)^2))
  print(paste("sqrtmse=",round(sqrtmse,2),sep=""))
}


