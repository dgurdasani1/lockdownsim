library("EnvStats", lib="/data/home/hhy494/R/x86_64-pc-linux-gnu-library/3.3/")
library(rstan)
library(data.table)
library(lubridate)
library(gdata)


 #rlnormAlt(30, mean = 10, cv = 2) - mean of SI is 4.7, and SD is 2.9 days - cv is sd/mean=0.62
 d<-read.table("data_ONS.txt", T)
 #d$X5to14[76]=4140
d=d[is.na(d$ONS_deaths)==0, ]
#identify on which day in the cumulative no of deaths is >=10
index1=which(cumsum(d$ONS_deaths)>=10)[1]
index2=1
N=nrow(d)
N2=N+90

library(arules)

#generate serial interval distribution - main model
x3=rlnormAlt(5e6,4.7,0.62)
#sensitivity analysis - longer serial interval - hongkong paper
#x3=rgammaAlt(5e6,6.5,0.72)
 g = ecdf(x3)

k = rep(0,N2)
#k[1] = (g(1.5) - g(0)) 
   for(i in 1:length(k)) {
      k[i] = (g(i) - g(i-1))
    }
serial.interval=k


mean1 = 5.1; cv1 = 0.86; # infection to onset
    mean2 = 18.8; cv2 = 0.45 # onset to death
    ## assume that CFR is probability of dying given infection
    x1 = rgammaAlt(5e6,mean1,cv1) # infection-to-onset ----> do all people who are infected get to onset?
    x2 = rgammaAlt(5e6,mean2,cv2) # onset-to-death
    f = ecdf(x1+x2)

    CFR=0.0088

    convolution = function(u) ((CFR) * f(u))
   
   h = rep(0,N2) 
  h[1] = (convolution(1.5) - convolution(0)) 
    for(i in 2:length(h)) {
      h[i] = (convolution(i+.5) - convolution(i-.5)) / (1-convolution(i-.5))
    }
s = rep(0,N2)
  s[1] = 1 
  for(i in 2:length(s)) {
    s[i] = s[i-1]*(1-h[i-1])
  }
  f = s * h

  
stan_data = list(N=nrow(d),M=48,N2=N+90,y=0,deaths=NULL,f=NULL,N0=6,SI=serial.interval[1:N2]) # N0 = 6 to make it consistent with Rayleigh

stan_data$f = f
  stan_data$N2=N2
  stan_data$x=1:N2
  stan_data$deaths=c(d$ONS_deaths, rep(-1, 90))


  m = stan_model('lockdownfinal.stan')
   fit = stan('lockdownfinal.stan',data=stan_data,iter=4000,warmup=1000,chains=4,thin=4)
   out = rstan::extract(fit)


library(bayesplot)

dates=rep(0, N2)
 dates[1]=as.character(as.Date(d$Date[1], "%d/%m/%Y"))
for (i in 2:N2) {
	dates[i]=as.character(as.Date(d$Date[1], "%d/%m/%Y") + i-1)
}
dfdeaths=NULL
for (i in seq(1:48)) {
x=as.matrix(out$E_deaths0[,,i])
colnames(x)=dates
g = (mcmc_intervals(x,prob = .95))
z=as.data.frame(g$data)
dfdeaths=cbind(dfdeaths, z$ll, z$m, z$hh)
}
 x=cbind(dfdeaths[1:111,2], d$ONS_deaths)
 x=as.data.frame(x)
 x$diffsq=abs(x$V2-x$V1)^2
 RMSE=sqrt(mean(x$diffsq))



library("loo")
 LLarray <- loo::extract_log_lik(stanfit = fit, parameter_name = "lp1",  merge_chains = FALSE)
 r_eff <- loo::relative_eff(x = exp(LLarray))
 loo::loo.array(LLarray[,,2:111],r_eff=r_eff[2:111])
  loor<-  loo::loo.array(LLarray[,,2:111],r_eff=r_eff[2:111])


fit2<-readRDS("fitfinal_altmodel4.rds")
LLarray <- loo::extract_log_lik(stanfit = fit2, parameter_name = "lp1",  merge_chains = FALSE)
 r_eff <- loo::relative_eff(x = exp(LLarray))
  loor2<-  loo::loo.array(LLarray[,,2:111],r_eff=r_eff[2:111])

  comp <- loo_compare(loor, loor2)


library(bayesplot)



Rt = (as.matrix(out$Rt))
dates=rep(0, N2)
 dates[1]=as.character(as.Date(d$Date[1], "%d/%m/%Y"))
for (i in 2:N2) {
	dates[i]=as.character(as.Date(d$Date[1], "%d/%m/%Y") + i-1)
}
colnames(Rt) = dates
q = (mcmc_intervals(Rt, prob = .95))
ggsave("Rt.pdf",q,width=4,height=6)
dev.off()
df=as.data.frame(q$data)
write.table(df, "predictionRt.txt", row.names=F, col.names=T, quote=F)


dfdeaths=NULL
for (i in seq(1:48)) {
x=as.matrix(out$E_deaths0[,,i])
colnames(x)=dates
g = (mcmc_intervals(x,prob = .95))
z=as.data.frame(g$data)
dfdeaths=cbind(dfdeaths, z$ll, z$m, z$hh)
}
write.table(dfdeaths, "predictiondeaths_longSI.txt", row.names=F, col.names=T, quote=F)


z=z$m[1:N]
jpeg("deathscomp.jpeg")
x=seq(1:nrow(d))
g <- ggplot(d, aes(x))
g<- g + geom_line(aes(y=d$ONS_deaths), colour="red")
g <- g + geom_line(aes(y=z), colour="blue")
g
dev.off()

dfprediction0=NULL
dfcases=NULL
for (i in seq(1:48)) {
x=as.matrix(out$prediction0[,,i])
colnames(x)=dates
g = (mcmc_intervals(x,prob = .95))
z=as.data.frame(g$data)
dfprediction0=cbind(dfprediction0, z$ll, z$m, z$hh)
}
write.table(dfprediction0, "prediction0cases.txt", row.names=F, col.names=T, quote=F)



dfdiffcsdeaths=NULL
for (i in seq(1:48)) {
x=as.matrix(out$diffcsdeaths81[,i])
colnames(x)=i
g = (mcmc_intervals(x,prob = .95))
z=as.data.frame(g$data)
dfdiffcsdeaths=cbind(dfdiffcsdeaths, z$ll, z$m, z$hh)
}
write.table(dfdiffcsdeaths, "diffcsdeaths81.txt", row.names=F, col.names=T, quote=F)


dfdiffcscases=NULL
for (i in seq(1:48)) {
x=as.matrix(out$diffcscases81[,i])
colnames(x)=i
g = (mcmc_intervals(x,prob = .95))
z=as.data.frame(g$data)
dfdiffcscases=cbind(dfdiffcscases, z$ll, z$m, z$hh)
}
write.table(dfdiffcscases, "diffcscases81.txt", row.names=F, col.names=T, quote=F)



dfcsdeaths=NULL
for (i in seq(1:48)) {
x=as.matrix(out$csdeaths[,i])
colnames(x)=i
g = (mcmc_intervals(x,prob = .95))
z=as.data.frame(g$data)
dfcsdeaths=cbind(dfcsdeaths, z$ll, z$m, z$hh)
}
write.table(dfcsdeaths, "csdeaths.txt", row.names=F, col.names=T, quote=F)



dfcscases=NULL
for (i in seq(1:48)) {
x=as.matrix(out$cscases[,i])
colnames(x)=i
g = (mcmc_intervals(x,prob = .95))
z=as.data.frame(g$data)
dfcscases=cbind(dfcscases, z$ll, z$m, z$hh)
}
write.table(dfcscases, "cscases.txt", row.names=F, col.names=T, quote=F)


dfRt0=NULL
for (i in seq(1:48)) {
x=as.matrix(out$Rt0[,,i])
colnames(x)=dates
g = (mcmc_intervals(x,prob = .95))
z=as.data.frame(g$data)
dfRt0=cbind(dfRt0, z$m)
}
write.table(dfRt0, "dfRt0.txt", row.names=F, col.names=T, quote=F)


 fit_summary<-summary(fit)
fit_summary$summary
write.table(as.array(fit), "fit.txt", row.names=F, col.names=T, quote=F)
 posterior<-as.array(fit)
  x=as.data.frame(fit_summary$summary)


jpeg("rhatdist.jpeg", res=300, height=1500, width=1500)
ggplot(x, aes(x=x$Rhat)) + 
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.0001,
                   colour="black", fill="white") +
    xlim(0.995,1.005)
    dev.off()


#prediction figures


dfdeaths=NULL
for (i in seq(1:48)) {
x=as.matrix(out$E_deaths0[,,i])
colnames(x)=dates
g = (mcmc_intervals(x,prob = .95))
z=as.data.frame(g$data)
dfdeaths=cbind(dfdeaths, z$m)
}



dfcases=NULL
for (i in seq(1:48)) {
x=as.matrix(out$prediction0[,,i])
colnames(x)=dates
g = (mcmc_intervals(x,prob = .95))
z=as.data.frame(g$data)
dfcases=cbind(dfcases, z$m)
}



jpeg("pred1June15June1July.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfdeaths)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=l[,24]), colour="forest green")
k <- k + geom_line(aes(y=l[,28]), colour="dark blue")
k <- k + geom_line(aes(y=l[,32]), colour="red")
k <- k + geom_line(aes(y=l[,36]), colour="purple")
k <- k + geom_line(aes(y=l[,40]), colour="brown")
k <- k + geom_line(aes(y=l[,45]), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Incident Deaths")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()

jpeg("predcase1June15June1July.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfcases)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=l[,24]), colour="forest green")
k <- k + geom_line(aes(y=l[,28]), colour="dark blue")
k <- k + geom_line(aes(y=l[,32]), colour="red")
k <- k + geom_line(aes(y=l[,36]), colour="purple")
k <- k + geom_line(aes(y=l[,40]), colour="brown")
k <- k + geom_line(aes(y=l[,45]), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Incident Cases")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()

jpeg("pred1June15June.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfdeaths)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=l[,23]), colour="forest green")
k <- k + geom_line(aes(y=l[,27]), colour="dark blue")
k <- k + geom_line(aes(y=l[,31]), colour="red")
k <- k + geom_line(aes(y=l[,35]), colour="purple")
k <- k + geom_line(aes(y=l[,39]), colour="brown")
k <- k + geom_line(aes(y=l[,45]), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Incident Deaths")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()


jpeg("predcase1June15June.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfcases)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=l[,23]), colour="forest green")
k <- k + geom_line(aes(y=l[,27]), colour="dark blue")
k <- k + geom_line(aes(y=l[,31]), colour="red")
k <- k + geom_line(aes(y=l[,35]), colour="purple")
k <- k + geom_line(aes(y=l[,39]), colour="brown")
k <- k + geom_line(aes(y=l[,45]), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Incident Cases")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()


jpeg("pred1June3July.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfdeaths)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=l[,22]), colour="forest green")
k <- k + geom_line(aes(y=l[,26]), colour="dark blue")
k <- k + geom_line(aes(y=l[,30]), colour="red")
k <- k + geom_line(aes(y=l[,34]), colour="purple")
k <- k + geom_line(aes(y=l[,38]), colour="brown")
k <- k + geom_line(aes(y=l[,45]), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Incident Deaths")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()


jpeg("predcase1June3July.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfcases)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=l[,22]), colour="forest green")
k <- k + geom_line(aes(y=l[,26]), colour="dark blue")
k <- k + geom_line(aes(y=l[,30]), colour="red")
k <- k + geom_line(aes(y=l[,34]), colour="purple")
k <- k + geom_line(aes(y=l[,38]), colour="brown")
k <- k + geom_line(aes(y=l[,45]), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Incident Cases")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()


jpeg("predcase1June.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfcases)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=l[,21]), colour="forest green")
k <- k + geom_line(aes(y=l[,25]), colour="dark blue")
k <- k + geom_line(aes(y=l[,29]), colour="red")
k <- k + geom_line(aes(y=l[,33]), colour="purple")
k <- k + geom_line(aes(y=l[,37]), colour="brown")
k <- k + geom_line(aes(y=l[,1]), colour="aquamarine")
k <- k + geom_line(aes(y=l[,9]), colour="darkorange1")
k <- k + geom_line(aes(y=l[,45]), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Incident Cases")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()

jpeg("pred1June.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfdeaths)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=l[,21]), colour="forest green")
k <- k + geom_line(aes(y=l[,25]), colour="dark blue")
k <- k + geom_line(aes(y=l[,29]), colour="red")
k <- k + geom_line(aes(y=l[,33]), colour="purple")
k <- k + geom_line(aes(y=l[,37]), colour="brown")
k <- k + geom_line(aes(y=l[,1]), colour="aquamarine")
k <- k + geom_line(aes(y=l[,9]), colour="darkorange1")
k <- k + geom_line(aes(y=l[,45]), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Incident Deaths")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()



jpeg("predcsdeaths1June.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfdeaths)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=cumsum(l[,21])), colour="forest green")
k <- k + geom_line(aes(y=cumsum(l[,25])), colour="dark blue")
k <- k + geom_line(aes(y=cumsum(l[,29])), colour="red")
k <- k + geom_line(aes(y=cumsum(l[,33])), colour="purple")
k <- k + geom_line(aes(y=cumsum(l[,37])), colour="brown")
k <- k + geom_line(aes(y=cumsum(l[,1])), colour="aquamarine")
k <- k + geom_line(aes(y=cumsum(l[,9])), colour="darkorange1")
k <- k + geom_line(aes(y=cumsum(l[,45])), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Cumulative Deaths")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()



jpeg("predcscases1June.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfcases)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=cumsum(l[,21])), colour="forest green")
k <- k + geom_line(aes(y=cumsum(l[,25])), colour="dark blue")
k <- k + geom_line(aes(y=cumsum(l[,29])), colour="red")
k <- k + geom_line(aes(y=cumsum(l[,33])), colour="purple")
k <- k + geom_line(aes(y=cumsum(l[,37])), colour="brown")
k <- k + geom_line(aes(y=cumsum(l[,1])), colour="aquamarine")
k <- k + geom_line(aes(y=cumsum(l[,9])), colour="darkorange1")
k <- k + geom_line(aes(y=cumsum(l[,45])), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Cumulative Cases")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()


jpeg("predcsdeaths1June3July.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfdeaths)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=cumsum(l[,22])), colour="forest green")
k <- k + geom_line(aes(y=cumsum(l[,26])), colour="dark blue")
k <- k + geom_line(aes(y=cumsum(l[,30])), colour="red")
k <- k + geom_line(aes(y=cumsum(l[,34])), colour="purple")
k <- k + geom_line(aes(y=cumsum(l[,38])), colour="brown")
k <- k + geom_line(aes(y=cumsum(l[,45])), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Cumulative Deaths")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()

jpeg("predcscases1June3July.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfcases)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=cumsum(l[,22])), colour="forest green")
k <- k + geom_line(aes(y=cumsum(l[,26])), colour="dark blue")
k <- k + geom_line(aes(y=cumsum(l[,30])), colour="red")
k <- k + geom_line(aes(y=cumsum(l[,34])), colour="purple")
k <- k + geom_line(aes(y=cumsum(l[,38])), colour="brown")
k <- k + geom_line(aes(y=cumsum(l[,45])), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Cumulative Cases")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()


jpeg("predcsdeaths1June15June.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfdeaths)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=cumsum(l[,23])), colour="forest green")
k <- k + geom_line(aes(y=cumsum(l[,27])), colour="dark blue")
k <- k + geom_line(aes(y=cumsum(l[,31])), colour="red")
k <- k + geom_line(aes(y=cumsum(l[,35])), colour="purple")
k <- k + geom_line(aes(y=cumsum(l[,39])), colour="brown")
k <- k + geom_line(aes(y=cumsum(l[,45])), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Cumulative Deaths")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()

jpeg("predcscases1June15June.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfcases)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=cumsum(l[,23])), colour="forest green")
k <- k + geom_line(aes(y=cumsum(l[,27])), colour="dark blue")
k <- k + geom_line(aes(y=cumsum(l[,31])), colour="red")
k <- k + geom_line(aes(y=cumsum(l[,35])), colour="purple")
k <- k + geom_line(aes(y=cumsum(l[,39])), colour="brown")
k <- k + geom_line(aes(y=cumsum(l[,45])), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Cumulative Cases")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()



jpeg("predcsdeaths1June15June3July.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfdeaths)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=cumsum(l[,24])), colour="forest green")
k <- k + geom_line(aes(y=cumsum(l[,28])), colour="dark blue")
k <- k + geom_line(aes(y=cumsum(l[,32])), colour="red")
k <- k + geom_line(aes(y=cumsum(l[,36])), colour="purple")
k <- k + geom_line(aes(y=cumsum(l[,40])), colour="brown")
k <- k + geom_line(aes(y=cumsum(l[,45])), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Cumulative Deaths")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()



jpeg("predcscases1June15June3July.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=as.data.frame(dfcases)
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=cumsum(l[,24])), colour="forest green")
k <- k + geom_line(aes(y=cumsum(l[,28])), colour="dark blue")
k <- k + geom_line(aes(y=cumsum(l[,32])), colour="red")
k <- k + geom_line(aes(y=cumsum(l[,36])), colour="purple")
k <- k + geom_line(aes(y=cumsum(l[,40])), colour="brown")
k <- k + geom_line(aes(y=cumsum(l[,45])), colour="black")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Cumulative Cases")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()






jpeg("Rt1June15June1July.jpeg", height=1500, width=1500, res=300 )
x=dates[95:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
q=dfRt0[95:nrow(dfRt0), ]
m<-ggplot(as.data.frame(q), aes(x, group=1))
m<-m+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
m<-m+ geom_line(aes(y=q[,24]), colour="forest green")
m<-m+ geom_line(aes(y=q[,28]), colour="dark blue")
m <-m+geom_line(aes(y=q[,32]), colour="red")
m <-m+geom_line(aes(y=q[,36]), colour="purple")
m <-m+geom_line(aes(y=q[,40]), colour="brown")
m <-m+geom_line(aes(y=rep(0.8122, 107)), colour="black")
#m<-m+geom_hline(aes(yintercept=0.81), colour="black", linetype="solid")
m <-m+geom_vline(aes(xintercept=x[20]), color="#999999", linetype="dashed")
m <-m+geom_vline(aes(xintercept=x[34]), color="#999999", linetype="dashed")
m <-m+geom_vline(aes(xintercept=x[52]), color="#999999", linetype="dashed")
m<-m+ylim(0.6, 1.2)
m <- m +theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
m<-m+labs(x = "Date")
m<-m+labs(y="Rt")
m
dev.off()



jpeg("predcs1June15June1July.jpeg", height=1500, width=1500, res=300)
x=dates[1:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
l=g[, c(19:88)]
k<- ggplot(l, aes(x))
k<-k+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
k <- k + geom_line(aes(y=g$csbasedeath), colour="forest green")
k <- k + geom_line(aes(y=g$csRt859095death), colour="dark blue")
k <- k + geom_line(aes(y=g$csRt90951death), colour="red")
k <- k + geom_line(aes(y=g$csRt951105death), colour="purple")
k <- k + geom_line(aes(y=g$csRt110511death), colour="coral1")
k <- k + geom_line(aes(y=g$csRt10511115death), colour="brown")
k <- k +  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
k<-k + labs(x = "")
k<-k+ labs(y="Cumulative Deaths")
k<-k+theme(axis.text.x = element_text(angle = 45))
k
dev.off()



jpeg("Rt1June15June1July.jpeg", height=1500, width=1500, res=300 )
x=dates[95:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
q=dfRt0[95:nrow(dfRt0), ]
m<-ggplot(as.data.frame(q), aes(x, group=1))
m<-m+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
m<-m+ geom_line(aes(y=q[,24]), colour="forest green")
m<-m+ geom_line(aes(y=q[,28]), colour="dark blue")
m <-m+geom_line(aes(y=q[,32]), colour="red")
m <-m+geom_line(aes(y=q[,36]), colour="purple")
m <-m+geom_line(aes(y=q[,36]), colour="coral1")
m <-m+geom_line(aes(y=q[,40]), colour="brown")
m <-m+geom_vline(aes(xintercept=x[20]), color="#999999", linetype="dashed")
m <-m+geom_vline(aes(xintercept=x[34]), color="#999999", linetype="dashed")
m <-m+geom_vline(aes(xintercept=x[52]), color="#999999", linetype="dashed")
m<-m+ylim(0.6, 1.2)
m <- m +theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
m<-m+labs(x = "Date")
m<-m+labs(y="Rt")
m
dev.off()



jpeg("Rt1June15June1July.jpeg", height=1500, width=1500, res=300 )
x=dates[95:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
q=dfRt0[95:nrow(dfRt0), ]
m<-ggplot(as.data.frame(q), aes(x, group=1))
m<-m+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
m<-m+ geom_line(aes(y=q[,24]), colour="forest green")
m<-m+ geom_line(aes(y=q[,28]), colour="dark blue")
m <-m+geom_line(aes(y=q[,32]), colour="red")
m <-m+geom_line(aes(y=q[,36]), colour="purple")
m <-m+geom_line(aes(y=q[,40]), colour="brown")
m <-m+geom_line(aes(y=rep(0.8122, 107)), colour="black")
#m<-m+geom_hline(aes(yintercept=0.81), colour="black", linetype="solid")
m <-m+geom_vline(aes(xintercept=x[20]), color="#999999", linetype="dashed")
m <-m+geom_vline(aes(xintercept=x[34]), color="#999999", linetype="dashed")
m <-m+geom_vline(aes(xintercept=x[52]), color="#999999", linetype="dashed")
m<-m+ylim(0.6, 1.2)
m <- m +theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
m<-m+labs(x = "Date")
m<-m+labs(y="Rt")
m
dev.off()



jpeg("Rt1June15June.jpeg", height=1500, width=1500, res=300 )
x=dates[95:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
q=dfRt0[95:nrow(dfRt0), ]
m<-ggplot(as.data.frame(q), aes(x, group=1))
m<-m+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
m<-m+ geom_line(aes(y=q[,23]), colour="forest green")
m<-m+ geom_line(aes(y=q[,27]), colour="dark blue")
m <-m+geom_line(aes(y=q[,31]), colour="red")
m <-m+geom_line(aes(y=q[,35]), colour="purple")
m <-m+geom_line(aes(y=q[,39]), colour="brown")
m <-m+geom_line(aes(y=rep(0.8122, 107)), colour="black")
#m<-m+geom_hline(aes(yintercept=0.81), colour="black", linetype="solid")
m <-m+geom_vline(aes(xintercept=x[20]), color="#999999", linetype="dashed")
m <-m+geom_vline(aes(xintercept=x[34]), color="#999999", linetype="dashed")
m <-m+geom_vline(aes(xintercept=x[52]), color="#999999", linetype="dashed")
m<-m+ylim(0.6, 1.2)
m <- m +theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
m<-m+labs(x = "Date")
m<-m+labs(y="Rt")
m
dev.off()



jpeg("Rt15June.jpeg", height=1500, width=1500, res=300 )
x=dates[95:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
q=dfRt0[95:nrow(dfRt0), ]
m<-ggplot(as.data.frame(q), aes(x, group=1))
m<-m+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
m<-m+ geom_line(aes(y=q[,22]), colour="forest green")
m<-m+ geom_line(aes(y=q[,26]), colour="dark blue")
m <-m+geom_line(aes(y=q[,30]), colour="red")
m <-m+geom_line(aes(y=q[,34]), colour="purple")
m <-m+geom_line(aes(y=q[,38]), colour="brown")
m <-m+geom_line(aes(y=rep(0.8122, 107)), colour="black")
#m<-m+geom_hline(aes(yintercept=0.81), colour="black", linetype="solid")
m <-m+geom_vline(aes(xintercept=x[20]), color="#999999", linetype="dashed")
m <-m+geom_vline(aes(xintercept=x[34]), color="#999999", linetype="dashed")
m <-m+geom_vline(aes(xintercept=x[52]), color="#999999", linetype="dashed")
m<-m+ylim(0.6, 1.2)
m <- m +theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
m<-m+labs(x = "Date")
m<-m+labs(y="Rt")
m
dev.off()


jpeg("Rt1June.jpeg", height=1500, width=1500, res=300 )
x=dates[95:length(dates)]
 x=format(as.Date(x))
 x=as.Date(x)
q=dfRt0[95:nrow(dfRt0), ]
m<-ggplot(as.data.frame(q), aes(x, group=1))
m<-m+scale_x_date(breaks = x[seq(1, length(x), by = 28)])
m<-m+ geom_line(aes(y=q[,21]), colour="forest green")
m<-m+ geom_line(aes(y=q[,25]), colour="dark blue")
m <-m+geom_line(aes(y=q[,29]), colour="red")
m <-m+geom_line(aes(y=q[,33]), colour="purple")
m <-m+geom_line(aes(y=q[,37]), colour="brown")
m <- m + geom_line(aes(y=q[,1]), colour="aquamarine")
m <- m + geom_line(aes(y=q[,9]), colour="darkorange1")
m <-m+geom_line(aes(y=rep(0.8122, 107)), colour="black")
#m<-m+geom_hline(aes(yintercept=0.81), colour="black", linetype="solid")
m <-m+geom_vline(aes(xintercept=x[20]), color="#999999", linetype="dashed")
m <-m+geom_vline(aes(xintercept=x[34]), color="#999999", linetype="dashed")
m <-m+geom_vline(aes(xintercept=x[52]), color="#999999", linetype="dashed")
m<-m+ylim(0.6, 1.2)
m <- m +theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
m<-m+labs(x = "Date")
m<-m+labs(y="Rt")
m
dev.off()



d<-read.table("predictioncases.txt", T)
e<-read.table("predictionRt.txt", T)
jpeg("predcases_ci.jpg", height=1500, width=1500, res=300)

d$parameter=as.Date(e$parameter)
d=d[1:111, ]
q<-ggplot(d)
 q<-q+geom_line(aes(y=m, x=parameter, colour="red")) +
     geom_ribbon(aes(ymin = l, ymax = h, x=parameter, fill="band"), alpha = .3) +
     scale_colour_manual("",values="blue")+
    scale_fill_manual("",values="grey12") +
     theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA)) +
     theme(axis.line = element_line(color="black", size = 0.5)) +
     labs(x = "") + labs(y="Predicted Incident Cases")
  #   geom_vline(aes(xintercept=d$parameter[44]), color="#999999", linetype="dashed") +
   #  geom_vline(aes(xintercept=d$parameter[37]), color="#999999", linetype="dashed")
 q
dev.off()

d<-read.table("predictiondeaths.txt", T)
e<-read.table("predictionRt.txt", T)
f<-read.table("data_ONS.txt", T)
d$parameter=as.Date(e$parameter)
jpeg("preddeaths_ci_longSI.jpg", height=1500, width=1500, res=300)
d=d[1:111, ]
q<-ggplot(d)
 q<-q+geom_line(aes(y=V2, x=parameter), colour="blue") +
     geom_ribbon(aes(ymin = V1, ymax = V3, x=parameter, fill="band"), alpha = .3) +
     scale_colour_manual("",values="blue")+
    scale_fill_manual("",values="grey12") +
     geom_line(aes(y=f$ONS_deaths[1:111], x=parameter), colour="red") +
     theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA)) +
     theme(axis.line = element_line(color="black", size = 0.5)) +
     labs(x = "") + labs(y="Predicted Incident Deaths")
  #   geom_vline(aes(xintercept=d$parameter[44]), color="#999999", linetype="dashed") +
   #  geom_vline(aes(xintercept=d$parameter[37]), color="#999999", linetype="dashed")
 q
dev.off()



d<-read.table("predictionRt.txt", T)
jpeg("Rt_ci_broadprior.jpg", height=1500, width=1500, res=300)
d=d[1:111, ]
d$parameter=as.Date(d$parameter)
q<-ggplot(d)
 q<-q+geom_line(aes(y=m, x=parameter, colour="red")) +
     geom_ribbon(aes(ymin = l, ymax = h, x=parameter, fill="band"), alpha = .3) +
     scale_colour_manual("",values="blue")+
    scale_fill_manual("",values="grey12") +
     theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA)) +
     theme(axis.line = element_line(color="black", size = 0.5)) +
     labs(x = "") + labs(y="Rt") + lims(y = c(0, 4.5)) +
     geom_vline(aes(xintercept=d$parameter[44]), color="#999999", linetype="dashed") +
     geom_vline(aes(xintercept=d$parameter[37]), color="#999999", linetype="dashed")
 q
dev.off()


