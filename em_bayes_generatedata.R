### generate toy data for inference exercises

setwd("C://Files") ## set working directory

### generate data

mu0_true = 12
sd0_true = 1

mu1_true = 8
sd1_true = 1

n0_true = 200
n1_true = 100

p_observed = 0.3

controls_true = rnorm(n0_true,mu0_true,sd0_true)
cases_true = rnorm(n1_true,mu1_true,sd1_true)

newsortorder = sample(c(1:n0_true,(n0_true+1):(n0_true+n1_true)))

newdata = rep(NA,n0_true+n1_true)
newdata[newsortorder<=n0_true] = controls_true
newdata[newsortorder>n0_true] = cases_true

disease_true = rep(NA,n0_true+n1_true)
disease_true[newsortorder<=n0_true] = 0
disease_true[newsortorder>n0_true] = 1

disease_observed = disease_true
disease_observed[sample(n0_true+n1_true,floor((n0_true+n1_true)*(1-p_observed)))] = NA

hist(newdata)

## true data distribution

hist0 = hist(newdata[disease_true==0],breaks=seq(from=floor(min(newdata)),to=ceiling(max(newdata)),by=0.5))
hist1 = hist(newdata[disease_true==1],breaks=seq(from=floor(min(newdata)),to=ceiling(max(newdata)),by=0.5))
a=barplot(rbind(hist0$count,hist1$count),beside=FALSE,col=c("blue","red"),xlab="Measure",ylab="Frequency")
legend("topright",fill=c("blue","red"),legend=c("Control","Case"))
axis(1,at=a+0.5,lab=hist0$breaks[-1])


## observed data distribution

hist0 = hist(newdata[disease_observed==0],breaks=seq(from=floor(min(newdata)),to=ceiling(max(newdata)),by=0.5))
hist1 = hist(newdata[disease_observed==1],breaks=seq(from=floor(min(newdata)),to=ceiling(max(newdata)),by=0.5))
histNA = hist(newdata[is.na(disease_observed)],breaks=seq(from=floor(min(newdata)),to=ceiling(max(newdata)),by=0.5))
a=barplot(rbind(hist0$count,hist1$count,histNA$count),beside=FALSE,col=c("blue","red","grey"),xlab="Measure",ylab="Frequency")
legend("topright",fill=c("blue","red","grey"),legend=c("Control","Case","Unknown"))
axis(1,at=a+0.5,lab=hist0$breaks[-1])


# write final generated dataset

finaldata_observed = cbind(status=disease_observed,measure = newdata)
finaldata_true = cbind(status=disease_true,measure = newdata)
sd(newdata[disease_true==1])
sd(newdata[disease_true==0])

write.csv(finaldata_observed,"finaldata_observed.csv")
write.csv(finaldata_true,"finaldata_true.csv")


