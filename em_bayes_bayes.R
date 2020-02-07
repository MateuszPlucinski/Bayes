## implement MCMC Bayesian algorithm to sample means and standard deviation in cases and controls

setwd("C://Files") ## set working directory

library(invgamma)

data_observed = read.csv("finaldata_observed.csv") # read in data


nmax_runs = 10000 # maximum number of iterations
burnin = 1000 # how many iterations to discard from beginning of markov chain
n_samples = dim(data_observed)[1] # number of samples (data points)

parameter_matrix = matrix(NA,4,nmax_runs) # matrix that stores parameter estimates (4 parameters X nmax_runs iterations)
state_matrix = matrix(NA,n_samples,nmax_runs) # matrix that stores state of hidden variables (n_samples X nmax_runs iterations)

### seed initial guesses for parameters

mu0_old = 11
sd0_old = 3 
mu1_old = 10 
sd1_old = 3
parameter_matrix[,1] = c(mu0_old,sd0_old,mu1_old,sd1_old)

## priors for means 
mu0_0 = 11
sigma0_0 = 100
mu0_1 = 10
sigma0_1 = 100

## priors for standard deviations

alpha_0 = 0.001
beta_0 = 0.001
alpha_1 = 0.001
beta_1 = 0.001

### randomly seed hidden observations with 0 or 1
state_old = sample(c(0,1),n_samples,replace=TRUE)
state_matrix[,1] = state_old 
state_matrix[data_observed$status==1,1] = 1
state_matrix[data_observed$status==0,1] = 0

for (i in 2:nmax_runs) {

	### load previous parameters
	state_old=state_matrix[,i-1]
	sd0_old = parameter_matrix[2,i-1]
	sd1_old = parameter_matrix[4,i-1]
	mu0_old = parameter_matrix[1,i-1]
	mu1_old = parameter_matrix[3,i-1]

	### sample from posterior sample for new estimates of mu0, sigma0, mu1, sigma1  (using conjugate priors)

	n0 = sum(state_old==0)
	mean0 = mean(data_observed$measure[state_old==0])
	mu0_new = rnorm(1,(1/sigma0_0^2+n0/sd0_old^2)^-1*(mu0_0/sigma0_0^2+mean0*n0/sd0_old^2)
					,sqrt((1/sigma0_0^2+1/sd0_old^2)^-1))
	sumsquares0 = sum((data_observed$measure[state_old==0]-mu0_new)^2)
	sd0_new = sqrt(rinvgamma(1, alpha_0+n0/2, beta_0+sumsquares0/2))

	n1 = sum(state_old==1)
	mean1 = mean(data_observed$measure[state_old==1])
	mu1_new = rnorm(1,(1/sigma0_1^2+n1/sd1_old^2)^-1*(mu0_1/sigma0_1^2+mean1*n1/sd1_old^2)
					,sqrt((1/sigma0_1^2+1/sd1_old^2)^-1))
	sumsquares1 = sum((data_observed$measure[state_old==1]-mu1_new)^2)
	sd1_new = sqrt(rinvgamma(1, alpha_1+n1/2, beta_1+sumsquares1/2))

	if (mu1_new > mu0_new) {  # constraint that mu1 < mu0
		mu0_new = mu0_old
		mu1_new = mu1_old
		sd0_new = sd0_old
		sd1_new = sd1_old
	}

	### sample from hidden variables 

	likelihoodratio = rep(NA,n_samples)
	state_new = state_old

	# calculate likelihood ratio P(D|Hnew) / P(D|Hold)
	likelihoodratio[state_old==0] = dnorm(data_observed$measure[state_old==0],mu1_new,sd1_new) / 
								dnorm(data_observed$measure[state_old==0],mu0_new,sd0_new)
	likelihoodratio[state_old==1] = dnorm(data_observed$measure[state_old==1],mu0_new,sd0_new) / 
								dnorm(data_observed$measure[state_old==1],mu1_new,sd1_new)

	z = runif(n_samples,0,1) # generate n_samples random variables between 0 and 1

	# if z < likelihood ratio, accept new state, reject otherwise
	state_new[z<likelihoodratio] = as.numeric(!(state_old[z<likelihoodratio]))

	# for points that are observed (ie not hidden), use observed values
	state_new[data_observed$status==1] = 1
	state_new[data_observed$status==0] = 0

	parameter_matrix[,i] = c(mu0_new,sd0_new,mu1_new,sd1_new)
	state_matrix[,i] = state_new
}

windows(16,8) # show markov chain for the mu0 parameter (to judge convergence)
plot(1:i,parameter_matrix[1,],typ="l",pch=16,ylab="mu_0",xlab="Iteration")
points(c(-100*i,i*100),c(median(parameter_matrix[1,]),median(parameter_matrix[1,])),typ="l",lty=2,col="red")

# remove burn in samples
parameter_matrix = parameter_matrix[,seq(from=burnin,to=nmax_runs,length.out=nmax_runs/5)]
state_matrix = state_matrix[,seq(from=burnin,to=nmax_runs,length.out=nmax_runs/5)]

windows(16,8)
par(mfrow=c(1,2))
### plot trajectory through parameter space
plot(parameter_matrix[1,],parameter_matrix[3,],typ="l",xlab="mu_0",ylab="mu_1",col=rgb(0,0,0,0.01),xlim=c(5,20),ylim=c(1,15))
points(parameter_matrix[1,],parameter_matrix[3,],pch=16,col=rgb(0,0,0,0.15))
points(median(parameter_matrix[1,]),median(parameter_matrix[3,]),pch=16,col=rgb(1,0,0,0.8),cex=2)
### plot evolution of hidden observations over time
image(t(state_matrix),xaxt="n",yaxt="n",col = hcl.colors(100, "YlOrRd", rev = TRUE),xlab="Iteration",ylab="Sample ID")
axis(1,at=seq(from=0,to=1,length.out=20),lab=burnin+floor(seq(from=1,to=dim(state_matrix)[1]*5,length.out=20)))
axis(2,at=seq(from=0,to=1,length.out=n_samples),lab=1:n_samples)

# create density plots of posterior estimates of parameters

windows()
par(mfrow=c(2,2))
plot(density(parameter_matrix[1,]),xlim=c(0,25),main="mu_0")
plot(density(parameter_matrix[2,]),xlim=c(0,25),main="sd_0")
plot(density(parameter_matrix[3,]),xlim=c(0,25),main="mu_1")
plot(density(parameter_matrix[4,]),xlim=c(0,25),main="sd_1")



quantile(parameter_matrix[1,],c(0.025,0.5,0.975))
quantile(parameter_matrix[2,],c(0.025,0.5,0.975))
quantile(parameter_matrix[3,],c(0.025,0.5,0.975))
quantile(parameter_matrix[4,],c(0.025,0.5,0.975))

rowMeans(parameter_matrix)

windows(16,8)
par(mfrow=c(1,2))

# plot posterior probability of hidden observations 

state_matrix_hidden = state_matrix
state_matrix_hidden[!is.na(data_observed$status),] = NA 
plot(1:n_samples,rowMeans(state_matrix),pch=19,xlab="Sample",ylab="Prob Case",col=rgb(0,0,0,0.8),ylim=c(0,1))

# plot final expectation of hidden observations, sorting by measure

sortorder = order(data_observed$measure)
plot(data_observed$measure[sortorder],rowMeans(state_matrix_hidden[sortorder,]),pch=19,xlab="Measure",ylab="Prob Case",col=rgb(0,0,0,0.8),ylim=c(0,1))

# print out final posterior probability of case
cbind("sample"=1:n_samples,"Prob Case"=rowMeans(state_matrix))

