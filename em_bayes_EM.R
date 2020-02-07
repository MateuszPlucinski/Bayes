## implement E-M algorithm to estimate means and standard deviation in cases and controls

setwd("C://Files") ## set working directory

data_observed = read.csv("finaldata_observed.csv") # read in data


nmax_runs = 1000 # maximum number of iteratios
epsilon = 10e-5 # convergence criterion
n_samples = dim(data_observed)[1] # number of samples (data points)

parameter_matrix = matrix(NA,4,nmax_runs) # matrix that stores parameter estimates (4 parameters X nmax_runs iterations)
state_matrix = matrix(NA,n_samples,nmax_runs) # matrix that stores expectation of variables (n_samples X nmax_runs iterations)

### seed initial guesses for parameters

mu0_old = 8
sd0_old = 1 
mu1_old = 10
sd1_old = 1
parameter_matrix[,1] = c(mu0_old,sd0_old,mu1_old,sd1_old)

### seed initial guesses for hidden observations

state_old = sample(c(0,1),n_samples,replace=TRUE)

state_matrix[,1] = state_old 

# for points that are observed (ie not hidden), use observed values
state_matrix[data_observed$status==1] = 1
state_matrix[data_observed$status==0] = 0

for (i in 2:nmax_runs) {

	### maximization step 
	### estimate mu0, sd0, mu1, sd1 using closed form solutions for mean and standard deviation of normal distributions
	state_old=state_matrix[,i-1]

	mu0_new = sum(data_observed$measure*(1-state_old)/sum(1-state_old))
	sd0_new = sqrt(sum((data_observed$measure-mu0_new)^2*(1-state_old)/sum(1-state_old)))
	mu1_new = sum(data_observed$measure*state_old/sum(state_old))
	sd1_new = sqrt(sum((data_observed$measure-mu1_new)^2*(state_old)/sum(state_old)))

	### expectation
	### calculate expected value (ie prob of case) for each hidden variable

	state_new=state_matrix[,i-1]
	state_new[is.na(data_observed$status)] = dnorm(data_observed$measure[is.na(data_observed$status)],mu1_new,sd1_new)/
			   (dnorm(data_observed$measure[is.na(data_observed$status)],mu0_new,sd0_new)+
				dnorm(data_observed$measure[is.na(data_observed$status)],mu1_new,sd1_new))

	### record new iteration
	parameter_matrix[,i] = c(mu0_new,sd0_new,mu1_new,sd1_new)
	state_matrix[,i] = state_new

	### check convergence

	if (abs(parameter_matrix[1,i] -parameter_matrix[1,i-1])< epsilon) {
		break
	}
}

parameter_matrix = parameter_matrix[,1:i]
state_matrix = state_matrix[,1:i]


windows(16,8)
par(mfrow=c(1,2))
### plot trajectory through parameter space
plot(parameter_matrix[1,],parameter_matrix[3,],typ="l",xlab="mu_0",ylab="mu_1",xlim=c(5,20),ylim=c(1,15))
points(parameter_matrix[1,],parameter_matrix[3,],pch=16)
points(parameter_matrix[1,1],parameter_matrix[3,1],pch=16,col="blue")
points(parameter_matrix[1,i],parameter_matrix[3,i],pch=16,col="red")
legend("topright",fill=c("blue","red"),legend=c("Start","End"))
### plot expectation of hidden observations over time

image(t(state_matrix),xaxt="n",yaxt="n",col = hcl.colors(100, "YlOrRd", rev = TRUE),xlab="Iteration",ylab="Sample ID")
axis(1,at=seq(from=0,to=1,length.out=i),lab=1:i)
axis(2,at=seq(from=0,to=1,length.out=n_samples),lab=1:n_samples)


windows(16,8)
par(mfrow=c(1,2))
# plot final expectation of hidden observations 
state_matrix_hidden = state_matrix
state_matrix_hidden[!is.na(data_observed$status),] = NA 
plot(1:n_samples,state_matrix_hidden[,i],pch=19,xlab="Sample",ylab="Prob Case",col=rgb(0,0,0,0.8))

# plot final expectation of hidden observations, sorting by measure

sortorder = order(data_observed$measure)
plot(data_observed$measure[sortorder],state_matrix_hidden[sortorder,i],pch=19,xlab="Measure",ylab="Prob Case",col=rgb(0,0,0,0.8))
