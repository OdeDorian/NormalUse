import numpy as np
import emcee
import scipy
from scipy.optimize import minimize
from scipy import special


class emo:

    nominator=np.sqrt(2) #*sigma
    logsqrt2pi=0.5*np.log(2*np.pi)
    agmin= 0 #np.min(dataAG)
    agmax= 10 #np.max(dataAG)
    #AGErrorSquare=AGError*AGError
    sqrtPi2i=2./np.sqrt(2*np.pi) 


    logsqrtPi2i=np.log(sqrtPi2i)



    def calProb(self,theta,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars):
        disCloud,  mu1,imu1sigma,mu2,imu2sigma= theta
        ##calculate the probablitydisSigma
        #trueError=np.sqrt(disError**2+disSigma**2)
        #trueError= #np.sqrt(disError**2+disSigma**2)

        weight=scipy.stats.norm.cdf(disCloud,dataDis,disError)
        
        #weight=1-scipy.stats.norm.cdf(dataDis,disCloud, trueError)
     
        #fore ground
        errorVar= 1./imu1sigma**2+AGErrorSquare
        errorVar_i=1./errorVar
        convolveSTDi=  np.sqrt(errorVar_i)
        a=  special.erf((self.agmax-mu1)/self.nominator*convolveSTDi)+ special.erf((mu1-self.agmin)/self.nominator*convolveSTDi)  
        
        #common factor,  np.exp(-0.5*(dataAG-mu1)**2/errorVar)*convolveSTDi*sqrtPi2i/a
        
        #pForeground= np.exp(-0.5*(dataAG-mu1)**2/errorVar)*convolveSTDi*sqrtPi2i/a
        pForeground= np.exp(-0.5*np.square(dataAG-mu1)*errorVar_i)*convolveSTDi

        
        #pForeground=     np.sum(-0.5*(dataAG-mu1)**2/errorVar - np.log(a)-0.5*np.log( errorVar)  )
     
        #fore ground
        errorVar2= 1./imu2sigma**2 +AGErrorSquare
        errorVar2_i=1./errorVar2
        convolveSTDi2=  np.sqrt(errorVar2_i)
        a2=  special.erf((self.agmax-mu2)/self.nominator*convolveSTDi2)+ special.erf((mu2-self.agmin)/self.nominator*convolveSTDi2)  
        
        
        pBackground= np.exp(-0.5*np.square(dataAG-mu2)*errorVar2_i)*convolveSTDi2
        
        totalP= pBackground+(pForeground-pBackground)*weight #true norm
        
        return np.sum(np.log(totalP) )+Nstars*self.logsqrtPi2i


    def logprior(self,theta):
        disCloud,  mu1,imu1sigma,mu2,imu2sigma= theta
        lp = 0
        
        lp = 0. if self.dismin < disCloud < self.dismax  and self.agmin<mu1<self.agmax and self.agmin<mu2<self.agmax and 0 < imu1sigma  and 0 < imu2sigma else -np.inf

        lp -= self.mu1c*mu1  +0
        lp -= self.mu2c*mu2 +0
        lp -= self.sigma1c*imu1sigma +0
        lp -= self.sigma2c*imu2sigma +0

        return lp

    def logposterior(self,theta,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars):
        lp=self.logprior(theta)

        if not np.isfinite(lp):
            return -np.inf

        return lp + self.calProb(theta,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars)



    def mcmc_sample(self,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars,priors,Nensembles=32,Nburns=1000,Nsamples=10000,thin=50):

        self.dismin,self.dismax,self.mu1c,self.mu2c,self.sigma1c,self.sigma2c=priors

        disinit=np.random.uniform(self.dismin,self.dismax,Nensembles)

#        mu1init=np.random.exponential(self.mu1c,Nensembles)
#        mu2init=np.random.exponential(self.mu2c,Nensembles)
#        imu1sigmainit=np.random.exponential(self.sigma1c,Nensembles)
#        imu2sigmainit=np.random.exponential(self.sigma2c,Nensembles)

        mu1init=self.mu1c+np.random.randn(Nensembles)*0.001
        mu2init=self.mu2c+np.random.randn(Nensembles)*0.001
        imu1sigmainit=self.sigma1c+np.random.randn(Nensembles)*0.001
        imu2sigmainit=self.sigma2c+np.random.randn(Nensembles)*0.001

        initsample=np.array([disinit,mu1init,imu1sigmainit,mu2init,imu2sigmainit]).T

        ndims=initsample.shape[1]

        argslist = (dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars)

        sampler = emcee.EnsembleSampler(Nensembles, ndims, self.logposterior, args=argslist)

        et=0
        print('Preparing MCMC sampling >>>>>>>')
  #      while et==0:
        try:
            sampler.run_mcmc(initsample, Nsamples+Nburns, progress=True, store=True)
            et+=1
            print(et)
        except:
            et=0

       # postsample = sampler.chain[:, Nburns:, :].reshape((-1, ndims))
        postsample = sampler.get_chain(discard=Nburns, thin=thin, flat=True)

        return postsample
