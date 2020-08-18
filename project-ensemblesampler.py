import numpy as np
import scipy as sp
import emcee
import pylab
from math import log
pylab.rcParams["font.family"] = "DejaVu Sans"
import corner

def log_likelihood(theta, x, y, yerr):
    m, b, log_f = theta
    model = m * x + b
    sigma2 = yerr ** 2 + model ** 2 * np.exp(2 * log_f)
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))

def log_prior(theta):
    '''m = gradient, b = intercept, log_f = fractional underestimation.
    function determining likely ranges of variables.'''
    m, b, log_f = theta
    if 0.0 < m < 1.5 and -15.0 < b < -5.0 and -20.0 < log_f < 0.0:
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)




def refine_parameters(x,y,yerr,f=0.454,minimize = '1',plotcorner='no',project='no'):
    '''
    uses emcee and returns a 4x3 array. column 1 = m, column 2 = c, column 3 = log_f.
    row 1 = first estimate by numpy polyfit, row 2 = scipy minimized. Can also use scipy minimize to further refine.
    '''
    
    V = np.polyfit(x,y,1) # get initial parameters
    m_initial = V[0]
    c_initial = V[1]
    f_initial = f
    
    log_f_initial = log(f_initial,10)
    
    soln = m_initial, c_initial, log_f_initial
    
    parameters = np.zeros((5,3))
    parameters[0,0] = m_initial
    parameters[0,1] = c_initial
    parameters[0,2] = log_f_initial
    
    
    if minimize == 2:
        from scipy.optimize import minimize
        
        n11 = lambda *args: -log_likelihood(*args)
        initial = np.array([m_initial,c_initial,log_f_initial])
        soln = minimize(n11,initial,args=(x,y,yerr))
        soln = soln.x
        
        parameters[1,0] = soln[0]
        parameters[1,1] = soln[1]
        parameters[1,2] = soln[2]
        
    if minimize ==1:
        np.delete(parameters,1,0) #this isn't working right now?
        
    m_m1, c_m1, log_f_m1 = soln
    
    pos = soln + 1e-4 * np.random.randn(32, 3) #5000
    print(pos.shape)
    nwalkers, ndim = pos.shape
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x, y, yerr))
    sampler.run_mcmc(pos, 5000, progress=True) #10000
    
    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
    
    labels = ["slope", "offset", 'Log'+ r'$_1$'+'$_0$'+u'\u0192']
    
    if plotcorner == 'yes': #plot the corner plot
        fig = corner.corner(flat_samples, labels=labels, truths=[m_m1, c_m1, log_f_m1])
        pylab.savefig('corner.png',dpi=1000)
        pylab.show()
    if project == 'yes':
        inds = np.random.randint(len(flat_samples), size=100)
        
        for ind in inds:
            x0 = np.linspace(min(x),max(x),2)
            sample = flat_samples[ind] #put this outside if statement
            pylab.plot(x0, np.dot(np.vander(x0, 2), sample[:2]), "C1", alpha=0.1)
            pylab.scatter(x[-2],y[-2],s=500,facecolors='none',edgecolors='red')
            pylab.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0)
            pylab.plot(x0, m_m1 * x0 + c_m1, "k", label="truth")
            pylab.xlabel('Log'+r'$_1$'+'$_0$'+' Luminosity',size=15)
            pylab.ylabel('Log'+r'$_1$'+'$_0$'+' Lag',size=15)
            pylab.grid(True)
        pylab.savefig('project.png',dpi=1000)
        pylab.show()
    
    
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:,i],[16,50,84])
        q = np.diff(mcmc)
        parameters[2,i] = mcmc[1]
        parameters[3,i] = -q[0]
        parameters[4,i] = q[1] 
    
    return parameters

#xi - luminosity
#yi - lag
#yerri - lag errors

x = np.log10(xi) #log x
y = np.log10(yi) #log y
yerr = 0.434/yi * yerri #log lag errors
        
p = refine_parameters(x,y,yerr,f=0.1,minimize=2,plotcorner='no',project='no') #change this to loop through different f?
#pylab.savefig('cornerwithakn120.png',dpi=1000)

sx = np.linspace(min(x),max(x),2)


pylab.scatter(x,y,color='black')    
pylab.plot(sx,sx*p[2,0] + p[2,1],linestyle='-.',label = 'ensemble',color='blue')

#pylab.legend(loc=2,prop={'size':15})

pylab.hlines(y=35.5*p[2,0] + p[2,1], xmin=35.5,xmax=36.5,linestyle='--',color='red')
pylab.vlines(x=36.5,ymin=35.5*p[2,0]+p[2,1],ymax=36.5*p[2,0]+p[2,1],linestyle='--',color='red')

pylab.text(35.75,1.3,'dX = 1',fontsize=15)
pylab.text(36.6,1.6,'dY = 0.324',fontsize = 15)

pylab.axhline(color='black')
pylab.xlabel('Log'+r'$_1$'+'$_0$'+' Luminosity',size=15)
pylab.ylabel('Log'+r'$_1$'+'$_0$'+' Lag',size=15)
pylab.ylim(min(y)-yerr[np.argmin(y)]*2)
pylab.grid(True)
#pylab.savefig('ex2.png',dpi=1000)
pylab.show()



        
        

    
