import pylab
import numpy as np
from numpy import pi
import time
import scipy as sp
from scipy.signal import find_peaks
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 
from scipy import interpolate
pylab.rcParams["font.family"] = "Times New Roman"

array = []

with open('AGN_K.txt', 'r') as f:       
   for line in f:
      array.append((line.split()))
AGN_K = np.asarray(array)
AGN_K = AGN_K[:,0:2].astype(float)
#NGC4051_K = NGC4051_K[:,3:-1].astype(float)


array = []

with open('AGN_V_Sim.txt', 'r') as f:       
   for line in f:
      array.append((line.split()))
AGN_VS = np.asarray(array)
AGN_VS = AGN_VS.astype(float)




'''
 In the data, the first column is Julian Date (JD) - 2450000, the second column
 is flux in mJy (“milli-Jansky”), and the third column is the error in the flux.
 For the simulated data, the first column is JD - 2450000, and column 2 to 11 
 are the fluxes for 10 interpolations. fc is the interpolated flux column currently being correlated.
'''


    
def crosscorrelate(a1,a2,lrange,step,fc):
    """introduces lag to a data set and cross-correlates that lagged set with another. could
    maybe lag in opposite directions to save time?"""
    lags = np.arange(-lrange,+lrange,step) #create range of lags
    ccc = np.zeros(len(lags)) #create input for ccc
    a1old = a1
    for i in range(0,len(lags)): #loop for each            
        a3 = a2 #could just replace a3 in line below with a2 and save the hassle of creating a whole new array
        a3[:,0] = a2[:,0] + lags[i] #lag array. for some reason after this a2 = a2+lags? must be a syntax thing
        a1,a3 = matcharr2(a1,a3,fc,0)
        #a1,a3 = subspline(a1),subspline(a3)
        a1,a3 = matcharr2(a1,a3,1,1)
        #a1,a3=submu(a1),submu(a3)
        ccc[i] = correlatecof(a1[:,1],a3[:,1]) #calculates ccc at this lag value
        a2[:,0] = a2[:,0] - lags[i] #set a2 back to original value
        a1 = a1old #set a1 back to original value
    return ccc,lags
    
def correlatecof(y1,y2):
    """calculates cross-correlation coefficient for two sets of data. see readme for details"""
    mu1 = np.mean(y1)
    sd1 = np.std(y1)
    mu2 = np.mean(y2)
    sd2 = np.std(y2)
    y1 = y1 - mu1
    y2 = y2 - mu2
    dn = sd1*sd2
    N = float(len(y1))
    dn = 1.0 / (dn*N)
    ccc = np.dot(y1,y2)*dn
    return ccc
    
    
def matcharr(y1,y2,fc):
    """matches two data sets by JD and discards non-matching data."""
    y3 = np.zeros((len(y1[:,0]),2)) #creates empty array the length of y1
    y2 = np.column_stack((y2[:,0],y2[:,fc])) #get JD of y2 and relevant flux column.
    y1d = []    #which rows to delete from y1
    y3d = []    #which rows to delete from y3
    for j in range(0,len(y1[:,0])):          #loop through y1
        match = 0                              #no match yet
        for k in range(0,len(y2[:,0])):     #loop through y2
            if int(y1[j,0]) == int(y2[k,0]):  #do the times match?
                y3[j,:] = y2[k,:]               #add to y3
                match = 1                      #match found
        if match == 0:                      #no match found
            y1d.append(j)                   #note rows to be
            y3d.append(j)                   #deleted
    y1 = np.delete(y1,y1d,0)            #and delete at       
    y3 = np.delete(y3,y3d,0)            #the end
    """an issue with this function is that data sets with a fixed offset, e.g. x=10,20,30 and y=11,21,31 will never be
    matched. it may be useful to introduce a parameter that looks for matches within a certain range and assign the
    most suitable match."""
    return y1,y3

def matcharr2(y1,y2,fc,test):
    '''matches two data sets by JD using linear interpolation.'''
    y2 = np.column_stack((y2[:,0],y2[:,fc])) #get JD of y2 and relevant flux column
    '''we need to truncate y1 to fit within the range of JD in y2'''
    '''we should find the first jd of y1 in the range'''
    a = np.where(y1[:,0]<=y2[0,0])
    y1 = np.delete(y1,a,0)
    '''and the last'''
    a = np.where(y1[:,0]>=y2[-1,0])
    y1 = np.delete(y1,a,0)
    y3 = np.zeros((len(y1[:,0]),2))
    for j in range(0,len(y1[:,0])): #loop through JD of y1
        idx = int(np.argmin(np.abs(y2[:,0]-y1[j,0]))) #find the jd in y2 closest to the jd in y1
        sign = int(np.sign(y2[idx,0]-y1[j,0])) #sign = -1 if jd2 is lower, +1 if higher     
        '''if sign = 0 both values are at the same JD and we don't need to interpolate'''
        if sign == 0.0:
            y3[j,:] = y2[j,:]
        else: #linearly interpolate
            xl = min(y2[idx-sign,0],y2[idx,0]) #lower x value
            xh = max(y2[idx-sign,0],y2[idx,0]) #higher x value
            yl = y2[np.argmin(np.abs(y2[:,0]-xl)),1]
            yh = y2[np.argmin(np.abs(y2[:,0]-xh)),1] #corresponding y values
            grad = (yh-yl)/(xh-xl)                                                      #XH-XL ARE ZERO - WHY????
            yint = yl - (grad*xl)
            yj = grad*y1[j,0] + yint #linearly interpolated flux  #YINT AND GRAD ARE PROBLEMS
            y3[j,0] = y1[j,0]
            y3[j,1] = yj
    return y1,y3

def subspline(x,epc=90):
    '''generates a spline from the averages at different epochs.
        then subtracts the spline from the dataset.'''
    '''x values are in JD, y are corresponding flux'''
    '''need to divide data into epochs - say 90JD break is enough'''
    h = [0] #the epoch preceding a break in observation;the last epoch in an observing period
    xmean = []
    ymean = []    
    for i in range(0,len(x[:,0])-1):
        if (x[i+1,0]-x[i,0]) >= epc:
            h.append(i+1)
    for j in range(0,len(h)-1):
        xtemp = x[h[j]:h[j+1],:]
        ymean.append(np.mean(xtemp[:,1])) #mean flux
        xmean.append(xtemp[np.argmin(np.abs(xtemp[:,1]-ymean[-1])),0])
    if len(h)<4:
        print("spline gen error!")
        return x
    #yspline = sp.interpolate.spline(xmean,ymean,x[:,0],order=3) #spline is deprecated
    
    yspline = sp.interpolate.CubicSpline(xmean,ymean)
    nuy = yspline(x[:,0])

    x[:,1] = np.subtract(x[:,1],nuy)
    return x
    
def submu(x,epc=90):
    '''subtracts the mean value of the flux from each value in an epoch.'''
    h = [0] #the epoch preceding a break in observation;the last epoch in an observing period    
    for i in range(0,len(x[:,0])-1):
        if (x[i+1,0]-x[i,0]) >= epc:
            h.append(i+1)
    for j in range(0,len(h)-1):
        ymean = np.mean(x[h[j]:h[j+1],1])
        x[h[j]:h[j+1],1] = np.subtract(x[h[j]:h[j+1],1],ymean)
    return x

def exlag(x,y,lim=180):
    '''extracts the lag from a 1d dataset x.'''
    for i in np.arange(0,1,0.05):
        pks = find_peaks(x,prominence=i)
        pks = pks[0]    #indices of peaks
        pkv = y[pks]
        pkv = pkv[ (pkv >= 0) & (pkv <= lim) ]
        if len(pkv) == 1:
            peak = pkv[0]
            pc = pks[0]
            return pc,peak
    print('unable to retrieve peak!')
    pkv = 0
    peak = 0
    return peak, pkv

def exlag2(x,y,lim=180):
    '''extracts the lag from a 1d dataset x.'''
    
    posx = np.asarray(x[(y > 0) & (y <= lim)]) #values of ccc for positive lag
    posy = np.asarray(y[(y > 0) & (y <= lim)]) #corresponding lag values
    
    pkxm = posx[np.argmax(posx)] #maximum ccc
    #print(pkxm)
    pkym = posy[np.argmax(posx)] #corresponding x
    
    for i in np.arange(0,1,0.05):
        pks = find_peaks(x,prominence=i)
        pks = pks[0]    #indices of peaks
        
        pky = y[pks] #peak values for lag
        pkx = x[pks]
        
        pk = pks[(pky > 0) & (pky <= lim)] #look for peaks in the given range
        
        if len(pk) == 1: #best peak selected
            if y[pk[0]] <= pkym:
                return y[pk[0]], x[pk[0]]
            else:
                return pkym, pkxm
            
    return pkym, pkxm

rangel = 500 #RANGE OF LAGS


pks = []
pkc = []


for i in range(1,11):
    ccc,lags = crosscorrelate(AGN_K,AGN_VS,rangel,1,i)
    pylab.plot(lags,ccc)
    a,b = exlag2(ccc,lags) #a is lag, b is the ccc
    pks.append(b) 
    pkc.append(a)
    

fn = 'AGN.png'     
pylab.axhline(color='black')
pylab.xlabel('Lag (JD)',size=15)
pylab.ylabel('Normalized CCC',size=15)
pylab.ylim(-1.0,1.0)
pylab.title('AGN',size=20,weight='bold')
#pylab.text(-500,0.6,'t = '+str(pk2)+'\nsd = '+str(std1)[:4],size=15)
pylab.grid(True)
#pylab.savefig(fn,dpi=1000)
pylab.show()

pks = [float("%.3g" % x) for x in pks]
pks = np.asarray(pks)
pkc = np.asarray(pkc)

Q3 = np.percentile(pkc,75)
Q1 = np.percentile(pkc,25)
IQR = Q3-Q1
Lowerlimit = Q1 - 1.5*IQR
Upperlimit = Q3 + 1.5*IQR

pks = pks[(pkc > Lowerlimit) & (pkc < Upperlimit)]
pkc = pkc[(pkc > Lowerlimit) & (pkc < Upperlimit)]

