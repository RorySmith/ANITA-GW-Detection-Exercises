import numpy as np


def htilde_of_f(m1,m2,D,flow,fhigh,deltaF):
    
    ## see https://arxiv.org/pdf/0907.0700.pdf ##
    fseries = np.linspace(flow,fhigh,int((fhigh-flow)/deltaF)+1)
    nu=(m1*m2)/((m1+m2)**2.)
    
    M = m1+m2
    mc = (m1*m2)**(5./3.) / ((m1+m2)**(1./5.))
    piM = M * 4.92549102554e-06 ## convert solar masses to seconds
    v2 = np.cbrt(piM*fseries)**2.
    v5 = np.cbrt(piM*fseries)**5.
    D *= 1e6*3.08567758149e+16 ## convert to MPC in SI units

    gwphase =  3./(128. * nu* v5) * ( (20./9.)*(743./336.+11./4.*nu)*v2   )
    amp = mc**(5./6) / D * fseries**(-7./6.)
    gwphase[0]=0
    amp[0]=0
    deltaF = fseries[1]-fseries[0]
    zeros = np.zeros(int(fseries[0]/deltaF))
    return np.hstack((zeros,amp*np.exp(-1j*gwphase)*(1+1j)))

