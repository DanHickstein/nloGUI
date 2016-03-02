import scipy.interpolate
import numpy as np
import os

polyOrder=15

def get_betas_from_dimensions(thickness, width, polyOrder=polyOrder):
    fileprefix = os.path.split(__file__)[0]
    data = np.load(fileprefix+'/beta_expansion_order=%i.npz'%polyOrder)
    t = data['t']
    w = data['w']
    b = data['b']
    
    betas = np.zeros(polyOrder+1)
    
    for i in range(polyOrder+1):
        betas[i] = scipy.interpolate.griddata((t, w), b[:,i], (thickness,width), method='nearest')
        # betas[i] = f(thickness, width)
        
    print 'thick:%i, width:%i' %(thickness, width),
    for num, b in enumerate(betas):
        print 'beta%i:%.1e\t'%(num+2,b),
    print ''
    
    return betas


def main():
    b = get_betas_from_dimensions(600, 1500)
    print b
    
if __name__=='__main__':
    main()