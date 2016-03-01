import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import re
import glob

folder = 'data_files_for_supercontinuum'

def calcGVD(wavelengths,RI):
    """ Calculates the GVD
    
    Parameters
    ----------
    wavelengths : array of float
        array of wavelengths (in nm)
        NOTE: This array must be evenly spaced!
    RI  : array of float
        array of refractive indices
        
    Returns
    -------
    waves : array of floats
        the array of wavelengths cropped to omit the 
        first two and last two values.
        the array is cropped because the second 
        derivative is not well defined at the edges.
    GVD : array of floats
        the group velocity dispersion (D, not beta2) in ps/nm/k    
    """

    c = 3e-7 # in km/ps
    
    # delta_l = np.mean(np.gradient(wavelengths))
    delta_l = wavelengths[1]-wavelengths[0]
    
    d2ndl2 = np.gradient( np.gradient(RI, delta_l), delta_l)
    D = -wavelengths/c * d2ndl2

    
    return wavelengths[2:-2], D[2:-2]




for thick in (350,500,600,660):#,950,1050):
    for filename in glob.glob(folder+'/t_%inm/*.mat'%thick):
        width = int(re.findall(r'\d+', filename)[1])
        print filename, width
        
        data = scipy.io.loadmat(filename)
        
        # ['param_vary', '__globals__', 'n_eff_aux', 'vg_aux', 'gamma_aux', '__header__', '__version__']
        waves = data['param_vary'][0]
        neff  = data['n_eff_aux'][:,4]
        
        # n_avg = np.mean(neff,axis=0)
        # print np.argmax(n_avg)
        waves = waves*1e9 # convert to nm
        

        waves, GVD = calcGVD(waves,neff)
        print thick,width
        if thick == 500 and width == 1000:
            print 'found it!'
            plt.plot(waves, GVD)
        

plt.show()        
        


# data = scipy.io.loadmat(filename)
#
# data2 = data['gvd']

# waves   = np.linspace(800,2500,40)[1:-1]
# heights = (0.3, 0.5, 0.6, 0.7)
# widths  = np.arange(0.8, 3, 0.2)
#
# print np.shape(data2)
#
# for height, datarow in zip(heights,data2.transpose()):
#     if height != 0.6: continue
#     for width, D in zip(widths, datarow):
#         alpha = 1
#         D = D[:,0]
#
#
#         relative_grad = np.gradient(D[waves>1200])/D[waves>1200]
#         if np.max(np.abs(relative_grad)) > 5:
#             ls = 'dotted'
#             alpha = 0.2
#         else:
#             ls = 'solid'
#             alpha = 1.0
#
#             print np.shape(waves), np.shape(D)
#
#         plt.plot(waves, D, alpha=alpha, ls=ls)
#
#
# plt.show()
