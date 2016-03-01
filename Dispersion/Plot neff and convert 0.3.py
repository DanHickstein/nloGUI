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
    group_index : 
        the group index. This doesn't need to be cropped, but it is anyway...
    GVD : array of floats
        the group velocity dispersion (D, not beta2) in ps/nm/k    
    """

    c = 3e-7 # in km/ps
    
    delta_l = wavelengths[1]-wavelengths[0]
    group_index = RI - wavelengths*np.gradient(RI, delta_l)
    d2ndl2 = np.gradient( np.gradient(RI, delta_l), delta_l)
    D = -wavelengths/c * d2ndl2

    return wavelengths[2:-2], group_index[2:-2], D[2:-2]


# thicknesses = (350,500,600,660,950,1050)
thicknesses = (500,600,660,950)

fig, axs = plt.subplots(len(thicknesses),3,figsize=(15,11))

# loop over each row of axes and each thickness
for axrow, thick in zip(axs, thicknesses):
    
    # iterate over each width 
    # for filename in glob.glob(folder+'/t_%inm/*.mat'%thick)     : # Plot them all
    for filename in glob.glob(folder+'/t_%inm/*.mat'%thick)[::2]: # skipping every second width
        
        width = int(re.findall(r'\d+', filename)[1])
        print filename, width
        
        data = scipy.io.loadmat(filename)
        
        # keys: ['param_vary', '__globals__', 'n_eff_aux', 'vg_aux', 'gamma_aux', '__header__', '__version__']
        waves = data['param_vary'][0]
        neff  = data['n_eff_aux'][:,4]
        neff  = np.abs(neff)
        
        waves = waves*1e9 # convert meters to nm
        
        # calculate GVD and group index:
        waves2, group_index, GVD = calcGVD(waves,neff)
        
        axrow[0].plot(waves, neff, label='w=%i'%width)
        axrow[1].plot(waves2, group_index, label='w=%i'%width)
        axrow[2].plot(waves2,GVD, label='w=%i'%width)
        
        axrow[0].set_ylim(1.4, 2.1)
        axrow[1].set_ylim(2.0, 2.2)
        axrow[2].set_ylim(-400, 400)
    
    # make the plots pretty:
    for ax in axrow:
        ax.legend(fontsize=10, frameon=False)
        ax.text(0.5, 0.96, 'height: %i nm'%thick, weight='bold', fontsize=11, 
                transform=ax.transAxes, ha='center', va='top')

        ax.set_xlim(800,3050)
        ax.axvline(1550,alpha=0.4,color='r')
        ax.axvline(1064,alpha=0.4,color='purple')
        ax.axhline(0,alpha=0.4, color='k')

axs[0,0].set_title('Effective index', weight='bold')
axs[0,1].set_title('Group index'    , weight='bold')
axs[0,2].set_title('D (ps/nm/km)'   , weight='bold')

plt.tight_layout() # make the axes closer together
plt.savefig('Kartik dispersion.png',dpi=100)
plt.show()