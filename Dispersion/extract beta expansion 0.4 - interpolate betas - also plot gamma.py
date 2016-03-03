import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import re
import glob
import scipy.interpolate


folder = 'data_files_for_supercontinuum'
polyOrder = 15


def get_beta_expansion(l, D, l0=1550, polyOrder=polyOrder, return_everything=True):
    """ Calculates the coefficients of the beta expansion given D.
    
    Parameters
    ----------
    l : array of floats
        array of wavelengths (in nm)
    D : array of float
        GVD in ps/nm/km (engineering definition, not physics definition)
    l0 : float
        central wavelengths of the Taylor expansion (nm)
        
    Returns
    -------
    
    """
    from scipy.misc import factorial
    
    c = 3e8 
    
    l = l * 1e-9 # convert lambda from nm to m
    l0 = l0 * 1e-9
    
    # Units of D are ps/nm/km
    # Convert to s/m/m 
    D = D * 1e-12 * 1e9 * 1e-3
    
    freq = c/l
    w_offset = 2*np.pi*c/l - 2*np.pi*c/l0
    
    w_offset = w_offset * 1e-12
   
    # Convert from D to beta via  beta2 = -D * lambda^2 / (2*pi*c) 
    beta2 = - D* l**2 / (2*np.pi*c)
    beta2 = beta2 * 1e24 * 1e3 # convert to ps^2 / km
    
    coefs = np.polyfit(w_offset, beta2, polyOrder)
    
    beta_coefs = coefs[::-1]
    
    polyFit = np.zeros((len(w_offset),))   

    for i in range(len(beta_coefs)):
        beta_coefs[i] = beta_coefs[i] * factorial(i)
        polyFit = polyFit + beta_coefs[i] / factorial(i)*w_offset**i
    
    if return_everything:
        return beta_coefs, w_offset, beta2, polyFit
    
def get_betas_from_dimensions(thickness, width, polyOrder=polyOrder):
    data = np.load('beta_expansion_order=%i.npz'%polyOrder)
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

def expand_betas(beta_coefs, w_offset):
    from scipy.misc import factorial
    
    polyFit = np.zeros((len(w_offset),))   
    
    for i in range(len(beta_coefs)):
        polyFit = polyFit + beta_coefs[i] / factorial(i)*w_offset**i
    
    return polyFit
        


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


thicknesses = (350,500,600,660,950,1050)
# thicknesses = (500,600,660,950)

fig, axs = plt.subplots(len(thicknesses),3,figsize=(16,10))

thick1D = []
width1D = []
b_list  = []

# loop over each row of axes and each thickness
for axrow, thick in zip(axs, thicknesses):
    
    # iterate over each width 
    for filename in glob.glob(folder+'/t_%inm/*.mat'%thick): 
        
        width = int(re.findall(r'\d+', filename)[1])
        # print filename, width
        
        data = scipy.io.loadmat(filename)
        
        # keys: ['param_vary', '__globals__', 'n_eff_aux', 'vg_aux', 'gamma_aux', '__header__', '__version__']
        waves = data['param_vary'][0]
        neff  = data['n_eff_aux'][:,4]
        gamma = data['gamma_aux'][:,4]
        neff  = np.abs(neff)
        
        waves = waves*1e9 # convert meters to nm
        
        # calculate GVD and group index:
        waves2, group_index, GVD = calcGVD(waves,neff)

        axrow[0].plot(waves2,GVD, label='w=%i'%width)
        axrow[0].set_ylim(-400, 400)
        
        beta_coefs, w, beta2, fit = get_beta_expansion(waves2, GVD)
        print 'thick:%i, width:%i' %(thick, width),
        for num, b in enumerate(beta_coefs):
            print 'beta%i:%.1e\t'%(num+2,b),
        print ''
        
        l1, = axrow[1].plot(w, beta2, label='w=%i'%width)
        color = l1.get_color()
        
        
        
        try:
            bcs = get_betas_from_dimensions(thick, width)
            fit = expand_betas(bcs, w)
        except:
            print 'oh no'
                
        axrow[1].plot(w,fit,color=color, ls='dashed')
        
        thick1D.append(thick)
        width1D.append(width)
        b_list.append(beta_coefs)
        
        axrow[2].plot(waves,gamma)
        
    
    # make the plots pretty:
    for ax in axrow:
        ax.legend(fontsize=10, frameon=False)
        ax.text(0.5, 0.96, 'height: %i nm'%thick, weight='bold', fontsize=11, 
                transform=ax.transAxes, ha='center', va='top')

        ax.axvline(1550,alpha=0.4,color='r')
        ax.axvline(1064,alpha=0.4,color='purple')
        ax.axhline(0,alpha=0.4, color='k')
        





for ax in axs[:,0]:
    ax.set_xlim(800,3050)

axs[0,0].set_title('D (ps/nm/km)', weight='bold')
# axs[0,1].set_title('Group index'    , weight='bold')
# axs[0,2].set_title('D (ps/nm/km)'   , weight='bold')

plt.tight_layout() # make the axes closer together
# plt.savefig('convert to beta.png',dpi=100)



thick1D = np.array(thick1D)
width1D = np.array(width1D)
b_list  = np.array(b_list)



np.savez('beta_expansion_order=%i.npz'%polyOrder, t=thick1D, w=width1D, b=b_list)

print thick1D[0],width1D[0]
# print f(350.,400.), b_list[0,0]






plt.show()