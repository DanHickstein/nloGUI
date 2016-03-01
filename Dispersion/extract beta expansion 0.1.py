import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import re
import glob


folder = 'data_files_for_supercontinuum'


def get_beta_expansion(l, D, l0=1550, polyOrder=4, return_everything=True):
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

def find_central_wavelength(l, D, low=1000, high=3000, polyOrder=4):
    
    guesses = np.linspace(low,high, 200)
    errors  = np.zeros(np.shape(guesses))
    
    for index, l0 in enumerate(guesses):
        beta_coefs, w, beta2, fit = get_beta_expansion(waves2, GVD, l0=l0)
        errors[index] = np.median(np.abs(beta2-fit))
    
    best_guess = guesses[np.argmin(errors)]
    beta_coefs, w, beta2, fit = get_beta_expansion(waves2, GVD, l0=best_guess)

    print 'Best guess: %.2f nm, B2: %.2f, B3: %.4f'%(best_guess, beta_coefs[0], beta_coefs[1])
    
    return beta_coefs, w, beta2, fit
    
    





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

fig, axs = plt.subplots(len(thicknesses),2,figsize=(12,11))

# loop over each row of axes and each thickness
for axrow, thick in zip(axs, thicknesses):
    
    # iterate over each width 
    # for filename in glob.glob(folder+'/t_%inm/*.mat'%thick)     : # Plot them all
    for filename in glob.glob(folder+'/t_%inm/*.mat'%thick)[::3]: # skipping every second width
        
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

        axrow[0].plot(waves2,GVD, label='w=%i'%width)
        axrow[0].set_ylim(-400, 400)
        
        beta_coefs, w, beta2, fit = find_central_wavelength(waves2, GVD)
        
        l1, = axrow[1].plot(w, beta2, label='w=%i'%width)
        color = l1.get_color()
        axrow[1].plot(w,fit,color=color, ls='dashed')
        
    
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
plt.show()