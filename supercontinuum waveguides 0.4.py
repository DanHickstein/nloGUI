import numpy as np
import matplotlib.pyplot as plt
from pynlo.interactions.FourWaveMixing import SSFM
from pynlo.media.fibers import fiber
from pynlo.light.DerivedPulses import SechPulse
import time
import scipy.integrate
from autoscale_magic import autoscale_y
from scipy.misc import factorial
from Dispersion.katrik_betas import get_betas_from_dimensions


Pulse   = 0.125  # pulse duration (ps)
pulseWL = 1550    # pulse central wavelength (nm)
EPP     = 40e-12     # Pulse energy (J)
GDD     = 0.0     # Group delay dispersion (ps^2)
TOD     = 0.0     # Third order dispersion (ps^3)

Window  = 10.0  # simulation window (ps)
Steps   = 50   # simulation steps
Points  = 2**13  # simulation points

Length  = 30.0    # fiber length (mm)
Alpha   = 0.4  # attentuation coefficient
Gamma   = 3000    # Gamma (1/(W km) -- 1400 is for Silicon Nitride

# thick   = 600
# widths  = (1200,1600,2000,2400)

thick   = 660
widths  = (1700,2100,2500,2900)

fibWL   = 1550   # Center WL (nm)

Raman = True 
Steep = True

def dB(num):
    return 10 * np.log10(np.abs(num)**2)


fig0, axs0 = plt.subplots(3,4, figsize=(13, 9), sharex=True, sharey='row') # this one is for the frequency domain
fig1, axs1 = plt.subplots(2,4, figsize=(13, 7), sharex=True, sharey='row') # this one is for the time domain
fig2, axs2 = plt.subplots(1,2, figsize=(14, 5))   # this is for 1D plot

axs0b = [ax.twinx() for ax in axs0[-1]] # set up the twin axes
for ax in axs0b[:-1]:
    for label in ax.get_yticklabels(): 
        label.set_visible(False)
        
for label in axs0b[-1].get_yticklabels():
    label.set_fontsize(8)

for ax in axs0[-1]: # axs0 bottom row
    ax.set_xlabel('Frequency (THz)')
    for label in ax.get_xticklabels():
        label.set_fontsize(8)

for ax in axs1[-1]: # axs1 bottom row
    ax.set_xlabel('Time (ps)')

for axs in (axs0, axs1):
    axs[0,0].set_ylabel('Intensity (dB)')
    axs[1,0].set_ylabel('propagation length (mm)')


axs2[0].set_ylabel('Intensity (dB)')
axs2[0].set_xlabel('Frequency (THz)')
axs2[1].set_xlabel('Time (ps)')

axs0[2,0].set_ylabel(r'$\beta_2$ (ps$^2$/km)', color='g')
axs0b[-1].set_ylabel(  'Integrated dispersion (1/m)', color='m')

axs0b[0].set_xlabel('Frequency (THz)')
axs0b[1].set_xlabel('Time (ps)')


for i, width in enumerate(widths):

    fiber_length = Length * 1e-3
    alpha = np.log((10**(Alpha * 0.1))) * 100  # convert from dB/cm to 1/m

    # set up the pulse parameters
    pulse = SechPulse(1, Pulse, pulseWL, time_window_ps=Window,
                      GDD=GDD, TOD=TOD, NPTS=Points, frep_MHz=100, power_is_avg=False)

    pulse.set_epp(EPP)
    
    betas = get_betas_from_dimensions(thick, width)
    
    fiber1 = fiber.FiberInstance()
    fiber1.generate_fiber(fiber_length, center_wl_nm=fibWL, betas=betas,
                                  gamma_W_m=Gamma * 1e-3, gvd_units='ps^n/km', gain=-alpha)
        

    laser_freq = 3e8 / (pulseWL * 1e-9) * 1e-12
    fiber_freq = 3e8 / (fibWL * 1e-9)   * 1e-12

    W = pulse.W_mks
    F = W / (2 * np.pi) * 1e-12 # convert to THz
    T = pulse.T_ps
    xW = W[W > 0]

    D = fiber1.Beta2_to_D(pulse)
    beta = fiber1.Beta2(pulse) * 1e3  # convert ps^2/m to ps^2/km

    # Plot the dispersion in the bottom panels
    axs0[2,i].plot(F[W > 0], beta[W > 0],       color='g')
    axs0b[ i].plot(F, fiber1.get_betas(pulse), color='m')
    
    axs0[2,i].axhline(0,color='g', alpha=0.5)
    axs0b[ i].axhline(0,color='m', alpha=0.5)
    

    # plot *before* in the top panels (blue)
    axs0[0,i].plot(F[W > 0], dB(pulse.AW[W > 0]), color='b') # freq
    axs1[0,i].plot(T,        dB(pulse.AT),        color='b') # time

    # plot some vertical lines
    for ax in (axs0[0,i], axs0[2,i]):
        # ax.axvline(fiber_freq,   color='g', alpha=0.1)
        ax.axvline(laser_freq,   color='b', alpha=0.3)
        ax.axvline(laser_freq*2, color='g', alpha=0.2)

    # Propagation
    t = time.time() # start a timer

    # set up the propagation parameters
    evol = SSFM.SSFM(local_error=0.001, USE_SIMPLE_RAMAN=True,
                     disable_Raman=np.logical_not(Raman), 
                     disable_self_steepening=np.logical_not(Steep))

    # propagate the pulse!
    y, AW, AT, pulse1 = evol.propagate( pulse_in=pulse, 
                                        fiber=fiber1, 
                                        n_steps=Steps)

    zW_in = np.transpose(AW)[:, (W > 0)]
    zT_in = np.transpose(AT)
    zW = dB(zW_in)
    zT = dB(zT_in)
    
    y = y * 1e3 # convert distance to mm

    axs0[0,i].plot(F[F > 0], zW[-1], color='r')
    axs1[0,i].plot(T,        zT[-1], color='r')
    
    axs2[0].plot(F[F > 0], zW[-1], label='%i'%width)
    axs2[1].plot(T,        zT[-1], label='%i'%width)
    

    extent = (np.min(F[F > 0]), np.max(F[F > 0]), np.min(y), np.max(y[:-1]))
    axs0[1,i].imshow(zW, extent=extent, vmin=np.max(zW) - 60.0,
               vmax=np.max(zW), aspect='auto', origin='lower')

    extent = (np.min(T), np.max(T), np.min(y), np.max(y[:-1]))
    axs1[1,i].imshow(zT, extent=extent, vmin=np.max(zT) - 60.0,
               vmax=np.max(zT), aspect='auto', origin='lower')

    print 'Total time: %.2f sec' % (time.time() - t)
    
    
    for ax in (axs0[0,i], axs1[0,i]):
        ax.set_title('Width: %i nm'%width)
    



for fig in (fig0, fig1, fig2):
    fig.suptitle('Thick:%inm, Pulse:%ifs, WL:%inm, Energy:%.1fpJ, Alpha:%.2fdB/cm'%(thick, 
                              Pulse*1e3,  pulseWL, EPP*1e12,      Alpha) )

for ax in axs2:
    ax.legend(fontsize=10, frameon=False)
    

for ax in axs0[0]: # top row of freq
    ax.set_xlim(0,400)
    ax.set_ylim(-80,0)

axs2[0].set_xlim(0,400)
axs2[0].set_ylim(-80,0)

    
for ax, axb in zip(axs0[2], axs0b):
    ax.set_ylim(-500,500)
    axb.set_ylim(-2e4,2e4)
    
for ax in axs1[0]: # top row of time
    ax.set_ylim(-50,30)
    
axs2[1].set_xlim(-Window/2,Window/2) # axs2 time panel
axs2[1].set_ylim(-50,30) # axs2 time panel

for i,fig in enumerate((fig0, fig1, fig2)):
    fig.savefig('Thick=%inm_Pulse=%ifs_WL=%inm_Energy=%.1fpJ_Alpha=%.2fdB_per_cm__Fig. %i.png'%(thick, 
                              Pulse*1e3,  pulseWL, EPP*1e12,      Alpha, i+1), dpi=150)

plt.show()


# def reinterpolate(x,y,Z,nx=500,ny=500):
#     # x and y are 1D arrays.
#     # z is a 2D array
#     # nx and ny are the number of points in the output image
#
#     print 'Reinterpolating....'
#     X,Y = np.meshgrid(x,y) # make 2D arrays
#
#     x_even = np.linspace(np.min(x),np.max(x),nx)
#     y_even = np.linspace(np.min(y),np.max(y),ny)
#     grid_x, grid_y = np.meshgrid(x_even,y_even)
#
#     grid_z = griddata((X.ravel(),Y.ravel()), Z.ravel(), (grid_x, grid_y), method='nearest')
#
#     return grid_z
