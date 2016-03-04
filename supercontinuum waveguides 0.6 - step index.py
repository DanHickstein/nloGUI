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

# Length1  = 10.0   # fiber length (mm)
# Length2  = 15.0    # fiber length (mm)
# deltaW   = 1100     # change in fiber width (nm)

Alpha   = 0.2 # attentuation coefficient
Gamma   = 2000    # Gamma (1/(W km) -- 1400 is for Silicon Nitride

# thick   = 600
# widths  = (1200,1600,2000,2400)

thick   = 660
# widths  = (1700,2100,2500,2900)
# widths  = (1700, 2600)

combos = []



# study of step versus regular

narrow = 1500.
middle = 2200
wide   = 2800.

combos.append((narrow, wide, # Width1, Width2 (nm)
               30.0,  0) ) # Length1, Length2 (mm)

combos.append((wide, wide, # Width1, Width2 (nm)
               30.0,  0) ) # Length1, Length2 (mm)

combos.append((narrow, wide, # Width1, Width2 (nm)
               15.0, 15.0) ) # Length1, Length2 (mm)

combos.append((middle, middle, # Width1, Width2 (nm)
                30.0, 0) ) # Length1, Length2 (mm)


# scan step positoin

# narrow = 1700.
# wide   = 2900.
#
# combos.append((narrow, wide, # Width1, Width2 (nm)
#                0.1,  29.9) ) # Length1, Length2 (mm)
#
# combos.append((narrow, wide, # Width1, Width2 (nm)
#                5,  25) ) # Length1, Length2 (mm)
#
# combos.append((narrow, wide, # Width1, Width2 (nm)
#                15, 15) ) # Length1, Length2 (mm)
#
# combos.append((narrow, wide, # Width1, Width2 (nm)
#                 25, 5) ) # Length1, Length2 (mm)




# opposite case

# wide = 1700.
# narrow   = 2900.
#
# combos.append((narrow, wide, # Width1, Width2 (nm)
#                0.1,  39.9) ) # Length1, Length2 (mm)
#
# combos.append((narrow, wide, # Width1, Width2 (nm)
#                5,  35) ) # Length1, Length2 (mm)
#
# combos.append((narrow, wide, # Width1, Width2 (nm)
#                15, 25) ) # Length1, Length2 (mm)
#
# combos.append((narrow, wide, # Width1, Width2 (nm)
#                 25, 15) ) # Length1, Length2 (mm)

# fine step scan

# narrow = 1700.
# wide   = 2900.
#
# combos.append((narrow, wide, # Width1, Width2 (nm)
#                10,  10) ) # Length1, Length2 (mm)
#
# combos.append((narrow, wide, # Width1, Width2 (nm)
#                12,  8) ) # Length1, Length2 (mm)
#
# combos.append((narrow, wide, # Width1, Width2 (nm)
#                14, 6) ) # Length1, Length2 (mm)
#
# combos.append((narrow, wide, # Width1, Width2 (nm)
#                16, 4) ) # Length1, Length2 (mm)






# # opposite case:
# narrow = 1700.
# wide   = 2900.
#
# combos.append((wide, narrow, # Width1, Width2 (nm)
#                1.,  24.) ) # Length1, Length2 (mm)
#
# combos.append((wide, narrow, # Width1, Width2 (nm)
#                5,  20) ) # Length1, Length2 (mm)
#
# combos.append((wide, narrow, # Width1, Width2 (nm)
#                10., 15.0) ) # Length1, Length2 (mm)
#
# combos.append((wide, narrow, # Width1, Width2 (nm)
#                 20., 5.) ) # Length1, Length2 (mm)



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


for i, (Width1, Width2, Length1, Length2) in enumerate(combos):
    
    Width1  = float(Width1)
    Width2  = float(Width2)
    Length1 = float(Length1)
    Length2 = float(Length2)
    
    Length_tot = Length1+Length2
    
    Steps1 = np.ceil(Steps * Length1 / Length_tot)
    Steps2 = np.ceil(Steps * Length2 / Length_tot)
    
    print Steps1, Steps2

    alpha = np.log((10**(Alpha * 0.1))) * 100  # convert from dB/cm to 1/m

    # set up the pulse parameters
    pulse = SechPulse(1, Pulse, pulseWL, time_window_ps=Window,
                      GDD=GDD, TOD=TOD, NPTS=Points, frep_MHz=100, power_is_avg=False)

    pulse.set_epp(EPP)
    
    betas1 = get_betas_from_dimensions(thick, Width1)
    betas2 = get_betas_from_dimensions(thick, Width2)
    
    fiber1 = fiber.FiberInstance()
    fiber1.generate_fiber(Length1 * 1e-3, center_wl_nm=fibWL, betas=betas1,
                                  gamma_W_m=Gamma * 1e-3, gvd_units='ps^n/km', gain=-alpha)
                                  
    fiber2 = fiber.FiberInstance()
    fiber2.generate_fiber(Length2 * 1e-3, center_wl_nm=fibWL, betas=betas2,
                                  gamma_W_m=Gamma * 1e-3, gvd_units='ps^n/km', gain=-alpha)
        

    laser_freq = 3e8 / (pulseWL * 1e-9) * 1e-12
    fiber_freq = 3e8 / (fibWL * 1e-9)   * 1e-12

    W = pulse.W_mks
    F = W / (2 * np.pi) * 1e-12 # convert to THz
    T = pulse.T_ps
    xW = W[W > 0]

    D = fiber1.Beta2_to_D(pulse)
    betaFib1 = fiber1.Beta2(pulse) * 1e3  # convert ps^2/m to ps^2/km
    betaFib2 = fiber2.Beta2(pulse) * 1e3  # convert ps^2/m to ps^2/km
    

    # Plot the dispersion in the bottom panels
    axs0[2,i].plot(F[W > 0], betaFib1[W > 0],  color='g')
    axs0[2,i].plot(F[W > 0], betaFib2[W > 0],  color='g', ls='dashed')
    
    # Weird beta
    axs0b[ i].plot(F, fiber1.get_betas(pulse), color='m')
    axs0b[ i].plot(F, fiber2.get_betas(pulse), color='m', ls='dashed')
    
    # zero-lines
    axs0[2,i].axhline(0,color='g', alpha=0.5)
    axs0b[ i].axhline(0,color='m', alpha=0.5)

    # plot *before* in the top panels (blue)
    axs0[0,i].plot(F[W > 0], dB(pulse.AW[W > 0]), color='b') # freq
    axs1[0,i].plot(T,        dB(pulse.AT),        color='b') # time

    # plot some vertical lines
    for ax in (axs0[0,i], axs0[2,i]):
        # ax.axvline(fiber_freq,   color='g', alpha=0.1)
        ax.axvline(laser_freq*0.5,   color='g', alpha=0.3)
        ax.axvline(laser_freq,   color='b', alpha=0.3)
        ax.axvline(laser_freq*2, color='g', alpha=0.2)
    

    # Propagation
    t = time.time() # start a timer

    # set up the propagation parameters
    evol = SSFM.SSFM(local_error=0.001, USE_SIMPLE_RAMAN=True,
                     disable_Raman=np.logical_not(Raman), 
                     disable_self_steepening=np.logical_not(Steep))


    y, AW, AT, pulse_out = evol.propagate(pulse_in=pulse, fiber=fiber1, n_steps=Steps1)
    pulse = pulse_out
    
    if Length2 > 0:
        y2, AW2, AT2, pulse_out = evol.propagate(pulse_in=pulse_out, fiber=fiber2, n_steps=Steps2)
        AW = np.hstack((AW, AW2))
        AT = np.hstack((AT, AT2))
 
    zW = dB( np.transpose(AW)[:, (W > 0)] )
    zT = dB( np.transpose(AT) )

    y = y * 1e3 # convert distance to mm

    axs0[0,i].plot(F[F > 0], zW[-1], color='r')
    axs1[0,i].plot(T,        zT[-1], color='r')

    axs2[0].plot(F[F > 0], zW[-1], label='%i nm for %i mm\n%inm for %i mm'%(Width1, Length1, Width2, Length2))
    axs2[1].plot(T,        zT[-1], label='%i nm for %i mm\n%inm for %i mm'%(Width1, Length1, Width2, Length2))
    

    # top = np.max(zW)
    top = -10
    extent = (np.min(F[F > 0]), np.max(F[F > 0]), 0, Length_tot)
    axs0[1,i].imshow(zW, extent=extent, vmin=top - 60.0, 
                     vmax=top, aspect='auto', origin='lower')
    
    extent = (np.min(T), np.max(T), np.min(y), Length_tot)
    axs1[1,i].imshow(zT, extent=extent, vmin=np.max(zT) - 60.0,
               vmax=np.max(zT), aspect='auto', origin='lower')
            
    axs0[1,i].axhline(Length1, color='w', alpha=1, ls='dashed')
    axs1[1,i].axhline(Length1, color='w', alpha=1, ls='dashed')
        

    print 'Total time: %.2f sec' % (time.time() - t)
    
    
    for ax in (axs0[0,i], axs1[0,i]):
        ax.set_title('%i nm for %i mm\n%inm for %i mm'%(Width1, Length1, Width2, Length2), fontsize=10)
    
        
    # break

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


axs2[0].axvline(laser_freq* 0.5,   color='b', alpha=0.3)
axs2[0].axvline(laser_freq,   color='b', alpha=0.3)
axs2[0].axvline(laser_freq*2, color='g', alpha=0.2)

    
for ax, axb in zip(axs0[2], axs0b):
    ax.set_ylim(-500,500)
    axb.set_ylim(-2e4,2e4)
    
for ax in axs1[0]: # top row of time
    ax.set_ylim(-50,30)
    
axs2[1].set_xlim(-Window/2,Window/2) # axs2 time panel
axs2[1].set_ylim(-50,30) # axs2 time panel

for i,fig in enumerate((fig0, fig1, fig2)):
    fig.savefig('STEP - Thick=%inm_Pulse=%ifs_WL=%inm_Energy=%.1fpJ_Alpha=%.2fdB_per_cm__Fig. %i.png'%(thick, 
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
