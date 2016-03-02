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


Pulse = 0.125  # pulse duration (ps)
Wave  = 1550    # pulse central wavelength (nm)
EPP   = 40e-12     # Pulse energy (J)
GDD   = 0.0     # Group delay dispersion (ps^2)
TOD   = 0.0     # Third order dispersion (ps^3)

Window = 10.0  # simulation window (ps)
Steps  = 50    # simulation steps
Points = 2**13  # simulation points

Length = 40    # fiber length (mm)
Alpha = 0.4  # attentuation coefficient
Gamma = 3000    # Gamma (1/(W km) -- 1400 is for Silicon Nitride
thick = 660
width = 2100 
FibWL = 1550   # Center WL (nm)

iRaman = False # disable raman
iSteep = False

def dB(num):
    return 10 * np.log10(np.abs(num)**2)


fig0, axs0 = plt.subplots(3,4, figsize=(13, 9), sharex=True, sharey='row')

plt.show()
axs

# axL1 = plt.subplot2grid((3, 4), (2, 0))
# axL2 = plt.subplot2grid((3, 4), (2, 1), sharex=axL1)
# axL3 = plt.subplot2grid((3, 4), (1, 1), sharex=axL1)
# axL4 = plt.subplot2grid((3, 4), (0, 1), sharex=axL1)
#
# ax1 = plt.subplot2grid((3, 4), (0, 2), sharex=axL1)
# ax2 = plt.subplot2grid((3, 4), (0, 3))
#
# ax3 = plt.subplot2grid((3, 4), (1, 2), rowspan=2, sharex=axL1)
# ax4 = plt.subplot2grid((3, 4), (1, 3), rowspan=2, sharey=ax3, sharex=ax2)


for label in ax4.get_yticklabels():
    label.set_visible(False)

plt.subplots_adjust(left=0.05, bottom=0.09, right=0.96,
                    top=0.96, wspace=0.3, hspace=0.29)

axL1.set_xlabel('Frequency (Hz)')
axL1.set_ylabel('D (ps/nm/km)')
axL2.set_xlabel('Frequency (Hz)')
axL2.set_ylabel(r'$\beta_2$ (ps$^2$/km)')
axL3.set_ylabel('Propagation const (beta) (1/m)')

ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('Intensity (dB)')
ax2.set_xlabel('Time (ps)')


lineL1,  = axL1.plot(0, 0)  # make some dummy lines to hold future data
lineL2,  = axL2.plot(0, 0)
lineL3a, = axL3.plot(0, 0)
lineL3b, = axL3.plot(0, 0)

lineL4, = axL4.plot(0, 0)

line1a, = ax1.plot(0, 0, color='b')
line1b, = ax1.plot(0, 0, color='r', lw=1.5)

line2a, = ax2.plot(0, 0, color='b')
line2b, = ax2.plot(0, 0, color='r', lw=1.5)

vline1a = axL1.axvline(0, color='m', alpha=0.5)
vline1b = axL1.axvline(0, color='g', lw=1.5, alpha=0.5)
vline2a = axL2.axvline(0, color='m', alpha=0.5)
vline2b = axL2.axvline(0, color='g', lw=1.5, alpha=0.5)

hline3 = axL3.axhline(0, color='k', lw=1.5, alpha=0.5)

axL1.axhline(0, alpha=0.5, color='k')
axL2.axhline(0, alpha=0.5, color='k')


def run_simulation():

    pump_pulse_length = Pulse
    centerwl = Wave  # in nm

    steps = Steps

    fiber_length = Length * 1e-3
    alpha = np.log((10**(Alpha * 0.1))) * 100  # convert from dB/cm to 1/m

    fibWL = FibWL

    # set up the pulse parameters
    pulse = SechPulse(1, pump_pulse_length, centerwl, time_window_ps=Window,
                      GDD=GDD, TOD=TOD, NPTS=Points, frep_MHz=100, power_is_avg=False)

    pulse.set_epp(EPP)
    
    betas = get_betas_from_dimensions(thick, width)
    
    fiber1 = fiber.FiberInstance()
    fiber1.generate_fiber(fiber_length, center_wl_nm=fibWL, betas=betas,
                                  gamma_W_m=Gamma * 1e-3, gvd_units='ps^n/km', gain=-alpha)
        

    laser_freq = 3e8 / (centerwl * 1e-9)
    fiber_freq = 3e8 / (fibWL * 1e-9)

    W = pulse.W_mks
    F = W / (2 * np.pi)
    T = pulse.T_ps
    xW = W[W > 0]

    D = fiber1.Beta2_to_D(pulse)
    beta = fiber1.Beta2(pulse) * 1e3  # convert ps^2/m to ps^2/km

    # Plot the dispersion in the left plots
    lineL1.set_data(F[W > 0], D[W > 0])
    lineL2.set_data(F[W > 0], beta[W > 0])
    fibWLindex = (np.abs(F[W > 0] - fiber_freq)).argmin()

    # plot the integral of beta
    # intBeta1 = scipy.integrate.cumtrapz(beta[W > 0])
    # intBeta1 = intBeta1 - intBeta1[fibWLindex]
    #
    # intBeta2 = scipy.integrate.cumtrapz(intBeta1)
    # intBeta2 = intBeta2 - intBeta2[fibWLindex]
    # lineL3a.set_data(F[W > 0][:-2], intBeta2)
    lineL3a.set_data(F,fiber1.get_betas(pulse))

    # plot the pulse in the top plots
    line1a.set_data(F[W > 0], dB(pulse.AW[W > 0]))
    line2a.set_data(T, dB(pulse.AT))

    vline1a.set_xdata((laser_freq, laser_freq))
    vline1b.set_xdata((fiber_freq, fiber_freq))
    vline2a.set_xdata((laser_freq, laser_freq))
    vline2b.set_xdata((fiber_freq, fiber_freq))

    # Propagation
    t = time.time()

    # set up the propagation parameters
    evol = SSFM.SSFM(local_error=0.001, USE_SIMPLE_RAMAN=True,
                     disable_Raman=iRaman, disable_self_steepening=iSteep)

    # propagate the pulse!
    y, AW, AT, pulse1 = evol.propagate(
        pulse_in=pulse, fiber=fiber1, n_steps=steps)

    zW_in = np.transpose(AW)[:, (W > 0)]
    zT_in = np.transpose(AT)
    zW = dB(zW_in)
    zT = dB(zT_in)

    line1b.set_data(F[F > 0], zW[-1])
    line2b.set_data(T, zT[-1])

    extent = (np.min(F[F > 0]), np.max(F[F > 0]), np.min(y), np.max(y[:-1]))
    ax3.imshow(zW, extent=extent, vmin=np.max(zW) - 60.0,
               vmax=np.max(zW), aspect='auto', origin='lower')
    ax3.set_xlabel('Frequency (Hz)')
    ax3.set_ylabel('Distance (m)')

    extent = (np.min(T), np.max(T), np.min(y), np.max(y[:-1]))
    ax4.imshow(zT, extent=extent, vmin=np.max(zT) - 60.0,
               vmax=np.max(zT), aspect='auto', origin='lower')
    ax4.set_xlabel('Delay (ps)')

    print 'Total time: %.2f sec' % (time.time() - t)

    reset_plots()


def reset_plots(caller=None):
    for ax in (axL1, axL2, axL3, axL4, ax1, ax2):  # rescale the x and ylims
        ax.relim()
        ax.autoscale_view()

    ax1.set_ylim(-40, 10)
    ax2.set_ylim(0, 50)

    axL1.set_xlim(100e12, 300e12)

    for ax in (axL1, axL2, axL3):  # rescale the x and ylims
        autoscale_y(ax, linenum=0)


run_simulation()
plt.suptitle('thickness=%i, width=%i'%(thick, width))
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
