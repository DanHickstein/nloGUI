import numpy as np
import matplotlib.pyplot as plt
from pynlo.interactions.FourWaveMixing import SSFM
from pynlo.media.fibers import fiber
from pynlo.light.DerivedPulses import SechPulse
import time
import scipy.integrate
from autoscale_magic import autoscale_y
from scipy.misc import factorial


Pulse = 0.200  # pulse duration (ps)
Wave  = 1550    # pulse central wavelength (nm)
EPP   = 40e-12     # Pulse energy (J)
GDD   = 0.0     # Group delay dispersion (ps^2)
TOD   = 0.0     # Third order dispersion (ps^3)

Window = 10.0  # simulation window (ps)
Steps  = 50    # simulation steps
Points = 2**13  # simulation points

Length = 50    # fiber length (mm)
Alpha = 0.5   # attentuation coefficient
Gamma = 3000    # Gamma (1/(W km) -- 1400 is for Silicon Nitride
Beta2 = -500 # Beta_2 (ps^2/km)
Beta3 = 0   # Beta_3 (ps^3/km)
Beta4 = 0.5   # Beta_4 (ps^4/km)
FibWL = 1550   # Center WL (nm)

iRaman = True
iSteep = True

def dB(num):
    return 10 * np.log10(np.abs(num)**2)

fig1 = plt.figure(figsize=(13, 8))

axL1 = plt.subplot2grid((3, 4), (2, 0))
axL2 = plt.subplot2grid((3, 4), (2, 1), sharex=axL1)
axL3 = plt.subplot2grid((3, 4), (1, 1), sharex=axL1)
axL4 = plt.subplot2grid((3, 4), (0, 1), sharex=axL1)

ax1 = plt.subplot2grid((3, 4), (0, 2), sharex=axL1)
ax2 = plt.subplot2grid((3, 4), (0, 3))

ax3 = plt.subplot2grid((3, 4), (1, 2), rowspan=2, sharex=axL1)
ax4 = plt.subplot2grid((3, 4), (1, 3), rowspan=2, sharey=ax3, sharex=ax2)


for label in ax4.get_yticklabels():
    label.set_visible(False)

plt.subplots_adjust(left=0.05, bottom=0.09, right=0.96,
                    top=0.96, wspace=0.3, hspace=0.29)

axL1.set_xlabel('Frequency (Hz)')
axL1.set_ylabel('D (ps/nm/km)')
axL2.set_xlabel('Frequency (Hz)')
axL2.set_ylabel(r'$\beta_2$ (ps$^2$/km)')
axL3.set_ylabel('FWM phase mismatch (1/m)')

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


class BraggFiber(fiber.FiberInstance):
    
    # def add_bragg_grating(self, lambdaB=1000e-9, delta_n=0.001, scaling=1):
    def add_bragg_grating(self, offset=50, amp=1e9, sigma=1):
        
        """ The allows the dispersion from a periodic bragg grating to be included
        
        Parameters
        ----------
        lambdaB : float
            the wavelength of the bragg grating (in meters)
            you can calculate this as lambdaB = 2 * (index of ref.) * spacing
        scaling : float
            this is the "radial overlap" in the Westbrook papers
            it just scales the effect of the grating
        delta_n : float
            the difference between the high and low index portions of the grating
            this is 0.003 in the Westbrook papers.
        
        """
        self.fiberspecs["bragg_grating"] = True
        
        # self.lambdaB = lambdaB
        # self.delta_n = delta_n
        # self.scaling = scaling
        
        self.amp = amp
        self.sigma = sigma
        self.offset = offset
    
    def get_betas(self,pulse):
        B = np.zeros((pulse.NPTS,))
        if self.fiberspecs["dispersion_format"] == "D":
            self.betas = DTabulationToBetas(pulse.center_wavelength_nm,
                                            np.transpose(np.vstack((self.x,self.y))),
                                            self.poly_order,
                                            DDataIsFile = False)
            for i in range(len(self.betas)):
                B = B + self.betas[i]/factorial(i+2)*pulse.V_THz**(i+2)

        elif self.fiberspecs["dispersion_format"] == "GVD":
            # calculate beta[n]/n! * (w-w0)^n
            # w0 is the center of the Taylor expansion, and is defined by the
            # fiber. the w's are from the optical spectrum
            fiber_omega0 =  2*np.pi*self.c / self.center_wavelength # THz
            betas = self.betas
            for i in range(len(betas)):
                betas[i] = betas[i]
                B = B + betas[i] / factorial(i + 2) * (pulse.W_THz-fiber_omega0)**(i + 2)
        else:
            return -1 

        if "bragg_grating" in self.fiberspecs and self.fiberspecs["bragg_grating"]:
            wl = pulse.wl_mks
            # print self.lambdaB
            # B_at_lambdaB = B[np.argmin(np.abs(wl-self.lambdaB))]
            #
            # plus_minus = np.ones(np.shape(B))
            # plus_minus[wl < self.lambdaB] = -1
            #
            # k = np.pi * self.delta_n * self.scaling / self.lambdaB
            # print k
            #
            # B = B_at_lambdaB + plus_minus*np.sqrt( (B - B_at_lambdaB)**2 - k**2)
            # # B = np.sqrt( (B - B_at_lambdaB)**2 - k**2)
            
            def gaussian(x, mu, a, sigma):
                return a * np.exp(-((x - mu) ** 2) / 2 / sigma ** 2)
            
            print pulse.center_frequency_THz
            
            B = B + gaussian(pulse.W_THz/(2*np.pi), pulse.center_frequency_THz - self.offset, self.amp, self.sigma)
            B = B + gaussian(pulse.W_THz/(2*np.pi), pulse.center_frequency_THz + self.offset, self.amp, self.sigma)
            
            
                
        return B
        


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
    
    fiber1 = BraggFiber()
    fiber1.generate_fiber(fiber_length, center_wl_nm=fibWL, betas=(
        Beta2, Beta3, Beta4), gamma_W_m=Gamma * 1e-3, gvd_units='ps^n/km', gain=-alpha)
    
    fiber1.add_bragg_grating()
    

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
    intBeta1 = scipy.integrate.cumtrapz(beta[W > 0])
    intBeta1 = intBeta1 - intBeta1[fibWLindex]

    intBeta2 = scipy.integrate.cumtrapz(intBeta1)
    intBeta2 = intBeta2 - intBeta2[fibWLindex]
    lineL3a.set_data(F[W > 0][:-2], intBeta2)

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
                     disable_Raman=False, disable_self_steepening=False)

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
