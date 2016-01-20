import numpy as np
import matplotlib.pyplot as plt
from pynlo.interactions.FourWaveMixing import SSFM
from pynlo.media.fibers import fiber
from pynlo.light.DerivedPulses import SechPulse
from scipy.interpolate import griddata

from matplotlib.widgets import Button, TextBox


def reinterpolate(x,y,Z,nx=500,ny=500):
    # x and y are 1D arrays. 
    # z is a 2D array
    # nx and ny are the number of points in the output image
    
    X,Y = np.meshgrid(x,y) # make 2D arrays
    
    x_even = np.linspace(np.min(x),np.max(x),nx)
    y_even = np.linspace(np.min(y),np.max(y),ny)
    grid_x, grid_y = np.meshgrid(x_even,y_even)
    
    print np.shape(X),np.shape(Y), np.shape(Z)
    grid_z = griddata((X.ravel(),Y.ravel()), Z.ravel(), (grid_x, grid_y), method='nearest')
    
    return grid_z
    
    

fig = plt.figure(figsize=(12,8))

ax1 = plt.subplot2grid((3,4), (0,2))
ax2 = plt.subplot2grid((3,4), (0,3))

ax3 = plt.subplot2grid((3,4), (1,2),rowspan=2)
ax4 = plt.subplot2grid((3,4), (1,3),rowspan=2)

plt.subplots_adjust(left=0.05,bottom=0.09,right=0.96,top=0.96,wspace=0.21,hspace=0.29)

axRun = plt.axes([0.4,0.05,0.09,0.04])
bRun = Button(axRun, 'Run Simulation')


###Parameter axes###
sl = 0.1
sw = 0.1
sh = 0.04
sv = 0.97 #start vertical
ss = -0.045 #spacing

ax01 = plt.axes([sl, sv+1*ss, sw, sh])
ax01 = plt.axes([sl, sv+2*ss, sw, sh])

ax02   = plt.axes([sl, sv+4*ss, sw, sh])
ax04    = plt.axes([sl, sv+5*ss, sw, sh])

ax05 = plt.axes([sl, sv+7*ss, sw, sh])
ax06 = plt.axes([sl, sv+8*ss, sw, sh])
ax07 = plt.axes([sl, sv+9*ss, sw, sh])

ax07 = plt.axes([sl, sv+11*ss, sw, sh])
ax08 = plt.axes([sl, sv+12*ss, sw, sh])

# axtnum = plt.axes([sl, sv+14*ss, sw, sh])
# axtsta = plt.axes([sl, sv+15*ss, sw, sh])
# axtend = plt.axes([sl, sv+16*ss, sw, sh])

# axcolor = plt.axes([sl, sv+18*ss, sw/2.5, sh*4])
# axrot   = plt.axes([sl+0.2, sv+18*ss, sw/2.5, sh*4])





def run_simulation(caller=None):
    print caller
    dz = 1e-4
    steps = 10
    range1 = np.arange(steps)

    centerwl = 835.0
    fiber_length = 0.01

    pump_power = 1.0e4
    pump_pulse_length = 28.4e-3
    npoints = 2**13


    init = SechPulse(pump_power, pump_pulse_length, centerwl, time_window = 10.0,
                        GDD = 0, TOD = 0.0, NPTS = npoints, frep_MHz = 100, power_is_avg = False)

    fiber1 = fiber.FiberInstance() 
    fiber1.load_from_db( fiber_length, 'dudley')

    evol = SSFM.SSFM(dz = 1e-6, local_error = 0.001, USE_SIMPLE_RAMAN = True)
    y = np.zeros(steps)
    AW = np.zeros((init.NPTS, steps))
    AT = np.copy(AW)

    y, AW, AT, pulse1 = evol.propagate(pulse_in = init, fiber = fiber1, 
                                       n_steps = steps)
                           
    wl = init.wl_nm

    loWL = 400
    hiWL = 1400
                         
    iis = np.logical_and(wl>loWL,wl<hiWL)

    iisT = np.logical_and(init.T_ps>-1,init.T_ps<5)

    xW = wl[iis]
    xT = init.T_ps[iisT]
    zW_in = np.transpose(AW)[:,iis]
    zT_in = np.transpose(AT)[:,iisT]
    zW = 10*np.log10(np.abs(zW_in)**2)
    zT = 10*np.log10(np.abs(zT_in)**2)
    mlIW = np.max(zW)
    mlIT = np.max(zT)

    D = fiber1.Beta2_to_D(init)
    beta = fiber1.Beta2(init)

    print np.shape(wl),np.shape(D)
    ax1.plot(wl,D)
    ax1.set_xlim(400,1600)
    ax1.set_ylim(-400,300)
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('D (ps/nm/km)')

    ax2.plot(wl,beta*1000)
    ax2.set_xlim(400,1600)
    ax2.set_ylim(-350,200)
    ax2.set_xlabel('Wavelength (nm)')
    ax2.set_ylabel(r'$\beta_2$ (ps$^2$/km)')

    zW_grid = reinterpolate(xW, y[:-1], zW)
    extent = (np.min(xW),np.max(xW),np.min(y),np.max(y[:-1]))
    ax3.imshow(zW_grid, extent=extent, vmin = mlIW - 40.0, vmax = mlIW,aspect='auto',origin='lower')

    ax3.autoscale(tight=True)
    ax3.set_xlim([loWL, hiWL])
    ax3.set_xlabel('Wavelength (nm)')
    ax3.set_ylabel('Distance (m)')

    # plt.pcolormesh(xT, y, zT, vmin = mlIT - 40.0, vmax = mlIT)
    
    zT_grid = reinterpolate(xT, y[:-1], zT)
    extent = (np.min(xT),np.max(xT),np.min(y),np.max(y[:-1]))
    ax4.imshow(zT_grid, extent=extent, vmin = mlIT - 40.0, vmax = mlIT,aspect='auto',origin='lower')
    ax4.autoscale(tight=True)
    ax4.set_xlabel('Delay (ps)')
    ax4.set_ylabel('Distance (m)')

bRun.on_clicked(run_simulation)

# run_simulation()

plt.show()