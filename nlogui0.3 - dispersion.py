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
    
    print 'Reinterpolating....'
    X,Y = np.meshgrid(x,y) # make 2D arrays
    
    x_even = np.linspace(np.min(x),np.max(x),nx)
    y_even = np.linspace(np.min(y),np.max(y),ny)
    grid_x, grid_y = np.meshgrid(x_even,y_even)
    
    grid_z = griddata((X.ravel(),Y.ravel()), Z.ravel(), (grid_x, grid_y), method='nearest')
    
    return grid_z
    

fig = plt.figure(figsize=(13,9))

axL1 = plt.subplot2grid((3,4), (2,0))
axL2 = plt.subplot2grid((3,4), (2,1))

ax1 = plt.subplot2grid((3,4), (0,2))
ax2 = plt.subplot2grid((3,4), (0,3))

ax3 = plt.subplot2grid((3,4), (1,2),rowspan=2)
ax4 = plt.subplot2grid((3,4), (1,3),rowspan=2)

plt.subplots_adjust(left=0.05,bottom=0.09,right=0.96,top=0.96,wspace=0.21,hspace=0.29)


axRun = plt.axes([0.4,0.90,0.09,0.04])
bRun = Button(axRun, 'Run Simulation')


###Parameter axes###
sl = 0.2
sw = 0.1
sh = 0.04
sv = 0.97 #start vertical
ss = -0.045 #spacing

ax01 = plt.axes([sl, sv+1*ss, sw, sh])
ax02 = plt.axes([sl, sv+2*ss, sw, sh])
ax03 = plt.axes([sl, sv+3*ss, sw, sh])

ax04 = plt.axes([sl, sv+5*ss, sw, sh])

ax05 = plt.axes([sl, sv+7*ss, sw, sh])
ax06 = plt.axes([sl, sv+8*ss, sw, sh])
ax07 = plt.axes([sl, sv+9*ss, sw, sh])

bPulse = TextBox(ax01,'Pulse duration (ps)',initial='0.0284')
bWave  = TextBox(ax02,'Center wavelength (nm)',initial='835')
bPower = TextBox(ax03,'Pulse peak power (kW)',initial='10') # should change to pulse energy

bLength = TextBox(ax04,'Fiber length (mm)',initial='10') 

bDz    = TextBox(ax05,'dz',initial='0.0001') 
bSteps = TextBox(ax06,'Steps',initial='100') 
bPoints = TextBox(ax07,'Simulations points',initial='2**13')

lineL1, = axL1.plot(0,0) # make some dummy lines to hold future data
lineL2, = axL2.plot(0,0)

line1a, = ax1.plot(0,0,color='b') # make some dummy lines to hold future data
line1b, = ax1.plot(0,0,color='r') # make some dummy lines to hold future data

line2a, = ax2.plot(0,0,color='b') # make some dummy lines to hold future data
line2b, = ax2.plot(0,0,color='r') # make some dummy lines to hold future data



def run_simulation(caller=None,prop=False):
    print 'prop: ',prop
    print caller
    pump_pulse_length = float(bPulse.text)
    centerwl          = float(bWave.text) # in nm
    pump_power        = float(bPower.text) * 1e3
    
    fiber_length = float(bLength.text) * 1e-3

    dz    = float(bDz.text)
    steps = int(bSteps.text)

    npoints = int(eval(bPoints.text))
    
    print 'Running simulation with'
    print 'Pulse duration (fs): %f'%pump_pulse_length
    print 'Center wavelength (nm): %f'%centerwl
    print 'Pump power (W): %f'%pump_power
    
    print 'Fiber length (m): %f'%fiber_length
    
    print 'dz: %f'%dz
    print 'Steps: %i'%steps
    print 'Number of points: %i'%npoints
    

    init = SechPulse(pump_power, pump_pulse_length, centerwl, time_window = 10.0,
                        GDD = 0, TOD = 0.0, NPTS = npoints, frep_MHz = 100, power_is_avg = False)

    fiber1 = fiber.FiberInstance() 
    fiber1.load_from_db( fiber_length, 'dudley')
    
    evol = SSFM.SSFM(dz = dz, local_error = 0.001, USE_SIMPLE_RAMAN = True)
    y = np.zeros(steps)
    AW = np.zeros((init.NPTS, steps))
    AT = np.copy(AW)
    
                           
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

    for ax in (ax1,ax2,ax3,ax4):
        ax.clear()
    
    print wl,D
    lineL1.set_data(wl[wl>0],D[wl>0])
    axL1.set_xlim(400,1600)
    axL1.set_ylim(-400,300)
    axL1.set_xlabel('Wavelength (nm)')
    axL1.set_ylabel('D (ps/nm/km)')

    lineL2.set_data(wl,beta*1000)
    axL2.set_xlim(400,1600)
    axL2.set_ylim(-350,200)
    axL2.set_xlabel('Wavelength (nm)')
    axL2.set_ylabel(r'$\beta_2$ (ps$^2$/km)')
    
     
    if prop==False: return
    
    y, AW, AT, pulse1 = evol.propagate(pulse_in = init, fiber = fiber1, 
                                       n_steps = steps)

    zW_grid = reinterpolate(xW, y[:-1], zW)
    extent = (np.min(xW),np.max(xW),np.min(y),np.max(y[:-1]))
    ax3.imshow(zW_grid, extent=extent, vmin = mlIW - 40.0, vmax = mlIW,aspect='auto',origin='lower')

    ax3.autoscale(tight=True)
    ax3.set_xlim([loWL, hiWL])
    ax3.set_xlabel('Wavelength (nm)')
    ax3.set_ylabel('Distance (m)')
        
    zT_grid = reinterpolate(xT, y[:-1], zT)
    extent = (np.min(xT),np.max(xT),np.min(y),np.max(y[:-1]))
    ax4.imshow(zT_grid, extent=extent, vmin = mlIT - 40.0, vmax = mlIT,aspect='auto',origin='lower')
    ax4.autoscale(tight=True)
    ax4.set_xlabel('Delay (ps)')

    for label in ax4.get_yticklabels():
        label.set_visible(False)
        
    ax1.plot(xW,zW[0],color='b',label='Before')
    ax1.plot(xW,zW[-1],color='r',label='After')
    ax1.set_ylim(-40,0)
    
    ax2.plot(xT,zT[0],color='b',label='Before')
    ax2.plot(xT,zT[-1],color='r',label='After')
    ax2.set_ylim(-10,40)

def run_simulation_full():
    run_simulation(prop=True)


for box in (bPulse, bWave, bPower, bLength, bDz, bSteps, bPoints):
    box.on_submit(run_simulation)

bRun.on_clicked(run_simulation_full)

plt.show()