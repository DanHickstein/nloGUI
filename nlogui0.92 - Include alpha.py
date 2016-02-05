import numpy as np
import matplotlib
# matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
from pynlo.interactions.FourWaveMixing import SSFM
from pynlo.media.fibers import fiber
from pynlo.light.DerivedPulses import SechPulse
from matplotlib.widgets import Button, TextBox, CheckButtons
import time
import scipy.integrate
from autoscale_magic import autoscale_y

iPulse = '0.200' # pulse duration (ps)
iWave  = '1550'    # pulse central wavelength (nm)
iEPP   = '500'     # Pulse energy (pJ)
iGDD   = '0.0'     # Group delay dispersion (ps^2)
iTOD   = '0.0'     # Third order dispersion (ps^3)

iWindow = '10.0'  # simulation window (ps)
iSteps  = '50'    # simulation steps
iPoints = '2**13' # simulation points

iLength = '10'    # fiber length (mm)
iAlpha  = '0.5'   # attentuation coefficient
iGamma =  '1000'    # Gamma (1/(W km) -- 1400 is for Silicon Nitride
iBeta2 =  '-200'  # Beta_2 (ps^2/km)
iBeta3 = '-0.1'   # Beta_3 (ps^3/km)
iBeta4 = '0.02'   # Beta_4 (ps^4/km)
iFibWL = '1550'   # Center WL (nm)

# conner
# iPulse = '0.200' # pulse duration (ps)
# iWave  = '1000'    # pulse central wavelength (nm)
# iEPP   = '10000'     # Pulse energy (pJ)
# iGDD   = '0.0'     # Group delay dispersion (ps^2)
# iTOD   = '0.0'     # Third order dispersion (ps^3)
#
# iWindow = '10.0'  # simulation window (ps)
# iSteps  = '50'    # simulation steps
# iPoints = '2**13' # simulation points
#
# iLength = '15'    # fiber length (mm)
# iAlpha  = '0.001'   # attentuation coefficient
# iGamma =  '20'    # Gamma (1/(W km) -- 1400 is for Silicon Nitride
# iBeta2 =  '-200'  # Beta_2 (ps^2/km)
# iBeta3 = '0.0'   # Beta_3 (ps^3/km)
# iBeta4 = '0.04'   # Beta_4 (ps^4/km)
# iFibWL = iWave   # Center WL (nm)


# # working supercontinuum with SiN
# iPulse = '0.200' # pulse duration (ps)
# iWave  = '1550'    # pulse central wavelength (nm)
# iEPP   = '100'     # Pulse energy (pJ)
# iGDD   = '0.0'     # Group delay dispersion (ps^2)
# iTOD   = '0.0'     # Third order dispersion (ps^3)
#
# iWindow = '10.0'  # simulation window (ps)
# iSteps  = '50'    # simulation steps
# iPoints = '2**13' # simulation points
#
# iLength = '20'    # fiber length (mm)
# iAlpha  = '0.8'   # attentuation coefficient
# iGamma =  '1400'    # Gamma (1/(W km) -- 1400 is for Silicon Nitride
# iBeta2 =  '-200'  # Beta_2 (ps^2/km)
# iBeta3 = '-1'   # Beta_3 (ps^3/km)
# iBeta4 = '0.02'   # Beta_4 (ps^4/km)
# iFibWL = '1550'   # Center WL (nm)

# # Underpowered SiN superC
# iPulse = '0.200' # pulse duration (ps)
# iWave  = '1550'    # pulse central wavelength (nm)
# iEPP   = '50'     # Pulse energy (pJ)
# iGDD   = '0.0'     # Group delay dispersion (ps^2)
# iTOD   = '0.0'     # Third order dispersion (ps^3)
#
# iWindow = '10.0'  # simulation window (ps)
# iSteps  = '50'    # simulation steps
# iPoints = '2**13' # simulation points
#
# iLength = '20'    # fiber length (mm)
# iAlpha  = '0.8'   # attentuation coefficient
# iGamma =  '1400'    # Gamma (1/(W km) -- 1400 is for Silicon Nitride
# iBeta2 =  '-200'  # Beta_2 (ps^2/km)
# iBeta3 = '-1'   # Beta_3 (ps^3/km)
# iBeta4 = '0.02'   # Beta_4 (ps^4/km)
# iFibWL = '1550'   # Center WL (nm)


# # Long, underpowered
# iPulse = '0.200' # pulse duration (ps)
# iWave  = '1550'    # pulse central wavelength (nm)
# iEPP   = '50'     # Pulse energy (pJ)
# iGDD   = '0.0'     # Group delay dispersion (ps^2)
# iTOD   = '0.0'     # Third order dispersion (ps^3)
#
# iWindow = '10.0'  # simulation window (ps)
# iSteps  = '50'    # simulation steps
# iPoints = '2**13' # simulation points
#
# iLength = '200'    # fiber length (mm)
# iAlpha  = '0.8'   # attentuation coefficient
# iGamma =  '1400'    # Gamma (1/(W km) -- 1400 is for Silicon Nitride
# iBeta2 =  '-200'  # Beta_2 (ps^2/km)
# iBeta3 = '-1'   # Beta_3 (ps^3/km)
# iBeta4 = '0.02'   # Beta_4 (ps^4/km)
# iFibWL = '1550'   # Center WL (nm)
#
# # Long, underpowered, more loss
# iPulse = '0.200' # pulse duration (ps)
# iWave  = '1550'    # pulse central wavelength (nm)
# iEPP   = '50'     # Pulse energy (pJ)
# iGDD   = '0.0'     # Group delay dispersion (ps^2)
# iTOD   = '0.0'     # Third order dispersion (ps^3)
#
# iWindow = '10.0'  # simulation window (ps)
# iSteps  = '50'    # simulation steps
# iPoints = '2**13' # simulation points
#
# iLength = '200'    # fiber length (mm)
# iAlpha  = '1.4'   # attentuation coefficient
# iGamma =  '1400'    # Gamma (1/(W km) -- 1400 is for Silicon Nitride
# iBeta2 =  '-200'  # Beta_2 (ps^2/km)
# iBeta3 = '-1'   # Beta_3 (ps^3/km)
# iBeta4 = '0.02'   # Beta_4 (ps^4/km)
# iFibWL = '1550'   # Center WL (nm)
#
# # 10 pJ doesn't work
# iPulse = '0.050' # pulse duration (ps)
# iWave  = '1550'    # pulse central wavelength (nm)
# iEPP   = '10'     # Pulse energy (pJ)
# iGDD   = '0.0'     # Group delay dispersion (ps^2)
# iTOD   = '0.0'     # Third order dispersion (ps^3)
#
# iWindow = '10.0'  # simulation window (ps)
# iSteps  = '50'    # simulation steps
# iPoints = '2**13' # simulation points
#
# iLength = '200'    # fiber length (mm)
# iAlpha  = '0.5'   # attentuation coefficient
# iGamma =  '1400'    # Gamma (1/(W km) -- 1400 is for Silicon Nitride
# iBeta2 =  '-200'  # Beta_2 (ps^2/km)
# iBeta3 = '-1'   # Beta_3 (ps^3/km)
# iBeta4 = '0.02'   # Beta_4 (ps^4/km)
# iFibWL = '1550'   # Center WL (nm)



iRaman = True
iSteep = True

auto_run = True  # determines if the simulation should be run upon running the script


def drawnow():
    fig1.canvas.draw()
    plt.pause(0.00001)

def dB(num):
    return 10*np.log10(np.abs(num)**2)
    
    
fig1 = plt.figure(figsize=(13,8))

axL1 = plt.subplot2grid((3,4), (2,0))
axL2 = plt.subplot2grid((3,4), (2,1),sharex=axL1)
axL3 = plt.subplot2grid((3,4), (1,1),sharex=axL1)
axL4 = plt.subplot2grid((3,4), (0,1),sharex=axL1)

ax1 = plt.subplot2grid((3,4), (0,2),sharex=axL1)
ax2 = plt.subplot2grid((3,4), (0,3))

ax3 = plt.subplot2grid((3,4), (1,2),rowspan=2,sharex=axL1)
ax4 = plt.subplot2grid((3,4), (1,3),rowspan=2,sharey=ax3,sharex=ax2)


for label in ax4.get_yticklabels():
    label.set_visible(False)



plt.subplots_adjust(left=0.05,bottom=0.09,right=0.96,top=0.96,wspace=0.3,hspace=0.29)

fig2 = plt.figure(figsize=(5,9))

# fig1.canvas.manager.window.move(500,200)
# fig2.canvas.manager.window.move(0,200)



###Parameter axes###
sl = 0.5
sw = 0.3
sh = 0.030
sv = 0.97 #start vertical
ss = -(sh+0.006) #spacing

ax01 = plt.axes([sl, sv+1*ss, sw, sh])
ax02 = plt.axes([sl, sv+2*ss, sw, sh])
ax03 = plt.axes([sl, sv+3*ss, sw, sh])
ax04 = plt.axes([sl, sv+4*ss, sw, sh])
ax05 = plt.axes([sl, sv+5*ss, sw, sh])

ax06 = plt.axes([sl, sv+7*ss, sw, sh])
ax07 = plt.axes([sl, sv+8*ss, sw, sh])
ax08 = plt.axes([sl, sv+9*ss, sw, sh])

ax09 = plt.axes([sl, sv+11*ss, sw, sh])
ax10 = plt.axes([sl, sv+12*ss, sw, sh])
ax11 = plt.axes([sl, sv+13*ss, sw, sh])
ax12 = plt.axes([sl, sv+14*ss, sw, sh])
ax13 = plt.axes([sl, sv+15*ss, sw, sh])
ax14 = plt.axes([sl, sv+16*ss, sw, sh])
ax15 = plt.axes([sl, sv+17*ss, sw, sh])


ax16 = plt.axes([sl-0.3, sv+19*ss, sw*1.8, sh*2])

ax17 = plt.axes([sl, sv+22*ss, sw, sh])
ax18 = plt.axes([sl, sv+23*ss, sw, sh])
ax19 = plt.axes([sl, sv+24*ss, sw, sh])

bPulse = TextBox(ax01,'Pulse duration (ps)',   initial=iPulse)
bWave  = TextBox(ax02,'Center wavelength (nm)',initial=iWave); bWave.label.set_color('m')
bEPP   = TextBox(ax03,'Energy per pulse (pJ)', initial=iEPP)
bGDD   = TextBox(ax04,'GDD', initial=iGDD)
bTOD   = TextBox(ax05,'TOD', initial=iTOD) 

bWindow= TextBox(ax06,'Time window (ps)',  initial=iWindow) 
bSteps = TextBox(ax07,'Steps',             initial=iSteps) 
bPoints= TextBox(ax08,'Simulation points', initial=iPoints)

bLength= TextBox(ax09,'Fiber length (mm)', initial=iLength)
bAlpha = TextBox(ax10,'Loss (alpha) (dB/cm)',initial=iAlpha)

bGamma = TextBox(ax11,'Gamma (1/(W km))',  initial=iGamma) 
bBeta2 = TextBox(ax12,'Beta_2 (ps^2/km)',  initial=iBeta2) 
bBeta3 = TextBox(ax13,'Beta_3 (ps^3/km)',  initial=iBeta3)
bBeta4 = TextBox(ax14,'Beta_4 (ps^4/km)',  initial=iBeta4) 
bFibWL = TextBox(ax15,'Center WL (nm)',    initial=iFibWL); bFibWL.label.set_color('g')

bCheck = CheckButtons(ax16,('Include Raman response', 'Include self steepening'),(iRaman,iSteep))

bScale = Button(ax17, 'Rescale plots')
bCalc  = Button(ax18, 'Calc Dispersion')
bRun   = Button(ax19, 'Run Simulation')

buttons = (bPulse, bWave, bEPP, bWindow, bSteps, bPoints, bLength, bGamma, bBeta2, bBeta3, bFibWL, bAlpha)


axL1.set_xlabel('Frequency (Hz)')
axL1.set_ylabel('D (ps/nm/km)')
axL2.set_xlabel('Frequency (Hz)')
axL2.set_ylabel(r'$\beta_2$ (ps$^2$/km)')
axL3.set_ylabel('FWM phase mismatch (1/m)')

ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('Intensity (dB)')
ax2.set_xlabel('Time (ps)')


lineL1,  = axL1.plot(0,0) # make some dummy lines to hold future data
lineL2,  = axL2.plot(0,0)
lineL3a, = axL3.plot(0,0)
lineL3b, = axL3.plot(0,0)

lineL4, = axL4.plot(0,0)

line1a, = ax1.plot(0,0,color='b') 
line1b, = ax1.plot(0,0,color='r',lw=1.5) 

line2a, = ax2.plot(0,0,color='b') 
line2b, = ax2.plot(0,0,color='r',lw=1.5) 

vline1a = axL1.axvline(0,color='m',alpha=0.5)
vline1b = axL1.axvline(0,color='g',lw=1.5,alpha=0.5)
vline2a = axL2.axvline(0,color='m',alpha=0.5)
vline2b = axL2.axvline(0,color='g',lw=1.5,alpha=0.5)

hline3 =  axL3.axhline(0,color='k',lw=1.5,alpha=0.5)

axL1.axhline(0,alpha=0.5,color='k')
axL2.axhline(0,alpha=0.5,color='k')



def run_simulation(caller=None,prop=False):
    bRun.label.set_text('Setting up...'); drawnow()
    
    # Grab the parameters from the text boxes: 
    pump_pulse_length = float(bPulse.text)
    centerwl          = float(bWave.text) # in nm
    pump_EPP          = float(bEPP.text) * 1e-12
    GDD               = float(bGDD.text)
    TOD               = float(bTOD.text)
    
    window  = float(bWindow.text)
    steps   = int(bSteps.text)
    npoints = int(eval(bPoints.text))
    
    fiber_length = float(bLength.text) * 1e-3
    alpha = np.log((10**(float(bAlpha.text)*0.1)))*100 # convert from dB/cm to 1/m
    gamma = float(bGamma.text)
    beta2 = float(bBeta2.text)
    beta3 = float(bBeta3.text)
    beta4 = float(bBeta4.text)
    fibWL = float(bFibWL.text)
    
    isNotRaman = np.logical_not(bCheck.lines[0][0].get_visible())
    isNotSteep = np.logical_not(bCheck.lines[1][0].get_visible())
    
    
    # set up the pulse parameters
    pulse = SechPulse(1, pump_pulse_length, centerwl, time_window=window,
                        GDD=GDD, TOD=TOD, NPTS=npoints, frep_MHz=100, power_is_avg=False)

    pulse.set_epp(pump_EPP)

    fiber1 = fiber.FiberInstance() 
    fiber1.generate_fiber(fiber_length,center_wl_nm=fibWL,betas=(beta2,beta3,beta4),gamma_W_m=gamma*1e-3,gvd_units='ps^n/km',gain=-alpha)
    
    
    # fiber1.load_from_db(fiber_length, 'dudley')
    
    print 'Parameters\n------------------'
    print '--Pulse--'
    print 'Duration (ps)    %f'%pump_pulse_length
    print 'Center wl (nm)   %g'%centerwl
    print 'Pulse energy(pJ) %g'%(pump_EPP * 1e12)
    print 'GDD (ps^2)       %f'%(GDD)
    print 'TOD (ps^3)       %f'%(TOD)
    print '--Simulation--'
    print 'Window (ps)    %f'%window
    print 'Steps          %i'%steps
    print 'Time points    %i'%npoints
    print '--Fiber--'
    print 'Length (mm)     %f'%(fiber_length*1e3)
    print 'Alpha (1/m)     %f'%(alpha)
    print 'Gamma (W-1km-1) %f'%(gamma)
    print 'Beta2 (ps^2/km) %f'%(beta2)
    print 'Beta3 (ps^3/km) %f'%(beta3)
    print 'Beta4 (ps^4/km) %f'%(beta4)
    
    print 'FiberWL (nm)    %f'%(fibWL)
    
    laser_freq = 3e8/(centerwl*1e-9)
    fiber_freq = 3e8/(fibWL*1e-9)

    W = pulse.W_mks
    F = W/(2*np.pi)
    T = pulse.T_ps
    xW = W[W>0]
    
    
    D = fiber1.Beta2_to_D(pulse)
    beta = fiber1.Beta2(pulse) * 1e3 # convert ps^2/m to ps^2/km
    
    # Plot the dispersion in the left plots
    lineL1.set_data(F[W>0],D[W>0])
    lineL2.set_data(F[W>0],beta[W>0])
    fibWLindex = (np.abs(F[W>0]-fiber_freq)).argmin()

    intBeta1 = scipy.integrate.cumtrapz(beta[W>0])
    intBeta1 = intBeta1 - intBeta1[fibWLindex]
    # lineL3a.set_data(F[W>0][:-1],intBeta1)
    
    intBeta2 = scipy.integrate.cumtrapz(intBeta1)
    intBeta2 = intBeta2 - intBeta2[fibWLindex]
    lineL3a.set_data(F[W>0][:-2],intBeta2)
    
    
    # dW = W-laser_freq*2*np.pi
    # beta2_SI = beta2*1e3*1e-24 # convert ps^2/km to s^2/m
    # beta4_SI = beta4*1e3*1e-48 # convert ps^4/km to s^4/m
    # phase_mismatch = beta2_SI*dW**2 + 1/12.*beta4_SI*dW**4
    # phase_mismatch2 = phase_mismatch + 2*gamma*pump_EPP*1e-3
    #
    # lineL3a.set_data(F[W>0],phase_mismatch)
    # lineL3b.set_data(F[W>0],phase_mismatch2)
    
    
    
    
    # plot the pulse in the top plots
    line1a.set_data(F[W>0],dB(pulse.AW[W>0]))
    line2a.set_data(T,dB(pulse.AT))
    

    vline1a.set_xdata((laser_freq,laser_freq))
    vline1b.set_xdata((fiber_freq,fiber_freq))
    vline2a.set_xdata((laser_freq,laser_freq))
    vline2b.set_xdata((fiber_freq,fiber_freq))
    
    
    if prop==False: 
        bRun.label.set_text('Run simulation'); drawnow()
        return # stop here unless we've actually pressed the "Run" button
    
   
    ######################## Propagation
    t = time.time()
    bRun.label.set_text('Running...'); drawnow()
    
    # set up the propagation parameters
    evol = SSFM.SSFM(local_error=0.001, USE_SIMPLE_RAMAN=True, disable_Raman = isNotRaman, disable_self_steepening = isNotSteep)
    
    # propagate the pulse!
    y, AW, AT, pulse1 = evol.propagate(pulse_in=pulse, fiber=fiber1, n_steps=steps)
          
    zW_in = np.transpose(AW)[:,(W>0)]
    zT_in = np.transpose(AT)
    zW = dB(zW_in)
    zT = dB(zT_in)

    bRun.label.set_text('Drawing plots...'); drawnow()
    
    line1b.set_data(F[F>0],zW[-1])
    line2b.set_data(T,zT[-1])
    
    
    for ax in (ax3,ax4):
        ax.clear()
    
    extent = (np.min(F[F>0]),np.max(F[F>0]),np.min(y),np.max(y[:-1]))
    ax3.imshow(zW, extent=extent, vmin = np.max(zW) - 60.0, vmax = np.max(zW),aspect='auto',origin='lower')
    ax3.set_xlabel('Frequency (Hz)')
    ax3.set_ylabel('Distance (m)')
    
    
    extent = (np.min(T),np.max(T),np.min(y),np.max(y[:-1]))
    ax4.imshow(zT, extent=extent, vmin = np.max(zT) - 60.0, vmax = np.max(zT),aspect='auto',origin='lower')
    ax4.set_xlabel('Delay (ps)')
    
    
    print 'Total time: %.2f sec'%(time.time()-t)
    
    reset_plots()
    bRun.label.set_text('Run simulation'); drawnow()


    
def reset_plots(caller=None):
    for ax in (axL1,axL2,axL3,axL4,ax1,ax2): # rescale the x and ylims
        ax.relim(); ax.autoscale_view()
        
    ax1.set_ylim(-40,10)
    ax2.set_ylim(0,50)
    
    axL1.set_xlim(100e12,300e12)
    
    for ax in (axL1,axL2,axL3): # rescale the x and ylims
        autoscale_y(ax,linenum=0)

    drawnow()
    


def run_simulation_full(caller=None):
    run_simulation(prop=True)


# for but in buttons:
#     but.on_submit(run_simulation)


bScale.on_clicked(reset_plots)
bCalc.on_clicked(run_simulation)
bRun.on_clicked(run_simulation_full)



run_simulation(prop=False)

reset_plots()

if auto_run==True:
    run_simulation_full()
    reset_plots()


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
