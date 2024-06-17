
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp

def g(r):

    '''
    returns components of the Swartzchild metric gt(r), gr(r), and gphi(r) at a given radial coordinate (r) in the frame of reference of the black hole

    Use as g = g(r) or gt, gr, gphi = g(r) with:
    r       : radial coordinate in Swartzchild radii    [RS]
    '''
    # the metric in geometric units in the theta = pi/2 plane
    gt = 1 - 1 / r
    gr = -1 / (1 - 1 / r)
    gphi = - r ** 2

    return gt, gr, gphi

def dUdtau(tau : float, X : np.ndarray):
    '''
    returns the derivatives of the temporal, radial, and azimuthal components of the four position (ct, r, phi) and velocity (Ut, Ur, Uphi) with respect to proper time (tau) in a plane inclined with polar angle theta = pi / 2.

    For use in a suitable integrator such as scipy.solve_ivp as X = solve_ivp(dUdtau, (0, tf), X0), where tf is the endtime and X0 is a tuple of initial conditions as outlined in traj:
    tau     : proper time           [RS / c]
    X       : (c * t, r, phi, Ut, Ur, Uphi)
        t       : temporal coordinate   [RS / c]
        r       : radial coordinate     [RS]
        phi     : azimuthal coordinate  [rad]
        Ut      : temporal velocity     [c]
        Ur      : radial velocity       [c]
        Uphi    : azimuthal velocity    [rad c / RS]
    '''
    r, Ut, Ur, Uphi = X[1], *X[3:6]
    
    # the geodesic equations in geometric units in the theta = pi/2 plane
    gt = 1 - 1 / r
    dUtdtau = - 1 / (r ** 2 * (1 - 1 / r)) * Ur * Ut
    dUrdtau = 1 / (2 * r ** 2 * (1 - 1 / r)) * Ur ** 2 - 1 / (2 * r ** 2) * (1 - 1 / r) * (Ut) ** 2 + (r - 1) * Uphi ** 2
    dUphidtau = - 2 / r * Ur * Uphi

    return Ut, Ur, Uphi, dUtdtau, dUrdtau, dUphidtau

def horizon(tau : float, X : np.ndarray):
    '''detects if the particle crosses the event horizon at r = RS to 2 d.p.. For use as an event in scipy.solve_ivp
    tau     : proper time           [RS / c]
    X       : (c * t, r, phi, Ut, Ur, Uphi)
    '''
    # events are evaluated as the zeros of functions, therefore returns the difference of the rounded radial coordinate and the schwarzschild radius in geometric units
    return abs(round(X[1], 2) - 1)

def trajectory(input_file : str, output_file : str):
    '''
    Sets up and graphs the trajectory of a particle with negligible mass in the Schwarzchild metric starting from an initial four position (R0) and velocity (U0) for a duration of proper time (tauf). returns arrays of proper time (tau), four position (R) and four velocity (U) at a given number of grid points (Ntau). Initial conditions should be setup in an input file. In addition to plotting the trajectory of the particle, the code stores the returned arrays in a given output file
    
    Use as trajectory(input_file, output_file) or tau, R, phi = trajectory(input_file, output_file)
    input_file      : name of text file to read initial conditions from
    output_file     : name of text file to write four positions and velocities to
    
    Input file must contain at least four lines. A line must be formatted as the given name of the quantity followed by an equals symbol and its value(s), comma delimited if necessary. Required quantites are as follows
    R0 = ct0, r0, phi0
        t0       : initial temporal position   [RS / c]
        r0       : initial radial position     [RS]
        phi0     : initial azimuthal position  [rad]
    U0 = Ur0, Uphi0 
        Ur0      : initial radial velocity       [c]
        Uphi0    : initial azimuthal velocity    [rad c / RS]
    Ntau = number of proper time points to evaluate R and U at [none]
    tauf = length of the proper time duration for which to model the trajectory    
    The input file may also contain any of the following on seperate lines. The order of the lines is immaterial. 
    method = integration method for scipy.solve_ivp (default is Radau)
    rtol = relative tolerance for scipy.solve_ivp (default is 1e-12)
    atol = absolute tolerance for scipy.solve_ivp (default is 1e-14)  
    rmax = radial limit of trajectory plot
    '''
    
    # define default values
    rtol = 1e-12
    atol = 1e-14
    method = 'Radau'
    component_names = ['t','r','phi','Ut','Ur','Uphi']
    
    # stop the particle if it reaches the event horizon
    horizon.terminal = True
    
    # reading in inputs as string and split into a list of seperate lines
    input_par = open(input_file, "r")
    inputs = input_par.read().replace(' ', '')
    input_par.close()
    inputlist = inputs.splitlines()
    
    # check each line in the input file for the initial conditions and return them as seperate variables as the required datatype
    for con in inputlist:
        con_type = con[:con.find('=')]
        match con_type:
            case 'R0': 
                con = con.replace('R0=','')
                R0 = np.array(con.split(',')).astype(float)
            case 'U0':
                con = con.replace('U0=','')
                U0 = np.array(con.split(',')).astype(float)
            case 'Ntau':
                Ntau = int(float(con.replace('Ntau=','')))
            case 'tauf':
                tauf = float(con.replace('tauf=',''))
            case 'method':
                method = con.replace('method=','')
            case 'rtol':
                rtol = float(con.replace('rtol=',''))
            case 'atol':
                rtol = float(con.replace('atol=',''))
            case 'rmax':
                rmax = float(con.replace('rmax=',''))

    # check that the neccesary input conditions are present in the file
    try: r0 = R0[1]
    except: raise Exception('Input file must contain a line with initial position \'R0 = t0, r0, phi0\'')
    try: Ur0, Uphi0 = U0
    except: raise Exception('Input file must contain a line with initial velocity \'U0 = Ur0, Uphi0\'')
    try: tauf
    except:  raise Exception('Input file must contain a line with proper time duration \'tauf = ')   
    try: tpoints = np.linspace(0, tauf, Ntau) # setup a grid of Ntau points over the proper time period [0, tauf]
    except:  raise Exception('Input file must contain a line with proper time duration \'Ntau = ')   
    
    # calculate the time component of the four velocity using the magnitude identity and the values of the schwarzschild metric at r0
    gt0, gr0, gphi0 = g(r0)
    Ut0 = np.sqrt((1 - gr0 * Ur0 ** 2 - gphi0 * Uphi0 ** 2) / gt0)
    U0 = np.insert(U0, 0, Ut0)
    
    # perform the integration using the initial conditions
    X = solve_ivp(dUdtau, (0, tauf), (*R0, *U0), method = method, t_eval = tpoints, rtol = rtol, atol = atol, events = horizon)
    
    # check if the integration ends due to the particle reaching the event horizon
    if X.status == 1: print('The test particle reaches the event horizon after proper time %.2f RS / c' % round(X.t[-1],2))

    # convert the returned arrays of proper time, four position and four velocity into strings to write to the output file
    out_string = 'tau' + ' = ' + ', '.join(X.t.astype(str)) + '\n'
    for y, f_type in zip(X.y, component_names):
        out_string += f_type + ' = ' + ', '.join(y.astype(str)) + '\n'

    # writes the proper time, four position and four velocity to the output file
    output_par = open(output_file,'w')
    output_par.write(out_string)

    # plots the azimuthal angular and radial coordinates on an 8*8 polar plot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize = (8,8))
    ax.plot(X.y[2] , X.y[1])

    # checks if rmax was given as an input and if so sets the maximum radius of the polar plot
    try: ax.set_rmax(rmax)
    except: pass

    # labels the radial direction on the polar plot
    label_position = ax.get_rlabel_position()
    ax.text(np.radians(label_position - 10), ax.get_rmax() / 2,'Orbital radius [$R_\\text{S}$]', rotation = label_position, ha = 'center',va = 'center')

    # returns arrays of proper time, four position and four velocity
    return X.t, X.y[0:3], X.y[3:6]

    
print('this is a program to simulate the trajectory of a particle with negligable mass in the vicinity of a black hole')

input_file = input('''to begin the simulation, create an input text file containing at least four lines. Each line must contain a seperate quantity and be be formatted as the given name of the quantity followed by an equals symbol and its value(s), comma delimited if necessary. Required quantites are as follows
    R0 = ct0, r0, phi0
        t0       : initial temporal position   [RS / c]
        r0       : initial radial position     [RS]
        phi0     : initial azimuthal position  [rad]
    U0 = Ur0, Uphi0 
        Ur0      : initial radial velocity       [c]
        Uphi0    : initial azimuthal velocity    [rad c / RS]
    Ntau = number of proper time points to evaluate R and U at [none]
    tauf = length of the proper time duration for which to model the trajectory    
    
    The input file may also contain any of the following on seperate lines. The order of the lines is immaterial. 
    method = integration method for scipy.solve_ivp (default is Radau)
    rtol = relative tolerance for scipy.solve_ivp (default is 1e-12)
    atol = absolute tolerance for scipy.solve_ivp (default is 1e-14)  
    rmax = radial limit of trajectory plot

    name of input file:\n''')

output_file = input('\nthe simulation also requires the name of an output text file to write the proper times, four positions and velocities to. If no file already has the given filename, then one will be created.\n name of output file:\n')

    



