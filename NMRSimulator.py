import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
from scipy import *
from scipy.integrate import quad
from matplotlib.widgets import TextBox, Button
from matplotlib.gridspec import GridSpec


import os.path
from os import path

from Variables import Variables as var
pi = math.pi


def generate_spectrum():
    """
    This will create a user interface that plots FID and NMR spectra
    ----------------------------------------------------------------
    Parameters:
        Lambda: proportionality constant between magnetic moment decay and signal
        M Naught: initial magnitude of magnetic moment of sample
        Alpha: angle of r.f. pulse
        Omega: Phase difference between r.f. pulse ans signal
        R2: spin-spin relaxation rate constant
    ----------------------------------------------------------------
    Function will create txt file with parameters of last run to use as initial 
    values of next run, and call plot function.
    """
    filename = "parameters.txt"
    #If the file exists already, we should read the parameters from memory. 
    #If the file does not exist, default variables are in Variables.py
    if path.exists(filename): 
        para= open(filename,"r")
        para_lines = para.readlines()
        var.LAMBDA = float(para_lines[0])
        var.M0 = float(para_lines[1])
        var.ALPHA = float(para_lines[2])
        var.OMEGA = float(para_lines[3])
        var.R2 = float(para_lines[4])
    
    plot()

    para = open(filename,"w") #Save the variables that the user input to memory
    para.write(str(var.LAMBDA)+"\n")
    para.write(str(var.M0)+"\n")
    para.write(str(var.ALPHA)+"\n")
    para.write(str(var.OMEGA)+"\n")
    para.write(str(var.R2)+"\n")
    para.close()

def signal_of_time(t):
    #Mathematical function plotted to create FID graph
    return var.LAMBDA*var.M0*np.sin(var.ALPHA)*np.exp(1j*var.OMEGA*t-var.R2*t) 

def signal_integrated(t):
    #Mathematical function integrated over time to create NMR spectra
    return var.LAMBDA*var.M0*np.sin(var.ALPHA)*np.exp(1j*var.OMEGA*t-var.R2*t)*np.exp(-1j*var.omega*t)

def complex_quadrature(func, a, b, **kwargs):
    #This function returns the integration over the real and complex number spaces
    def real_func(x):
        return np.real(func(x))
    def imag_func(x):
        return np.imag(func(x))
    real_integral = quad(real_func, a, b, **kwargs)
    imag_integral = quad(imag_func, a, b, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])

def integrate(omega_vals):
    #Integrates signal_integrated for all values of omega_vals 
    #Returns an array of numbers that signify signal strength 
    integrals = np.empty(len(omega_vals),dtype=complex)
    count = 0
    for omega in omega_vals:
        var.omega = omega
        integral = complex_quadrature((signal_integrated),0,np.inf)
        integrals[count] = integral[0]
        count+=1
    return integrals

#In order for buttons to call plot, plot must take in an event. We do not use the event though.
def plot(event = nan):
    #Plot two graphs, one with time and one with omega_vals
    t = np.arange(0, 5, 0.1)
    omega_vals = np.arange(-5, 5, 0.1)

    signal1 = signal_of_time(t) #Get array of y values for FID graph
    signal2 = integrate(omega_vals) #Get array of y values function for NMR spectra

    #Change the graph titles depending on whether it's on or off resonance
    res = ": Off Resonance"
    if var.OMEGA == 0:
        res = ": On Resonance"

    #Graph the FID on the top left of the window
    plt.subplot2grid((2, 3), (0, 0), colspan=2)
    plt.plot(t, np.real(signal1), label = 'real')
    plt.plot(t, np.imag(signal1), "r", label = "imaginary" )
    plt.title('FID' + res)
    plt.xlabel('time (s)')
    plt.legend()

    #Graph the NMR spectra on the bottom left of the window
    plt.subplot2grid((2, 3), (1, 0),colspan=2)
    plt.plot(omega_vals, np.real(signal2), label = 'real' )
    plt.plot(omega_vals, np.imag(signal2), 'r', label = "imaginary" )
    plt.title('NMR Signal' + res)
    plt.xlabel('frequency (Hz)')
    plt.legend()

    #Create the user interface on the right side of the window
    plt.subplot2grid((2, 3), (0, 2), rowspan=2)
    plt.axis("off")
    plt.text(0.2,0.95,"Input Values:")
    button = Button(plt.axes([0.72, 0.07, 0.23, 0.1]), 'Generate Spectrum', color="lightgoldenrodyellow", hovercolor='c')
    button.on_clicked(plot)

    #Create a user text widget for each parameter and set the parameter values to the inputted values
    var_list = [str(var.LAMBDA),str(var.M0),str(round(np.degrees(var.ALPHA),3)),str(round(np.degrees(var.OMEGA),3)),str(var.R2)]
    var_names = ["Lambda: ","M Naught: ","Alpha: ","Omega: ","R2: "]
    separation = 0.0
    text_boxes = []
    for i in range (len(var_list)):
        text_boxes.append(TextBox(plt.axes([.80, 0.7-separation, 0.1, 0.07]),var_names[i], initial = var_list[i]))
        separation += 0.12
    func_list = [set_lambda,set_M0,set_aplha,set_omega,set_R2]
    for i in range (len(var_list)):
        text_boxes[i].on_submit(func_list[i])
    
    plt.subplots_adjust(hspace = 0.6)
    plt.show() 

#Functions to set the values of each parameter to the inputed value of that parameter
#Called on line 132
def set_lambda(text):
    var.LAMBDA = float(text)
def set_M0(text):
    var.M0 = float(text)
def set_aplha(text):
    var.ALPHA = np.radians(float(text))
def set_omega(text):
    var.OMEGA = np.radians(float(text))
def set_R2(text):
    var.R2 = float(text)


generate_spectrum()








