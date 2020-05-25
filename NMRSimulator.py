import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import *
from scipy.integrate import quad

from decimal import Decimal

import os.path
from os import path

from Variables import Variables as var
pi = math.pi



def user_var_input():
    filename = "parameters.txt"
    change_val = "*"
    if path.exists(filename):
        para= open(filename,"r")
        para_lines = para.readlines()
        var.LAMBDA = int(para_lines[0])
        var.M0 = int(para_lines[1])
        var.ALPHA = int(para_lines[2])
        var.OMEGA = int(para_lines[3])
        var.R2 = int(para_lines[4])
        print("------Previous parameters------ \nLambda = " + para_lines[0] + "M naught = " + para_lines[1] + "Alpha = " + para_lines[2] + "Omega = " + para_lines[3] + "R2 = " + para_lines[4] + "")
        change_val = input("---------To change enter-------- \nL for lambda \nM for M naught \nA for alpha \nO for Omega \nR for R2 \n* for all values \nExample: LMO \nenter for no changes: \n").lower()
    if "l" in change_val or "*" in change_val:     
        var.LAMBDA = int(input("Input value for lambda: "))  
    if "m" in change_val or "*" in change_val:    
        var.M0 = int(input("Input value for M naught: "))
    if "a" in change_val or "*" in change_val:       
        var.ALPHA = int(input("Input value for alpha (degrees): ")) 
    if "o" in change_val or "*" in change_val:    
        var.OMEGA = int(input("Input value for (capital) omega (degrees): "))
    if "r" in change_val or "*" in change_val:    
        var.R2 = int(input("Input value for R2: "))
    
    para = open(filename,"w")
    para.write(str(var.LAMBDA)+"\n")
    para.write(str(var.M0)+"\n")
    para.write(str(var.ALPHA)+"\n")
    para.write(str(var.OMEGA)+"\n")
    para.write(str(var.R2)+"\n")
    para.close()
    var.ALPHA = np.radians(var.ALPHA)
    var.OMEGA = np.radians(var.OMEGA)
    
def signal_of_time(t):
    return var.LAMBDA*var.M0*np.sin(var.ALPHA)*np.exp(1j*var.OMEGA*t-var.R2*t) 

def signal_integrated(t):
    return var.LAMBDA*var.M0*np.sin(var.ALPHA)*np.exp(1j*var.OMEGA*t-var.R2*t)*np.exp(1j*var.omega*t)*np.exp(1j*var.omega*t)

def complex_quadrature(func, a, b, **kwargs):
    def real_func(x):
        return np.real(func(x))
    def imag_func(x):
        return np.imag(func(x))
    real_integral = quad(real_func, a, b, **kwargs)
    imag_integral = quad(imag_func, a, b, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])

def integrate(omega_vals):
    integrals = np.empty(len(omega_vals),dtype=complex)
    count = 0
    for omega in omega_vals:
        var.omega = omega
        integral = complex_quadrature((signal_integrated),0,np.inf)
        integrals[count] = integral[0]
        count+=1
    return integrals

def plot_frequency():
   
    omega_vals = np.arange(-5, 5, 0.01)
    signal = integrate(omega_vals)
    fig, ax = plt.subplots()
    ax.plot(omega_vals, np.real(signal), label = 'real')
    ax.plot(omega_vals, np.imag(signal,),'r', label = 'imaginary')
    res = ": Off Resonace"
    if var.OMEGA == 0:
        res = ": On Resonance"
    ax.set(xlabel='Frequency (Hz)', ylabel='', title='NMR Signal'+res)
    ax.legend()

    ax.grid()

    fig.savefig("test.png")
    plt.show()

def plot_time():
   
    t = np.arange(0, 5, 0.01)
    signal = signal_of_time(t)
    fig, ax = plt.subplots()
    ax.plot(t, np.real(signal), label = 'real')
    ax.plot(t, np.imag(signal),'r', label = 'imaginary')
    res = ": Off Resonace"
    if var.OMEGA == 0:
        res = ": On Resonance"
    ax.set(xlabel='Time (s)', ylabel='', title='FID'+res)
    ax.legend()

    ax.grid()

    fig.savefig("test.png")
    plt.show()

user_var_input()
plot_time()
plot_frequency()


