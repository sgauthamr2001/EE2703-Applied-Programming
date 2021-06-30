""" 
EE2703 Applied Programming Lab - Jan-May 2021 
Purpose : Assignment-8 - Submission, 
Date    : 26th April 2021, 
Author  : Sai Gautham Ravipati (EE19B053), 
Input   : None, 
Output  : Plots of Q1,Q2,Q3,Q4,Q5 
Usage   : python3 EE2703_ASSIGN8_EE19B053.py 
"""
"""
Importing necessary modules
"""
import numpy as np 
import scipy.signal as sp
import sympy as sy
import matplotlib.pyplot as plt
import os  
import warnings
warnings.filterwarnings("ignore")

os.makedirs('EE2703_ASN8', exist_ok = True)          #Creating a directory to store plots and data

"""
Set of function calls used elsewhere in the code including Low pass, High pass transfer function, 
conversion from sympy to scipy, step response of filters.
""" 
def lp_tf(R1,R2,C1,C2,G,Vi):

    """
    Solving matrix form of nodal analysis equation using sympy and the transfer function of 
    low pass filter is returned. 
    Inputs - Resistances  -  R1,R2 
             Capacitances -  C1,C2
             Op-amp Gain  -  G 
             Imput        -  Vi 
    Output - Transfer function when Vi is 1 
    """ 
    s = sy.symbols('s')    
    A = sy.Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])    #Conductance matrix
    b = sy.Matrix([0,0,0,-Vi/R1])                                                                    #Source Vector 
    V = A.inv()*b                                                                                    #Inversion gives node voltages 
    return V[3]                                                                                      
    
def hp_tf(R1,R2,C1,C2,G,Vi):
    
    """
    Solving matrix form of nodal analysis equation using sympy and the transfer function of 
    high pass filter is returned. 
    Inputs - Resistances  -  R1,R2 (R2 is R3 given in Circuit)
             Capacitances -  C1,C2
             Op-amp Gain  -  G 
             Imput        -  Vi 
    Output - Transfer function when Vi is 1 
    """ 
    s = sy.symbols('s')
    A = sy.Matrix([[0,0,1,-1/G],[-s*R2*C2/(1+ s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-s*C1-s*C2,s*C2,0,1/R1]])    #Conductance matrix 
    b = sy.Matrix([0,0,0,-Vi*s*C1])                                                                         #Source vector  
    V = A.inv()*b                                                                                           #Inversion gives node voltages 
    return V[3]

def sy_to_sp(H_sy, s = sy.symbols('s')):

    """
    Converts a polynomial transfer function in sympy to numerator and denominator which can 
    be fed into sp.lti.
    """
    n, d = sy.simplify(H_sy).as_numer_denom()            #Expressions of numerator and denominator   
    p_nd = sy.Poly(n, s), sy.Poly(d, s)                  #Converting these expressions to polynomials 
    c_nd = [p.all_coeffs() for p in p_nd]                #Reading coefficients from polynomials 
    Hl, Hd = [sy.lambdify((),c)() for c in c_nd]         #Changing these to float quantities  
    return Hl,Hd
    
def step_response(n,d,t):

    """
    Given the coefficient vector of numerator, denominator of the transfer function and time
    range, returns the step response in time domain. 
    """ 
    H_sp = sp.lti(n,d)                        #Transfer function
    F_sp = sp.lti(n, np.polymul(d,[1,0]))     #Step Response in laplace domain
    _,x  = sp.impulse(F_sp,T = t)             #Step response in time domain 
    return x  

def plot_graph(x, y, xl, yl, pt, fn, flag) : 

    """
    Inputs - x  - datapoints of x-axis 
             y  - datapoints of y-axis 
             xl - x_label
             yl - y_label
             pt - plot_title
             fn - file name to be stored
             flag - 0 : linear plot 
                    1 : semilog plot
    """
    plt.figure()
    if flag == 0 : 
        plt.plot(x,y)
    elif flag == 1 : 
        plt.semilogx(x,y) 
    else :
        print('Only 0,1 are allowed for flag')
    plt.xlabel(xl,fontsize = 12)                
    plt.ylabel(yl, fontsize = 12)
    plt.title(pt, fontsize = 14)
    plt.grid('True')
    plt.savefig(fn)  
    return 0
    
print('Generating plots, please wait.') 

"""
Lines of code corresponding to Q1,2 (Low pass filter) 
"""
Vo_lp = lp_tf(1e4,1e4,1e-9,1e-9,1.586,1)     #Low pass transfer function in terms of symbolic variable s 
n_lp,d_lp = sy_to_sp(Vo_lp)                  #Obtaining the coefficients of numerator and denominator to be fed into sp.lti 
H_lp  = sp.lti(n_lp,d_lp)                    #Transfer function obtained using scipy.lti 

#---------------Q1---------------
t_1  = np.arange(0,1e-3,1e-6)                                                                    #Time range for which step response is computed 
x_1  = step_response(n_lp,d_lp,t_1)                                                              #Obtaining step response using function defined previously  
plot_graph(t_1,x_1,'t (sec)','Vo','Step Response of Low Pass Filter','EE2703_ASN8/Q1.jpg',0)     #Plotting step response vs time   

#---------------Q2---------------   
t_2  = np.arange(0,1e-2,1e-7)                                                                                       #Time range for which response to sum of sinusoids is computed 
Vi_2 = (np.sin(2e3*np.pi*t_2) + np.cos(2e6*np.pi*t_2))*(t_2>=0)                                                     #Variable holding the sinusodial function computed at values of t_2
t_2,x_2,_ = sp.lsim(H_lp,Vi_2,T = t_2)                                                                              #Using sp.lism to obtain the reponse to sum of sinusoids 
plot_graph(t_2,x_2,'t (sec)','Vo','Response of LPF to sum of 2 kHz and 2 MHz sinusoids','EE2703_ASN8/Q2.jpg',0)     #Plotting the function response vs time                            

"""
Lines of code corresponding to Q3,4,5 (High pass filter)
"""
Vo_hp = hp_tf(1e4,1e4,1e-9,1e-9,1.586,1)     #High pass transfer function in term of symbolic variable s  
n_hp,d_hp = sy_to_sp(Vo_hp)                  #Obtaining the coefficients of numerator and denominator to be fed into sp.lti 
H_hp  = sp.lti(n_hp,d_hp)                    #Transfer function obtained using scipy.lti 

#---------------Q3---------------
wrange  = np.logspace(0,8,801)                                                                                  #Range of w for which bode plot is plotted 
w,S,phi = H_hp.bode(w = wrange)                                                                                 #Using .bode to get the phase and the magnitude     
plot_graph(w,S,'w','|H(jw)|','Magnitude Response of High pass filter','EE2703_ASN8/Q3_mag.jpg',1)               #Plotting the bode magnitude plot 
plot_graph(w,phi,'w','$\phi$ in degrees','Phase Response of High pass filter','EE2703_ASN8/Q3_phase.jpg',1)     #Plotting the phase plot  

#---------------Q4---------------
t_4a  = np.arange(0,1e-3,1e-8)                                                                                       #Time range for which response is obtained 
Vi_4a = np.cos(2*np.pi*1e7*t_4a)*np.exp(-1000*t_4a)*(t_4a>=0)                                                        #Vi = cos(2*pi*1e7*t)*exp(-1000*t)*u(t)
t_4a,x_4a,_ = sp.lsim(H_hp,Vi_4a,T = t_4a)                                                                           #Using sp.lsim to get the response 
plot_graph(t_4a,x_4a,'t (sec)','Vo','Response of HPF to $cos(2\pi{t} x 10^7)e^{-1000t}$','EE2703_ASN8/Q4a.jpg',0)    #Plotting the response of decayed sinusoid 

t_4b  = np.arange(0,1e-3,1e-6)                                                                                       #Time range for which response is obtained 
Vi_4b = np.cos(2*np.pi*1e3*t_4b)*np.exp(-0.1*t_4b)*(t_4b>=0)                                                         #Vi = cos(2*pi*1e3*t)*exp(-0.1*t)*u(t)
t_4b,x_4b,_ = sp.lsim(H_hp,Vi_4b,T = t_4b)                                                                           #Using sp.lsim to get the response   
plot_graph(t_4b,x_4b,'t (sec)','Vo','Response of HPF to $cos(2\pi{t} x 10^3)e^{-0.1t}$','EE2703_ASN8/Q4b.jpg',0)     #Plotting the response of decayed sinusoid 

#---------------Q5---------------
t_5  = np.arange(0,1e-3,1e-6)                                                                     #Time range for which step response is computed
x_5  = step_response(n_hp,d_hp,t_5)                                                               #Obtaining step response using function defined previously 
plot_graph(t_5,x_5,'t (sec)','Vo','Step Response of High Pass Filter','EE2703_ASN8/Q5.jpg',0)     #Plotting step response vs time

print('All plots have been saven in EE2703_ASN8 directory')

"""
---------------------------------------End-Of-The-Code-----------------------------------------
"""
