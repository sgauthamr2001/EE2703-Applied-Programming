""" 
EE2703 Applied Programming Lab - Jan-May 2021 
Purpose : Assignment-7 - Submission, 
Date    : 16th April 2021, 
Author  : Sai Gautham Ravipati (EE19B053), 
Input   : None, 
Output  : Plots of Time Domain response for Q1,2,3 and 6, solved equations of Q4, Bode plots for Q5, 
Usage   : python3 EE2703_ASSIGN7_EE19B053.py 
"""
"""
Importing necessary modules
"""
import numpy as np 
import scipy.signal as sp
import matplotlib.pyplot as plt
import os  

os.makedirs('EE2703_ASN7', exist_ok = True)          #Creating a directory to store plots and data
 
"""
Code corresponding to Q1 - decay 0.5
"""
s1 = np.polymul([1,0.5],[1,0.5])                     #Forming the polynomial (s + 0.5)^2
s2 = np.polyadd(s1,[2.25])                           #Forming the polynomial (s + 0.5)^2 + 2.25
d1 = np.polymul(s2,[1,0,2.25])                       #Forming the polynomial ((s + 0.5)^2 + 2.25)(s^2 + 2.25)
H = sp.lti([1,0.5],d1)                               #Constructing H as (s + 0.5)/[((s + 0.5)^2 + 2.25)(s^2 + 2.25)]      
t,x = sp.impulse(H,T = np.linspace(0,50,1250))       #Solving for x, for 0 < t < 50
plt.figure()
plt.plot(t,x)                                        #Plotting x vs t 
plt.xlabel('t (sec)',fontsize = 12)                
plt.ylabel('x',fontsize = 12)
plt.title('Time Response of spring system for f = 1.5, d = 0.5',fontsize = 14)
plt.grid('True')                                  
plt.savefig('EE2703_ASN7/Q1.jpg')                    #Saving the plot in directory 

"""
Code corresponding to Q2 - decay 0.05
"""
s1 = np.polymul([1,0.05],[1,0.05])                   #Forming the polynomial (s + 0.05)^2
s2 = np.polyadd(s1,[2.25])                           #Forming the polynomial (s + 0.05)^2 + 2.25
d1 = np.polymul(s2,[1,0,2.25])                       #Forming the polynomial ((s + 0.05)^2 + 2.25)(s^2 + 2.25)
H = sp.lti([1,0.05],d1)                              #Constructing H as (s + 0.05)/[((s + 0.05)^2 + 2.25)(s^2 + 2.25)]      
t,x = sp.impulse(H,T = np.linspace(0,50,1250))       #Solving for x, for 0 < t < 50
plt.figure()
plt.plot(t,x)                                        #Plotting x vs t 
plt.xlabel('t (sec)',fontsize = 12)                    
plt.ylabel('x',fontsize = 12)
plt.title('Time Response of spring system for f = 1.5, d = 0.05',fontsize = 14)
plt.grid('True')
plt.savefig('EE2703_ASN7/Q2.jpg')                    #Saving the plot in directory 

"""
Code corresponding to Q3 - decay 0.05, f 1.4-1.6
"""
for f in np.arange(1.4,1.65,0.05):
    H = sp.lti(1,[1,0,2.25])                         #Transfer function considering the system to be LTI, 1/(s^2 + 0.05)
    t = np.linspace(0,100,2500)                      #Considering 0 < t < 100
    f_t = np.cos(f*t)*np.exp(-0.05*t)*(t>0)          #Constructing f(t) of the form given in question  
    t,x,svec = sp.lsim(H,f_t,t)                      #Performing CT convoluton between impusle response and function
    plt.figure()
    plt.plot(t,x)                                    #Plotting x vs t 
    plt.xlabel('t (sec)',fontsize = 12)              
    plt.ylabel('x',fontsize = 12)
    plt.grid(True)                        
    plt.title('Time Response of spring system for f = {}, d = 0.05'.format(f),fontsize = 14)
    plt.savefig('EE2703_ASN7/Q3_' + str(f) + '.jpg')  #Saving the plot in directory 

"""
Code corresponding to Q4
""" 
Hx = sp.lti([1,0,2],[1,0,3,0])                       #Hx is (s^2 + 2)/(s^3 + 3s) 
Hy = sp.lti([2],[1,0,3,0])                           #Hy is 2/(s^3 + 3s) 
t,x = sp.impulse(Hx,T = np.linspace(0,20,1000))      #Solving for x      
t,y = sp.impulse(Hy,T = np.linspace(0,20,1000))      #Solving for y 
plt.figure()
plt.plot(t,x,label = 'x')                            #Plotting x vs t
plt.plot(t,y,label = 'y')                            #Plotting y vs t 
plt.xlabel('t',fontsize = 12)                        
plt.ylabel('x and y',fontsize = 12)
plt.title('Solution of coupled spring problem', fontsize = 14)
plt.legend(loc = 'upper right')
plt.grid(True)
plt.savefig('EE2703_ASN7/Q4.jpg')

"""
Code corresponding to Q5 
"""
R = 100                                              #Value of Resistance 
L = 1e-6                                             #Value of Inductance 
C = 1e-6                                             #Value of Capacitance 
H = sp.lti([1],[L*C,R*C,1])                          #Tranfer function is 1/(LCs^2 + RCs + 1)
w,S,phi = H.bode()                                   #Using .bode to obtain magnitude and phase 
plt.figure()
plt.semilogx(w,S)                                    #Plotting magnitude plot 
plt.xlabel('w',fontsize = 12)
plt.ylabel('|H(jw)|',fontsize = 12)
plt.title('Magnitude Response of $H(jw)$', fontsize = 14)
plt.grid(True)
plt.savefig('EE2703_ASN7/Q5_mag.jpg')                #Saving the plot in directory
plt.figure() 
plt.semilogx(w,phi)                                  #Plotting the phase plot 
plt.xlabel('w',fontsize = 12)
plt.ylabel('$\phi$ in degrees',fontsize = 12)
plt.title('Phase Response of $H(jw)$', fontsize = 14)
plt.grid(True)
plt.savefig('EE2703_ASN7/Q5_phase.jpg')              #Saving the plot in directory 

"""
Code corresponding to Q6
"""
vi = lambda t : (np.cos(1e3*t)-np.cos(1e6*t))*(t>0)  #Input to the RLC circuit 
t1  = np.arange(0,1e-2,1e-5)                         #0 < t1 < 10ms 
t2  = np.arange(0,3e-5,3e-8)                         #0 < t2 < 30us
t1,y1,_ = sp.lsim(H,vi(t1),t1)                       #Performing convolution for impulse response for both the cases
t2,y2,_ = sp.lsim(H,vi(t2),t2)
plt.figure()
plt.plot(t1,y1)                                      #Plotting y1 vs t1                   
plt.xlabel('t (sec)',fontsize = 12)
plt.ylabel('Vc (in V)',fontsize = 12)
plt.title('Response of RLC circuit for $v_1(t)$ for 0 < t < 10ms')
plt.grid(True)
plt.savefig('EE2703_ASN7/Q6_T1.jpg')                 #Saving the plot in directory 
plt.figure()
plt.plot(t2,y2)                                      #Plotting y2 vs t2
plt.xlabel('t (sec)',fontsize = 12)              
plt.ylabel('Vc (in V)',fontsize = 12)
plt.grid(True)
plt.title('Response of RLC circuit for $v_1(t)$ for 0 < t < 30$\mu$s ')
plt.savefig('EE2703_ASN7/Q6_T2.jpg')                 #Saving the plot in directory 

print('All plots have been saven in EE2703_ASN7 directory')

"""
---------------------------------------End-Of-The-Code-----------------------------------------
"""
