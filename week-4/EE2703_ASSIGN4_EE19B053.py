""" 
EE2703 Applied Programming Lab - Jan-May 2021 
Purpose : Assignment-4 - Submission, 
Date    : 11th March 2021, 
Author  : Sai Gautham Ravipati (EE19B053), 
Input   : None, 
Output  : Different plots for fourier coefficeints, periodic extension, fourier approximations, and
          deviation in coefficients obtained using least squares approach, 
Usage   : python3 EE2703_ASSIGN4_EE19B053.py.
"""

"""
Importing the neccasary modules.
"""
import numpy as np 
from math import pi 
import matplotlib.pyplot as plt
from scipy.linalg import lstsq
from scipy import integrate as integrate

"""
The following lines are function definitons used in the code.
"""
def exp(x):
    """
    Takes either a vector/scalar as input and returns exp(x) applied element wise. 
    Ex : x = [1,2,3] , function returns [exp(1),exp(2),exp(3)]
    """ 
    return np.exp(x)
    
def cos_cos(x) : 
    """
    Takes either a vector/scalar as input and returns cos(cos(x)) applied element wise.
    """
    return np.cos(np.cos(x))

def periodic_extension(x,lower,upper) : 
    """
    Takes in a vector/scalar x, lower and upper limits of one period and translates every element i in x
    to the range [lower,upper) such that f(i) = f(i+a*p) , where a is an integer and p is the period. This
    works for any other range also instead of hard coding [-2*pi,4*pi)
    """ 
    if not(upper>lower):              #Checking if limits of the interval are valid. 
        print('Please give proper range of limits')
        exit()
    p = upper - lower                 #Period of the function 
    x_out = 2*pi*np.ceil(-x/p) + x    #Translating every element to the range [lower,upper)
    return x_out 
        
def coeff_intg(f,n):
    """
    Takes in the name of the function, exp or cos_cos which are defined previously, and the highest 
    index of fourier coefficients to be generated. 
    Ex : 
    coeff_intg(exp,25) returns the coefficient vector (a0,a1,b1,...,a25,b25) for exp(x)
    """
    if n != int(n) or n < 0 :                           #Checking if n is a positive integer. 
        print('Please enter a positive integer value for n')
        exit()
    u = lambda x,k : f(x)*np.cos(k*x)                   #Function definition for u(x,k) = f(x)*cos(kx) 
    v = lambda x,k : f(x)*np.sin(k*x)                   #Function definition for v(x,k) = f(x)*sin(kx) 
    rtnval = np.zeros(2*n+1)                            #Initialising the array to be returned
    rtnval[0] = (integrate.quad(f,0,2*pi)[0])/(2*pi)    #Computing the coefficients using integrate.quad 
    for k in range(1,n+1) :
        rtnval[2*k-1] = integrate.quad(u,0,2*pi,args=(k))[0]/pi
        rtnval[2*k]   = integrate.quad(v,0,2*pi,args=(k))[0]/pi
    return rtnval
    
def vec_lstsq(f,n,x):
    """
    Takes in the name of the function, exp or cos_cos which are defined previously, and the highest 
    index of fourier coeffiecients to be generated and returns the matrices A,b.  
    """
    if n != int(n) or n < 0 :     #Checking if n is a positive integer 
        print('Please enter a integer value for n')
        exit()
    pt  = x.shape[0]              #Number of x datapoints generated  
    b = f(x)                      #Generating the matrix b
    A = np.zeros((pt,2*n+1))      #Initialisng A with zeros 
    A[:,0] = 1                    #Assinging a constant value of n to the first column of A 
    for k in range(1,n+1) :       #Filling the columns of A with alternate cos(kx) and sin(kx) columns
        A[:,2*k-1] = np.cos(k*x)
        A[:,2*k]   = np.sin(k*x)
    return A,b  
    
"""
The following lines of code implements of each of the sub-question in assignment 
"""                      
x_range = np.linspace(-2*pi,4*pi,401)             #Generating 401 points in the range [-2*pi,4*pi]
x_range = x_range[:-1]                            #Removing 4*pi to make x_range as [-2*pi,4*pi)
xp_range = periodic_extension(x_range,0,2*pi)     #Translating all the points to the range [0,2*pi)

"""
The following lines of code semi-log plot exp(x) and its periodic extension for x in [-2*pi,4*pi). For the 
periodic extension exp is taken for xp_range and is plotted vs x_range. 
""" 
plt.figure(1)                                 
plt.semilogy(x_range,exp(x_range), 'r', label = 'Actual value')              #Plotting original curve in red 
plt.semilogy(x_range,exp(xp_range),'--g', label = 'Periodic Extension')      #Plotting periodic extension with green dotted lines 
plt.grid(True)
plt.ylabel(r'$e^{x}$ (log scale) $\rightarrow$',fontsize=12)
plt.xlabel(r'x $\rightarrow$',fontsize=12)
plt.title('Semi-log plot of $e^{x}$ in the interval [$-2\pi$,$4\pi$)',fontsize=14)
plt.legend(loc='upper right')

"""
The following lines of code plot cos_cos(x) and its periodic extension for x in [-2*pi,4*pi). For the 
periodic extension cos_cos is taken for xp_range and is plotted vs x_range. 
""" 
plt.figure(2)                                                                 
plt.plot(x_range,cos_cos(x_range), 'r',label = 'Actual value')               #Plotting original curve in red             
plt.plot(x_range,cos_cos(xp_range),'--g', label = 'Periodic Extension')      #Plotting periodic extension witg green dotted lines 
plt.grid(True)
plt.ylabel(r'$cos(cos(x))\rightarrow$',fontsize=12)
plt.xlabel(r'x$\rightarrow$',fontsize=12)
plt.title('Plot of $cos(cos(x))$ in the interval [$-2\pi$,$4\pi$)',fontsize=14)
plt.legend(loc='upper right')
 
"""
Generating first 51 fourier coefficients using both integration.quad as well as least squares approach.
"""
n_max = 25                                       #Highest index of a,b i.e. in this case {a0,a1,b1,...,a25,b25} are generated. 
pts   = 400                                      #Number of points to be generated in the interval [0,2*pi)
coeff_intg_exp = coeff_intg(exp,n_max)           #Generating exp(x) coefficients using integrate.quad
coeff_intg_ccx = coeff_intg(cos_cos,n_max)       #Generating cos(cos(x)) coefficients using integrate.quad 
x = np.linspace(0,2*pi,pts+1)                    #Generating 401 points in range [0,2*pi]    
x = x[:-1]                                       #Removing the last point, to avoid overlap while periodic extension 
A_exp,b_exp = vec_lstsq(exp,n_max,x)             #Unpacking A_exp,b_exp to be passed to lstsq
A_ccx,b_ccx = vec_lstsq(cos_cos,n_max,x)         #Unpacking A_ccx,b_ccx to be passed to lstsq
coeff_lstq_exp = lstsq(A_exp,b_exp)[0]           #Generating exp(x) coefficients using least squares approach
coeff_lstq_ccx = lstsq(A_ccx,b_ccx)[0]           #Generating cos(cos(x)) coefficients using least squares approach 

"""
Semi-log plot of magnitude of coefficients of exp(x) 
"""
plt.figure(3)
plt.semilogy(np.abs(coeff_intg_exp), 'ro', label = 'Coefficients - Integrated')                        #coeff_intg_exp are plotted in red 
plt.semilogy(np.abs(coeff_lstq_exp), 'go', label = 'Coefficients - Least Squares',markersize = 4.5)    #coeff_lstq_exp are plotted in green and of lower size 
plt.grid(True)  
plt.ylabel(r'Magnitude $\rightarrow$',fontsize=12)
plt.xlabel(r'$n\rightarrow$',fontsize=12)
plt.title('semi-log plot of fourier coeffcients of $e^{x}$',fontsize=14)
plt.legend(loc='upper right')

"""
log-log plot of magnitude of coefficients of exp(x) 
"""
plt.figure(4)
plt.loglog(np.abs(coeff_intg_exp), 'ro', label = 'Coefficients - Integrated')                        #coeff_intg_exp are plotted in red 
plt.loglog(np.abs(coeff_lstq_exp), 'go', label = 'Coefficients - Least Squares',markersize = 4.5)    #coeff_lstq_exp are plotted in green and of lower size 
plt.grid(True)
plt.ylabel(r'Magnitude $\rightarrow$',fontsize=12)
plt.xlabel(r'$n\rightarrow$',fontsize=12)
plt.title('log-log plot of fourier coeffcients of $e^{x}$',fontsize=14)
plt.legend(loc='upper right')

"""
semi-log plot of magnitude of coefficients of cos(cos(x)) 
"""
plt.figure(5)
plt.semilogy(np.abs(coeff_intg_ccx), 'ro', label = 'Coefficients - Integrated')                        #coeff_intg_exp are plotted in red 
plt.semilogy(np.abs(coeff_lstq_ccx), 'go', label = 'Coefficients - Least Squares',markersize = 4.5)    #coeff_lstq_exp are plotted in green and of lower size     
plt.grid(True)
plt.ylabel(r'Magnitude $\rightarrow$',fontsize=12)
plt.xlabel(r'$n\rightarrow$',fontsize=12)
plt.title('semi-log plot of fourier coeffcients of $cos(cos(x))$',fontsize=14)
plt.legend(loc='upper right')

"""
log-log plot of magnitude of coefficients of cos(cos(x)) 
"""
plt.figure(6)
plt.loglog(np.abs(coeff_intg_ccx), 'ro', label = 'Coefficients - Integrated')                        #coeff_intg_exp are plotted in red 
plt.loglog(np.abs(coeff_lstq_ccx), 'go', label = 'Coefficients - Least Squares',markersize = 4.5)    #coeff_lstq_exp are plotted in green and of lower size     
plt.grid(True)
plt.ylabel(r'Magnitude $\rightarrow$',fontsize=12)
plt.xlabel(r'$n\rightarrow$',fontsize=12)
plt.title('log-log plot of fourier coeffcients of $cos(cos(x))$',fontsize=14)
plt.legend(loc='upper right')

"""
The following lines of code computes the maximum deviation of coefficients generated using least squares
from that of integrated value for both the functions 
"""
delta_exp = np.max(abs(coeff_intg_exp - coeff_lstq_exp))
delta_ccx = np.max(abs(coeff_intg_ccx - coeff_lstq_ccx)) 

print('The deviation in case of exp(x) is %e' % delta_exp)
print('The deviation in case of cos(cos(x)) is %e' % delta_ccx)

"""
The following lines of code computes the product of A, coeff_lstq_f which is a function of x, 
plotting these points using green dots and true value of the function in red, where x has 400
points in range [0,2*pi). 
"""
Ac_exp = np.dot(A_exp,coeff_lstq_exp)
Ac_ccx = np.dot(A_ccx,coeff_lstq_ccx)

plt.figure(7)
plt.semilogy(x,Ac_exp,'go',label = 'Using Least squares',markersize = 5)    #x vs Ac_exp using green dots 
plt.semilogy(x,exp(x),'r', label = 'Actual value')                          #x vs exp(x) in red
plt.grid(True)
plt.ylabel(r'$e^{x}\rightarrow$',fontsize=12)
plt.xlabel(r'x$\rightarrow$',fontsize=12)
plt.title('Semilog plot of $e^{x}$',fontsize=14)
plt.legend(loc='upper right')

plt.figure(8)
plt.plot(x,Ac_ccx,'go',label = 'Using Least squares',markersize = 5)    #x vs Ac_ccx using green dots 
plt.plot(x,cos_cos(x),'r', label = 'Actual value')                      #x vs cos(cos(x)) in red 
plt.grid(True)
plt.ylabel(r'$cos(cos(x))\rightarrow$',fontsize=12)
plt.xlabel(r'x$\rightarrow$',fontsize=12)
plt.title('Plot of $cos(cos(x))$',fontsize=14)
plt.legend(loc='upper right')

"""
Plots linear plot of approximated function, uncomment to plot
"""
"""
plt.figure(9)
plt.plot(x,Ac_exp,'g',label = 'Using Least squares',markersize = 5)     #x vs Ac_exp using green 
plt.plot(x,exp(x),'r', label = 'Actual value')                          #x vs exp(x) in red
plt.grid(True)
plt.ylabel(r'$e^{x}\rightarrow$',fontsize=12)
plt.xlabel(r'x$\rightarrow$',fontsize=12)
plt.title('Linear plot of $e^{x}$',fontsize=14)
plt.legend(loc='upper right')

"""
plt.show()    #Displaying all the plots 

"""
---------------------------------------End-Of-The-Code-----------------------------------------
"""
 
