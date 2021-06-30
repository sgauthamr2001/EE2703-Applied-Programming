""" 
EE2703 Applied Programming Lab - Jan-May 2021 
Purpose : Assignment-5 - Submission, 
Date : 25th March 2021, 
Author : Sai Gautham Ravipati (EE19B053), 
Input : Arguments like, Size along x,y, radius of lead and number of iterations, 
Output : Plots of current density, potential contour, errors, temperture contour, 
Usage : python3 EE2703_ASSIGN5_EE19B053.py --Nx x --Ny y --R r --n Niter.  
"""

"""
Importing necessary modules
"""
import numpy as np  
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import argparse
import os 

"""
Parsing the command line arguments, if the user doesn't pass any default values are taken, further execution is stopped if 
radius of the lead is more then 0.5 cm. 
"""
parser = argparse.ArgumentParser()
parser.add_argument("--Nx", metavar = 'x', type = int, default = 25, help='Size along x')                            #Default value of size along x is 25 
parser.add_argument("--Ny", metavar = 'y', type = int, default = 25, help='Size along y')                            #Default value of size along y is 25 
parser.add_argument("--R", metavar = 'r', type = float, default = 0.35, help ='Radius of central lead')              #Default Radius of central lead is 0.35 cm 
parser.add_argument("--n", metavar = 'Niter', type = int, default = 1500, help = 'Number of iterations to perform')  #Default number of iterations is 1500
args = parser.parse_args()
[Nx, Ny, r, Niter] = [args.Nx, args.Ny, args.R, args.n]                                                              #Storing the parsed arguments into [Nx,Ny,radius,Niter] array 

print('Current values (Nx, Ny, Radius, Niter) - (%d, %d, %.3f, %d)' % (Nx, Ny, r, Niter))                            #Being informed of the arguments passed or the default once   
print('To change the values, please update arguments with flag, for say, --Nx 10')                                       

if r > 0.5 :                                                                                                         #Asserting that the radius of the lead is less than 0.5 cm 
    print('Please check the radius of the lead, it should not exceed 0.5 cm, as dimensions of plate are 1 cm x 1 cm')
    exit()
    
os.makedirs('EE2703_plots', exist_ok = True)                                                                         #Creating a directory to store plots    

"""
Function definitions used elsewhere in the code 
""" 
def fit(x,y):
    """
    The following function takes in x data-points and y data-points, fits them to y = a.exp(bx), this is done by taking 
    log base, log y = log a + bx, which can be represented as a vector product of (log a, b) and (1, x), thus this is of 
    the form m = Kn and least squares regression is performed to obtain (log a, b). 
    Output -  np.exp(a) and b
    """
    x_in = np.vstack([x, np.ones(len(x))]).T            #Stacking the array x with ones to form the array (log x, 1) 
    y_in = np.log(y).T
    b,a  = np.linalg.lstsq(x_in,y_in,rcond=None)[0]     #Performing least squares regression 
    return np.exp(a),b  
                                    
def error_estimate(a,b,N) : 
    """
    Input  - a,b and N (Iteration number) 
    Outbut - Bound of the error approximated by cumulative sum 
    """
    return  -a*np.exp(b*(N+0.5))/b    
    
"""
Initializing and filling the arrays
"""
phi = np.zeros((Ny,Nx))                               #Initializing the array phi of dimensions (Ny,Nx) 
a = np.linspace(-0.5,0.5,Ny)                          #Generating Ny points between (-0.5,0.5) assuming the center of the plate is at (0,0)
b = np.linspace(-0.5,0.5,Nx)                          #Generating Nx points between (-0.5,0.5) assuming the center of the plate is at (0,0)
Y,X = np.meshgrid(b,-a)                               #Forming a mesh grid of these points 
ii  = np.where((X*X + Y*Y) <= (r)**2)                 #Asserting the locations of the nodes which lie within the radius of the central lead 
phi[ii] = 1.0                                         #Setting the V of nodes in this region to 1, and the rest shall remain zero
x_ii,y_ii = ii[1]/(Nx-1) - 0.5, ii[0]/(Ny-1) - 0.5    #Translating from meshgrid indices to the range (-0.5,0.5) to plot as red dots 
    
"""
Code block to contour plot the potential and mark V = 1 with red nodes
"""   
plt.figure(1)                                     
plt.contourf(Y, X, phi, cmap = plt.cm.pink)                   #Plots the contour plot using mesh grid. 
plt.plot(x_ii,y_ii,'ro',markersize = 3, label = 'V = 1')      #Plots the points corresponding to V = 1 in red. 
plt.xlabel(r'$\leftarrow$ x $\rightarrow$',fontsize = 12)     
plt.ylabel(r'$\leftarrow$ y $\rightarrow$',fontsize = 12)
plt.colorbar()
plt.legend(loc = 'upper right')
plt.title('Potential Contour')
plt.axis('square')
plt.savefig('EE2703_plots/fig1.png')                          #Saving the plot 
                                      
"""
Performing iteration, Updating potential, asserting boundary conditions 
"""                                                        
phi_old = phi.copy()                                                                    #Saving the previous copy of potential into phi_old
error = np.zeros(Niter)                                                                 #Initializing error vector 
N     = np.arange(Niter) + 1                                                            #N is an array [1,2,...,Niter]
for i in N:                                                                             #Performing iterartions and asserting boundary conditions at each stage 
    phi_old = phi.copy()                                  
    phi[1:-1,1:-1] = (phi[1:-1,:-2] + phi[1:-1,2:] + phi[:-2,1:-1] + phi[2:,1:-1])/4    #Updating the potential at each stage 
    phi[1:-1,0]    = phi[1:-1,1]                                                        #Boundary condition - left
    phi[1:-1,-1]   = phi[1:-1,-2]                                                       #Boundary condition - right
    phi[0,:]       = phi[1,:]                                                           #Boundary condition - top 
    phi[-1,1:-1]   = 0                                                                  #Ground condition 
    phi[ii]  = 1.0
    error[i-1] = (np.abs(phi - phi_old).max())                                          #Populating the error vector

"""
Code block to plot error vs iteration in different scales
"""

plt.figure(2)                             #Plotting error on a semilog scale
plt.semilogy(error,'g')                                         
plt.title('Error - semi log scale')                      
plt.grid(True)
plt.xlabel(r'Iterations  $\rightarrow$',fontsize = 12)
plt.ylabel(r'Error (log) $\rightarrow$',fontsize = 12)
plt.savefig('EE2703_plots/fig2.png')                         

plt.figure(3)                             #Plotting error on a log-log scale as well every 50th points for better representation
plt.loglog(N,error, 'b', label = 'Error')                      
plt.loglog(N[::50],error[::50],'ro',label = '50th point sampled')
plt.title('Error - log log scale')
plt.grid(True)
plt.legend(loc = 'upper right') 
plt.xlabel(r'Iterations (log)  $\rightarrow$',fontsize = 12)
plt.ylabel(r'Error (log)  $\rightarrow$',fontsize = 12) 
plt.savefig('EE2703_plots/fig3.png')

"""
Code block to fit the error to exponential and plot the fits, further to plot the cumulative error  
"""
a_500,b_500   = fit(N[500:],error[500:])                       #Fitting using only points after 500
a_1500,b_1500 = fit(N,error)                                   #Fitting using all the points 
error_fit_1 = a_500*np.exp(b_500*N)                            #Generating the data-points after fitting through points after 500
error_fit_2 = a_1500*np.exp(b_1500*N)                          #Generating the data-points after fitting the entire data    

plt.figure(4)                                                  #Plotting every 50th point on a semi-log scale       
plt.semilogy(N,error_fit_1,'r',label = 'fit_1')                
plt.semilogy(N,error_fit_2,'--g',label = 'fit_2')
plt.semilogy(N[::50],error[::50],'bo',label = 'True Error',markersize = 3) 
plt.title('Error fit - semilog scale')
plt.grid(True)
plt.legend(loc = 'upper right')
plt.xlabel(r'Iterations $\rightarrow$',fontsize = 12)
plt.ylabel(r'Error (log) $\rightarrow$',fontsize = 12)
plt.savefig('EE2703_plots/fig4.png')                 

plt.figure(5)                                                  #Plotting every 50th point on a log-log scale       
plt.loglog(N,error,'b',label = 'True Error') 
plt.loglog(N[::50],error_fit_1[::50],'ro',label = 'fit_1')     #Considering every 50th point 
plt.loglog(N[::50],error_fit_2[::50],'g*',label = 'fit_2')
plt.title('Error fit - log log scale')
plt.grid(True)
plt.legend(loc = 'upper right')
plt.xlabel(r'Iterations (log) $\rightarrow$',fontsize = 12)
plt.ylabel(r'Error (log) $\rightarrow$',fontsize = 12)
plt.savefig('EE2703_plots/fig5.png') 

error_bound = error_estimate(a_1500,b_1500,N)                  #Finding the maximum posible error uing the function defined previously

plt.figure(6)                                                  #Plot of maximum bound error in log-log scale 
plt.loglog(N[::50],error_bound[::50],'ro',label = 'Error')     #Considering every 50th point 
plt.title('Accumulated Error - log log scale')    
plt.grid(True)
plt.legend(loc = 'upper right')
plt.xlabel(r'Iteration (log) $\rightarrow$',fontsize = 12)
plt.ylabel(r'Error Bound (log) $\rightarrow$',fontsize = 12)
plt.savefig('EE2703_plots/fig6.png') 

"""
Plotting the surface plot of potential
""" 
fig = plt.figure(7)                                            
ax =  plt.axes(projection ='3d')                              
ax.set_title('The 3D surface plot of the potential')
ax.set_xlabel(r'$\leftarrow$ x $\rightarrow$',fontsize = 12)     
ax.set_ylabel(r'$\leftarrow$ y $\rightarrow$',fontsize = 12)
surf = ax.plot_surface(Y, X, phi, rstride=1, cstride=1, cmap = plt.cm.jet)    #Plots 3D surface plot of potential after iterating
fig.colorbar(surf, ax = ax,shrink = 0.5, aspect = 5)
plt.savefig('EE2703_plots/fig7.png') 

"""
Plotting the potential contour after performing Niter iterations 
"""
plt.figure(8)                                                    
plt.contourf(Y,X,phi,cmap = plt.cm.pink)                          #Plot contour plot after iterations
plt.plot(x_ii,y_ii,'ro',markersize = 3, label = 'V = 1')          #Marks the points V = 1 with red dots 
plt.xlabel(r'$\leftarrow$ x $\rightarrow$',fontsize = 12)          
plt.ylabel(r'$\leftarrow$ y $\rightarrow$',fontsize = 12)
plt.title('Potential Contour after %d Iterations' % Niter)
plt.legend(loc = 'upper right')
plt.colorbar()
plt.axis('square')
plt.savefig('EE2703_plots/fig8.png') 

"""
Plotting the quiver plot of current density 
""" 
Jx = (phi[1:-1,0:-2] - phi[1:-1,2:])/2    #Evaluating the density with x-gradient                  
Jy = (phi[0:-2,1:-1] - phi[2:,1:-1])/2    #Evaluating the density with y-gradient 

plt.figure(9)
plt.quiver(Y[1:-1,1:-1], X[1:-1,1:-1], Jx, -Jy, scale = 6)
plt.plot(x_ii,y_ii,'ro',markersize = 3, label = 'V = 1')          #Marks the points V = 1 with red dots 
plt.legend(loc = 'upper right')
plt.title('Vector plot of current flow')
plt.xlabel(r'$\leftarrow$ x $\rightarrow$',fontsize = 12)          
plt.ylabel(r'$\leftarrow$ y $\rightarrow$',fontsize = 12)
plt.axis('square')
plt.savefig('EE2703_plots/fig9.png') 

"""
To find the contour of of Temperature of iteration, a similar procedure is followed as done in the 
case of potential to solve the equation and to do poisson update
"""
T     = np.zeros((Ny,Nx)) + 300                                                                   #Initializing a vector with all entries as 300 
T_old = T.copy()                                                                                  #Saving the previous copy of T 
for i in N:                                                                                       #Performing iterartions and asserting boundary conditions at each stage 
    T_old = T.copy()                                  
    T[1:-1,1:-1] = (T[1:-1,:-2] + T[1:-1,2:] + T[:-2,1:-1] + T[2:,1:-1] + (Jx)**2 + (Jy)**2)/4    #Updating T at each stage 
    T[1:-1,0]    = T[1:-1,1]                                                                      #Boundary condition - Left 
    T[1:-1,-1]   = T[1:-1,-2]                                                                     #Boundary condition - Right
    T[0,:]       = T[1,:]                                                                         #Boundary condition - Top
    T[-1,1:-1]   = 300                                                                            #Asserting nodes of the wire and ground are at 300K
    T[ii]        = 300

plt.figure(10)                                                    
plt.contourf(Y,X,T,100,cmap = plt.cm.hot)                             #Plots contour plot after iterations
plt.plot(x_ii,y_ii,'ro',markersize = 3, label = 'V = 1')              #Marks the points V = 1 with red dots 
plt.xlabel(r'$\leftarrow$ x $\rightarrow$',fontsize = 12)          
plt.ylabel(r'$\leftarrow$ y $\rightarrow$',fontsize = 12)
plt.title('Temperature Contour after %d Iterations' % Niter)
plt.legend(loc = 'upper right')
plt.colorbar()
plt.axis('square')
plt.savefig('EE2703_plots/fig10.png') 

print('\nAll the plots have been saved in EE2703_plots/ directory') 
"""
---------------------------------------End-Of-The-Code-----------------------------------------
"""

