""" 
EE2703 Applied Programming Lab - Jan-May 2021 
Purpose : Final Assesment - Submission
Date    : 31th May 2021
Author  : Sai Gautham Ravipati (EE19B053)  
Usage   : python3 ee19b053.py  [Please refer Psuedo-Code at the end before execution] 
Input   : None 
Version : Python 3.8.3 [four spaces - one indentation level]
Output  : Plots of Vector Current, Magnetic Field along z-axis: Linear and LogLog scale, Fit of Magnetic Field 
          using Least Squares - stored in 'EE2703_plots' directory. Value of coefficients for fit-Cz^b, Max. &
          Min. and Mean value of field - printed in Terminal. [Please view the plots on full-screen] 
"""

"""
---------------------------Start-Of-The-Code------------------------------
"""

"""
Importing necessary modules and creating new directory. PyLab is  not used 
instead numpy and matplotlib.pyplot  have been  imported as recommended by 
the official documentation. 
"""
import os                                       #os module has been used to create directory to store plots 
import numpy as np                              #numpy module has been used to solve all vector equations
import matplotlib.pyplot as plt                 #matplotlib.pyplot is been used to create & store the plots

os.makedirs('EE2703_plots', exist_ok = True)    #Creating new directory to store the plots using os module

"""
Below two lines of codes define the style of annotaion used  in the plots, 
and have nothing to do with the algorithm used to solve the equation. 
"""
bbox = dict(boxstyle = "round",fc = "0.8")                                     
arrowprops = dict(arrowstyle = "->",connectionstyle = "angle,angleA = 0,angleB = 90,rad = 10")

"""
The following function has been used to fit the Bz to Cz^b, by  performing
Least Squares Minimisation. np.linalg.lstsq perfors the same. By taking l-
-og on both sides we obtain log(Bz) = log(C) + b*log(z). By expressing the
same in Matrix form yields [b,log(C)][log(z),1].T = log(Bz). This  is very 
much similar to a linear matrix equation M.p = g(t), and  least squares is 
performed to obtain b and log(C). Coefficients b and C are returned by the 
function. [Q10]
"""
def fit(B,z1):
    """
    Inputs - B : Field Vector 
             z1: z coordinates to obatin Cz^b
    Ouputs - C,b : Coefficients of the fit obtained using least squares 
    """    
    x_in  = np.c_[np.log(z1), np.ones(len(z1))]         #Stacking the array z1 with ones to form the array (log z1, 1)         
    y_in  = np.log(B).T                                 #Taking y_in a B array 
    b,C_l = np.linalg.lstsq(x_in,y_in,rcond=None)[0]    #Performing least squares regression to get b and log(C)
    
    return np.exp(C_l),b                                                                   

print('\nGenerating plots and data. Please Wait.')  
"""
Firstly radius, a is assigned a  value of 10 and number of sections n  are
assigned a value of 100.d1, d2, d3 generate points along x,y,z  directions 
which when operated by  np.meshgrid break the volume into (3,3,1000) mesh. 
The polar coordinate vector phi contains n points in the interval [0,2pi).
The length of each segment s is  obtained by  dividing circumference by n. 
The current is scales by 4*pi/u0 for better scale in the plot, and is giv-
-en by cos(phi)*[-sin(phi)*i + cos(phi)*j]. Similarly the  vectors r' & dl 
are given by a[cos(phi)*i + sin(phi)*j] & -s[sin(phi)*i + cos(phi)*j] res-
-pectively. np.c_ has been used to store the vectors. Once we have the ve-
-ctors plt.quiver has been used to represent the vector currents & plt.pl-
-ot is used to plot the center of each segment. The plot is saved as 'vec-
-tor_current.png' in 'EE2703_plots' directory. [Q2,Q3,Q4]
"""
a  = 10                                                         #Setting the radius to 10                          
n  = 100                                                        #Setting the number of sections to 100

d1 = np.linspace(0,2,3)                                         #Generating array [0,1,2] which act as points along x
d2 = np.linspace(0,2,3)                                         #Generating array [0,1,2] which act as points along y
d3 = np.linspace(1,1000,1000)                                   #Generating array [1,2,..,1000] which act as points along z
X,Y,Z = np.meshgrid(d1,d2,d3)                                   #Generating meshgrid using these points, which shall break the space into 3x3x1000 grid  

phi = np.linspace(0, 2*np.pi, n, endpoint = False)              #Generating the polar vector phi, where phi is in [0,2pi) with n points 
s   = 2*np.pi*a/n                                               #Finding the length of each segment s 

I  = np.c_[-np.sin(phi)*np.cos(phi),np.cos(phi)*np.cos(phi)]    #Finding the scaled crrent vector, scaled by 4*pi/u0 
r1 = a*np.c_[np.cos(phi),np.sin(phi)]                           #Finding the vector r'[l]
dl = s*np.c_[-np.sin(phi),np.cos(phi)]                          #Finding the vector dl[l] 
 
plt.figure(0,figsize = (8,8))                                                          #Initialising new figure to plot vector currents and segments 
plt.quiver(r1[:,0], r1[:,1], I[:,0], I[:,1], scale = 20, label = 'Vector Currents')    #Plotting the quiver plot of vector currents 
plt.plot(r1[:,0],r1[:,1],'r.',label = 'Centre Points of Elements',markersize = 5)      #Plotting the segments represented by their centers 
plt.title('Current Elements in x-y plane, split into {} sections'.format(n),fontsize = 14)    
plt.legend(loc = 'upper left')                                                              
plt.xlabel(r'$\leftarrow$ x (in cm) $\rightarrow$',fontsize = 12)          
plt.ylabel(r'$\leftarrow$ y (in cm) $\rightarrow$',fontsize = 12)
plt.grid(True) 
plt.axis('square')
plt.annotate('Max. Amplitude',(a,0),xytext = (-135,0),textcoords = 'offset points',bbox = bbox,arrowprops = arrowprops)    #Annotating the Max Amplitude current vector at (a,0)
plt.annotate('Min. Amplitude',(0,-a),xytext = (0,45), textcoords = 'offset points',bbox = bbox,arrowprops = arrowprops)    #Annotating the Min Amplitude current vector at (0,-a)
plt.annotate('Max. Amplitude',(-a,0),xytext = (57,0), textcoords = 'offset points',bbox = bbox,arrowprops = arrowprops)    #Annotating the Max Amplitude current vector at (-a,0)
plt.annotate('Min. Amplitude',(0,a),xytext = (0,-50), textcoords = 'offset points',bbox = bbox,arrowprops = arrowprops)    #Annotating the Min Amplitude current vector at (0,a) 
plt.savefig('EE2703_plots/vector_current.png')                                                                             #Saving the plot in directory as 'vector_current.png'

"""
The following function calc(l), for a  given index l, relates to r'(l) and
dl(l), which corresponds to the location of a segment of the loop antenna.
For a given segment, firstly the R_ijkl = |r_ijk - r'(l)| is to be calcul-
-ated for every point in the mesh grid. This can be done in a straight way
using for loops. But numpy vectorisation saves time & the same is  used in 
the function call below. This is done by forming an array of meshgrid poi-
-ts, after r'x, r'y is removed from X,Y points, & Z remains same as r'z is
0. Taking the norm translates to [(X-r'x)^2 + (Y-r'y)^2 + Z^2]^0.5 for all
the mesh grid points. The norm is taken along axis 0, and the  array is R, 
which is of the shape 3 x 3 x 1000, and holds R_ijkl, for all ijk. Then it
is extended to calculate the terms A_l that are to be summed up to give A. 
This is done by element wise operation  for R elements, and multiplication 
with dl[l] which is of the dimension (1,2). Returned value of the function
is of the shape-(3,3,1000,2). A for-loop is used to sum over all values of
l, which finally gives the A vector at all the locations in mesh. Here, we
can use for-loop for two reasons. The index l is assumed to be a scalar in 
the function calc(l) and and point-wise subtraction is done. To solve this 
calc has to be modified to accept vector inputs for l. Once the modificat-
-ion is done, altough we can get the vectorised norm for all values of l &
the dimension of the returned array is of the form - (100,3,3,1000,2). We
should sum over the first axis which has 100 elements, to get the A at all
points in space. Here again np.sum uses a for-loop internally by looking a
C program, making it a bit faster then the ordinary for-loop, but the time 
shall be of same order. Further the extra matrix operation in calc(l) adds 
a finite computation time which shall compensate the time saved by np.sum,
& at the end both the approaches end up almost at equal speeds, or atleast
at the same order, for say ms. So to avoid complexities in  using calc(l), 
a single for-loop is used to compute A, the input l to calc(l)  is kept as
a scalar. At the end we have A at locations in space. Note that both A_x &
A_y are taken as single vector A, and incremented  at the same time in the 
loop which reduces the computation time when compared to the case of doing
seperately. [Q5,Q6,Q7]
"""
def calc(l):
    """
    Inputs - l : The index of the segment  
             (X,Y,Z,r1,dl taken from the global namespace)
    Ouputs - A_l of the dimension (3,3,1000,2) 
    """ 
    R1  = np.array((X-r1[l,0], Y-r1[l,1], Z))         #Forming an array of meshgrid points after point-wise removal of r1_x,r1_y from X,Y,Z
    R   = np.linalg.norm(R1, axis = 0)                #Taking the norm to compute R_ijk 
    R   = R.reshape((3,3,1000,1))                     #Reshaping from (3,3,1000) to (3,3,1000,1) to enable the multiplication with dl - (1,2) 
    A_l = (np.cos(phi[l])*np.exp(-0.1j*R)/R)*dl[l]    #Computing A_l as per the equation, when added up gives A in space 
    
    return A_l

A = calc(0)              #Initialising A with calc(0)
for i in range(1,n):
    A += calc(i)         #Accumulating with calc(i) where i :[1,n) to get A in space 

"""
Vector field B is given as curl(A) in space. But we have only A_x & A_y in
the vector A. These compinents are sufficient to calculate the z-component 
of magnetic field, which is given as Bz = partial(A_y,x)-partial(A_x,y).To
vectorise the same we have to consider small variation in x & y, so taking
dx= dy= 2, partial(A_y,x) = {A_y[y:1,x:2,:] - A_y[y:1,x:0,:]}/dx and other 
term partial(A_x,y) = {A_x[y:2,x:1,:] - A_y[y:0,x:1,:]}/dy. Then  the term 
Bz = {A_y[y:1,x:2,:] - A_y[y:1,x:0,:] - A_x[y:2,x:1,:] + A_y[y:0,x:1,:]}/2.
Here the field is obatined on a line parallel to z-axis located at d1[1] &
d2[1] which correspond to (1, 1) in terms of x-y coordinates. The field is 
being considered at a location using off-grid beacuse, by  symmetry of the 
distribution, the field on the z-axis is 0. Due to precision of the CPU,if
we try to find the field at origin shall give very small values of the or-
-der 1e-17 and it has no theoritical interest. So field is being consider-
-ed at (1,1) on a line parallel to z-axis & this gives us non-zero values.
Secondly while computing the  vectorised version of  Bz, it is to be noted 
that sub-indices i, j are interchanged since we have  used a mesh-grid, so 
the first sub-index corresponds to y & the second sub-index corresponds to 
x. Once we get the Bz field the same is plotted in both linear and log-log
plots using 'plt.plot' and 'plt.loglog'. The plots are saved as 'field-li-
-near.png' and 'field_loglog.png' in the  directory. Note that the unit of 
field is T, as B is independent of length. [Q8,Q9]
"""
Bz = 0.5*(A[1,2,:,1]-A[2,1,:,0]-A[1,0,:,1]+A[0,1,:,0])                        #Computing the field using the vectorised implementation and printing the data

print('\nData of magnetic field parallel to z-axis at (1,1,0)')
print('1.The mean value of field for z : 1 cm to 1000 cm is {:.3e} T.'.format(np.mean(abs(Bz))))
print('2.The max. value of field for z : 1 cm to 1000 cm is {:.3e} T.'.format(np.max(abs(Bz))))
print('3.The min. value of field for z : 1 cm to 1000 cm is {:.3e} T.'.format(np.min(abs(Bz))))

plt.figure(1, figsize = (9,6))                                                #Initializing new figure to plot field in linear scale on z-axis at (1,1,0)
plt.plot(d3,np.abs(Bz),'-r.',label = 'Field Magnitude',markersize = 7.5)      #Plotting the points at [1,1000] joined by lines  
plt.title('Magnetic field along z-axis at (1,1) - Linear',fontsize = 14)   
plt.grid(True) 
plt.legend(loc = 'upper right')
plt.xlabel(r'z (in cm) $\rightarrow$',fontsize = 12)
plt.ylabel(r'$B_z$ (in T) $\rightarrow$',fontsize = 12)
plt.savefig('EE2703_plots/field_linear.png')                                  #Saving the plot in directory as 'field_linear.png'

plt.figure(2, figsize = (9,6))                                                #Initializing new figure to plot field in log-log scale on z-axis at (1,1,0)        
plt.loglog(d3,np.abs(Bz),'-r.',label = 'Field Magnitude',markersize = 7.5)    #Plotting the points at [1,1000] joined by lines  in log scale
plt.title('Magnetic field along z-axis at (1,1) - Loglog',fontsize = 14)
plt.grid(True) 
plt.legend(loc = 'upper right')
plt.xlabel(r'z (in cm) (log) $\rightarrow$',fontsize = 12)
plt.ylabel(r'$B_z$ (in T) (log) $\rightarrow$',fontsize = 12)
plt.savefig('EE2703_plots/field_loglog.png')                                  #Saving the plot in directory as 'field_loglog.png'

"""
Once we have the field, it can be verified from the loglog plot that after
around 30 cm, the field decays linearly. So assuming a fit-Cz^b, and using 
the function defined  previously, the data is fit using all the points and
points only after 30 cm being considered. Both the fits  have been plotted 
and the plot is saved as 'field_fit.png' in the directory. The C,b obtain-
-ed has been printed out. Note that the fit with all the points is plotted 
in red with solid line and the other fit is plotted in green with a dashed
line. [Q10] 
"""
C0,b0 = fit(np.abs(Bz),d3)              #Fitting all the data using the function fit and same is printed out
C1,b1 = fit(np.abs(Bz[30:]),d3[30:])    #Fitting the data only after 30 cm using the function fit and same is printed out 

print('\nThe magnetic along line parallel to z-axis at (1,1,0) is fit to Cz^b')    
print('The Coefficients with all points fit :')
print('1. C = {:.4f}'.format(C0))
print('2. b = {:.4f}'.format(b0))
print('The Coefficients with points after 30 cm fit :')
print('1. C = {:.4f}'.format(C1))
print('2. b = {:.4f}'.format(b1))

plt.figure(3, figsize = (9,6))           
plt.loglog(d3,np.abs(Bz),'k.',label = 'Original Field Magnitude',markersize = 7.5)          #Plotting the original points on a log-log plot
plt.loglog(d3,C1*d3**b1,label = 'Fit - Points from 30 cm', color = 'g',linestyle = '--')    #Plotting the fit with points after 30 cm being considered 
plt.loglog(d3,C0*d3**b0,label = 'Fit - All points',color = 'r')                             #Plotting the fit with all points being considered 
plt.title('Fit of Magnetic field along z-axis at (1,1) - Loglog',fontsize = 14)
plt.grid(True) 
plt.legend(loc = 'upper right')
plt.xlabel(r'z (in cm) (log) $\rightarrow$',fontsize = 12)
plt.ylabel(r'$B_z$ (in T) (log) $\rightarrow$',fontsize = 12)
plt.savefig('EE2703_plots/field_fit.png')                                                   #Saving the plot in directory as 'field_fit.png'

print('\nAll the plots have been saved in EE2703_plots/ directory') 

"""
-----------------------------End-Of-The-Code------------------------------
"""

"""
Pseudo-Code [Q1] 

Start
Set radius, a to 10
Set number of sections, n to 100
Select the range of points in space separated by 1 cm
    x - [0,2];
    y - [0,2]; 
    z - [1,1000]
Form a mesh of points, r_ijk in space in the selected range

Break the loop to n segments
    polar coordinate, p - [0,2*pi) with n points  
Set segment length, s to 2*pi*a/n
Compute scaled current, I as [-cos(p)*sin(p), cos(p)*cos(p)] 
Compute location of segments, r' as [a*cos(p), a*sin(p)] 
Compute segment vectors, dl' as [-s*sin(p), s*cos(p)]  
Plot the vector currents, I and mid-points of segments, r' 

Function calc
    Pass In: Integer corresponding to index of the segment
    Compute the distance of the all points of the mesh, R from the segment as |r_ijk - r'(segment)|
    Compute the individual terms in vector A corresponding to a given segment as {cos(p)*exp(-0.1j*R)/R}*dl' 
    Pass Out: A_l, Sub-vector to be added vector A, corresponding to a given segment
Endfunction

Init vector A, call calc function return A_l, for the first segment 

For each segment (l : 1 to n - 1) 
   call calc function, return A_l, for segment l 
   Accumulate vector A with A_l
Endfor
     
Compute z-component of field at (1,1) using the curl of vector A given as [Ax_component,Ay_component]
Plot the magnetic field, Bz and the z values 

Function fit 
    Pass In: Vector B, holding the magnetic field and vector z, holding the z points 
    Compute C,b for the fit B = Cz^b, translates to log B = [b,log(C)][log(z),1].T using least squares minimisation 
    Pass Out: Coefficients of the fit C,b  
Endfunction

Fit Bz to Cz^b points using call fit function, return C,b 
Plot the orginal field points Bz, the fit Cz^b and the z points on the same plot 
End
""" 
