""" 
EE2703 Applied Programming Lab - Jan-May 2021 
Purpose : Assignment-10 - Submission, 
Date    : 24th May 2021, 
Author  : Sai Gautham Ravipati (EE19B053), 
Input   : None 
Output  : Phase and Magnitude spectrums of Non-Periodic Signals 
Usage   : python3 EE2703_ASSIGN10_EE19B053.py   
"""
"""
Importing neccesary modules
"""
import numpy as np 
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt 
import os 
import warnings
warnings.filterwarnings("ignore")

os.makedirs('EE2703_ASN10/', exist_ok = True)    #Creating a directory to store plots and data

"""
Different Functions which find the spectrum, plot the spectrum and function, estimate w and delta.
"""
def spectrum_gen(t_range, N, f, wnd_flag, fft_flag, noise):

    """
    This function generates Y,w of the spectrum for various input functions and returns them. The
    inputs to the function are :
    t_range   : [low_t,up_t] where low_t is the lower limit and up_t is the upper limit of t.  
    N         : Number of points for which the function is calculated.
    f         : Name of the function which has t as argument.
    wnd_flag  : Boolean value, 1 if hamming window is to be used orwise 0 is passed.
    fft_flag  : Boolean value, 1 to perform 'y[0] = 0' and fftshift(y), orwise 0 is passed.
    noise     : Boolean value, 1 is passes to add Gaussian noise of A - 0.1 is to be added, 
                orwise 0 is passed. 
  
    Output - Y : DFT Spectrum points which is of the dimension N. 
             w : The corresponding values of w, which map to values of Y. 
    """
    t  = np.linspace(t_range[0], t_range[-1], N+1)                 #Generating N+1 t values between low_t, up_t      
    t  = t[:-1]                                                    #Removing the last value to avoid repetition, for say [0,2pi) we count eiter 0 or 2pi
    dt = t[1]-t[0]                                                 #Finding the difference between adjacent time points, this is the sampling rate
    w  = np.linspace(-np.pi, np.pi, N+1)/dt                        #Generating corresponding range of w, for which Y is mapped, w_max is given as pi*N/T           
    w  = w[:-1]                                                    #Removing the last value to avoid repetition similar to time scale
   
    wnd = 1                                                        #If wnd_flag is 0, asserts f(t) remains the same
    if wnd_flag :                                                  #If wnd_flag is 1, window is used and wnd (n) is sampled from hamming window function
        n   = np.arange(N)
        wnd = np.fft.fftshift(0.54+0.46*np.cos(2*np.pi*n/(N-1)))
        
    y = f(t)*wnd                                                   #y is given by the product of f and wnd 
    
    if noise : 
      y = y + 0.1*np.random.randn(N)                               #If noise is not 0, random noise of amplitude 0.1 is added to y
      
    if fft_flag :                                                  #If fft_flag is not 0, then y[0] is set to zero and fftshift makes y to start with 0 
        y[0] = 0                                                
        y = np.fft.fftshift(y)   
      
    Y = np.fft.fftshift(np.fft.fft(y))/N                           #Computes fft, fftshift is done to get Y(0) to center, scaled by N to use as spectrum
   
    return Y,w 
    
def spectrum_plot(Y, w, y_title, x_range, plt_flag, file_name, thresh): 

    """
    This function takes in Y,w and plots the spectrum various input functions and then saves them
    in the directory EE2703_ASN10. The inputs to the function are :
    Y         : DFT Spectrum points which is of the dimension N.
    w         : The corresponding values of w, which map to values of Y.
    y_title   : Title to be added to the plot, displayed as 'Spectrum of <y_ttl>'
    x_range   : Range of x axis data-points to be displayed in the plot. 
    plt_flag  : 1 -    Plots only those points having mag > thresh in phase spectrum.
                else - Plots all phase points in red with point mentioned in 1 as green dots.
    file_name : Output file name of the image
    thresh    : Phase of points with mag < thresh is taken as zero if plt_flag is 1.           
                
    Output - None
    The plots shall be saved as EE2703_ASN10/<file_name>
    """ 
    plt.figure()                                                   #Generates new figure for every instance the function is called              
    plt.subplot(2,1,1)                                             #Sub-plot to plot the magnitude spectrum                 
    plt.plot(w, abs(Y), lw=2)                                      #Plots mag spectrum obtained using DFT
    plt.xlim(x_range)                                              #Limiting the range of values to be displayed to x_range 
    plt.ylabel(r"$|Y|$", size=12)                                  #Setting y-label                    
    plt.title(r"Spectrum of {}".format(y_title), size=14)          #Setting the title using <y_title>   
    plt.grid(True)                                                 #Setting the grid 
    
    plt.subplot(2,1,2)                                             #Sub-plot to plot the phase spectrum  
    
    if plt_flag == 1 :                                             #If flag is 1, points with mag > thresh are plotted in green
        ii = np.where(abs(Y) > thresh)                             #Checking the points with mag > thresh  
        plt.plot(w[ii], np.angle(Y[ii]), 'go', lw=2, markersize=5)      
       
    else :                                                         #If flag is not 1 all the points are plotted in red, and points with mag > thresh in green
        plt.plot(w,np.angle(Y),'ro',lw=2, label = 'All Points')       
        ii = np.where(abs(Y) > thresh)                                  
        plt.plot(w[ii], np.angle(Y[ii]), 'go', lw=2, markersize=5, label = 'Mag > {}'.format(thresh))   
        plt.legend(loc = 'upper right')                               
        
    plt.xlim(x_range)                                              #Limiting the range of values to be displayed to x_range 
    plt.ylim([-4,4])                                               #Setting y-limits     
    plt.ylabel(r"Phase of $Y$", size=12)                           #Setting y-label
    plt.xlabel(r"$\omega$", size=12)                               #Setting x-label
    plt.grid(True)                                                 #Setting grid 
    
    plt.savefig('EE2703_ASN10/' + file_name)                       #Saving the plot in the directory 
    
    return 0
    
def func_plot(t_range, N, f, y_title, wnd_flag, noise, file_name):

    """
    This function plots the sampled points, in the time period being considered in blue and extends
    it to two more periods after the given period and one before and plots those points in red. The
    inputs are :
    t_range   : [low_t,up_t] where low_t is the lower limit and up_t is the upper limit of t. 
    N         : Number of points for which the function is calculated.
    f         : Name of the function which has t as argument.
    wnd_flag  : Boolean value, 1 if hamming window is to be used orwise 0 is passed.
    y_title   : Title to be added to the plot, displayed as 'Plot of <y_title>'.
    file_name : Output file name of the image.
    noise     : Boolean value, 1 is passes to add Gaussian noise of A - 0.1 is to be added, 
                orwise 0 is passed. 
    
    Output - None
    The plots shall be saved as EE2703_ASN10/<file_name>
    """
    t1 = np.linspace(t_range[0], t_range[-1], N+1)                 #Generating N+1 t1 values between low_t, up_t                                                   
    t1 = t1[:-1]                                                   #Removing the last value to avoid repetition, for say [0,2pi) we count eiter 0 or 2pi
    T  = t_range[1] - t_range[0]                                   #Finding the sampling periods 
    t2 = np.linspace(t_range[0]-T,t_range[-1]-T, N+1)              #Considering the period t1 - T
    t2 = t2[:-1]
    t3 = np.linspace(t_range[0] + T,t_range[-1] + T, N+1)          #Considering the period t1 + T
    t3 = t3[:-1]
    
    wnd = 1                                                        #If wnd_flag is 0, asserts f(t) remains the same
    if wnd_flag :                                                  #If wnd_flag is 1, window is used and wnd (n) is sampled from hamming window function
        n   = np.arange(N)
        wnd = np.fft.fftshift(0.54+0.46*np.cos(2*np.pi*n/(N-1)))
        
    y = f(t1)*wnd                                                  #y is given by the product of f and wnd 
    
    if noise : 
      y = y + 0.1*np.random.randn(N)                               #If noise is not 0, random noise of amplitude 0.1 is added to y
    
    plt.figure()                                                                 #Initializing new figure 
    plt.plot(t1,y,'bo', lw=2, label = 'Sampled Points', markersize=3)            #Plotting the sampled points                            
    plt.plot(t2,y,'ro', lw=2, label = 'Extension', markersize=3)                 #Plotiing the points in adjacent periods 
    plt.plot(t3,y,'ro', lw=2, markersize=3) 
    plt.legend(loc = 'upper right')                                              #Displays legend                                           
    plt.ylabel(r"$y$", size=12)                                                  #Setting ylabel                              
    plt.xlabel(r"$t$", size=12)                                                  #Setting xlabel
    plt.title(r"Plot of sampled {}, wrapped".format(y_title), size=14)           #Setting the title                   
    plt.grid(True)                                                               #Setting the grid 
    plt.savefig('EE2703_ASN10/' + file_name)                                     #Saving the plot in directory 
   
    return 0 
    
def wd_estm(w0, d, N, t_range, thresh, noise, dw=4, dd=2, p=1):

    """
    This function estimates w0 and delta for the case pf sinusoid from the DFT spectrum. The input 
    arguments are : 
    w0 : The true value of w0 used to generate the function points. 
    d  : The true value of d  used to generate the function points. 
    N         : Number of points for which the function is calculated.
    t_range   : [low_t,up_t] where low_t is the lower limit and up_t is the upper limit of t. 
    thresh    : Used in estimation of delta, where only points greater than thresh are used. 
    dw        : Radius of the Y points from 0, used in estimation w0.
    dd        : Radius of the Y points from 1, used in estimation d. 
    p         : Weight of the norm used in estimation of w0. 
    noise     : Boolean value, 1 is passes to add Gaussian noise of A - 0.1 is to be added, 
                orwise 0 is passed. 
    """
    f = lambda t : np.cos(w0*t + d)                                    #Defining the sinusoidal function  
    Y1,w1 = spectrum_gen(t_range, N, f, 1, 1, noise)                   #Getting the spectrum points Y and w using spectrum_gen()
    ii  = np.where(w1>=0)                                              #Finding indices of points with w >= 0
    Y,w = Y1[ii],w1[ii]                                                #Considering only one side of the spectrum  
    w_estm = np.sum(abs(Y[:dw]**p)*w[:dw])/np.sum(abs(Y[:dw]**p))      #Estimating w0 as Y weighted norm of dw number points of w.  
    jj = np.where(abs(Y) > thresh)[0]                                  #Finding inices of points with mag > thresh
    d_estm = np.mean(np.angle(Y[jj[1:dd]]))                            #Estimating d_estm as mean phase of points with indices in dd length
    
    return w_estm, d_estm, Y1, w1
    
print('Generating spectrums, Please Wait.\n') 

"""
Question - 2
""" 
f2 = lambda t : np.cos(0.86*t)**3                                                  #Function for Q2  
q2 = '$cos^3(0.86t)$ '                                                             #Title of the function

func_plot([-4*np.pi,4*np.pi], 256, f2, q2 + 'without wnd', 0, 0, 'q2.jpg')         #Plot of sampled function without window 
func_plot([-4*np.pi,4*np.pi], 256, f2, q2 + 'with wnd', 1, 0, 'q2_wd.jpg')         #Plot of sampled function with window                    

Y2a,w2a = spectrum_gen([-4*np.pi,4*np.pi], 256, f2, 0, 1, 0)                       #Generating the spectrum without hamming window     
Y2b,w2b = spectrum_gen([-4*np.pi,4*np.pi], 256, f2, 1, 1, 0)                       #Generating the spectrum with hamming window 

spectrum_plot(Y2a, w2a, q2 + 'without wnd', [-8,8], 1, 'q2_dft.jpg', 1e-3)         #Plotting the spectrum without hamming window being considered 
spectrum_plot(Y2b, w2b, q2 + 'with wnd', [-8,8], 1, 'q2_dft_wd.jpg', 1e-3)         #Plotting the spectrum with hamming window being considered 

"""
Question - 3,4
"""
w0 = 0.8
d  = 1

f34 = lambda t : np.cos(w0*t + d)                                                          #Function for Question 3,4
q34 = '$cos({}t + {})$ '.format(w0,d)                                                      #Title for the function

w0_em3, d_em3, Y3, w3 =  wd_estm(w0, d, 128, [-np.pi,np.pi], 1e-4, 0)                      #Estimating w0,d and generating spectrum with hamming window      
spectrum_plot(Y3, w3, q34 + 'without noise', [-8,8], 1, 'q3_dft.jpg', 1e-4)                #Plotting the spectrum with hamming window being considered 

print('Without noise :')                                                                   #Printing the error in estimates 
print('Error in estimate of w0 is {:.3f}'.format(abs(w0 - w0_em3)))
print('Error in estimate of delta is {:.3f}'.format(abs(d - d_em3)))

w0_em4, d_em4, Y4, w4 =  wd_estm(w0, d, 128, [-np.pi,np.pi], 1e-4, 1)                      #Estimating w0,d and generating spectrum with hamming window, noise      
spectrum_plot(Y4, w4, q34 + 'with noise', [-8,8], 1, 'q3_dft_wn.jpg', 1e-4)                #Plotting the spectrum with hamming window, noise being considered 

print('With noise :')                                                                      #Printing the error in estimates 
print('Error in estimate of w0 is {:.3f}'.format(abs(w0 - w0_em4)))
print('Error in estimate of delta is {:.3f}'.format(abs(d - d_em4)))

"""
Question - 5,6
"""
f56 = lambda t : np.cos(16*(1.5 + t/(2*np.pi))*t)                                          #Function for Question 5,6  
q56 = r'$\cos(16(1.5+\frac{t}{2\pi})t)$ '                                                  #Title for the function    

Y5a,w5 = spectrum_gen([-np.pi,np.pi], 1024, f56, 1, 0, 0)                                  #Generating the spectrum with hamming window    
spectrum_plot(Y5a, w5, q56 + 'with wnd', [-60,60], 0, 'q5_dft_wd.jpg', 1e-3)               #Plotting the spectrum with hamming window being considered 

Y5b,w5 = spectrum_gen([-np.pi,np.pi], 1024, f56, 0, 0, 0)                                  #Generating the spectrum without hamming window    
spectrum_plot(Y5b, w5, q56 + 'without wnd', [-60,60], 0, 'q5_dft.jpg', 1e-3)               #Plotting the spectrum without hamming window being considered 

t6 = np.linspace(-np.pi,np.pi, 1025)                                                       #Generating t in [-pi,pi] and 1025 points                             
Y6a = np.zeros((16,64),dtype = 'complex_')                                                 #Initializing 2D array of spectrum data-points for with window case 
Y6b = np.zeros((16,64),dtype = 'complex_')                                                 #Initializing 2D array of spectrum data-points for without window case 

for i in range(16):                                                               
  Yia,w6 = spectrum_gen(t6[64*i:i*64 + 65], 64, f56, 1, 0, 0)                              #Breaking the 1024 vector into pieces that are 64 samples wide, extracting the DFT
  Yib,w6 = spectrum_gen(t6[64*i:i*64 + 65], 64, f56, 0, 0, 0)        
  Y6a[i][:] = Yia                                                                          #Storing the DFTs as a column in a 2D array
  Y6b[i][:] = Yib

t6 = t6[:-1]                                                                               #Making T as [-pi,pi) with 1024 points 
t6 = t6[::64]                                                                              #Considering every 64th point
t6,w6 = np.meshgrid(t6,w6)                                                                 #Generating mesh grid with w and t

fig_p1 = plt.figure()                                                                      #Initialising new figure for surface plot of phase with window
ax1 =  plt.axes(projection ='3d')                                                          #Setting axes as 3d
ax1.set_title('The 3D surface plot of Phase of Y with wnd', fontsize = 14)                 #Setting title 
ax1.set_xlabel(r'$\leftarrow$ w $\rightarrow$', fontsize = 12)                             #Setting xlabel
ax1.set_ylabel(r'$\leftarrow$ t $\rightarrow$', fontsize = 12)                             #Setting ylabel
surf1 = ax1.plot_surface(w6, t6 ,np.angle(Y6a.T), cmap = plt.cm.GnBu)                      #Plotting surface plot                     
fig_p1.colorbar(surf1, ax = ax1, shrink = 0.5, aspect = 5)                                 #Displaying colorbar 
plt.savefig('EE2703_ASN10/q6_surf_phase_wd.jpg')                                           #Saving the plot 


plt.figure()                                                                               #Initialising new figure for contour plot of phase with window
plt.contourf(t6, w6, np.angle(Y6a.T), cmap = plt.cm.GnBu)                                  #Plotting the contour plot    
plt.xlabel(r'$\leftarrow$ t $\rightarrow$',fontsize = 12)                                  #Setting xlabel
plt.ylabel(r'$\leftarrow$ w $\rightarrow$',fontsize = 12)                                  #Setting ylabel
plt.title('Contour plot of Phase of Y with wnd')                                           #Setting title 
plt.colorbar()                                                                             #Setting the colour bar 
plt.savefig('EE2703_ASN10/q6_cont_phase_wd.jpg')                                           #Saving the plot

plt.figure()                                                    
plt.contourf(t6, w6, abs(Y6a.T), cmap = plt.cm.GnBu)                                       #Initialising new figure for contour plot of mag with window
plt.xlabel(r'$\leftarrow$ t $\rightarrow$',fontsize = 12)                                  #Setting x-label
plt.ylabel(r'$\leftarrow$ w $\rightarrow$',fontsize = 12)                                  #Setting y-label
plt.title('Contour plot of $|Y|$ with wnd')                                                #Setting title 
plt.colorbar()                                                                             #Setting colourbar 
plt.ylim([-120,120])                                                                       #Limiting y-axis 
plt.savefig('EE2703_ASN10/q6_cont_mag_wd.jpg')                                             #Saving the plot 
                                                      
inds = np.where(abs(w6) > 120)                                                             #Removing points to make the plot readable
Y6a[:,inds] = np.nan

fig_m1 = plt.figure()                                                                      #Initialising new figure for surface plot of mag with window
ax2 =  plt.axes(projection ='3d')                                                          #Setting axes as 3d
ax2.set_title('The 3D surface plot of $|Y|$ with wnd', fontsize = 14)                      #Setting title                                     
ax2.set_xlabel(r'$\leftarrow$ w $\rightarrow$', fontsize = 12)                             #Setting xlabel
ax2.set_ylabel(r'$\leftarrow$ t $\rightarrow$', fontsize = 12)                             #Setting ylabel 
surf = ax2.plot_surface(w6, t6, abs(Y6a.T), cmap = plt.cm.GnBu, rstride=1, cstride=1)      #Plotting surface plot
fig_m1.colorbar(surf, ax = ax2, shrink = 0.5, aspect = 5)                                  #Displaying colorbar 
ax2.view_init(55,107)                                                                      #Setting the angle 
ax2.set_xlim(-120,120)                                                                     #Setting x_limit 
plt.savefig('EE2703_ASN10/q6_surf_mag_wd.jpg')                                             #Saving the plot 

fig_p2 = plt.figure()                                                                      #Initialising new figure for surface plot of phase without window
ax3 =  plt.axes(projection ='3d')                                                          #Setting axes as 3d
ax3.set_title('The 3D surface plot of Phase of Y without wnd', fontsize = 14)              #Setting title 
ax3.set_xlabel(r'$\leftarrow$ w $\rightarrow$', fontsize = 12)                             #Setting xlabel
ax3.set_ylabel(r'$\leftarrow$ t $\rightarrow$', fontsize = 12)                             #Setting ylabel
surf1 = ax3.plot_surface(w6, t6 ,np.angle(Y6b.T), cmap = plt.cm.GnBu)                      #Plotting surface plot                     
fig_p2.colorbar(surf1, ax = ax3, shrink = 0.5, aspect = 5)                                 #Displaying colorbar 
plt.savefig('EE2703_ASN10/q6_surf_phase.jpg')                                              #Saving the plot 


plt.figure()                                                                               #Initialising new figure for contour plot of phase without window
plt.contourf(t6, w6, np.angle(Y6b.T), cmap = plt.cm.GnBu)                                  #Plotting the contour plot    
plt.xlabel(r'$\leftarrow$ t $\rightarrow$',fontsize = 12)                                  #Setting xlabel
plt.ylabel(r'$\leftarrow$ w $\rightarrow$',fontsize = 12)                                  #Setting ylabel
plt.title('Contour plot of Phase of Y without wnd')                                        #Setting title 
plt.colorbar()                                                                             #Setting the colour bar 
plt.savefig('EE2703_ASN10/q6_cont_phase.jpg')                                              #Saving the plot

plt.figure()                                                    
plt.contourf(t6, w6, abs(Y6b.T), cmap = plt.cm.GnBu)                                       #Initialising new figure for contour plot of mag without window
plt.xlabel(r'$\leftarrow$ t $\rightarrow$',fontsize = 12)                                  #Setting x-label
plt.ylabel(r'$\leftarrow$ w $\rightarrow$',fontsize = 12)                                  #Setting y-label
plt.title('Contour plot of $|Y|$ without wnd')                                             #Setting title 
plt.colorbar()                                                                             #Setting colourbar 
plt.ylim([-120,120])                                                                       #Limiting y-axis 
plt.savefig('EE2703_ASN10/q6_cont_mag.jpg')                                                #Saving the plot 
                                                     
inds = np.where(abs(w6) > 120)                                                             #Removing points to make the plot readable
Y6b[:,inds] = np.nan

fig_m2 = plt.figure()                                                                      #Initialising new figure for surface plot of mag without window
ax4 =  plt.axes(projection ='3d')                                                          #Setting axes as 3d
ax4.set_title('The 3D surface plot of $|Y|$ without wnd', fontsize = 14)                   #Setting title                                     
ax4.set_xlabel(r'$\leftarrow$ w $\rightarrow$', fontsize = 12)                             #Setting xlabel
ax4.set_ylabel(r'$\leftarrow$ t $\rightarrow$', fontsize = 12)                             #Setting ylabel 
surf = ax4.plot_surface(w6, t6, abs(Y6b.T), cmap = plt.cm.GnBu, rstride=1, cstride=1)      #Plotting surface plot
fig_m2.colorbar(surf, ax = ax4, shrink = 0.5, aspect = 5)                                  #Displaying colorbar 
ax4.view_init(55,107)                                                                      #Setting the angle 
ax4.set_xlim(-120,120)                                                                     #Setting x_limit 
plt.savefig('EE2703_ASN10/q6_surf_mag.jpg')                                                #Saving the plot 

"""
Example - 1
"""
f1 = lambda t : np.sin(2**0.5*t)
q1 = '$sin(\sqrt{2}t)$ '

func_plot([-np.pi, np.pi], 64, f1, q1 + 'without wnd', 0, 0, 'q1.jpg')             #Plot of sampled function without window 
func_plot([-np.pi, np.pi], 64, f1, q1 + 'with wnd', 1, 0, 'q1_wd.jpg')             #Plot of sampled function with window                    

Y1a,w1a = spectrum_gen([-np.pi,np.pi], 64, f1, 0, 1, 0)                            #Generating the spectrum without hamming window     
Y1b,w1b = spectrum_gen([-np.pi,np.pi], 64, f1, 1, 1, 0)                            #Generating the spectrum with hamming window 

spectrum_plot(Y1a, w1a, q1 + 'without wnd', [-8,8], 1, 'q1_dft.jpg', 1e-3)         #Plotting the spectrum without hamming window being considered 
spectrum_plot(Y1b, w1b, q1 + 'with wnd', [-8,8], 1, 'q1_dft_wd.jpg', 1e-3)         #Plotting the spectrum with hamming window being considered 

"""
Plot of Sin(2^0.5*t) with sampling period. Uncomment to plot 
"""
"""
t1 = np.linspace(-np.pi,np.pi,65);t1=t1[:-1]                       #Generating 65 t1 values between -pi to pi and removing last point
t2 = np.linspace(-3*np.pi,-np.pi,65);t2=t2[:-1]                    #Extending the interval 
t3 = np.linspace(np.pi,3*np.pi,65);t3=t3[:-1]

plt.figure()                                                       #Initializing new figure
plt.plot(t1,np.sin(np.sqrt(2)*t1),'b',lw=2)                        #Plotting the points to be sampled 
plt.plot(t2,np.sin(np.sqrt(2)*t2),'r',lw=2)                        #Plotting the extension
plt.plot(t3,np.sin(np.sqrt(2)*t3),'r',lw=2)                     
plt.ylabel(r"$y$",size=12)                                         #Setting y label
plt.xlabel(r"$t$",size=12)                                         #Setting x label
plt.title(r"$\sin\left(\sqrt{2}t\right)$",size = 14)               #Setting title 
plt.grid(True)                                                     #Setting grid 
plt.savefig("EE2703_ASN10/q1_unsampled.jpg")                       #Saving the plot 
"""

"""
Plot of Spectrum of brick wall function. Uncomment to plot 
"""
"""
t = np.linspace(-np.pi,np.pi,65);t=t[:-1]                           #Generating 65 t1 values between -pi to pi and removing last point
dt= t[1]-t[0];fmax=1/dt
y = t
y[0]=0                                                              #the sample corresponding to -tmax should be set zeroo
y = np.fft.fftshift(y)                                              #make y start with y(t=0)
Y = np.fft.fftshift(np.fft.fft(y))/64.0                             #Computing DFT
w = np.linspace(-np.pi*fmax,np.pi*fmax,65);w=w[:-1]                 #Finding w     
plt.figure()                                                        #Creating new figure 
plt.semilogx(abs(w),20*np.log10(abs(Y)))                            #PLotting the spectrum
plt.xlim([1,10])                                                    #Limiting x axis 
plt.ylim([-20,0])                                                   #LImiting y axis 
plt.xticks([1,2,3,4,6,10])                                          #Setting xticks
plt.ylabel(r"$|Y|$ (dB)",size=12)                                   #Setting ylabel
plt.title(r"Spectrum of a digital ramp", size =14)                  #Setting title 
plt.xlabel(r"$\omega$",size=12)                                     #Setting xlabel
plt.grid()                                                          #Setting grid 
plt.savefig("EE2703_ASN10/ramp.jpg")                                #Saving the plot 
"""

"""
Estimator of optimal dw,dd,p. Uncomment to run (takes longer to run)
"""
"""
dd_range = np.arange(2,5)                                                                     #Range over which dd is iterated 
dw_range = np.arange(2,10)                                                                    #Range over which dw is iterated 
p_range  = np.arange(1,3,0.1)                                                                 #Range over which p  is iterated 

we_dict = {}                                                                                  #Stores the error in w estimated with dw and p                                                           
de_dict = {}                                                                                  #Stores the error in d estimated with dd 

noise  = 0
 
for dw in dw_range:                                                                           #Initialising the dictionary we with zeros 
    for p in p_range:                        
        we_dict[dw,p] = 0

for dd in dd_range:                                                                                #Iterating over dd_range 
    de = []                                                                                        #Stores de in each iteration 
    for dw in dw_range: 
      for p in p_range: 
          w_error = []                                                                             #Stores error in w in each iteration 
          d_error = []                                                                             #Stores error in d in each iteration  
          for i in range(500):            
              w0 = np.random.rand() + 0.5                                                          #Sampling w randomly between 0.5 and 1.5
              d  = np.random.rand()*2*np.pi - np.pi                                                #Sampling d randomly between -pi and pi 
              w_estm, d_estm,_,_ = wd_estm(w0, d, 128, [-np.pi,np.pi], 1e-5, noise, dw, dd, p)     #Using the function call to get the estimated values 
              w_error.append(abs(w_estm - w0))                                                     #Appending the error in w
              d_error.append(abs(d_estm - d))                                                      #Appending the error in d 
          we = np.mean(w_error)                                                                    #Taking the mean error in w 
          de.append(np.mean(d_error))                                                              #Appending the error to de 
          we_dict[dw,p] += we                                                                      #Accumulating error after iterations 
    de_dict[dd] = np.mean(de)                                                                      #Storing the de error in de dict 
 
dw,p  = min(we_dict, key=we_dict.get)                                                         #Finidng the dw and p which gives the minimum error in w
dd    = min(de_dict, key=de_dict.get)                                                         #Finding the dd which gives the minimum error in d 

print('The values of dw and p which give minimum deviation is {},{}'.format(dw,p))
print('The value of dd which give minimum deviation is {}'.format(dd))
"""

print('\nAll plots have been saved to EE2703_ASN10 directory')
"""
---------------------------------------End-Of-The-Code-----------------------------------------
"""
