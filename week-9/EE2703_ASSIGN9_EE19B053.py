""" 
EE2703 Applied Programming Lab - Jan-May 2021 
Purpose : Assignment-9 - Submission, 
Date    : 24th May 2021, 
Author  : Sai Gautham Ravipati (EE19B053), 
Input   : None 
Output  : Phase and Mgnitude spectrums of different functions,
Usage   : python3 EE2703_ASSIGN9_EE19B053.py   
"""
"""
Importing neccesary modules
"""
from pylab import * 
import os  

os.makedirs('EE2703_ASN9/', exist_ok = True)    #Creating a directory to store plots and data

"""
Function call which plots the spectrum for a given function 
"""
def func_spectrum(low_t, up_t, N, y, y_ttl, file_name, x_range, flag) :

    """
    This function generates spectrum plots for various input functions and then it saves them in the
    directory EE2703_ASN9. The inputs to the function are :
    low_t     : Lower limit of the time-interval for which function data-points are calculated
    up_t      : Upper limit of the iime-interval for which function data-points are calculated
    N         : Number of data-points for which the function is calculated, note fft is of same size 
    y         : Name of the function which has t as argument 
    y_ttl     : Title to be added to the function, displayed as 'Spectrum of <y_ttl>'
    file_name : Output file name of the image
    x_range   : Range of x datapoints to be displayed in the plot 
    flag      : 1 - Plots only those points having magnitude > 1e-3 in phase spectrum 
                2 - Plots all phase points in red with point mentioned in 1 as green dots
                else - Plots all the phases points in red
                
    Output - None
    The plots shall be saved as EE2703_ASN9/<file_name>
    """
    t = linspace(low_t,up_t, N + 1)    #Generating N+1 t values between low_t, up_t  
    t = t[:-1]                         #Removing the last value to avoid repetition, for say [0,2pi) we count eiter 0 or 2pi
    T = up_t - low_t                   #Finding the length of sampling period
    w = linspace(-pi,pi, N + 1)*N/T    #Generating N+1 w values between low_w, up_w
    w = w[:-1]                         #Removing the last value of w makes the size of the vector N
    
    figure()                           #Generates new figure for every instance the function is called 
    subplot(2,1,1)                     #Sub-plot to plot the magnitude spectrum  
    
    if (y_ttl == '$e^{-t^2/2}$') : 
                
        """
        Checking if the input function is Gaussian, as N,t are to be monitered 
        since it is not band-limited, and printing the max absolute error. 
        """ 
        Y   = fft(ifftshift(y(t)))               #Using fft to generate the spectrum, y and Y are of same length
        Y   = fftshift(Y)/N                      #Shifting the values such that Y[0] is in middle
        Y   = Y*T                                #Multiplying with T 
        Y_t = exp(-w**2/2)*sqrt(2*pi)            #True transform obtained analytically, is evaluated at w
        
        print("The maximum error in case of bell curve for N : {} and t : [{:.2f} , {:.2f}] is {:.3e}".format(N,low_t,up_t,max(abs(Y_t - Y))))
        
        """
        Plotting both true and the one obained using DFT
        """
        plot(w, abs(Y), lw = 2, label = '\n DFT -  N: {} & \n t : [{:.2f},{:.2f})'.format(N,low_t,up_t))     #Plots one obtained using DFT
        plot(w, abs(Y_t), '--',lw = 2, label = 'Analytical spectrum')                                        #Plots the one obtained analytically 
        legend(loc = 'upper right')                                                           
        
    else :
      
        Y = fftshift(fft(y(t)))/N    #If the function is not gaussian, the spectrum is obtained using fft and shifting the same. 
        plot(w, abs(Y), lw = 2)      #Plotting the magnitude spectrum 
     
    xlim(x_range)                                    #Limiting the range of values to be displayed to x_range 
    ylabel(r"$|Y|$", size = 12)                      #Setting y-label
    title(r"Spectrum of {}".format(y_ttl))           #Setting the title using <y_ttl>
    grid(True)                                       #Setting the grid 
    
    subplot(2,1,2)    #Sub-plot to plot the phase spectrum  
    
    if flag == 1 :                                               #If flag is 1, only those points whose mag > 1e-3 are shown in green 
        ii = where(abs(Y) > 1e-3)                                #Checking the points with mag > 1e-3   
        plot(w[ii], angle(Y[ii]),'go',lw=2,markersize = 5)       #Plotting the points in phase spectrum 
        
    elif flag == 2 :                                             #If the flag is 2, previous points in green and all the points are plotted in red    
        plot(w,angle(Y),'ro',lw=2)                               #Plotting all the points in phase spectrum  
        ii = where(abs(Y) > 1e-3)
        plot(w[ii], angle(Y[ii]),'go',lw=2,markersize = 5)   
        
    else :                                                       #If flag is not 1 or 2 all the points are plotted in red 
        plot(w,angle(Y),'ro',lw=2)
        
    xlim(x_range)                           #Limiting the range of values to be displayed to x_range  
    ylim([-4,4])                            #Setting y-limits                 
    ylabel(r"Phase of $Y$", size = 12)      #Setting y-label
    xlabel(r"$\omega$", size = 12)          #Setting x-label
    grid(True)                              #Setting the grid 
    savefig('EE2703_ASN9/' + file_name)     #Saving the plot at file location mentioned 
    
    return 0 

"""
Uncomment to test the functions fft and iftt
"""
"""    
x = rand(128)
X = fft(x)
y = ifft(X)
c_[x,y]
print('The error after inverting back is {:.2e}'.format(abs(x-y).max()))
print('Original and Inverted samples:')
print(c_[x,y])
"""
"""
Function definitions of mathematical functions to be fed to func_spectrum 
"""
y_1a = lambda t: sin(5*t)                      #Example 1 
y_1b = lambda t: (1+0.1*cos(t))*cos(10*t)      #Example 2 
y_2a = lambda t: (cos(t))**3                   #Question 2a
y_2b = lambda t: (sin(t))**3                   #Question 2b 
y_3  = lambda t: cos(20*t + 5*cos(t))          #Question 3
y_4  = lambda t: exp(-0.5*t**2)                #Question 4
  
print('Generating spectrums, Please Wait.\n') 

"""
Generating and storing the spectrum for different input values and functions 
"""                             
#Spectrum for Example -1  
t_1a = linspace(0,2*pi,128)           #Generating 128 points in 0 to 2*pi 
Y_1a = fft(y_1a(t_1a))                #Using fft to find the fourier transform 

figure()                              #Generates new figure 

subplot(2,1,1)                        #Sub-plot to plot the magnitude spectrum  
plot(abs(Y_1a),lw=2)                  #Magnitude plot
grid(True)                            #Setting grid 
title(r'Spectrum of $sin(5t)$')       #Phase plot 
ylabel(r'$|Y|$', size = 12)           #Setting ylabel

subplot(2,1,2)                        #Sub-plot to plot the phase spectrum  
ylabel(r'Phase of Y', size = 12)      #Setting ylabel
xlabel(r'$k$',size = 12)              #Setting xlabel
plot(unwrap(angle(Y_1a)),lw=2)        #Phase plot 
grid(True)                            #Setting grid 
savefig('EE2703_ASN9/ex1.jpg')        #Saving the plot at file_location 

func_spectrum(0, 2*pi, 128, y_1a, '$Sin(5t)$ improved', 'ex1_improved.jpg', [-10,10], 2)                          #Improved spectrum for Example-1
func_spectrum(0, 2*pi, 128, y_1b,'$cos(10t) + 0.1cos(t)cos(10t)$','ex2.jpg', [-15,15], 1)                         #Spectrum for Example-2
func_spectrum(-4*pi, 4*pi, 512, y_1b,'$cos(10t) + 0.1cos(t)cos(10t)$ improved', 'ex2_improved.jpg', [-15,15], 1)  #Improved sepctrum for Example-2
func_spectrum(-4*pi, 4*pi, 512, y_2a,'$cos^3(t)$', 'q2a.jpg', [-15,15], 1)                                        #Spectrum for Question 2a 
func_spectrum(-4*pi, 4*pi, 512, y_2b,'$sin^3(t)$', 'q2b.jpg', [-15,15], 1)                                        #Spectrum for Question 2b
func_spectrum(-4*pi, 4*pi, 512, y_3,'$cos(20*t+5*cos(t))$','q3.jpg',[-45,45], 1)                                  #Spectrum for Question 3
func_spectrum(-2*pi, 2*pi, 256, y_4,'$e^{-t^2/2}$','q4a.jpg',[-15,15], 1)                                         #Spectrum for Question 4 t : [-2pi,2pi], N: 256
func_spectrum(-3*pi, 3*pi, 512, y_4,'$e^{-t^2/2}$','q4b.jpg',[-15,15], 1)                                         #Spectrum for Question 4 t : [-3pi,3pi], N: 512

print('\nAll plots have been saved to EE2703_ASN9 directory')

"""
The following code has been used to estimate [-T/2,T/2) for a given value of N, 
such that error is less than thresh. Uncomment to run.
"""
"""
N = 256                                            #Setting the value of N                                         
T = 2*pi                                           #Setting the starting value of T, which is increased pi every iteration  
thresh = 1e-6                                      #Checks the precision in estimate                    
err = 100                                          #Checks for convergence 
     
for i in range(15):                                #Running the iteration for 15 turns 
    t = np.linspace(-T/2,T/2,N+1)                  #Generating N+1 t points in  [-T/2,T/2]              
    t = t[:-1]                                     #Removing the last point 
    w = np.linspace(-N*pi/T,N*pi/T,N+1)            #Generating the corresponding N+1 w points 
    w = w[:-1]                                     #Removing the last point 
    Y   = fft(ifftshift(y_4(t)))                   #Using fft to generate the spectrum, y and Y are of same length
    Y   = fftshift(Y)/N                            #Shifting the values such that Y[0] is in middle
    Y   = Y*T                                      #Multiplying with T 
    Y_t = exp(-w**2/2)*sqrt(2*pi)                  #True transform obtained analytically, is evaluated at w
    error = max(abs(Y - Y_t))                      #Finding the error 
    if error < thresh :                            #Checking if its less than thresh
        T_est = T 
        err   = error
        break
    T += 2*pi                                      #Incrementing T 

if err != 100 :                                    #Printing the results 
    print('\nThe estimeted value of T is {:.3f}'.format(T_est))
    print('The value of error between true and obtained spectrum for N:{}, T:{:.3f} is {:.3e}'.format(N,T_est,err))
else : 
    print("Failed to Converge")
""" 
"""
---------------------------------------End-Of-The-Code-----------------------------------------
"""
