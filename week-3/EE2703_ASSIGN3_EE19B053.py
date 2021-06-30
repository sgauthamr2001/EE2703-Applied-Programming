""" 
EE2703 Applied Programming Lab - Jan-May 2021 
Purpose : Assignment-3 - Submission, 
Date    : 3rd March 2021, 
Author  : Sai Gautham Ravipati (EE19B053), 
Input   : None 
Output  : Different plots based on various aspects, Location of minimum in the contour, Equality of Matrices, 
Usage   : python3 EE2703_ASSIGN3_EE19B053.py

Note : 
Please run the script 'generate_data.py' before running this script.   
"""

"""
Importing the neccasary modules.
"""
import numpy as np 
import matplotlib.pyplot as plt
import scipy.special as sp
import scipy
from scipy.linalg import lstsq

"""
The following function takes time series, A, B as inputs and returns the value computed by the equation.
sp.jn(2,t) is a bessel function of first kind and of order 2. 
""" 
def g(t,A = 1.05 ,B = -0.105) : 
    return A*sp.jn(2,t) + B*t

"""
The data in the file 'fitting.dat' is being loaded into the variable 'data', user is informed if no such 
file is generated. 
"""   
try :
    data = np.loadtxt('fitting.dat')
except IOError :
    print('File not found. Please check if data file is generated.') 
    exit()

"""
The first column of the data is being assigned to a new variable t, which contains time series values, 
while rest of the 9 columns are stored into variable f, f_true is a column vector containing data of 
ground truth values generated using the previously defined function.
"""
t = data[:,0]
f = data[:,1:]
f_true = g(t)

"""
The following block of code plots figure_0 which includes the following graphs :
1) f(t) + noise versus time for different values of sigma. 
2) True value versus time. 
"""

#Plot for Q3
plt.figure(0)
plt.plot(t,f)
sigma = np.logspace(-1,-3,9)
[plt.plot(t,f[:,i],label = r'$\sigma_%d$ =%.4f' % (i,sigma[i])) for i in range(9)]
plt.plot(t,f_true,'k',label = r'True Value')
plt.title(r'Variation with different amounts of noise added (Q3)')
plt.xlabel(r'$t \rightarrow$',size = 10)
plt.ylabel(r'$f(t) + noise \rightarrow$',size=10)
plt.legend(loc = 'upper right')
plt.grid()

"""
The following block of code plots figure_1 :
1) Plots the error for every 5th point with noise produced using sigma_0 and ground truth.
2) Plots the original curve. 
"""
#Plot for Q4
f_dev = f[:,0]
stdev = np.std(f_dev - f_true)
plt.figure(1)
plt.title('Datapoints for $\sigma_i$ =' + str(sigma[0]) + ' with exact function (Q5)')
plt.xlabel(r'$t \rightarrow$',size=10)
plt.errorbar(t[::5],f[::5,0],stdev,fmt='ro',label = 'errorbar')
plt.plot(t,f_true, label = 'original')
plt.legend(loc = 'upper right')
plt.grid()

"""
The following block of code checks if the vector generated using the vector product of M,p
is same as that of the one generated using the function defined previously using bessel 
function.   
"""

#Code for Q6
m = t.shape[0]
M = np.zeros((m,2))
A = 1.05; B = -0.105
p = np.array([A,B])
M[:,0] = sp.jn(2,t)
M[:,1] = t
f_Mp = np.matmul(M,p)

if np.array_equal(f_Mp,f_true) : 
    print('Q.6: Both the matrices are equal') 
else : 
    print('Q.6: Both matrices are not equal')

"""
This part of code computes the mean squared error for A : {0,0.1,...,2} and B : {-0.2,0.19,...,0},
using the column data loaded previously, and value generated using the previous function.      
"""

#Code for Q7
A = np.arange(0,2.1,0.1)
B = np.arange(-0.2,0.01,0.01) 
eps = np.zeros((21,21,9))
for k in range(9):
    for i in range(21):
        for j in range(21):
            eps[i][j][k] = np.sum(np.square(f[:,k] - g(t,A[i],B[j])))/101 

"""
The following block of code plots figure_2 which is the contour plot of MSE for different values 
of A,B using the first column of the previous data.

"""

#Plot for Q8
plt.figure(2)
plot = plt.contour(A,B,eps[:,:,0],15)
plt.xlabel(r'$A \rightarrow$',size=10)
plt.ylabel(r'$B \rightarrow$',size=10)
plt.title('Contour plot of $\epsilon_{ij}$ (Q8)')
plt.clabel(plot,inline=1,fontsize=10)
i,j = np.unravel_index(np.argmin(eps[:,:,0], axis=None), eps[:,:,0].shape)
print('Location of the minimum in contour : (%f,%f)' % (A[i],B[j]))
plt.plot(A[i],B[j],'bo',label = 'Minimum')
plt.plot(1.05,-0.105,'ro',label = 'Actual Value')
plt.legend(loc='upper right')
plt.grid()


"""
Values of A,B are estimated using least square approximation by taking A = 1.05 and B = -0.105
and error in estimated values with respect to A,B is plotted in figure_3 and same is plotted in 
log-log scale in figure_4. First column of error reflects error in A and second column reflects 
error in B.

"""
estimate = []
[estimate.append(lstsq(M,f[:,i])[0]) for i in range(9)] 
error   = (estimate - p)**2 

#Plot for Q10
plt.figure(3)
plt.plot(sigma,error[:,0],'--ro',label = '$A_{err}$')
plt.plot(sigma,error[:,1],'--bo',label = '$B_{err}$')
plt.title('Variation of error with noise (Q10)')
plt.ylabel(r'MS Error in estimation of A,B $\rightarrow$',fontsize=10)
plt.xlabel(r'$\sigma_{n}\rightarrow$',fontsize=10)
plt.legend(loc='upper left')
plt.grid()

#Plot for Q11
plt.figure(4)
plt.stem(sigma,error[:,0],'--r')
plt.stem(sigma,error[:,1],'--g')
plt.loglog(sigma,error[:,0],'ro',label = '$A_{err}$')
plt.loglog(sigma,error[:,1],'go',label = '$B_{err}$')
plt.title("Variation of error with Noise on loglog scale (Q11)")
plt.xlabel(r'$\sigma_{n} \rightarrow$',fontsize=10)
plt.ylabel(r'MS Error in estimation of A,B $\rightarrow$',fontsize=10)
plt.legend(loc='upper left')
plt.grid()
plt.show()

"""
---------------------------------------End-Of-The-Code-----------------------------------------
"""



