""" 
EE2703 Applied Programming Lab - Jan-May 2021 
Purpose : Assignment-6 - Submission, 
Date    : 10th April 2021, 
Author  : Sai Gautham Ravipati (EE19B053), 
Input   : Arguments like n,M,nk,u0,p,Msig, 
Output  : Plots of electron density, intensity, phase space and intensity data, 
Usage   : python3 EE2703_ASSIGN6_EE19B053.py --n n_val --M M_val --p p_val --nk nk_val --u0 v_thsh --Msig sig_val
"""
"""
Importing necessary modules
Requirements : 
tabulate module - <pip install tabulate> 
"""
import numpy as np                       #Numpy Module  
import matplotlib.pyplot as plt          #Matplotlib for plottimg 
import argparse                          #Argparse to parse the input arguments 
from tabulate import tabulate            #Helps in displaying the intensity data in a structured format 
import os                                #To create directory to save plots and intensity data 
from tqdm import tqdm                    #Displays speed of the update step 

os.makedirs('EE2703_ASN6', exist_ok = True)    #Creating a directory to store plots and data 
	
"""
Function to fill the entries of I,X,V.
"""
def TL_sim(n,M,nk,u0,p,Msig) : 
	"""
	The following piece of code runs the simulation, by taking in the input arguments as
	those parsed previously and returns the populated arrayS I,V,X along with xx,dx,u. 
	"""
	xx = np.zeros(n*M)    #Initializing e- postion array with zeros of size nM
	u  = np.zeros(n*M)    #Electron Velocity 
	dx = np.zeros(n*M)    #Electron Displacement in current turn 

	I = []    #Holds Intensity of Emitted Light 
	X = []    #Holds Electron Position 
	V = []    #Holds Electron Velocity 

	ii = np.where(xx > 0)    #Finding the locations of e- present in the chamber 

	for k in tqdm(range(0, nk), desc = "Update step in Progress") : 	
	
		"""
		The following block of code runs iteratively and updates the variables in 
		each turn.
		"""
		dx[ii] = u[ii] + 0.5    #Updating the displacement in each turn, of those e- in the chamber 
		xx[ii] += dx[ii]        #Advancing e- position in each turn 
		u[ii]  += 1             #Advancing e- velocity in each turn  
		
		jj = np.where(xx >= n)    #Determing the positions of e- hitting the anode
		xx[jj] = 0                #Setting position to zero 
		u[jj]  = 0                #Setting velocity to zero 
		dx[jj] = 0                #Setting displacement to zero   
		
		kk = np.where(u >= u0)[0]                      #Indices correspoding to the ionised e- 
		ll = np.where(np.random.rand(len(kk))<=p)        
		kl = kk[ll]                                    #Indices of energetic e- that suffer a collision  

		u[kl] = 0                                      #Setting the velocity of these e- to zero assuming elastic collision       
		xx[kl] -= dx[kl]*np.random.rand(len(kl))       #Updating xx array, assuming location of collision randomly between previous and current locations  
		                               
		I.extend(xx[kl].tolist())    #Adding photon intenities to the I vector 
		
		m  = int(np.random.randn()*Msig+M)    #Injecting new e- 
		mm = np.where(xx == 0)                #Finding the unsued indices 
		
		min_idx = min(len(mm[0]),m)           #Finding the lower among number of e- and number of slots availbale  
		xx[mm[0][:min_idx]] = 1               #Setting the position to 1 for these e- 
		u[mm[0][:min_idx]]  = 0               #Setting the velocity of these e- to 0 
		
		ii = np.where(xx > 0)                 #Finding the indices of e- in the chamber 
		X.extend(xx[ii].tolist())             #Populating X,V in each turn 
		V.extend(u[ii].tolist())
		
	return I,X,V,xx,u,dx

"""
The following functions help while parsing the input arguments, and check for the following : 
1) The value of p entered lies between 0 and 1. 
2) Only positive integers are given as input for several arguments.  
3) The value of sigma entered is positive.
"""
def pf(x) :

	"""
	Checks if the probability that ionization will occur is an integer 
	"""
	try :
		p = float(x)
	except ValueError :
		raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))         #Error is returned if conversion to float is void 
	if p < 0.0 or p > 1.0 :                                                                    #Checking if p is [0,1]
		raise argparse.ArgumentTypeError("Probability %r not in range [0.0, 1.0]"%(x,))    
	return p
    
def nat(x) :

	"""
	Checks if the input argument is a natural number 
	"""
	try :  
		ix = int(x)
	except ValueError :
		raise argparse.ArgumentTypeError("%r not an integer literal" % (x,))               #Error is returned if conversion to int is void
	if ix <= 0 :                                                                               #Checking if a positive integer is entered
		raise argparse.ArgumentTypeError("%s is an invalid positive integer value, please enter a postive integer value" % x)
	return ix
	
def sig(x):

	"""
	Checks if the variance value entered is positive
	"""
	try : 
		s = float(x)                                                                                                                              
	except ValueError :
		raise argparse.ArgumentTypeError("%r not an floating-point literal" % (x,))       #Error is returned if conversion to int is void       
	if s <= 0 :                                                                               #Checking if variance value is positive 
		raise argparse.ArgumentTypeError("%s is an invalid positive value, please enter a postive integer value" % x)
	return s
	
"""
Parsing the input arguments, and value is parsed using previously defined functions  
"""
parser = argparse.ArgumentParser()
parser.add_argument("--n", metavar = 'n_val', type = nat, default = 100, help = 'Spatial grid size')                           
parser.add_argument("--M", metavar = 'M_val', type = nat, default = 5,   help = 'Number of electrons injected per turn') 
parser.add_argument("--p", metavar = 'p_val', type = pf,  default = 0.25, help = 'Probability that ionization will occur')                           
parser.add_argument("--nk", metavar = 'nk_val', type = nat, default = 500, help ='Number of turns to simulate')             
parser.add_argument("--u0", metavar = 'v_thsh', type = sig, default = 5,   help = 'Threshold velocity') 
parser.add_argument("--Msig", metavar = 'sig_val', type = sig, default = 2, help = 'Standard deviation of Electron distribution') 

args = parser.parse_args()
[n,M,p,nk,u0,Msig] = [args.n,args.M,args.p,args.nk,args.u0,args.Msig]    #Storing the parsed values in their corresponding holders  

print('Current values (n, M, p, nk, u0, Msig) - (%d, %d, %.2f, %d, %.2f, %.2f)' % (n,M,p,nk,u0,Msig))    #Being informed of the arguments passed or the default once   
print('To check usage, pass -h while execution')

I,X,V,xx,u,dx = TL_sim(n,M,nk,u0,p,Msig)    #Storing the populated vectors returned by the function call into the holders 
print('Please Wait, plots and data being generated.')

"""
Code block to plot electron density histogram, plot is saved in directory EE2703_ASN6
"""
plt.figure(0)
plt.hist(X,n,[0,n],edgecolor = 'black')
plt.title("Electron Density",fontsize = 14)
plt.ylabel("Electron Density", fontsize = 12)
plt.xlabel("position", fontsize = 12)
plt.savefig('EE2703_ASN6/fig1.png') 

"""
Code block to plot Intensity histogram
"""
plt.figure(1)
count,bins,_ = plt.hist(I,n,[0,n],edgecolor = 'black')    #Bin locations are extracted to store intensity data 
plt.title("Intensity of emitted light",fontsize = 14)
plt.ylabel("I", fontsize = 12)
plt.xlabel("position", fontsize = 12)
plt.savefig('EE2703_ASN6/fig2.png')

"""
Code block to plot e- phase space
"""
plt.figure(2)
plt.plot(X,V,'b.')
plt.title("Electron Phase space",fontsize = 14)
plt.xlabel("position", fontsize = 12)
plt.ylabel("velocity", fontsize = 12)
plt.savefig('EE2703_ASN6/fig3.png')

xpos = 0.5*(bins[0:-1]+bins[1:])         #Finding the midpoint values of the bins          

f = open('EE2703_ASN6/data.txt', 'w+')   #Intensity Data is written into this file 
data = []                                      
[data.append([x,y]) for x, y in zip(xpos, count)]
print("Intensity Data : \n", file = f)
print(tabulate((data), headers=['xpos', 'count'],tablefmt='orgtbl'), file = f)    #Using tabulate module to write the data as table
print('All the plots and intensity data have been saved in EE2703_ASN6/ directory.') 

"""
---------------------------------------End-Of-The-Code-----------------------------------------
"""
