""" 
EE2703 Applied Programming Lab - Jan-May 2021 
Purpose : Assignment-2 - Submission, 
Date    : 3rd March 2021, 
Author  : Sai Gautham Ravipati (EE19B053), 
Input   : Filepath of the LTspice net-list given as command line input, if netlist file is the same directory as that of code
          name of the file shall work. 
Output  : Circuit is solved for the variables and output is printed. 
Usage   : python3 EE2703_ASSIGN2_EE19B053.py <path_to_netlist>  

Please go through the following lines before execution : 
1. (V.. node_i node_j dc value) : node_j of a voltage source is at a lower potential than that of node_i. 
2. (I.. node_i node_j dc value) : Flow is from node_i to node_j for a current source .
3. All values are to be entered in scientific format, i.e. if 1K Ohm resistor value is to be entered as 1e3 in netlist. 
4. Please check if there is 'GND' node, if there isn't none execution is stopped. 
5. Please give the phase of AC source in degrees.  
6. In case of AC sources amplitudes are given as output, not rms.
"""

import sys              #Importing system module
from sys import exit    #Importing sys.exit() as exit()
from math import pi     #Importing pi
import cmath            #Importing cmath 
import numpy as np      #Importing numpy as np

"""
Constant variables used at several places in rest of the code, instead of hard-coding at every location.
"""
CKT = '.circuit'    #Defines the start of a circuit block  
END = '.end'        #Defines the end of a circuit block
AC  = '.ac'         #Keyword which checks for the pressence of AC sources in the circuit

"""
Checking if command is of the form "python3 <arg_0> <arg_1>" so that user doesn't pass anything apart from netlist file_path
as argument.  
"""
if len(sys.argv) != 2 : 
    print('Please check the number of arguments passed.\nSample usage: "python3 %s <path-to-netlist>"' % sys.argv[0])
    exit()
        
"""
Creating a class for passive elements RLC
""" 
class RLC:
    def __init__(self,name,node_i,node_j,value):
    
        """
        Inputs -
        name   : Name of the element - Ex : R1,R2,C1 etc. 
        node_i : One of the node between which the element is present 
        node_j : Other node between which the element is present 
        value  : Value of the element parsed from the token of the netlist 
        """ 
        self.name = name           
        self.node_i = node_i       
        self.node_j = node_j
        self.value = float(value) #Gives value of the element after converting string input to float  
        
"""
Creating a class for independent sources
"""         
class XS:
    def __init__(self, name, node_i, node_j,value, phase = 0):
    
        """
        Inputs - 
        name   : Name of the element - Ex : V1,I2 etc. 
        node_i : One of the node between which the element is present 
        node_j : Other node between which the element is present 
        value  : Value of the element parsed from the token of the netlist
        phase  : The phase of the source and is an attribute corresponding to AC sources.
                 Its default value is set to zero unless passed by the user. 
        """
        self.name = name
        self.value = value
        self.node_i = node_i
        self.node_j = node_j
        self.value = cmath.rect(float(value),phase*pi/180)  #Gives complex value of the element after accomodating both magnitude and phase 

"""
The following lines contain three functions which are being used in the main function of 
the program
"""
def file_process() :

    """
    1) Netlist file path is taken through command line argument and if there are issues related to file input, user is alerted
       about the same. 
    2) It is verified if '.netlist' file is given as input or some other file. 
    3) Lines in the file are read using f.readlines() 
    4) start_idx, end_idx store all the indices of the lines whose first charecters match with CKT and END respectively.  
    4) ac_idx store all the indices of the lines whose first charecters match with AC.  
    5) For a proper netlist file len(start_idx) and len(end_idx) should be one, so if there are multiple or no '.circuit' or 
       '.end',or '.' is in the middle of the line, then execution is stopped, further if '.end' occurs before '.circuit' 
       execution is stopped. 
    6) After extrating the lines between '.circuit' and '.end', comments are removed using line.split("#")[0] and .strip is used
       to remove white spaces as well as line breaks, after this is done only non-empty lines are stored in the array ckt_block.  
    7) If there is no string '.ac' present in the netlist after '.end' the circuit is considered as a DC circuit and frequency f 
       is set to zero. 
    8) If there '.ac' string present then frequency is obtained after parsing the line to obtain frequency.
    9) If multiple lines with '.ac' are present then the execution is terminated.
    """ 
    try :                                                                                                             #Checking for file input error
        file_path = sys.argv[1]                                                                                       #Considering file path as command line argument                               
        if(file_path.endswith('.netlist')) :                                                                          #Checking if it's a .netlist file 
            with open(file_path) as f:        
                lines = f.readlines()   
                start_idx = []; end_idx = []; ac_idx = [] 
                start_idx = [lines.index(line) for line in lines if line[0:len(CKT)] == CKT]                          #Storing indices of all lines that start with '.circuit'
                end_idx   = [lines.index(line) for line in lines if line[0:len(END)] == END]                          #Storing indices of all lines that end with '.end'
                ac_idx    = [lines.index(line) for line in lines if line[0:len(AC)] == AC]   
                if(len(start_idx) != 1 or len(end_idx) != 1) :                                                        #Checking that there are only one '.circuit' and '.end' in the file 
                    print('The number of lines containing "%s" and "%s" at the start of the line are not equal to one.' % (CKT,END))
                    print('This may be due to stray whitespaces before "%s" & "%s" or absence, repitition of these words.'% (CKT,END))
                    exit()
                if(start_idx[0] >= end_idx[0]) :                                                                      #Checking that '.circuit' occurs before '.end' in the file 
                    print('Please check the netlist file, %s should occur before %s.' % (CKT,END))
                    exit()
                if(start_idx[0] >= end_idx[0]) :                                                                      #Checking that '.circuit' occurs before '.end' in the file 
                    print('Please check the netlist file, %s should occur before %s.' % (CKT,END))
                    exit() 
                block = lines[start_idx[0]+1:end_idx[0]]                                                              #Extracting the circuit block between '.circuit' and '.end' 
                ckt_block = [line.split("#")[0].strip() for line in block if line.split("#")[0].strip() != '']        #Removing comments and empty lines   
                if(len(ac_idx) == 0) :                                                                                #Checking if there is '.ac' present or not in the netlist 
                    f = 0                                                                                             #If there isn't none setting frequency to zero 
                elif (len(ac_idx) == 1) :                                                                             #Checking if '.ac' is present after '.end' in the circuit 
                    if(end_idx[0] >= ac_idx[0]):  
                        print('"%s" should occur after "%s" in the circuit.' % (AC,END))
                        exit()
                    ac_token = lines[ac_idx[0]].split("#")[0].strip()                                                 #Extracting frequency from the last word of line starting with '.ac'
                    f = float(ac_token.split()[-1])                           
                else : 
                    print('This program supports only single frequency operation, so check that only one "%s" is present' % (AC))
                    exit()   
            return ckt_block,f                                                                                        #ckt_block with lines between '.start' and '.end' and f is returned          
        else :                                    
            print('Please give ".netlist" file as input')
            exit()
    except IOError:                         
        print('Invalid File')
        exit() 

def parse(line) :

    """
    1) A line in ckt_block is given as input.
    2) Line by line parsing is done and nodes are stored into nodes list if not already present in the list.
    3) Based on the first word, corresponding object is stored in a dictionary with key corresponding to that 
       element with properties defined in the classes mentioned previously. 
    4) If elements are not either of R,L,C,V,I execution is stopped.
    """ 
    ckt_token = line.split()             #Splitiing the lines into separate words 
    if ckt_token[1] not in nodes :       #Adding unique nodes to the list of nodes 
        nodes.append(ckt_token[1])
    if ckt_token[2] not in nodes :
        nodes.append(ckt_token[2])
    element_type = ckt_token[0][0]       #First letter of the first word is considered as element type, and the objects with properties are stored in a dictionary 
    if element_type == 'R' :             
        compt_dict['R'].append(RLC(ckt_token[0], ckt_token[1], ckt_token[2], ckt_token[3]))
    elif element_type == 'L' :
        compt_dict['L'].append(RLC(ckt_token[0], ckt_token[1], ckt_token[2], ckt_token[3])) 
    elif element_type == 'C' :
        compt_dict['C'].append(RLC(ckt_token[0], ckt_token[1], ckt_token[2], ckt_token[3])) 
    elif element_type == 'V' : 
            if ckt_token[3] == "ac":     #If 4th word of a line is 'ac' then Vp-p is converted to amplitudde to pass in as property to the object 
                compt_dict['V'].append(XS(ckt_token[0],ckt_token[1], ckt_token[2], 0.5*float(ckt_token[4]), float(ckt_token[5])))   
            else :                       #For DC case the last element is taken as value
                compt_dict['V'].append(XS(ckt_token[0],ckt_token[1], ckt_token[2], ckt_token[-1]))  
    elif element_type == 'I' : 
            if ckt_token[3] == "ac":
                compt_dict['I'].append(XS(ckt_token[0],ckt_token[1], ckt_token[2], 0.5*float(ckt_token[4]), float(ckt_token[5])))   
            else :
                compt_dict['I'].append(XS(ckt_token[0],ckt_token[1], ckt_token[2], ckt_token[-1]))    
    else :                               #If elements other than these is encountered, execution is stopped   
        print('Elements - R,L,C,V,I are only accepted')  
        exit()
    return 0         
               
def M_b_generator(compt_dict,node_dict,f) :

    """
    Fills M,b matrices in accordance  to the stamps of Modified Nodal Analysis
    Input  : compt_dict : Dictionary of component objects
    node_dict  : Dictionary of nodes where each node as a key is assigned a number 
    f          : The frequency of operation of the circuit  
    output : Numpy arrays M,b   
    """
    w = 2*pi*f                                                          #Conversion of frequency to angular frequency 
    N_n = len(node_dict.keys())                                         #Gives the number of nodes in the circuit 
    V_n = len(compt_dict['V'])                                          #Gives the number of Voltage sources in the circuit 
    L_n = len(compt_dict['L'])                                          #Gives the numner of inductors in the circuit 
    M_shape = (N_n + V_n, N_n + V_n)                                    #Setting appropriate dimensions of the matrics 
    b_shape = (N_n + V_n)
    if w == 0:                                                          #In case of DC operation in presence of inductors, current thorugh them is also taken as a variable 
        M_shape = (N_n + V_n + L_n, N_n + V_n + L_n)                    #Dimesions are set appropriately considering the same 
        b_shape = (N_n + V_n + L_n)               

    
    M = np.zeros(M_shape, np.complex)                                   #Initializing M,b with zeros and dtype of entries as complex 
    b = np.zeros(b_shape, np.complex)
    M[0][0] = 1                                                         #Equation corresponding to the ground node 
    for R in compt_dict['R']:                                           #For every resistor in the circuit, entries in M are made after checking if any of node is ground 
        if R.node_i != 'GND':
            M[node_dict[R.node_i],node_dict[R.node_i]] += 1/R.value
            M[node_dict[R.node_i],node_dict[R.node_j]] -= 1/R.value
        if R.node_j != 'GND':
            M[node_dict[R.node_j],node_dict[R.node_j]] += 1/R.value
            M[node_dict[R.node_j],node_dict[R.node_i]] -= 1/R.value
    rel_pos = 0                                                         #Relative position of an L related entry for DC steady state after nodes, voltage sources
    for L in compt_dict['L']:                                       
        try :                                                           #For every Inductor in an AC circuit, entries in M are made after checking if any of node is ground 
            if L.node_i != 'GND':
                M[node_dict[L.node_i],node_dict[L.node_i]] -= (1j)/(w*L.value)
                M[node_dict[L.node_i],node_dict[L.node_j]] += (1j)/(w*L.value)
            if L.node_j != 'GND':
                M[node_dict[L.node_j],node_dict[L.node_j]] -= (1j)/(w*L.value)
                M[node_dict[L.node_j],node_dict[L.node_i]] += (1j)/(w*L.value)
        except ZeroDivisionError:                                       #If the circuit is DC corresponding changes are made in the matrices M,b as inductor acts as a short at steady state. 
                if L.node_i != 'GND':                                    
                    M[node_dict[L.node_i],N_n + V_n + rel_pos] += 1 
                    M[N_n + V_n + rel_pos,node_dict[L.node_i]] -= 1
                    b[N_n + V_n + rel_pos] = 0
                if L.node_j != 'GND':
                    M[node_dict[L.node_j],N_n + V_n + rel_pos] -= 1 
                    M[N_n + V_n + rel_pos,node_dict[L.node_j]] += 1
                    b[N_n + V_n + rel_pos] = 0 
    for C in compt_dict['C']:                                           #For every Capacitor in the circuit, entries in M are made after checking if any of node is ground 
        if C.node_i != 'GND': 
            M[node_dict[C.node_i],node_dict[C.node_i]] += (1j*w*C.value)
            M[node_dict[C.node_i],node_dict[C.node_j]] -= (1j*w*C.value)
        if C.node_j != 'GND':
            M[node_dict[C.node_j],node_dict[C.node_j]] += (1j*w*C.value)
            M[node_dict[C.node_j],node_dict[C.node_i]] -= (1j*w*C.value)
    for I in compt_dict['I'] :                                          #For every current sorce in the circuit, entries in b are made after checking if any of node is ground                      
        if I.node_i != 'GND':
            b[node_dict[I.node_i]] = -I.value
        if I.node_j != 'GND':
            b[node_dict[I.node_j]] = +I.value
    rel_pos = 0                                                         #Relative position of an V related entry after nodes                                     
    for V in compt_dict['V'] :                                          #For every voltage source in the circuit, entries in M,b are made after checking if any of node is ground 
        if V.node_i != 'GND':
           M[node_dict[V.node_i],N_n + rel_pos] += 1
        if V.node_j != 'GND':
           M[node_dict[V.node_j],N_n + rel_pos] -= 1   
        M[N_n + rel_pos,node_dict[V.node_i]] += 1                    
        M[N_n + rel_pos,node_dict[V.node_j]] -= 1
        b[N_n + rel_pos] = V.value 
        rel_pos += 1
    return M,b                                                          #Matrices M,b are returned 
    
"""
The main functions calls each of the previous functions in a sequential manner and finally M,b are filled.
Using linalg.solve variable vector x is obtained. Finally required data is printed.  
"""                
if __name__ == "__main__":
    ckt_block,f = file_process()                                        #Unpacking returned values from file_process() to ckt_block and f  
    nodes = []                                                          #Initializing node list 
    compt_dict = { 'R': [], 'L': [], 'C': [], 'V': [], 'I': []}         #Initialising component dictionary 
    [parse(line) for line in ckt_block]                                 #Parsing line by line and line and filling up the dictionary as well as list 
    try:                       
        nodes.remove('GND')                                             #Checking for 'GND' node in the circuit and execution is stopped in absence 
        nodes = ['GND'] + sorted(nodes)
    except ValueError :                                              
        print('Please check if there is "GND" node in the netlist file')
        exit()
    node_dict = { node: nodes.index(node) for node in nodes }           #Assigning number to each entry in nodes list with key as node names 
    M,b = M_b_generator(compt_dict,node_dict,f)                         #Unpacking returned values to M,b    
    try:                                                                 
        x = np.linalg.solve(M,b)                                        #Solving for x 
    except Exception:
        print('Please check the netlist file, singular matrix encountered')
        exit()
    N_n = len(node_dict.keys())  
    V_n = len(compt_dict['V'])
    L_n = len(compt_dict['L'])
    if not(f) :                                                         #Printing the information in structured manner 
        for i in range(1,N_n):
            node_name = [key for key in node_dict.keys() if node_dict[key] == i]
            print("The voltage at node {} is {:.2e} V".format(node_name[0],x[i].real))
        for i in range(V_n):
            print('The current through source {} is {:.2e} A'.format(compt_dict['V'][i].name,x[N_n+i].real))
        for i in range(L_n):
            print("The current through inductor {} is {:.2e} A".format(compt_dict['L'][i].name,x[N_n + V_n + i].real))
    else:
        for i in range(1,N_n):
            node_name = [key for key in node_dict.keys() if node_dict[key] == i]
            mag,phase_rad = cmath.polar(x[i]) 
            phase         = phase_rad*180/pi
            print("The voltage at node {} is magnitude : {:.2e} V  phase: {:.4f}".format(node_name[0],mag,phase))
        for i in range(V_n):
            mag,phase_rad = cmath.polar(x[N_n + i]) 
            phase         = phase_rad*180/pi
            print('The current through source {} is magnitude {:.2e} A phase: {:.4f}'.format(compt_dict['V'][i].name,mag,phase))  
     
"""
----------------------------------------------------End-Of-The-Code---------------------------------------------------------
"""  
    
