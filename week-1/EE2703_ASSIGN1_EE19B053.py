""" 
EE2703 Applied Programming Lab - Jan-May 2021 
Purpose : Assignment-1 - Submission, 
Date : 17th February 2021, 
Author : Sai Gautham Ravipati (EE19B053), 
Input : Filepath of the LTspice net-list given as command line input, if netlist file is the same directory as that of code
        name of the file shall work. 
Output : Circuit definition traversed from the last element to first and each line with words in reversed order is printed. 
Usage : python3 asn1.py <path_to_netlist>  
"""

import sys    #Importing system module
from sys import exit

"""
Constant variables used at several places in rest of the code, instead of hard-coding at every location.
"""

CKT = '.circuit'    #Defines the start of a circuit block  
END = '.end'        #Defines the end of a circuit block

"""
Checking if command is of the form "python3 <arg_0> <arg_1>" so that user doesn't pass anything apart from netlist file_path
as argument.  
"""

if len(sys.argv) != 2 : 
    print('Please check the number of arguments passed.\nSample usage: "python3 %s <path-to-netlist>"' % sys.argv[0])
    exit()
    
"""
1) The following function takes one particular line of the netlist free of all comments, between the lines '.circuit' and 
   '.end' as input. 
2) After splitting the line into an array the first and last words of the array are considered to be type of the  element 
   and value of the element respectively. If the length of the array is not 4,5,6 the exceution is stopped pointing for an  
   invalid format in the netlist for elements being considered. 
3) If the length of the array is either of 4,5,6 corresponing node locations are checked if the names are alphanumeric and 
   if they are so, the reversed array of words is returned.
4) If the node names are not alphanumeric then execution is stopped and the user is informed of the same. 
""" 
    
def reverse_gen(ckt_line) : 
    ckt_line = ckt_line.split(' ')                                          #Splitting the input line into word array
    ckt_element = ckt_line[0]                                               #Taking the first word as element type 
    value = ckt_line[-1]                                                    #Taking the last word as element value
    if(len(ckt_line) == 4) :                                                #Checking if it's an R,L,C, Independent V,I sources 
        if(not(ckt_line[1].isalnum() and ckt_line[2].isalnum())) :          #Verifying whether the node names are alphanumeric and assigning them to variables
            print('Node names should be alphanumeric!')                     #If node names are not alphaumeric, execution is stopped, alerting the user 
            exit()
        node_i = ckt_line[1]
        node_j = ckt_line[2]
        return [value,node_j,node_i,ckt_element]                            #Reversed array of words is returned
    elif(len(ckt_line) == 5) :                                              #Checking if it's a CCVS or CCCS 
        if((not(ckt_line[1].isalnum() and ckt_line[2].isalnum()))) :        #Verifying whether the node names are alphanumeric and assigning nodes and voltage source to variables 
            print('Node names should be alphanumeric!')                     #If node names are not alphaumeric, execution is stopped, alerting the user 
            exit()
        node_i = ckt_line[1]
        node_j = ckt_line[2]
        V_m = ckt_line[3]
        return [value,V_m,node_j,node_i,ckt_element]                        #Reversed array of words is returned
    elif (len(ckt_line) == 6) :                                             #Checking if it's a VCVS or VCCS 
        if(not(all(node_name.isalnum() for node_name in ckt_line[1:5]))) :  #Verifying whether the node names are alphanumeric and assigning them to variables
            print('Node names should be alphanumeric!')                     #If node names are not alphaumeric, execution is stopped, alerting the user 
            exit()    
        node_i = ckt_line[1]
        node_j = ckt_line[2]
        node_m = ckt_line[3]
        node_n = ckt_line[4]
        return [value,node_n,node_m,node_j,node_i,ckt_element]              #Reversed array of words is returned 
    else :                                                                  #If the number of words is not 4,5 or 6 execution is stopped indicating the user
         print('Netlist may have elements other than R,L,C,V,I,E,G,H,F or invalid entries, Please check!')
         exit()
                 
"""
1) The following is the implementation of main function. 
2) Netlist file path is taken through command line argument and if there are issues related to file input, user is alerted
   about the same. 
3) It is verified if '.netlist' file is given as input or some other file. 
4) Lines in the file are read using f.readlines() 
5) start_idx, end_idx store all the indices of the lines whose first charecters match with CKT and END respectively.  
6) For a proper netlist file len(start_idx) and len(end_idx) should be one, so if there are multiple or no '.circuit' or 
   '.end',or '.' is in the middle of the line, then execution is stopped, further if '.end' occurs before '.circuit' 
   execution is stopped. 
7) After extrating the lines between '.circuit' and '.end', comments are removed using line.split("#")[0] and .strip is used
   to remove white spaces as well as line breaks, after this is done only non-empty lines are stored in the array ckt_block.  
8) Array is reversed and each of the line in the array is passed through the function reverse_gen, the returned values are 
   stored in the array ckt_reversed. 
9) ckt_reversed is printed with lines joined using line break between them.  
  
"""               
                 
if __name__ == "__main__":
    try :                                                                                                             #Checking for file input error
        file_path = sys.argv[1]                                                                                       #Considering file path as command line argument                               
        if(file_path.endswith('.netlist')) :                                                                          #Checking if it's a .netlist file 
            with open(file_path) as f:        
                lines = f.readlines()   
                start_idx = []; end_idx = [] 
                start_idx = [lines.index(line) for line in lines if line[0:len(CKT)] == CKT]                          #Storing indices of all lines that start with '.circuit'
                end_idx   = [lines.index(line) for line in lines if line[0:len(END)] == END]                          #Storing indices of all lines that end with '.end'
                if(len(start_idx) != 1 or len(end_idx) != 1) :                                                        #Checking that there are only one '.circuit' and '.end' in the file 
                    print('The number of lines containing "%s" and "%s" at the start of the line are not equal to one.' % (CKT,END))
                    print('This may be due to stray whitespaces before "%s" & "%s" or absence, repitition of these words.'% (CKT,END))
                    exit()
                if(start_idx[0] >= end_idx[0]) :                                                                      #Checking that '.circuit' occurs before '.end' in the file 
                    print('Please check the netlist file, %s should occur before %s.' % (CKT,END))
                    exit()
                block = lines[start_idx[0]+1:end_idx[0]]                                                              #Extracting the circuit block between '.circuit' and '.end' 
                ckt_block = [line.split("#")[0].strip() for line in block if line.split("#")[0].strip() != '']        #Removing comments and empty lines   
                ckt_reversed = [ ' '.join(reverse_gen(element)) for element in reversed(ckt_block)]                   #Using the function to generate the lines based on reversed definition 
                print("The circuit definition traversed from the last element to first and each line with words in reversed order is :")
                print("\n".join(ckt_reversed))
        else :                                    
            print('Please give ".netlist" file as input')
    except IOError:                         
        print('Invalid File')
        exit()
     
"""
----------------------------------------------------End-Of-The-Code---------------------------------------------------------
"""
