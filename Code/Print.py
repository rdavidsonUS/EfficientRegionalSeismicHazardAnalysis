# Input binary files generated using SolveMe code (i.e.,  Bianry files in which we have Pj's, Zj's, e+ir, e-ir, and objective function) 
# Output text file readable by MATLAB for Pj's, Zj's, e+ir, e-ir, and objective function 

##############################################################################################
# The body of code works as follows:   
       # i) Read binary files 
       # ii) post-processing for P and Z
       # iii) Save P,Z and E in txt

# User Input: Please define the user input 
# Notice: A user need to modify part iii based on the numebr of candidates and Jred.  
##############################################################################################

## LIBRARIES  
import pickle                     
import numpy as np
import scipy.io as sio

## User Input
Num_candidate=[2003, 507]                                                                                     # Number of candidates coresponding to finename [2003, 507, 4002, 5001, 1003, 132, 3005]     
Jred_values= [300, 500]                                                                                       # Jred
Num_sites=35                                                                                                  # Please enter the number of sites

# i) read binary files 
pkl_file = open('P', 'rb')
P = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('Z', 'rb')
Z = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('e_plus', 'rb')
E_plus = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('e_minus', 'rb')
E_minus = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('Objective_func', 'rb')
Obj= pickle.load(pkl_file)
pkl_file.close()

# ii) post-processing for P and Z. If Z is 0, p should be zerio and integeration of E
for count in list(xrange(len(Num_candidate))):
          for J_Red in Jred_values:
                    if J_Red<=Num_candidate[count]:                    
                              for j in list(xrange(max(Num_candidate))):
                                        if (J_Red<=Num_candidate[count] and j<Num_candidate[count]):                                                                                          
                                                  if Z[Num_candidate[count],J_Red,j]<=0.5:
                                                            P[Num_candidate[count],J_Red,j]=0 
                                                  if P[Num_candidate[count],J_Red,j]==0:
                                                            Z[Num_candidate[count],J_Red,j]=0 
                                                  if Z[Num_candidate[count],J_Red,j]>=0.5:
                                                            Z[Num_candidate[count],J_Red,j]=1
                                        if j>=Num_candidate[count]:
                                                  Z[Num_candidate[count],J_Red,j]=0
                                                  P[Num_candidate[count],J_Red,j]=0
                                                  
E={}
for count in list(xrange(len(Num_candidate))):
          for J_Red in Jred_values:                                          
                    for i in list(xrange(Num_sites)):
                              if J_Red<=Num_candidate[count]:                                                  
                                        E[Num_candidate[count],J_Red,i,0]=E_plus[Num_candidate[count],J_Red,i,0]+E_minus[Num_candidate[count],J_Red,i,0]
                                        E[Num_candidate[count],J_Red,i,1]=E_plus[Num_candidate[count],J_Red,i,1]+E_minus[Num_candidate[count],J_Red,i,1]
                                        E[Num_candidate[count],J_Red,i,2]=E_plus[Num_candidate[count],J_Red,i,2]+E_minus[Num_candidate[count],J_Red,i,2]
                                        E[Num_candidate[count],J_Red,i,3]=E_plus[Num_candidate[count],J_Red,i,3]+E_minus[Num_candidate[count],J_Red,i,3]
                                         
#iii) save P, Z and E in txt                             
f = open("Pjs.txt", "w")
f.write('%s %s %s %s\n' % (2003, 2003,507,507))                                                                          
f.write('%s %s %s %s\n' % (300,500, 300, 500))                                                                          

for j in list(xrange(max(Num_candidate))): 
          f.write('%s %s %s %s\n' % (P[2003, 300,j],P[2003, 500,j], P[507, 300,j],P[507, 500,j]))                                                          
f.close()

f = open("Zjs.txt", "w")
f.write('%s %s %s %s\n' % (2003, 2003, 507, 507))                                                                          
f.write('%s %s %s %s\n' % (300, 500, 300, 500))                                                                          

for j in list(xrange(max(Num_candidate))): 
          f.write('%s %s %s %s \n' % (Z[2003, 300,j], Z[2003, 500,j], Z[507, 300,j],Z[507, 500,j]))                                                          
f.close()                                                      

f = open("Obj.txt", "w")
f.write('%s %s %s %s\n' % (2003, 2003, 507,507))                                                                          
f.write('%s %s %s %s\n' % (300,500, 300, 500))                                                                          
f.write('%s %s %s %s\n' % (Obj[2003,300], Obj[2003,500],Obj[507,300],Obj[507,500]))   

print  "**********************" 

for count in list(xrange(len(Num_candidate))):
          for J_Red in Jred_values:  
                    if J_Red<=Num_candidate[count]:                                                                      
                              f = open("E%s_%s.txt" % (Num_candidate[count], J_Red), "w")
                              for i in list(xrange(Num_sites)):                    
                                        f.write('%s %s %s %s\n' % (E[Num_candidate[count],J_Red,i,0],E[Num_candidate[count],J_Red,i,1],E[Num_candidate[count],J_Red,i,2],E[Num_candidate[count],J_Red,i,3]))
                              f.close()


