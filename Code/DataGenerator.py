# Input .mat file from matlab 
# Output binary file compatible for python  

##############################################################################################
# User Input: Please define the user input parameters 

# The body of code works as follows:   
       # i)   This part of the code reads .mat file from Matlab
       # ii)  Creates compatible data in the form of a dictionary for Python 
       # iii) Data will be save into a binary file 
# We just need to run this code once since we save the data into a binary file which can be used over and over again 
##############################################################################################
## Libraries 
import pickle                                                                                                   # This is a library to convert data to binary file
import scipy.io                                                                                                 # This is a scientific computing library 

## User input                                   
R=[250,500,1000,2500]                                                                                           # Return periods
filename=['Pnew9','Pnew19']                                                                                     # Name(s) for input files 
Num_candidate=[2003,507]                                                                                        # Number of candidates coresponding to finename
Num_sites=35                                                                                                    # Please enter the number of sites
path='Pnew'                                                                                                     # Please enter the name of .mat file (.mat file and .py files in the same dic) 
## End of User Input

## The Body of code
# i)  This part of the code reads .mat file from Matlab
Return_periods=list(xrange(len(R)))                                                                                
count=0
for Earthquake_NumScneario in Num_candidate:                                                                   
    count=count+1            
    ## Indices                                                               
    erthquake_candidates= list(xrange(Earthquake_NumScneario))                                                   # Earthquake candidates,  j index 
    sites=list(xrange(Num_sites))                                                                                # Sites                ,  i index
    Coeff=scipy.io.loadmat(path)                                                                                 # This is P(y_ij>=y_ir)
    A= Coeff

# ii) Creates compatible data in the form of a dictionary for Python  
    Coeff={}                                                                                                     # Save data into a dictionary 
    for r in Return_periods:
        for j in erthquake_candidates:
            for i in sites:
                    Coeff[j,i,r]= A[filename[count-1]][0,j][0][i][r]
   
# iii) Data will be saved into a binary file  
    Namefile= filename[count-1]                                                                                  # Save the dictionary into a binary file 
    print Namefile
    output = open(Namefile, 'wb')                                                                                # This is done since we do not want to call this code over and over again
    pickle.dump(Coeff, output)
    output.close()
## The end of the body of code
##############################################################################################

    