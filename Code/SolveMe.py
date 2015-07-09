# Input binary files generated using DataGenerator code 
# Output Bianry files in which we have Pj's, Zj's, e+ir, e-ir, and objective functions 

##############################################################################################
# User Input: Please define the user input parameters 

# The body of code works as follows:   
       # i) Create the model in PULLP 
       # ii) Solve the model using Gurobi or COIN 
       # iii) Save data into a binary file. 

## LIBRARIES  
from pulp import *                                                                                              # This is a library to create the model
import gurobipy                                                                                                 # This is a library to solve the model (i.e., solver)
import pickle                                                                                                   # This is a library to convert data into binary file

## User input                                   
R=[250,500,1000,2500]                                                                                           # Return periods
Filename=['Pnew9','Pnew19']                                                                                     # Here is all names for input files 
Num_candidate=[2003,507]                                                                                        # Number of candidates coresponding to finename [2003, 507, 4002, 5001, 1003, 132, 3005]     
Jred_values= [300,500]                                                                                          # Jred
Num_sites=35                                                                                                    # Please enter the number of sites
TimeLimit=3600*3                                                                                                # This is the time you allow the optimization model runs (in seconds)
## End of User Input

# i) Create the model 
P_solution={}                                                                                                   # Solutions
Z_solution={}
e_minus_solution={}
e_plus_solution={}
Objective_solution={}  
            ################  Parameters and indices ####################
count=0                                                         
for ii in Num_candidate:                                                                                        # Num candidates 
    count=count+1    
    for jj in Jred_values:                                                                                      # Jred
        if jj<=ii:
            print "This is the solution for %s %s" % (ii,jj)
            Num_candidate=ii                                                                                   
            Return_periods=list(xrange(len(R)))                                                                      
            J_red=jj                                          

            ## Indices
            candidates = list(xrange(Num_candidate))                                                            # Earthquake candidates,  j index 
            sites=list(xrange(Num_sites))                                                                       # Sites                ,  i index
                
            W={}                                                                                                # Weights of objective function W_i_r
            for i in sites:
                for r in Return_periods:
                    W[i,r]=R[r]
            
            Namefile=Filename[count-1]
            
            pickle_file = open(Namefile, 'rb')                                                                  # Reading P(y(i,j)>=Y(i,r)), we call it Coeff (i,j,r)
            Coeff = pickle.load(pickle_file)
            pickle_file.close()
            
            ################ Mathematical Model #################### 
            ## Model
            prob = LpProblem("MIP", LpMinimize)                                                                # Create model with PULP
            
            ## Variables
            e_plus={}                                                                                          # e+(i,r) variable
            e_minus={}                                                                                         # e-(i,r) variable
            P={}                                                                                               # P(j) variable
            Z={}                                                                                               # Z(j) variable
            
            e_plus_name={}                                                                                     # We consider a name associated with each variable for the sake of clarity 
            e_minus_name={}
            for i in sites:
                for r in Return_periods:
                    e_plus_name[i,r]    =   'ep %s %s' %(i,r)
                    e_minus_name[i,r]   =   'en %s %s' %(i,r)
                    e_plus[i,r]         =   LpVariable(e_plus_name[i,r], 0, None)
                    e_minus[i,r]        =   LpVariable(e_minus_name[i,r], 0, None)
            
            P_name={}                                                                                           # We consider a name associated with each variable for the sake of clarity 
            Z_name={}
            for j in candidates:
                P_name[j] =  'P %s' % j
                Z_name[j] =  'z %s' % j
                P[j]                 =   LpVariable(P_name[j], 0, 1)
                Z[j]                 =   LpVariable(Z_name[j], 0, 1, LpBinary)
                
            ## Objective function
            Set_i_j=list(itertools.product(sites, Return_periods))                                              # here we create a list of all possible combinations of i, r          
            prob+=lpSum(W[(i,r)]*e_plus[(i,r)]+W[(i,r)]*e_minus[(i,r)] for (i,r) in Set_i_j)
            
            ## Constraints 
            for i in sites:                                                                                     # Constraint 3 in the paper 
                for r in Return_periods:
                    prob += lpSum(P[j]*Coeff[j,i,r] for j in candidates) +e_minus[i,r]-e_plus[i,r] == 1.0/R[r]
        
            for j in candidates:                                                                                # Constraint 4 in the paper
                prob += P[j]-Z[j] <=0
            
            prob += lpSum(Z[j] for j in candidates) <= J_red                                                    # Constraint 5 in the paper
            
# ii) Solve the model using Gurobi or COIN
            prob.solve(PULP_CBC_CMD(maxSeconds=TimeLimit))                                                      # Can replace PULP_CBC_CMD(maxSeconds=TimeLimit) with GUROBI(). 
            
            print "Status:", LpStatus[prob.status]
            
            ###################### Save solution in dictionary ######################## 
            
            for j in candidates:
                P_solution[ii,jj,j]=P[j].varValue 
                Z_solution[ii,jj,j]=Z[j].varValue
                
            for i in sites:
                for r in Return_periods:
                    e_plus_solution[ii,jj,i,r]=e_plus[i,r].varValue 
                    e_minus_solution[ii,jj,i,r]= e_minus[i,r].varValue
            
            Objective_solution[ii,jj]= value(prob.objective)
        
###################### Report solution ######################## 

#iii) Save solution into binary files 
Namefile_output= 'P' 
output = open(Namefile_output, 'wb')
pickle.dump(P_solution, output)
output.close()

Namefile_output= 'Z' 
output = open(Namefile_output, 'wb')
pickle.dump(Z_solution, output)
output.close() 

Namefile_output= 'e_plus' 
output = open(Namefile_output, 'wb')
pickle.dump(e_plus_solution, output)
output.close()    

Namefile_output= 'e_minus' 
output = open(Namefile_output, 'wb')
pickle.dump(e_minus_solution, output)
output.close()    

Namefile_output= 'Objective_func' 
output = open(Namefile_output, 'wb')
pickle.dump(Objective_solution, output)
output.close()      

