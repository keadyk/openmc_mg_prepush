import numpy as npy
import sys
import math
import cmath

#first things first: get relevant data (M, sigT, h_j & h_k, etc) from the file
execfile(sys.argv[1])

#fill weights and mus
if M > 0:
    execfile("weights_mus.inp")
else:
    print "Uh-oh, we didn't get a valid value for M from " + sys.argv[1] + "!"
    sys.exit(1)

#print out relevant data:
print "----------Relevant data:----------"
print "     Slab size: " + str(X)
print "  Total x-sect: " + str(sigT) 
print "   Scat. ratio: " + str(s_rat)
print "    Ang. order: " + str(M)
print "Fine grid size: " + str(h_k)
print "Crs. grid size: " + str(h_j)
print "  Fine per crs: " + str(fperc)
print "   # inner its: " + str(n_in)
print "-----------------------------------"

print "---------Mus and weights:----------"
for i in range(M/2):
    print "Mu: " + str(mus.item(i)) + " Wt: " + str(wts.item(i))
print "-----------------------------------"

#Calculated alphas for
alfas = [] 
for i in range(M):
    this_mu = mus.item(math.floor(float(i)/2))*math.pow(-1.0,i+1)
    #alfas.append(0)
    #if i==0:
    #   print "USING DIAMOND DIFFERENCE ALPHA!"
    alfas.append((1+math.exp(-sigT*h_k/this_mu))/(1-math.exp(-sigT*h_k/this_mu)) - 2.0*this_mu/(sigT*h_k))
    if i==0:
        print "USING STEP CHARACTERISTIC ALPHA!"
    #print "Mu: " + str(this_mu) + ", alpha: " + str(alfas[-1]) 

#initialize max eig, f-value and corresponding lambda:
max_lambda = 0.0
max_lambda_f = 0.0
max_f = 0
max_eig = 0.0
max_eig_real = 0.0
max_eig_imag = 0.0 

#SET A RANGE OF LAMBDA VALUES:
lambda_range = npy.arange((2*math.pi)/(sigT*X), (X/h_j - 1)*(2*math.pi)/(sigT*X)+0.001, (2*math.pi)/(sigT*X))

#lambda_range = npy.arange(0.01, 2*math.pi, 0.005)
#start with just a simple fine lambda guess:

for this_lambda in lambda_range:
    #lambda_range_f = npy.arange((2*math.pi)/(sigT*h_j), (fperc - 1)*(2*math.pi)/(sigT*h_j), (2*math.pi)/(sigT*h_j))
    if(fperc > 1):
        lambda_range_f = npy.arange((2*math.pi)/(sigT*h_j), (fperc-1)*(2*math.pi)/(sigT*h_j)+0.001, (2*math.pi)/(sigT*h_j))
    else:
        lambda_range_f = npy.arange((2*math.pi)/(sigT*h_j), (2*math.pi)/(sigT*h_j)+0.001, (2*math.pi)/(sigT*h_j))

    for this_lambda_f in lambda_range_f:
        #print "\r ****Running ({0:5.3f}, {1:5.3f}); current max = {2:7.5f} @ ({3:7.5f},{4:7.5f})****".format(this_lambda, this_lambda_f, max_eig, max_lambda, max_lambda_f),
        sys.stdout.flush()
        
        for g in range(int(X/h_j)):
            this_cell = (g+0.5)*h_j
            this_le   = g*h_j
            this_exp = complex(math.cos(sigT*this_lambda_f*this_cell), math.sin(sigT*this_lambda_f*this_cell))
            this_le_exp = complex(math.cos(sigT*this_lambda_f*this_le), math.sin(sigT*this_lambda_f*this_le))
            this_c_exp = complex(math.cos(sigT*this_lambda*this_cell), math.sin(sigT*this_lambda*this_cell))
            this_lec_exp = complex(math.cos(sigT*this_lambda*this_le), math.sin(sigT*this_lambda*this_le))
            #print "left-edge: " + '\t{0:7.5f}'.format(this_le_exp) + '\t{0:7.5f}'.format(this_lec_exp)
            #print "center: " + '\t{0:7.5f}'.format(this_exp) + '\t{0:7.5f}'.format(this_c_exp)
        #print "\nFINE CELLs...."
        for h in range(int(X/h_k)):
            this_cell = (h+0.5)*h_k
            this_le   = h*h_k
            this_exp = complex(math.cos(sigT*this_lambda_f*this_cell), math.sin(sigT*this_lambda_f*this_cell))
            this_le_exp = complex(math.cos(sigT*this_lambda_f*this_le), math.sin(sigT*this_lambda_f*this_le))
            #print "left-edge f: " + '\t{0:7.5f}'.format(this_le_exp)
            #print "center f: " + '\t{0:7.5f}'.format(this_exp)            
       
        #Allocate the eigenvalue vector for later
        eigs = 0.0
          
        #Calculate coefficients of M_1
        #allocate arrays:
        M_1 = npy.zeros((M, M), dtype=complex)
        M_2 = npy.zeros((M, M), dtype=complex)
        M_1_inv_sum = npy.zeros(M, dtype=complex)
        #Just save the ordered weights bc I'm lazy
        ord_wt = npy.zeros(M, dtype=complex)
        for i in range(M):
            this_mu = mus.item(math.floor(i/2.0))*math.pow(-1.0,(i+1))
            this_wt = wts.item(math.floor(i/2.0))  
            alfa = alfas[int(i)]
            #print "\nmu= " + str(this_mu) + " wt= " + str(this_wt) + " alfa= " + str(alfa) 
            M_1[i][i] +=  2.0*complex(math.cos(sigT*(this_lambda_f+this_lambda)*0.5*h_k), 0) + (2.0*this_mu/(sigT*h_k) + alfa)*complex(0, math.sin(sigT*(this_lambda_f+this_lambda)*0.5*h_k))
            
            M_2[i][i] += (2.0*complex(math.cos(sigT*(this_lambda_f+this_lambda)*0.5*h_k), 0) + alfa*complex(0, math.sin(sigT*(this_lambda_f+this_lambda)*0.5*h_k))) * this_wt
                
            #else:
                #M_1[i][i] += (this_mu/(sigT*h_k)+0.5*(1+alfa))*\
                #complex(math.cos(sigT*this_lambda*h_j), math.sin(sigT*this_lambda*h_j))*\
                #complex(math.cos(sigT*this_lambda_f*(h_j-h_k)), math.sin(sigT*this_lambda_f*(h_j-h_k)))            
                #M_1[i][i] += -this_mu/(sigT*h_k) + 0.5*(1-alfa)
                
                #M_2[i][i] += 0.5*(1+alfa)*complex(math.cos(sigT*this_lambda*h_j), math.sin(sigT*this_lambda*h_j)*complex(math.cos(sigT*this_lambda_f*(h_j-h_k)), math.sin(sigT*this_lambda_f*(h_j-h_k)))) * this_wt
                #M_2[i][i] += 0.5*(1-alfa)*this_wt
                #print " crs term: " + str(math.cos(sigT*this_lambda*h_j)) + " " + str(math.sin(sigT*this_lambda*h_j))
                #print " fine term: " + str(math.cos(sigT*this_lambda_f*(h_j-h_k))) + " " + str(math.sin(sigT*this_lambda_f*(h_j-h_k)))
        #Invert M_1 matrix   
        M_1_inv = npy.matrix(M_1).I

        #print "...Matrix M_1 inv..."
        #for i in range(M): 
        #   line_string = ""        
        #   for c in range(M):
        #        line_string += '{0:7.5f}\t'.format(M_1_inv.item(i,c))     
        #   print line_string 
        
        #Sum rows... (M_1_inv is diagonal)
        for l in range(M):
            M_1_inv_sum[l] = M_1_inv.item(l, l)

        #print "...Matrix M_1 ..."
        #for i in range(M): 
        #   line_string = ""        
        #   for c in range(M):
        #        line_string += '{0:7.5f}\t'.format(M_1.item(i,c))     
        #   print line_string 
           
        #print "...Matrix M_2 ..."
        #for i in range(M): 
        #   line_string = ""        
        #   for c in range(M):
        #        line_string += '{0:7.5f}\t'.format(M_2.item(i,c))     
        #   print line_string 
           
        M2_M1i = M_2 * M_1_inv
           
        #print "...Matrix M2_M1i..."
        #for i in range(M): 
        #   line_string = ""        
        #   for c in range(M):
        #        line_string += '{0:7.5f}\t'.format(M_1_inv.item(i,c))     
        #   print line_string 

        #Now, calculate H and Q (sans E_r):
        H = M2_M1i * 0.5*s_rat
        Q = M2_M1i * 0.5*(1-s_rat)

        #print "...Matrix H..."
        #for i in range(M): 
        #   line_string = ""        
        #   line_string += '{0:7.5f}\t'.format(H.item(i,i))     
        #   print line_string
        #print "...Matrix Q..."
        #for i in range(M): 
        #   line_string = ""        
        #   line_string += '{0:7.5f}\t'.format(Q.item(i,i))     
        #   print line_string
           
        #Sum up the H and Q entries:
        H_sum = 0
        Q_sum = 0
        for m in range(M):
            H_sum += H.item(m,m)
            Q_sum += Q.item(m,m)
        #print "\nH and Q sums: " + '{0:7.5f}\t'.format(H_sum) + " " + '{0:7.5f}\t'.format(Q_sum)  

        H_k = (H_sum**n_in)
        H_k1 = (H_sum**(n_in-1))
        #print "raised to power... " + str(H_k) + " " + str(H_k1)  

        #Use summation formula to calculate
        H_sum_k = (H_k - 1.0)/(H_sum - 1.0)
        H_sum_k1 = (H_k1 - 1.0)/(H_sum - 1.0)     
        
        #Compute and sum J-coeff matrix
        J_matrix = M_1_inv_sum * (0.5*s_rat * (H_k1 + Q_sum*H_sum_k1) + 0.5*(1-s_rat))
        
        Final_J_coeff = 0.0
        for m in range(M):
            this_mu = mus.item(math.floor(m/2.0))*math.pow(-1.0,(m+1))
            this_wt = wts.item(math.floor(m/2.0))  
            #print "\nmu= " + str(this_mu) + " wt= " + str(this_wt)
            Final_J_coeff += this_mu * this_wt * J_matrix.item(m)
            #print '{0:7.5f}\t'.format(J_matrix.item(m))  
        #Multiply H sum terms and Q_sum term to obtain scalar flux error ratio
        Final_coeff = H_k + Q_sum*H_sum_k
        
        #print '{0:7.5f}\t'.format(Final_J_coeff)      
        #print "...Final scalar flux error matrix..."
        #for i in range(fperc): 
        #   line_string = ""        
        #   for c in range(fperc):
        #        line_string += '{0:7.5f}\t'.format(Final_coeffs.item(i,c))    
        #   print line_string  
        
        #print "\n...Final scalar flux error ratios for coarse lambda " + str(this_lambda) + ", fine lambda " + str(this_lambda_f) + ", r=" + str(r) + "..."
        #print '{0:7.5f}\t'.format(Final_coeff)
        #print  '{0:7.5f}\t'.format(complex(math.cos(sigT*this_lambda_f*0.5*(h_j - (2*r+1)*h_k)), math.sin(sigT*this_lambda_f*0.5*(h_j - (2*r+1)*h_k))))          
        #print "...Final current error matrix..."
        #for i in range(fperc): 
        #   line_string = ""        
        #   for c in range(fperc):
        #        line_string += '{0:7.5f}\t'.format(Final_J_coeff)    
        #   print line_string 
        #line_string = ""        
        #for c in range(fperc):
        #    line_string += '{0:7.5f}\t'.format(Final_J_coeffs_pp1.item(c))    
        #print line_string      
        
        #NOW! We have the error matrices that relate the current and scalar flux error to the fission source error.
        
        #Calculate weird sum term:
        for a in range(int(X/h_j)):
            this_cell = (a+0.5)*float(h_j)
            weird_sum = complex(0.0, 0.0)
            for f in range(fperc):
                this_offset = 0.5*h_j - 0.5*(2*f+1)*h_k
                #print "Offset amount: " + str(f+1) + " " + str(this_offset)
                weird_sum += complex(math.cos(-sigT*(this_lambda_f + this_lambda)*this_offset), math.sin(-sigT*(this_lambda_f + this_lambda)*this_offset))
                #Try leaving the weird sum alone instead   
            weird_sum *= complex(math.cos(sigT*this_lambda_f*this_cell), math.sin(sigT*this_lambda_f*this_cell))
            #print "weird sum" + str(weird_sum)
        
        #calculate the MAIN big constant, const:
        #I TOOK A FACTOR OF FPERC OUT OF THIS! Was this a bad thing?!
        const = 3.0*sigT*h_j*complex(0.0, math.sin(sigT*this_lambda*0.5*h_j))/(complex(math.cos(sigT*this_lambda*h_j), 0.0) - 1)
        #print "constant " + str(const)
        C = const*(fperc/weird_sum)
    
        #Slightly different from notes due to zero-indexing
        #print "...Final scalar flux error matrix..."   
        #print '{0:7.5f}'.format(Final_coeff)     
        #print "...Final current error matrix..."
        #print '{0:7.5f}'.format(Final_J_coeff)
        #print "...Constant...\n" + '{0:7.5f}'.format(C) 
        eigs = C*Final_J_coeff + Final_coeff 
        #print "eigs " + str(eigs[r])
        
        #print '\n'
        if npy.abs(eigs) > max_eig:
            max_eig = npy.abs(eigs)
            max_eig_real = eigs.real
            max_eig_imag = eigs.imag
            max_lambda = this_lambda 
            max_lambda_f = this_lambda_f
        print '{0:7.5f}, {1:7.5f}: {2:7.5f}'.format(this_lambda, this_lambda_f, eigs)


print "\n"        
print "Spectral radius mag = " + str(max_eig)
print "Real part           = " + str(max_eig_real)
print "Imag part           = " + str(max_eig_imag)
print "   ...at lambda_c   = " + str(max_lambda)
print "   ...and lambda_f  = " + str(max_lambda_f)
print "    ...and cell     = " + str(max_f)      
            
        
    
    
    
    
