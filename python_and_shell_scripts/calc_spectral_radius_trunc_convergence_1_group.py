import numpy as npy
import sys
import math

#first things first: get relevant data (M, sigT, h_j & h_k, etc) from the file
#execfile(sys.argv[1])
exec(open(sys.argv[1]).read())

#fill weights and mus
if M > 0:
    #execfile("weights_mus.inp")
    exec(open("weights_mus.inp").read())
else:
    print( "Uh-oh, we didn't get a valid value for M from " + sys.argv[1] + "!")
    sys.exit(1)

#print(out relevant data:
print( "----------Relevant data:----------")
print("     Slab size: " + str(X))
print("  Total x-sect: " + str(sigT)) 
print("   Scat. ratio: " + str(s_rat))
print("    Ang. order: " + str(M))
print("Fine grid size: " + str(h_k))
print("Crs. grid size: " + str(h_j))
print("  Fine per crs: " + str(fperc))
print("   # inner its: " + str(n_in))
print("-----------------------------------")

print("---------Mus and weights:----------")
for i in range(M/2):
    print("Mu: " + str(mus.item(i)) + " Wt: " + str(wts.item(i)))
print("-----------------------------------")

#Calculated alphas for
alfas = [] 
for i in range(M):
    this_mu = mus.item(math.floor(float(i)/2))*math.pow(-1.0,i+1)
    #alfas.append(0)
    #if i==0:
    #   print("USING DIAMOND DIFFERENCE ALPHA!"
    alfas.append((1+math.exp(-sigT*h_k/this_mu))/(1-math.exp(-sigT*h_k/this_mu)) - 2.0*this_mu/(sigT*h_k))
    if i==0:
        print("USING STEP CHARACTERISTIC ALPHA!")
    #print("Mu: " + str(this_mu) + ", alpha: " + str(alfas[-1]) 

#initialize max eig, f-value and corresponding lambda:
max_lambda = 0.0
max_f = 0
max_eig = 0.0

#SET A RANGE OF LAMBDA VALUES:
lambda_range = npy.arange((2*math.pi)/(sigT*X), (X/h_j - 1)*(2*math.pi)/(sigT*X)+0.001, (2*math.pi)/(sigT*X))
#lambda_range = npy.arange(0, (X/h_j - 1)*(2*math.pi)/(sigT*X)+0.001, (2*math.pi)/(sigT*X))

#lambda_range = npy.arange(0.01, 2*math.pi, 0.005)

for this_lambda in lambda_range:
    #print("\r ****Running lambda =  {0:5.3f}; current max = {1:7.5f} @ {2:7.5f}****".format(this_lambda, max_eig, max_lambda),
    sys.stdout.flush()
    
    #allocate arrays:
    M_1 = npy.zeros((M*fperc, M*fperc), dtype=complex)
    M_2 = npy.zeros((M*fperc, M*fperc), dtype=complex)
    M2_M1i = npy.zeros((M*fperc, M*fperc), dtype=complex)
    M2_M1i_sum = npy.zeros((fperc, fperc), dtype=complex)
    mu_wt_M1i_sum = npy.zeros((fperc, fperc), dtype=complex)
    M3 = npy.zeros((M*fperc, fperc), dtype=complex)
    mu_wt_M3 = npy.zeros((M*fperc, fperc), dtype=complex)
    Final_J_coeffs_pp1 = npy.zeros((fperc), dtype=complex)
    ID = npy.identity(fperc, dtype=complex)
    
    #Calculate coefficients of M_1 and M_2
    for i in range(M*fperc):
        this_m = math.ceil(float(i+1)/fperc)
        this_f = (i+1) - ((this_m-1)*fperc)
        this_mu = mus.item(math.floor(float(this_m-1)/2))*math.pow(-1.0,this_m)
        this_wt = wts.item(math.floor((this_m-1)/2))  
        alfa = alfas[int(this_m-1)]
        #print("\nmu= " + str(this_mu) + " wt= " + str(this_wt) + " alfa= " + str(alfa) 
        #Using the periodicity condition from Kelley and Larsen:
        if(this_f == fperc):
            M_1[i][i] += -this_mu/(sigT*h_k) + 0.5*(1 - alfa)
            M_1[i][(this_m-1)*fperc] += (this_mu/(sigT*h_k) + 0.5*(1 + alfa))*complex(math.cos(sigT*this_lambda*h_j), math.sin(sigT*this_lambda*h_j))
            
            M_2[i][i] += 0.5*(1 - alfa)*this_wt
            M_2[i][(this_m-1)*fperc] += 0.5*(1 + alfa)*this_wt*complex(math.cos(sigT*this_lambda*h_j), math.sin(sigT*this_lambda*h_j)) 
            
            #print(" crs term: " + str(math.cos(sigT*this_lambda*h_j)) + " " + str(math.sin(sigT*this_lambda*h_j))
            
        else:
            M_1[i][i] += -this_mu/(sigT*h_k) + 0.5*(1 - alfa)
            M_1[i][i+1] += this_mu/(sigT*h_k) + 0.5*(1 + alfa) 
            
            M_2[i][i] += 0.5*(1 - alfa)*this_wt
            M_2[i][i+1] += 0.5*(1 + alfa)*this_wt  
           
    #Invert M_1 matrix   
    M_1_inv = npy.matrix(M_1).I
       
    #print("...Matrix M_1 inv..."
    #for i in range(M*fperc): 
    #   line_string = ""        
    #   for c in range(M*fperc):
    #        line_string += '{0:7.5f}\t'.format(M_1_inv.item(i,c))     
    #   print(line_string 
       
    #print("...Matrix M_1...")
    #for i in range(M*fperc): 
    #   line_string = ""        
    #   for c in range(M*fperc):
    #        line_string += '{0:7.5f} + {1:7.5f}i\t'.format(M_1.item(i,c).real, M_1.item(i,c).imag)     
    #   print(line_string) 
    #print("...Matrix M_2...")
    #for i in range(M*fperc): 
    #   line_string = ""        
    #   for c in range(M*fperc):
    #        line_string += '{0:7.5f} + {1:7.5f}i\t'.format(M_2.item(i,c).real, M_2.item(i,c).imag)     
    #   print(line_string  )  
    M2_M1i = npy.matrix(M_2)*M_1_inv

    #SUMS OVER M::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   
    #Sum over columns in a row
    for i in range(M*fperc):
        for j in range(fperc):
            for m in range(M):
                #print("*mu and wt: " + str(this_mu) + " " + str(this_wt)
                M3[i][j] += M2_M1i.item(i, j+m*fperc)
                mu_wt_M3[i][j] += M_1_inv.item(i, j+m*fperc)

    #Sum over rows in a column
    for j in range(fperc):
        for r in range(fperc):
            for m in range(M):
                this_mu = mus.item(math.floor(float(m)/2))*math.pow(-1.0,m+1)
                this_wt = wts.item(math.floor(float(m)/2)) 
            #Grab every "M"th value:
                M2_M1i_sum[j][r] += M3.item((j+m*fperc, r))
                mu_wt_M1i_sum[j][r] += this_mu * this_wt * mu_wt_M3.item((j+m*fperc, r))
    #SUMS OVER M:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

    #print("\n...Matrix M2*M1 inverse, summed over m...")
    #for i in range(fperc): 
    #   line_string = ""        
    #   for c in range(fperc):
    #        line_string += '{0:7.5f} + {1:7.5f}i\t'.format(M2_M1i_sum.item(i,c).real, M2_M1i_sum.item(i,c).imag)     
    #   print(line_string)
       
    #Now, calculate H and Q (sans E_r):
    H = npy.matrix(M2_M1i_sum) * 0.5

    #print("\n...Matrix H...")
    #for i in range(fperc): 
    #   line_string = ""        
    #   for c in range(fperc):
    #        line_string += '{0:7.5f} + {1:7.5f}i\t'.format(H.item(i,c).real, H.item(i,c).imag)     
    #   print(line_string)

    H_k1 = npy.identity(fperc, dtype=complex)
    
    for k in range(n_in):
        #reset H_k
        #H_k = npy.zeros((fperc, fperc), dtype=complex)       
        #calculate h^k
        H_k = npy.matrix(H) * npy.matrix(H_k1)
        
        #Store "old" H_k, unless it's the last time through
        if(k < (n_in-1)):
            H_k1 = npy.matrix(H_k)                 
               
    #Finally, calculate coefficients
    Final_coeffs = npy.matrix(H_k)
    Final_J_coeffs = npy.matrix(mu_wt_M1i_sum)*0.5*npy.matrix(H_k1)
     
    #print("...Final scalar flux error matrix..."
    #for i in range(fperc): 
    #   line_string = ""        
    #   for c in range(fperc):
    #        line_string += '{0:7.5f}\t'.format(Final_coeffs.item(i,c))    
    #   print(line_string          
    #print("...Final scalar flux error matrix..."
    #for i in range(fperc): 
    #   line_string = ""        
    #   for c in range(fperc):
    #        line_string += '{0:7.5f}\t'.format(Final_coeffs.item(i,c))    
    #   print(line_string  
    
    #print("\n...Final scalar flux error ratios for " + str(this_lambda) + "..."
    #line_string = "" 
    #for i in range(fperc): 
    #   sum = 0       
    #   for c in range(fperc):
    #        sum += (Final_coeffs.item(i,c))    
    #   line_string += '{0:7.5f}\t'.format(sum)
    #print(line_string
      
    #print("...Final current error matrix..."
    #for i in range(fperc): 
    #   line_string = ""        
    #   for c in range(fperc):
    #        line_string += '{0:7.5f}\t'.format(Final_J_coeffs.item(i,c))    
    #   print(line_string 
    #line_string = ""        
    #for c in range(fperc):
    #    line_string += '{0:7.5f}\t'.format(Final_J_coeffs_pp1.item(c))    
    #print(line_string      
    
    #NOW! We have the error matrices that relate the current and scalar flux error to the fission source error.
    
    #calculate the big constant, G:
    const = 1.5*sigT*h_j
    G = complex(const, const*math.sin(sigT*this_lambda*h_j)/(math.cos(sigT*this_lambda*h_j) - 1))
    #print("...Constant...\n" + '{0:7.5f} + {1:7.5f}i'.format(G.real, G.imag) 
    #Allocate the K-matrix
    K = npy.zeros((fperc, fperc), dtype=complex)
    
    for a in range(fperc):
        for b in range(fperc):
            K[a][b] += G*Final_J_coeffs.item(0, b) + Final_coeffs.item(a, b) 
    
    #print("...K matrix..."
    #for i in range(fperc): 
    #   line_string = ""        
    #   for c in range(fperc):
    #        line_string += '{0:7.5f}\t'.format(K.item(i,c))    
    #   print(line_string     
    
    #Now, FINALLY, calculate the eigenvalues!
    eigs = npy.linalg.eigvals(npy.matrix(K))
    
    #print('\n'
    for h in range(fperc):
        if math.fabs(eigs.item(h).real) > max_eig:
            max_eig = math.fabs(eigs.item(h).real)
            max_lambda = this_lambda 
            max_f = h+1
        #print('{0:7.5f}: {1:7.5f}'.format(this_lambda, eigs.item(h).real))
    #exit(1)

print("\n")        
print("Spectral radius = " + str(max_eig))
print("   ...at lambda = " + str(max_lambda))
print("    ...and cell = " + str(max_f))      
            
        
    
    
    
    
