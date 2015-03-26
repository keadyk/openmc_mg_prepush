import numpy as npy
import sys
import math

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
print " Total x-sects: " 
print "             1: " + str(sigT[0])
print "             2: " + str(sigT[1])
print " Scat. x-sects: " 
print "        1 -> 1: " + str(sigS[0][0])
print "        1 -> 2: " + str(sigS[0][1])
print "        2 -> 1: " + str(sigS[1][0])
print "        2 -> 2: " + str(sigS[1][1])
print "    Ang. order: " + str(M)
print "Fine grid size: " + str(h_k)
print "Crs. grid size: " + str(h_j)
print "  Fine per crs: " + str(fperc)
print "    # inner its: " + str(n_in)
print "-----------------------------------"

print "---------Mus and weights:----------"
for i in range(M/2):
    print "Mu: " + str(mus.item(i)) + " Wt: " + str(wts.item(i))
print "-----------------------------------"

if(sigS[1][0] != 0.0):
    print "We assume no upscattering. SigS 2 -> 1 must = 0.0!"
    exit(1)

#Calculated alphas for
alfas = [[] for i in range(2)]
for g in range(2): 
    for i in range(M):
        this_mu = mus.item(math.floor(float(i)/2))*math.pow(-1.0,i+1)
        #alfas.append(0)
        #if i==0:
        #   print "USING DIAMOND DIFFERENCE ALPHA!"
        alfas[g].append((1+math.exp(-sigT[g]*h_k/this_mu))/(1-math.exp(-sigT[g]*h_k/this_mu)) - 2.0*this_mu/(sigT[g]*h_k))
        if i==0 and g==0:
            print "USING STEP CHARACTERISTIC ALPHA!"
        print "Mu: " + str(this_mu) + ", alpha[" + str(g+1) + "]: " + str(alfas[g][-1]) 

#initialize max eig, f-value and corresponding lambda:
max_lambda = 0.0
max_f = 0
max_eig = 0.0
max_eig_imag = 0.0

#SET A RANGE OF LAMBDA VALUES:
lambda_range = npy.arange((2*math.pi)/X, (X/h_j - 1)*(2*math.pi)/X+0.001, (2*math.pi)/X)
#lambda_range = npy.arange(0, (X/h_j - 1)*(2*math.pi)/(sigT*X)+0.001, (2*math.pi)/(sigT*X))

#lambda_range = npy.arange(0.01, 2*math.pi, 0.005)

for this_lambda in lambda_range:
    #print "\r ****Running lambda =  {0:5.3f}; current max = {1:7.5f} @ {2:7.5f}****".format(this_lambda, max_eig, max_lambda),
    sys.stdout.flush()
    
    #allocate arrays:
    M_1 = npy.zeros((2*M*fperc, 2*M*fperc), dtype=complex)
    M_2 = npy.zeros((2*M*fperc, 2*M*fperc), dtype=complex)
    M2_M1i = npy.zeros((2*M*fperc, 2*M*fperc), dtype=complex)
    M2_M1i_sum = npy.zeros((2*fperc, 2*fperc), dtype=complex)
    mu_wt_M1i_sum = npy.zeros((2*fperc, 2*fperc), dtype=complex)
    M3 = npy.zeros((2*M*fperc, 2*fperc), dtype=complex)
    mu_wt_M3 = npy.zeros((2*M*fperc, 2*fperc), dtype=complex)
    S = npy.zeros((2*fperc, 2*fperc), dtype=complex)
    ID = npy.identity(fperc, dtype=complex)
    
    #Calculate coefficients of M_1 and M_2
    for i in range(M*fperc):
        this_m = math.ceil(float(i+1)/fperc)
        this_f = (i+1) - ((this_m-1)*fperc)
        this_mu = mus.item(math.floor(float(this_m-1)/2))*math.pow(-1.0,this_m)
        this_wt = wts.item(math.floor((this_m-1)/2))  
        alfa0 = alfas[0][int(this_m-1)]
        alfa1 = alfas[1][int(this_m-1)]
        #print "\nmu= " + str(this_mu) + " wt= " + str(this_wt) + " alfa= " + str(alfa) 
        #Using the periodicity condition from Kelley and Larsen:
        if(this_f == fperc):
            M_1[i][i] += -this_mu/(sigT[0]*h_k) + 0.5*(1 - alfa0)
            M_1[i][(this_m-1)*fperc] += (this_mu/(sigT[0]*h_k) + 0.5*(1 + alfa0))*complex(math.cos(this_lambda*h_j), math.sin(this_lambda*h_j))            
            M_1[M*fperc+i][M*fperc+i] += -this_mu/(sigT[1]*h_k) + 0.5*(1 - alfa1)
            M_1[M*fperc+i][M*fperc+(this_m-1)*fperc] += (this_mu/(sigT[1]*h_k) + 0.5*(1 + alfa1))*complex(math.cos(this_lambda*h_j), math.sin(this_lambda*h_j))

            M_2[i][i] += 0.5*(1 - alfa0)*this_wt
            M_2[i][(this_m-1)*fperc] += 0.5*(1 + alfa0)*this_wt*complex(math.cos(this_lambda*h_j), math.sin(this_lambda*h_j))   
            
            M_2[M*fperc+i][M*fperc+i] += 0.5*(1 - alfa1)*this_wt
            M_2[M*fperc+i][M*fperc+(this_m-1)*fperc] += 0.5*(1 + alfa1)*this_wt*complex(math.cos(this_lambda*h_j), math.sin(this_lambda*h_j))           
            #print " crs term: " + str(math.cos(sigT*this_lambda*h_j)) + " " + str(math.sin(sigT*this_lambda*h_j))
            
        else:
            M_1[i][i] += -this_mu/(sigT[0]*h_k) + 0.5*(1 - alfa0)
            M_1[i][i+1] += this_mu/(sigT[0]*h_k) + 0.5*(1 + alfa0) 

            M_1[M*fperc+i][M*fperc+i] += -this_mu/(sigT[1]*h_k) + 0.5*(1 - alfa1)
            M_1[M*fperc+i][M*fperc+i+1] += this_mu/(sigT[1]*h_k) + 0.5*(1 + alfa1) 
            
            M_2[i][i] += 0.5*(1 - alfa0)*this_wt
            M_2[i][i+1] += 0.5*(1 + alfa0)*this_wt  

            M_2[M*fperc+i][M*fperc+i] += 0.5*(1 - alfa1)*this_wt
            M_2[M*fperc+i][M*fperc+i+1] += 0.5*(1 + alfa1)*this_wt  
           
    #Invert M_1 matrix   
    M_1_inv = npy.matrix(M_1).I
       
    #print "...Matrix M_1 inv..."
    #for i in range(M*fperc): 
    #   line_string = ""        
    #   for c in range(M*fperc):
    #        line_string += '{0:7.5f}\t'.format(M_1_inv.item(i,c))     
    #   print line_string 
       
    print "...Matrix M_1..."
    for i in range(2*M*fperc): 
       line_string = ""        
       for c in range(2*M*fperc):
            line_string += '{0:7.5f} + {1:7.5f}i\t'.format(M_1.item(i,c).real, M_1.item(i,c).imag)     
       print line_string 
    print "...Matrix M_2..."
    for i in range(2*M*fperc): 
       line_string = ""        
       for c in range(2*M*fperc):
            line_string += '{0:7.5f} + {1:7.5f}i\t'.format(M_2.item(i,c).real, M_2.item(i,c).imag)    
       print line_string    
    M2_M1i = npy.matrix(M_2)*M_1_inv

    #SUMS OVER M::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   
    #Sum over columns in a row
    for i in range(M*fperc):
        for j in range(fperc):
            for m in range(M):
                #print "*mu and wt: " + str(this_mu) + " " + str(this_wt)
                M3[i][j] += M2_M1i.item(i, j+m*fperc)
                mu_wt_M3[i][j] += M_1_inv.item(i, j+m*fperc)

                M3[M*fperc+i][fperc+j] += M2_M1i.item(M*fperc+i, M*fperc+j+m*fperc)
                mu_wt_M3[M*fperc+i][fperc+j] += M_1_inv.item(M*fperc+i, M*fperc+j+m*fperc)

    #Sum over rows in a column
    for j in range(fperc):
        for r in range(fperc):
            for m in range(M):
                this_mu = mus.item(math.floor(float(m)/2))*math.pow(-1.0,m+1)
                this_wt = wts.item(math.floor(float(m)/2)) 
            #Grab every "M"th value:
                M2_M1i_sum[j][r] += M3.item((j+m*fperc, r))
                mu_wt_M1i_sum[j][r] += this_mu * this_wt *  mu_wt_M3.item((j+m*fperc, r))

                M2_M1i_sum[fperc+j][fperc+r] += M3.item((M*fperc+j+m*fperc, fperc+r))
                mu_wt_M1i_sum[fperc+j][fperc+r] += this_mu * this_wt *  mu_wt_M3.item((M*fperc+j+m*fperc, fperc+r))
    #SUMS OVER M:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

    print "\n...Matrix M2*M1 inverse, summed over m..."
    for i in range(2*fperc): 
       line_string = ""        
       for c in range(2*fperc):
            line_string += '{0:7.5f} + {1:7.5f}i\t'.format(M2_M1i_sum.item(i,c).real, M2_M1i_sum.item(i,c).imag)     
       print line_string
    #print "\n...Matrix mu_wt_M1 inverse, summed over m..."
    for i in range(2*fperc): 
       line_string = ""        
       for c in range(2*fperc):
            line_string += '{0:7.5f} + {1:7.5f}i\t'.format(mu_wt_M1i_sum.item(i,c).real, mu_wt_M1i_sum.item(i,c).imag)     
    #   print line_string
    #Calculate S matrix
    for r in range(fperc):
        #Group 1 eqns for this r:
        #1->1 coupling
        S[r][r] += sigS[0][0]/sigT[0]
        #1->2 coupling
        S[r][fperc+r] += ((sigT[0]-sigS[0][0])*(sigT[1]-sigS[1][1]))/(sigS[0][1]*sigT[0])

        #Group 2 eqns for this r:
        #2->1 coupling
        S[fperc+r][r] += sigS[0][1]/sigT[1]
        #2->2 coupling
        S[fperc+r][fperc+r] += sigS[1][1]/sigT[1]

    print "\n...Matrix S..."
    for i in range(2*fperc): 
       line_string = ""        
       for c in range(2*fperc):
            line_string += '{0:7.5f} + {1:7.5f}i\t'.format(S.item(i,c).real, S.item(i,c).imag)     
       print line_string

    #Now, calculate H and Q (sans E_r):
    #DOUBLE CHECK THIS!!!
    H = npy.matrix(M2_M1i_sum)*0.5*npy.matrix(S)

    #print "\n...Matrix H..."
    for i in range(2*fperc): 
       line_string = ""        
       for c in range(2*fperc):
            line_string += '{0:7.5f} + {1:7.5f}i\t'.format(H.item(i,c).real, H.item(i,c).imag)     
    #   print line_string

    H_k1 = npy.identity(2*fperc, dtype=complex)
    
    for k in range(n_in):       
        #calculate h^k
        H_k = npy.matrix(H) * npy.matrix(H_k1)
        
        #Store "old" H_k, unless it's the last time through
        if(k < (n_in-1)):
            H_k1 = npy.matrix(H_k)                 
               
    #Finally, calculate coefficients
    Final_coeffs = npy.matrix(H_k)
    Final_J_coeffs = 0.5*npy.matrix(mu_wt_M1i_sum)*npy.matrix(S)*H_k1

    #print "...H_k1..."
    for i in range(2*fperc): 
       line_string = ""        
       for c in range(2*fperc):
            line_string += '{0:7.5f} + {1:7.5f}i\t'.format(H_k1.item(i,c).real, H_k1.item(i,c).imag)    
    #   print line_string 
     
    #print "...Final scalar flux error matrix..."
    for i in range(2*fperc): 
       line_string = ""        
       for c in range(2*fperc):
            line_string += '{0:7.5f} + {1:7.5f}i\t'.format(Final_coeffs.item(i,c).real, Final_coeffs.item(i,c).imag)    
    #   print line_string          
    #print "...Final current error matrix..."
    for i in range(2*fperc): 
       line_string = ""        
       for c in range(2*fperc):
            line_string += '{0:7.5f} + {1:7.5f}i\t'.format(Final_J_coeffs.item(i,c).real, Final_J_coeffs.item(i,c).imag)    
    #   print line_string          
    
    #NOW! We have the error matrices that relate the current and scalar flux error to the fission source error.
    
    #calculate the big constant, G:
    const0 = 1.5*sigT[0]*h_j
    const1 = 1.5*sigT[1]*h_j

    Z = [const0/(math.cos(this_lambda*h_j)-1), const1/(math.cos(this_lambda*h_j)-1)]
    G = [Z[0]*(complex(math.cos(this_lambda*h_j), math.sin(this_lambda*h_j))-1), Z[1]*(complex(math.cos(this_lambda*h_j), math.sin(this_lambda*h_j))-1)]
    #print "...Constants...\n" + '{0:7.5f} + {1:7.5f}i, {2:7.5f} + {3:7.5f}i, {4:7.5f} + {5:7.5f}i, {6:7.5f} + {7:7.5f}i'.format(Z[0].real, Z[0].imag, Z[1].real, Z[1].imag, G[0].real, G[0].imag, G[1].real, G[1].imag) 

    #Allocate and fill the L-matrix (coefficients of I_g from low-order equations)
    L = npy.zeros((2, 2), dtype=complex)
    L[0][0] = (1.0-Z[0]*(sigT[0] - sigS[0][0])*h_j)
    L[0][1] = (Z[0]*(sigT[0] - sigS[0][0])*(sigT[1] - sigS[1][1])*h_j)/sigS[0][1]
    L[1][0] = Z[1]*sigS[0][1]*h_j
    L[1][1] = (1.0-Z[1]*(sigT[1] - sigS[1][1])*h_j)

    #Allocate and fill the T-matrix (coefficients of E_r,g from low-order equations)
    T = npy.zeros((2, 2*fperc), dtype=complex)
    for s in range(2*fperc):    
        #Every coefficient gets a potential contribution from 'r=1' J coeff:
        T[0][s] += G[0]*Final_J_coeffs.item(0,s)
        T[1][s] += G[1]*Final_J_coeffs.item(fperc,s)
    for a in range(fperc):
        for b in range(2*fperc):
            T[0][b] += 1.0/fperc * Final_coeffs.item(a,b)
            T[1][b] += 1.0/fperc * Final_coeffs.item(fperc+a,b)

    print "...L matrix..."
    for i in range(2): 
       line_string = ""        
       for c in range(2):
            line_string += '{0:7.5f} + {1:7.5f}i\t'.format(L.item(i,c).real, L.item(i,c).imag)    
       print line_string  
    #print "...T matrix..."
    for i in range(2): 
       line_string = ""        
       for c in range(2*fperc):
            line_string += '{0:7.5f} + {1:7.5f}i\t'.format(T.item(i,c).real, T.item(i,c).imag)    
    #   print line_string  

    L_inv_T = npy.matrix(L).I * npy.matrix(T)

    #print "...L_inv * T matrix..."
    for i in range(2): 
       line_string = ""        
       for c in range(2*fperc):
            line_string += '{0:7.5f} + {1:7.5f}i\t'.format(L_inv_T.item(i,c).real, L_inv_T.item(i,c).imag)    
       #print line_string     

    #Allocate the K-matrix
    K = npy.zeros((2*fperc, 2*fperc), dtype=complex)
    
    for c in range(fperc): 
        for d in range(2*fperc):
            #Coefficients relating F_r,g to E_r,g 
            K[c][d] += Final_coeffs.item(c,d)
            K[fperc+c][d] += Final_coeffs.item(fperc+c,d) 
            #Coeffs relating I_g to E_r,g
            K[c][d] += L_inv_T.item(0,d)
            K[fperc+c][d] += L_inv_T.item(1,d)
        for a in range(fperc):
            for b in range(2*fperc):
                #Coefficients relating (1/p * sum_r'=1^fperc F_r',g) to E_r,g
                K[c][b] += -1.0/fperc * Final_coeffs.item(a,b)
                K[fperc+c][b] += -1.0/fperc * Final_coeffs.item(fperc+a,b)
    #print "...K matrix..."
    for i in range(2*fperc): 
       line_string = ""        
       for c in range(2*fperc):
            line_string += '{0:7.5f} + {1:7.5f}i\t'.format(K.item(i,c).real, K.item(i,c).imag)    
       #print line_string     
  
    #Now, FINALLY, calculate the eigenvalues!
    eigs = npy.linalg.eigvals(K)
    
    #print '\n'
    for h in range(fperc):
        if math.fabs(eigs.item(h).real) > max_eig:
            max_eig = math.fabs(eigs.item(h).real)
            max_eig_imag = math.fabs(eigs.item(h).imag)
            max_lambda = this_lambda 
            max_f = h+1
        #print '{0:7.5f}: {1:7.5f} + {2:7.5f}i'.format(this_lambda, eigs.item(h).real, eigs.item(h).imag)
    exit(1)  

print "\n"        
print "Spectral radius = " + str(max_eig)
if(max_eig_imag > 1e-7):
    print "...uh oh, there's an imaginary part: " + str(max_eig_imag)
print "   ...at lambda = " + str(max_lambda)
print "    ...and cell = " + str(max_f)      
            
        
    
    
    
    
