import numpy as npy
import scipy as spy
import scipy.linalg as spyla
import sys
import math

#Make this whole shit a FUNCTION so we can put the result into a shell variable:
def calc_spectral_radius_implicit():
    #first things first: get relevant data (M, sigT, h_j & h_k, etc) from the file
    #execfile(sys.argv[1])
    #Let's just go ahead and hard-code the input file name
    try:
        #Load into the global namespace
        execfile("2g_test_input.inp", globals())
    except:
        print "Boo!"
    #fill weights and mus
    if M > 0:
        execfile("weights_mus.inp", globals())
    else:
        print "Uh-oh, we didn't get a valid value for M from " + sys.argv[1] + "!"
        sys.exit(1)

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

    #initialize max eig, f-value and corresponding lambda:
    max_lambda = 0.0
    max_f = 0
    max_g = 0
    max_eig = 0.0

    #SET A RANGE OF LAMBDA VALUES:
    #Try stopping just a ways past 2 pi...
    #print "ONLY SEARCHING FROM 6 to 6.5! BE CAREFUL!"
    lambda_range = npy.arange((2*math.pi)/X, (X/h_j - 1)*(2*math.pi)/X, (2*math.pi)/X)
    #lambda_range = npy.arange((2*math.pi)/X, (2*math.pi)/X+0.01, (2*math.pi)/X)
    for this_lambda in lambda_range:
        #print "\r ****Running lambda =  {0:5.3f}; current max = {1:7.5f} @ {2:7.5f}****".format(this_lambda, max_eig, max_lambda),
        #sys.stdout.flush()


        #allocate arrays for A, B, C, D:
        A = npy.zeros(((2*M*fperc + 2), (2*M*fperc + 2)), dtype=complex)
        B = npy.zeros(((2*M*fperc + 2), 2*fperc), dtype=complex)
        C = npy.zeros((2*fperc, (2*M*fperc + 2)), dtype=complex)
        #Not really needed... D = npy.zeros((2*fperc, 2*fperc), dtype=complex)

        #fill matrix A:
        #print " ...Matrix A... "

        #First M*fperc equations are coefficients of A, I from 3a for group 1 (I coeffs=0): 
        for i in range(M*fperc):
            this_m = math.ceil(float(i+1)/fperc)
            this_f = (i+1) - ((this_m-1)*fperc)
            this_mu = mus.item(math.floor((this_m-1)/2.0))*math.pow(-1.0,this_m)
            alfa = alfas[0][int(this_m-1)]
            #print "m and f " + str(this_m) + " " + str(this_f)
            #print "Mu: " + str(this_mu) + ", alpha: " + str(alfa) 
            
            #For each row, there will be 2 single entries: 
            #A_(r+1, m, 1), A_(r, m, 1)
            if this_f==fperc:
                value1 = (this_mu/h_k + 0.5*sigT[0]*(1.0+alfa))
                value2 = this_lambda*h_j
                #print "value 1: " + str(value1) + " & 2: " + str(value2) 
                A[i][(this_m-1)*fperc] += value1*complex(math.cos(value2), math.sin(value2))
                #print "cos: "
            else:
                A[i][i+1] += (this_mu/h_k + 0.5*sigT[0]*(1.0+alfa))
            A[i][i] += (-this_mu/h_k + 0.5*sigT[0]*(1.0-alfa))
            
            #AND two entries that are summed over M:
            #sum_1^M A_(r+1, m, 1), sum_1^M A_(r, m, 1)
            if this_f==fperc:
                value2 = this_lambda*h_j
                for p in range(M):
                    p_alfa = alfas[0][p]
                    A[i][p*fperc] += -0.25*sigS[0][0]*(1+p_alfa)*wts.item(math.floor(p/2.0))*complex(math.cos(value2), math.sin(value2))           
            else:
                for p in range(M):
                    p_alfa = alfas[0][p]
                    A[i][(this_f-1) + p*fperc + 1] += -0.25*sigS[0][0]*(1+p_alfa)*wts.item(math.floor(p/2.0))
            for p in range(M):
                p_alfa = alfas[0][p]
                A[i][(this_f-1) + p*fperc] += -0.25*sigS[0][0]*(1-p_alfa)*wts.item(math.floor(p/2.0))
            
            #line_string = ""
            #print '\n'
            #for c in range(2*M*fperc):
            #    line_string += str(A.item(i, c)).format('f') + "\t" 
            #    line_string += '{0:7.5f}\t'.format(A.item(i,c))     
            #print line_string          
        #Second M*fperc are coefficients of A, I, from 3a for group 2 (I coeffs=0):
        for m in range(M*fperc):
            this_m = math.ceil(float(m+1)/fperc)
            this_f = (m+1) - ((this_m-1)*fperc)
            this_mu = mus.item(math.floor(float(this_m-1)/2))*math.pow(-1.0,this_m)
            alfa = alfas[1][int(this_m-1)]
            #print "m and f " + str(this_m) + " " + str(this_f) + " " + str(this_mu)

            #For each row, there will be 2 single entries: 
            #A_(r+1, m, 2), A_(r, m, 2)
            if this_f==fperc:
                value1 = (this_mu/h_k + 0.5*sigT[1]*(1.0+alfa))
                value2 = this_lambda*h_j
                #print "value 1: " + str(value1) + " & 2: " + str(value2) 
                A[M*fperc+m][M*fperc+(this_m-1)*fperc] += value1*complex(math.cos(value2), math.sin(value2))
                #print "cos: "
            else:
                A[M*fperc+m][M*fperc+m+1] += (this_mu/h_k + 0.5*sigT[1]*(1.0+alfa))
            A[M*fperc+m][M*fperc+m] += (-this_mu/h_k + 0.5*sigT[1]*(1.0-alfa))
            
            #AND FOUR entries that are summed over M:
            #sum_1^M A_(r+1, m, 2), sum_1^M A_(r, m, 2)
            if this_f==fperc:
                value2 = this_lambda*h_j
                for p in range(M):
                    p_alfa = alfas[1][p]
                    A[M*fperc+m][M*fperc+p*fperc] += -0.25*sigS[1][1]*(1+p_alfa)*wts.item(math.floor(p/2.0))*complex(math.cos(value2), math.sin(value2))           
            else:
                for p in range(M):
                    p_alfa = alfas[1][p]
                    A[M*fperc+m][M*fperc+(this_f-1) + p*fperc + 1] += -0.25*sigS[1][1]*(1+p_alfa)*wts.item(math.floor(p/2.0))
            for p in range(M):
                p_alfa = alfas[1][p]
                A[M*fperc+m][M*fperc+(this_f-1) + p*fperc] += -0.25*sigS[1][1]*(1-p_alfa)*wts.item(math.floor(p/2.0))
          
            #sum_1^M A_(r+1, m, 1), sum_1^M A_(r, m, 1)      
            if this_f==fperc:
                value2 = this_lambda*h_j
                for p in range(M):
                    p_alfa = alfas[0][p]
                    A[M*fperc+m][p*fperc] += -0.25*sigS[0][1]*(1+p_alfa)*wts.item(math.floor(p/2.0))*complex(math.cos(value2), math.sin(value2))           
            else:
                for p in range(M):
                    p_alfa = alfas[0][p]
                    A[M*fperc+m][(this_f-1) + p*fperc + 1] += -0.25*sigS[0][1]*(1+p_alfa)*wts.item(math.floor(p/2.0))
            for p in range(M):
                p_alfa = alfas[0][p]
                A[M*fperc+m][(this_f-1) + p*fperc] += -0.25*sigS[0][1]*(1-p_alfa)*wts.item(math.floor(p/2.0))

            
            #line_string = ""
            #print '\n'
            #for c in range(2*M*fperc):
            #    line_string += str(A.item(i, c)).format('f') + "\t" 
            #    line_string += '{0:7.5f}\t'.format(A.item(M*fperc+m,c))     
            #print line_string        
        
        #Final 2 equations are coefficients of A, I from 3i/j/g combo, group 1 & 2 
        #calculate two constants G_1 and G_2:
        G1 = 2.0/(3.0 * sigT[0] * h_j) * (complex(math.cos(this_lambda*h_j), 0.0) - 1.0) 
        G2 = 2.0/(3.0 * sigT[1] * h_j) * (complex(math.cos(this_lambda*h_j), 0.0) - 1.0) 
        
        #Group 1 equation has two single coefficients (of I_1 and I_2)
        A[2*M*fperc][2*M*fperc] = -G1 + (sigT[0] - sigS[0][0])*h_j
        A[2*M*fperc][2*M*fperc + 1] = -(sigT[0] - sigS[0][0])*(sigT[1] - sigS[1][1])/sigS[0][1]*h_j
        #AND three summed entries sum_1^M mu_m A_(1,m,1), sum_1^M A_(1,m,1),
        for p in range(M):
            this_mu = mus.item(math.floor(p/2.0))*math.pow(-1.0,p+1)
            #A (1,m,1) sum terms
            A[2*M*fperc][p*fperc] += wts.item(math.floor(p/2.0)) * ((this_mu + G1/fperc)*complex(math.cos(this_lambda*h_j), math.sin(this_lambda*h_j)) + G1/fperc - this_mu)
            #and finally sum_1^M mu_m sum_2^(p-1) A_(r,m,1), 
            for r in range(fperc-1):
                A[2*M*fperc][p*fperc + r+1] += 2.0 * G1/fperc * wts.item(math.floor(p/2.0))
     
                
        #Group 2 equation has two single coefficients (of I_1 and I_2)
        A[2*M*fperc + 1][2*M*fperc + 1] = -G2 + (sigT[1] - sigS[1][1])*h_j
        A[2*M*fperc + 1][2*M*fperc] = -sigS[0][1]*h_j        
        #AND two summed entries sum_1^M A_(1,m,2),
        for p in range(M):
            this_mu = mus.item(math.floor(p/2.0))*math.pow(-1.0,p+1)
            #A (1,m,1) sum terms
            A[2*M*fperc + 1][M*fperc + p*fperc] += wts.item(math.floor(p/2.0)) * ((this_mu + G2/fperc)*complex(math.cos(this_lambda*h_j), math.sin(this_lambda*h_j)) + G2/fperc - this_mu)
            #and finally sum_1^M mu_m sum_2^(p-1) A_(r,m,2),
            for r in range(fperc-1):
                A[2*M*fperc + 1][M*fperc + p*fperc + r+1] += 2.0 * G2/fperc * wts.item(math.floor(p/2.0))  

        #for i in range(2*M*fperc + 2):
        #    line_string = ""
        #    for c in range(2*M*fperc + 2):
            #    line_string += str(A.item(i, c)).format('f') + "\t" 
        #        line_string += '{0:7.5f}\t'.format(A.item(i,c))     
            #print line_string
        #A_det = npy.linalg.det(A)
        #A_cn = npy.linalg.cond(A)
        #print "Determinant of A: " + str(A_det) + ", condit. # " + str(A_cn) + " digit loss, ~" + str(math.log(A_cn)) 
         
        #fill matrix B:
        #First M*fperc equations have E coefficients from 3a, group 1 eq 
        for k in range(M*fperc):        
            this_m = math.ceil(float(k+1)/fperc)
            this_f = (k+1) - ((this_m-1)*fperc)
            
            B[k][fperc + this_f-1] += -0.5*(sigT[0] - sigS[0][0])*(sigT[1] - sigS[1][1])/sigS[0][1]
        #Second M*fperc equations have E coefficients from 3a, group 2 eq (all = 0)
        #Final 2 equations have E coefficients from 3i/j/g, group 1 & 2 (all = 0)

        #print " ...Matrix B... "

        #for j in range(2*M*fperc + 2):
        #    line_string = ""
        #    for c in range(2*fperc):
        #        line_string += '{0:7.5f}\t'.format(B.item(j, c)) 
        #    print line_string
            
            #line_string = ""
            #for c in range(fperc):
            #    line_string += '{0:7.5f}\t'.format(B.item(M*fperc+k, c)) 
            #print line_string  
              

        #fill matrix C:
        #print "...Matrix C..."
        
        #First fperc are coeffs of I & A from 3k, group 1
        for l in range(fperc):
            #Coefficient of Ig:
            C[l][2*M*fperc] += rel_fact
            #ALL A entries are summed over all M
            for o in range(M):
                this_wt = wts.item(math.floor(o/2)) 
                alfa = alfas[0][o]
                if(l+1 == fperc):
                    C[l][o*fperc] += this_wt*0.5*(1+alfa)*complex(math.cos(this_lambda*h_j), math.sin(this_lambda*h_j))
                else:
                    C[l][l+1 + o*fperc] += this_wt*0.5*(1+alfa)
                C[l][l + o*fperc] += this_wt*0.5*(1-alfa)
                C[l][o*fperc] += -1.0*rel_fact/fperc*this_wt*(complex(math.cos(this_lambda*h_j), math.sin(this_lambda*h_j)) + 1)  
                
                for r in range(fperc-1):
                    C[l][o*fperc + r+1] += -2.0*rel_fact/fperc * this_wt   
                
        #Second fperc are coeffs of I & A from 3k, group 2 
        for k in range(fperc):
            #Coefficient of Ig:
            C[fperc + k][2*M*fperc + 1] += rel_fact
            #ALL A entries are summed over all M
            for o in range(M):
                this_wt = wts.item(math.floor(o/2)) 
                alfa = alfas[0][o]
                if(k+1 == fperc):
                    C[fperc + k][M*fperc + o*fperc] += this_wt*0.5*(1+alfa)*complex(math.cos(this_lambda*h_j), math.sin(this_lambda*h_j))
                else:
                    C[fperc + k][M*fperc + k+1 + o*fperc] += this_wt*0.5*(1+alfa)
                C[fperc + k][M*fperc + k + o*fperc] += this_wt*0.5*(1-alfa)
                C[fperc + k][M*fperc + o*fperc] += -1.0*rel_fact/fperc*this_wt*(complex(math.cos(this_lambda*h_j), math.sin(this_lambda*h_j)) + 1)  
                
                for r in range(fperc-1):
                    C[fperc + k][M*fperc + o*fperc + r+1] += -2.0*rel_fact/fperc * this_wt  
               
        #for n in range(2*fperc):        
        #    line_string = ""
        #    for c in range(2*M*fperc+2):
        #        line_string += '{0:7.5f}\t'.format(C.item(n, c)) 
        #    print line_string  

        #Convert array A to matrix:
        A_m = npy.matrix(A)
        A_m_inv = A_m.I

        #print " ...Matrix A^-1... "
        #for i in range(2*M*fperc + 2):
        #    line_string = ""
        #    for c in range(2*M*fperc + 2):
            #    line_string += str(A.item(i, c)).format('f') + "\t" 
        #        line_string += '{0:7.5f}\t'.format(A_m_inv.item(i,c))     
        #    print line_string


        B_m = npy.matrix(B)
        C_m = npy.matrix(C)
         
        #print out inverse of A:
        #print "...Matrix A^-1..."

        #for a in range(2*M*fperc):
        #    line_string = ""
        #    for b in range(2*M*fperc):
        #        line_string += '{0:7.5f}\t'.format(A_m_inv.item(a, b))
        #    print line_string 
            
        #Find product -C*A^-1*B
        E = npy.zeros((2*fperc, 2*fperc), dtype=complex)
        E = -C_m*(A_m_inv*B_m)

        #print "...Matrix C*A^-1*B..."

        #for d in range(2*fperc):
        #    line_string = ""
        #    for e in range(2*fperc):
        #        line_string += '{0:7.5f}\t'.format(E.item(d, e))
        #    print line_string
            
        eigs = npy.array((fperc), dtype=complex)
        eigs = npy.linalg.eigvals(E)

        for h in range(2*fperc):
            if math.fabs(eigs.item(h).real) > max_eig:
                max_eig = math.fabs(eigs.item(h).real)
                max_lambda = this_lambda 
                max_f = math.ceil((h+1)/2)
                max_g = math.ceil((h+1)/fperc)
            #print '{0:7.5f}: {1:7.5f}'.format(this_lambda, eigs.item(h))    
    return max_eig

#This is the "main-like" section-- it just calls the defined function and PRINTS the result so the shell can grab it
print calc_spectral_radius_implicit()


