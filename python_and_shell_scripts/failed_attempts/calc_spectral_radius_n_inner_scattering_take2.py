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
max_f = 0
max_eig = 0.0

#SET A RANGE OF LAMBDA VALUES:
#Try stopping just a ways past 2 pi...
#print "ONLY SEARCHING FROM 6 to 6.5! BE CAREFUL!"
lambda_range = npy.arange(0.01, 7.01, 0.01)
for this_lambda in lambda_range:
    print "\r ****Running lambda =  {0:5.3f}; current max = {1:7.5f} @ {2:7.5f}****".format(this_lambda, max_eig, max_lambda),
    sys.stdout.flush()


    #allocate arrays for A, B, C, D:
    A = npy.zeros((2*M*fperc, 2*M*fperc), dtype=complex)
    B = npy.zeros((2*M*fperc, fperc), dtype=complex)
    C = npy.zeros((fperc, 2*M*fperc), dtype=complex)
    D = npy.zeros((fperc, fperc), dtype=complex)

    #fill matrix A:
    #print " ...Matrix A... "

    #First M*fperc equations are coefficients of a, b from 4a: 
    for i in range(M*fperc):
        this_m = math.ceil(float(i+1)/fperc)
        this_f = (i+1) - ((this_m-1)*fperc)
        #print "m and f " + str(this_m) + " " + str(this_f)
        this_mu = mus.item(math.floor(float(this_m-1)/2))*math.pow(-1.0,this_m)
        alfa = alfas[int(this_m-1)]
        #print "Mu: " + str(this_mu) + ", alpha: " + str(alfa) 
        
        #For each row, there will be 3 non-zero entries: 
        #a_(m,n), a_(m,n+1), b_(m,n)
        if this_f==fperc:
            value1 = -(0.5*(1+alfa))
            value2 = sigT*this_lambda*fperc*h_k
            #print "value 1: " + str(value1) + " & 2: " + str(value2) 
            A[i][(this_m-1)*fperc] += complex((value1*math.cos(value2)), (value1*math.sin(value2)))
            #print "cos: "
        else:
            A[i][i+1] += -0.5*(1+alfa)
        A[i][i] += -0.5*(1-alfa)
        A[i][M*fperc+i] += 1.0
        #line_string = ""
        #for c in range(2*M*fperc):
            #line_string += str(A.item(i, c)).format('f') + "\t" 
        #    line_string += '{0:7.5f}\t'.format(A.item(i,c))     
        #print line_string   
                
    #Second M*fperc are coefficients of a, b, from 4b
    for m in range(M*fperc):
        this_m = math.ceil(float(m+1)/fperc)
        this_f = (m+1) - ((this_m-1)*fperc)
        this_mu = mus.item(math.floor(float(this_m-1)/2))*math.pow(-1.0,this_m)
        #print "m and f " + str(this_m) + " " + str(this_f) + " " + str(this_mu)
        
        #For each row, there will be M+3 non-zero entries: 
        #a_(m,n), a_(m,n+1), and b_(m',n) for ALL m' in M
        if this_f==fperc:
            value1 = this_mu/(h_k*sigT)
            value2 = sigT*this_lambda*fperc*h_k
            A[M*fperc+m][(this_m-1)*fperc] += complex((value1*math.cos(value2)), (value1*math.sin(value2)))
        else:
            A[M*fperc+m][m+1] += (this_mu/(h_k*sigT))
        A[M*fperc+m][m] += -(this_mu/(h_k*sigT))
        A[M*fperc+m][M*fperc+m] += 1.0
        
        #These are the coefficients of b (summed over all m for this n)
        for p in range(M):
            A[M*fperc+m][M*fperc + (this_f-1) + p*fperc] += -0.5*s_rat*((pow(s_rat, n_in-1) + 1)/(pow(s_rat, n_in) + 1))*wts.item(math.floor(p/2))
        #line_string = ""
        #for c in range(2*M*fperc):
            #line_string += str(A.item(M*fperc+m, c)).format('5.2f') + "\t"
            #line_string += '{0:7.5f}\t'.format(A.item(M*fperc+m, c))      
        #print line_string
        
    #fill matrix B:
    #print " ...Matrix B... "

    #for j in range(M*fperc):
        #line_string = ""
        #for c in range(fperc):
        #    line_string += '{0:7.5f}\t'.format(B.item(j, c)) 
        #print line_string
    #First M*fperc equations have no D_(n) coefficients    
    for k in range(M*fperc):        
        this_m = math.ceil(float(k+1)/fperc)
        this_f = (k+1) - ((this_m-1)*fperc)
        
        B[M*fperc+k][this_f-1] += -(0.5*(1-s_rat))
        
        #line_string = ""
        #for c in range(fperc):
        #    line_string += '{0:7.5f}\t'.format(B.item(M*fperc+k, c)) 
        #print line_string  
          

    #fill matrix C:
    #print "...Matrix C..."

    #Calc big ugly const:
    constant = 1.5 * (sigT * h_j) 
    big_const_real = constant

    #Need a case here to prevent python from getting mad:
    if this_lambda != 0: 
        big_const_imag = constant * ((math.sin(sigT*this_lambda*fperc*h_k))/(math.cos(sigT*this_lambda*fperc*h_k) - 1))
    else:
        big_const_imag = 0.0
        
    for n in range(fperc):
        #Coefficients of a_(m,1) and b_(m,n)-- all m's for this n
        for o in range(M):
            this_mu = mus.item(math.floor(o/2))*math.pow(-1.0, (o+1))
            this_wt = wts.item(math.floor(o/2))  
            #print "this mu and this weight: " + str(this_mu) + " " + str(this_wt)   
            C[n][0 + o*fperc] += complex((this_wt * big_const_real * this_mu), (this_wt * big_const_imag * this_mu)) 
            C[n][M*fperc + n + o*fperc] += this_wt
            
        #line_string = ""
        #for c in range(2*M*fperc):
        #    line_string += '{0:7.5f}\t'.format(C.item(n, c)) 
        #print line_string  

    #Convert array A to matrix:
    A_m = npy.matrix(A)
    A_m_inv = A_m.I

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
    E = npy.zeros((fperc, fperc), dtype=complex)
    E = -C_m*(A_m_inv*B_m)

    #print "...Matrix C*A^-1*B..."

    #for d in range(fperc):
    #    line_string = ""
    #    for e in range(fperc):
    #        line_string += '{0:7.5f}\t'.format(E.item(d, e))
    #    print line_string
        
    eigs = npy.array((fperc), dtype=complex)
    eigs = npy.linalg.eigvals(E)

    for h in range(fperc):
        if math.fabs(eigs.item(h).real) > max_eig:
            max_eig = math.fabs(eigs.item(h).real)
            max_lambda = this_lambda 
            max_f = h+1
        #print '{0:7.5f}: {1:7.5f}'.format(this_lambda, eigs.item(h))
 
print "\n"        
print "Spectral radius = " + str(max_eig)
print "   ...at lambda = " + str(max_lambda)
print "    ...and cell = " + str(max_f)

