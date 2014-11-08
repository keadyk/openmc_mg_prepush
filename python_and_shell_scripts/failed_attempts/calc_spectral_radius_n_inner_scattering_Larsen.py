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
    print "Mu: " + str(this_mu) + ", alpha: " + str(alfas[-1]) 

#initialize max eig, f-value and corresponding lambda:
max_lambda = 0.0
max_f = 0
max_eig = 0.0

#SET A RANGE OF LAMBDA VALUES:
#Try stopping just a ways past 2 pi...
#print "ONLY SEARCHING FROM 6 to 6.5! BE CAREFUL!"
lambda_range = npy.arange(0.0, 0.5*(2*math.pi)/10.0, (2*math.pi)/10.0)
for this_lambda in lambda_range:
    #print "\r ****Running lambda =  {0:5.3f}; current max = {1:7.5f} @ {2:7.5f}****".format(this_lambda, max_eig, max_lambda),
    sys.stdout.flush()
    
    #allocate arrays:
    M_1 = npy.zeros((M*fperc, M*fperc), dtype=complex)
    M_2 = npy.zeros((M*fperc, M*fperc), dtype=complex)
    M2_M1i = npy.zeros((M*fperc, M*fperc), dtype=complex)
    M2_M1i_sum = npy.zeros((fperc, fperc), dtype=complex)
    M3 = npy.zeros((M*fperc, fperc), dtype=complex)
    H = npy.zeros((fperc, fperc), dtype=complex)
    H_k = npy.zeros((fperc, fperc), dtype=complex)
    H_k1 = npy.zeros((fperc, fperc), dtype=complex)
    H_k1Q = npy.zeros((fperc, fperc), dtype=complex)
    H_k_sum = npy.zeros((fperc, fperc), dtype=complex)
    H_k_H_k1Qsum = npy.zeros((fperc, fperc), dtype=complex)
    #M_1_inv = npy.zeros((fperc, 2*M*fperc), dtype=complex)
    Q = npy.zeros((fperc, fperc), dtype=complex)
    F_rE_r = npy.zeros((fperc), dtype=complex)
    
    #Coefficients of M_2
    #print "...Matrix M_2..."
    for i in range(M*fperc):
        this_m = math.ceil(float(i+1)/fperc)
        this_f = (i+1) - ((this_m-1)*fperc)
        this_mu = mus.item(math.floor(float(this_m-1)/2))*math.pow(-1.0,this_m)
        this_wt = wts.item(math.floor((this_m-1)/2))  
        alfa = alfas[int(this_m-1)]
        
        print "m, f, mu, wt, alfa: " + str(this_m) + " " + str(this_f) + " " + str(this_mu) + " " + str(this_wt) + " " + str(alfa)
        #if this_mu > 0.0:
        #    M_1[i][i] = this_mu/(sigT*h_k) + 0.5*(1 + alfa)
        #    M_1[i][i-1] = -this_mu/(sigT*h_k) + 0.5*(1 - alfa)
        
            #NOTE!! We're going to "fold" the weight in here
            #because it's easiest to do it now!
        #    M_2[i][i] = 0.5*(1 + alfa)*this_wt
        #    M_2[i][i-1] = 0.5*(1 - alfa)*this_wt
        #else:
        #    M_1[i][i] = -this_mu/(sigT*h_k) + 0.5*(1 - alfa)
        #    M_1[i][i+1] = this_mu/(sigT*h_k) + 0.5*(1 + alfa)
        
            #NOTE!! We're going to "fold" the weight in here
            #because it's easiest to do it now!
        #    M_2[i][i] = 0.5*(1 - alfa)*this_wt
        #    M_2[i][i+1] = 0.5*(1 + alfa)*this_wt
        #Try changing the periodicity condition:
        if(this_f == fperc):
            M_1[i][i] += -this_mu/(sigT*h_k) + 0.5*(1 - alfa)
            #REAL PART ONLY!!
            M_1[i][(this_m-1)*fperc] += (this_mu/(sigT*h_k) + 0.5*(1 + alfa))*complex(math.cos(sigT*this_lambda*h_j), math.sin(sigT*this_lambda*h_j))
            
            M_2[i][i] += 0.5*(1 - alfa)*this_wt
            #REAL PART ONLY!!!
            M_2[i][(this_m-1)*fperc] += 0.5*(1 + alfa)*this_wt*complex(math.cos(sigT*this_lambda*h_j), math.sin(sigT*this_lambda*h_j))  
        else:
            M_1[i][i] += -this_mu/(sigT*h_k) + 0.5*(1 - alfa)
            M_1[i][i+1] += this_mu/(sigT*h_k) + 0.5*(1 + alfa) 
            
            M_2[i][i] += 0.5*(1 - alfa)*this_wt
            M_2[i][i+1] += 0.5*(1 + alfa)*this_wt  
        

        #line_string = ""
        #for c in range(M*fperc):
        #    line_string += '{0:7.5f}\t'.format(M_2.item(i,c))     
        #print line_string
        
    print "...Matrix M_1..."
    for i in range(M*fperc): 
       line_string = ""        
       for c in range(M*fperc):
            line_string += '{0:7.5f}\t'.format(M_1.item(i,c))     
       print line_string 
       
    print "...Matrix M_2..."
    for i in range(M*fperc): 
       line_string = ""        
       for c in range(M*fperc):
            line_string += '{0:7.5f}\t'.format(M_2.item(i,c))     
       print line_string 
       
    M_1_inv = npy.matrix(M_1).I
    
    print "...Matrix M_1_inv..."    
    for i in range(M*fperc): 
       line_string = ""        
       for c in range(M*fperc):
            line_string += '{0:7.5f}\t'.format(M_1_inv.item(i,c))     
       print line_string 
       
    M2_M1i = npy.matrix(M_2)*M_1_inv
    
    print "...Matrix M_2*M_1i..."
    for i in range(M*fperc): 
       line_string = ""        
       for c in range(M*fperc):
            line_string += '{0:7.5f}\t'.format(M2_M1i.item(i,c))     
       print line_string 
       
    #NEXT, we need to sum over M (angle):
    #for j in range(M*fperc):
    #    for r in range(fperc):
    #        for m in range(M):
            #Grab every "M"th value:
    #            M2_M1i_sum[r][r] += M2_M1i.item(j, r+fperc*m)
    for i in range(M*fperc):
        for j in range(fperc):
            for m in range(M):
                M3[i][j] += M2_M1i.item(i, j+m*fperc)
    #print "...Matrix M_3 summed..."
    #for i in range(M*fperc): 
    #   line_string = ""        
    #   for c in range(fperc):
    #        line_string += '{0:7.5f}\t'.format(M3.item(i,c))     
    #   print line_string 
    for j in range(fperc):
        for r in range(fperc):
            for m in range(M):
            #Grab every "M"th value:
                M2_M1i_sum[j][r] += M3.item((j+m*fperc, r))
    print "...Matrix M_2*M_1i summed..."
    for i in range(fperc): 
       line_string = ""        
       for c in range(fperc):
            line_string += '{0:7.5f}\t'.format(M2_M1i_sum.item(i,c))     
       print line_string 
       
    #Now, calculate H and Q (sans E_r):
    for i in range(fperc):
        for j in range(fperc):
            H[i][j] = M2_M1i_sum[i][j] * 0.5*s_rat
            Q[i][j] = M2_M1i_sum[i][j] * 0.5*(1-s_rat)

    print "...Matrix H..."
    for i in range(fperc): 
       line_string = ""        
       for c in range(fperc):
            line_string += '{0:7.5f}\t'.format(H.item(i,c))     
       print line_string
       
    print "...Matrix Q (sans E_r)..."
    for i in range(fperc): 
       line_string = ""        
       for c in range(fperc):
            line_string += '{0:7.5f}\t'.format(Q.item(i,c))    
       print line_string  
    
    #Now, calculate the error ratio F_r/E_r for each fine cell in the coarse cell!!!!
    #for j in range(fperc):
    #    this_sum = 0.0
    #    for k in range(n_in):
    #        #Calculate the matrix H^k
    #        this_sum += math.pow(H[j][j], k)*Q[j][j]
    #        print "added " + str(math.pow(H[j][j], k)*Q[j][j]) + " to sum " + str(this_sum) 
    #    print "...then added this to " + str(math.pow(H[j][j], n_in))
    #    F_rE_r[j] = math.pow(H[j][j], n_in) + this_sum

    #Set H_k-1 to H initially
    H_k1 = npy.identity(fperc, dtype=complex)
    #Also set H_k to H, really just for the case n_in = 1 
    #H_k = H
    H_k = npy.identity(fperc, dtype=complex)
    H_k_sum = npy.identity(fperc, dtype=complex)
    
    print "...Matrix H^0..." 
    for i in range(fperc): 
       line_string = ""        
       for c in range(fperc):
            line_string += '{0:7.5f}\t'.format(H_k1.item(i,c))    
       print line_string 
    
    for k in range(n_in-1):
        H_k = npy.zeros((fperc, fperc), dtype=complex)
        #Calculate the matrix H^k
        for a in range(fperc):
            for b in range(fperc):
                for c in range(fperc):
                    H_k[a][b] += H[a][c] * H_k1[c][b]
                    print "adding H(" + str(a) + ")(" + str(c) + ")= " + str(H.item(a,c)) + " x " + "Hk1(" + str(c) + ")(" + str(b) + ")= " + str(H_k1.item(a,c)) + "to H(" + str(a) + ")(" + str(b)
                H_k_sum[a][b] += H_k[a][b]
        H_k1= H_k
        print "...Matrix H^" + str(k+1) + "..." 
        for i in range(fperc): 
           line_string = ""        
           for c in range(fperc):
                line_string += '{0:7.5f}\t'.format(H_k.item(i,c))    
           print line_string  
        print "...Matrix H^" + str(k+1) + " sum..." 
        for i in range(fperc): 
           line_string = ""        
           for c in range(fperc):
                line_string += '{0:7.5f}\t'.format(H_k_sum.item(i,c))    
           print line_string  

    #Finally, calc H^n_in, mult. H^n' sum by Q, and add the two:
    #Reset H_k matrix first:
    H_k = npy.zeros((fperc, fperc), dtype=complex)
    for a in range(fperc):
        for b in range(fperc):
            for c in range(fperc):
                print "H entries: " + str(H[a][c]) + " " + str(Q[a][c]) + " " + str(H_k1[c][b])
                H_k[a][b] += H[a][c] * H_k1[c][b]
                H_k1Q[a][b] += Q[a][c] * H_k1[c][b]   
            print "H_k entry:" + str(H_k[a][b]) + ", H_k1Q: " + str(H_k1Q[a][b])
            H_k_H_k1Qsum[a][b] = H_k[a][b] + H_k1Q[a][b] 
    print "...Matrix H^(n-1)sum*Q..." 
    for i in range(fperc): 
       line_string = ""        
       for c in range(fperc):
            line_string += '{0:7.5f}\t'.format(H_k1Q.item(i,c))    
       print line_string 
    
    print "...H^k + Matrix H^(n-1)sum*Q..." 
    for i in range(fperc): 
       line_string = ""        
       for c in range(fperc):
            line_string += '{0:7.5f}\t'.format(H_k_H_k1Qsum.item(i,c))    
       print line_string 
               
    #Finally, calculate error ratio for each fine cell:  
    for r in range(fperc):         
        for p in range(fperc):
            print "Adding: " + str(H_k_H_k1Qsum[r][p])  
            F_rE_r[r] += H_k_H_k1Qsum[r][p]
        
    print "...ERROR RATIO after inner iterations..."
    line_string = ""        
    for c in range(fperc):
        line_string += '{0:11.9f}\t'.format(F_rE_r.item(c))    
    print line_string    
            
        
    
    
    
    
