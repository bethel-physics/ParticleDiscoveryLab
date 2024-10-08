# Created by Jack Sisson on 6/30/2020
# Email: jsis2000@outlook.com
# Translation of the matlab function pollsf.m https://github.com/AlejGarcia/NM4P

import numpy as np
import math

def pollsf(x, y, sigma, M):

    # Function to fit a polynomial to data
    # Inputs 
    #   x       Independent variable
    #   y       Dependent variable
    #   sigma   Estimate error in y
    #   M       Number of parameters used to fit data
    # Outputs:
    #   a_fit   Fit parameters; a(1) is intercept, a(2) is slope
    #   sig_a   Estimated error in the parameters a()
    #   yy      Curve fit to the data
    #   chisqr  Chi squared statistic
    
        
    #* Form the vector b and design matrix A
    b =np.zeros(len(y))
    for i in range(len(y)):
        b[i]= y[i]/sigma[i]

    N = len(x)
    A=np.zeros((N,M))
    
    for i in range(N):
        for j in range(M):
            A[i,j]= x[i]**(j)/sigma[i]
    
    #* Compute the correlation matrix C 
    trnspA = np.transpose(A)
    trnspA2=trnspA
    trnspA =trnspA.dot(A)
    C = np.linalg.inv(trnspA)
    
    #* Compute the least squares polynomial coefficients a_fit
    trnspB =np.transpose(b)        
    a_fit=trnspA2.dot(trnspB)
    a_fit=a_fit.dot(C)

    #* Compute the estimated error bars for the coefficients
    sig_a=np.zeros(M)
    for j in range(M):
        sig_a[j]= math.sqrt(abs(C[j,j]))
                

    #* Evaluate curve fit at each data point and compute Chi^2
    yy = np.zeros(len(x))

    for j in range(M):
        for i in range(len(x)):
            yy[i]=yy[i]+(x[i]**(j))*a_fit[j]

    chisqr = 0 
    for i in range(len(y)):
        chisqr=chisqr + ((y[i]-yy[i])/sigma[i])

    return a_fit, sig_a, yy, abs(chisqr)
