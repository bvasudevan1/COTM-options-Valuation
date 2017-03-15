# -*- coding: utf-8 -*-
"""
@author: BVasudevan

WIP: to calculate the price of a call to the max option.
NYI: Add module to pull in data from Quandl/Bloomberg, calc corr_coef, SDs and evaluate price. 
"""
import numpy as np
import matplotlib.pyplot as plt
import pdb

class value_options():
    global r
    global n
    global rho
    global expectedPayoff
    global payoffs
    global x
    r=0.05
    payoffs=[]    
    x=0
    
    def normalPairs(n):
        """Uses the marsaglia polar method to generate pairs of normals"""
        x=[]
        y=[]
        i=0
        while i<n:
            u=np.random.uniform(0,1,2)
            v=[2*u[0]-1,2*u[1]-1]
            s=np.dot(v,v)
            if s<=1:
                f=np.sqrt(-2*np.log(s)/s)
                x.append(f*v[0])
                y.append(f*v[1])
                i=i+1
        return x,y
    
    def txBVnormal(iidStd,rho):
        X_1 = []
        X_2 = []
        iidStd[0]
        r=0.05
        sig=[0.2,0.2]
        for i in range(n):
            X_1.append(r+sig[0]*(iidStd[0][i]-0.5*sig[0]))
            X_2.append(r+sig[1]*(-0.5*sig[1]+rho*iidStd[0][i]+np.sqrt(1-rho**2)*iidStd[1][i]))
        return X_1,X_2
    
    def logNormal(X_i):
        logX_ = []
        init = 100
        for i in range(n):
            logX_.append(init*np.exp(X_i[i]))
        return logX_
    
    def COTM(maxPair,k):
        dp=0
        for i in range(len(maxPair)):
            a=maxPair[i]-k
            dp=dp+(np.exp(-r)*max(a,0))
        V=dp/len(maxPair)
        
        return V
    
    def ValVsCor():
        global x
        rho = [-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        while x<19:
            return rho[x]
        
    def val_estimates():
        for i in range(0,100):
            value_options.main()
        
        
    def main():
        global rho
        global r
        global n
        global x
        global payoffs
        print("Enter the desired number of normal pairs")
        n = int(input())
        print("Enter the Correlation Coefficient")
        rho = float(input())
    
        PH=value_options.normalPairs(n) 
        """Note: while the normals are stored in 2 series
           they were generated as a pair using MS-polar"""
        X_1 = value_options.txBVnormal(PH,rho)[0]
        X_2 = value_options.txBVnormal(PH,rho)[1]
                
        logX_1 = value_options.logNormal(X_1)
        logX_2 = value_options.logNormal(X_2)
        
        plt.figure(0)
        plt.scatter(logX_1,logX_2,s=1)
        
        
        a=np.array([logX_1,logX_2]).T
        MaxPair=[]
        for i in range(len(a)):
            MaxPair.append(max(a[i,0],a[i,0]))
        plt.figure(1)
        plt.hist(MaxPair,bins=100)
        
                    
        print("Enter a value for strike price")
        k= float(input())
        
        expectedPayoff = value_options.COTM(MaxPair,k)
        payoffs.append(expectedPayoff)
        print('Payoff of this COTM option with the given parameters is ', expectedPayoff)
        #payoffCurve_Corr.append(expectedPayoff)
        
        
    

        #plt.plot(payoffCurve_Corr)
       # print(X)
       
        
        
#        r=0.05
#        sig=[0.2,0.2]
#        init_cond=[100,100]
#        rho=0

        
        


                
            
            