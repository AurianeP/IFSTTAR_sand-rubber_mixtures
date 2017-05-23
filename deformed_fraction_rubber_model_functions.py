# -*- coding: utf-8 -*-
"""
Created on Mon May 22 17:55:58 2017

@author: Auriane
"""
import numpy as np
from scipy.optimize import curve_fit,minimize,leastsq

#==============================================================================
# Error function for polyfit
#==============================================================================
def error(a,x,data): # return difference between data and p(x), polynom of x with coefficients a(k)
    poly=np.poly1d(a)
    return poly(x)-data

#==============================================================================
# Polynomial fit via logarithm derivation 
#==============================================================================
def Polyfit(pnew_cut,f_data,j0,i0,fid,fid1):
        ### Polynom over f: Choose the right degree
        print "Choosing right order of the polynom to fit experimental data"
        #*** initialize parameters
#        epsilon_de=0.05 #1%
#        epsilon_exp=0.005 #experimental uncertainty
#        epsilon_e1=size(pnew_cut)*epsilon_exp**2 #maximum error in the polynomial fit
        d=0
        dmax=3 # max order of tested polynom
        if (j0==0 and i0==4):
            dmax=3
        if (j0==1 and i0==0):
            dmax=2
        random_choice=-3*np.random.rand(d+1)
        init_values=np.power(10,random_choice)
        a=init_values
        e0=10E-10 # e0 very small : de will be very big at the end of the first iteration
#        e1=2*epsilon_e1
#        de=2*epsilon_de # this way (e1-e0)/e0=2*epsilon>epsilon always
        while(d<=dmax-1): # (e1>epsilon_e1) and (de>epsilon_de) and 
            #*** update parameters
            d=d+1
#            print a
#            print(np.poly1d(a))
#            init_values=r_[0.0,a]
            init_values=np.insert(a,0,0.)        
#            print init_values            
#            print(np.poly1d(init_values))
            # initial guess for polynom P(d+1) is given by optimized coefficients of P(d),
            # constant coefficient is randomly initialized
            #*** minimization of least squared error:
            res = leastsq(error,init_values,args=(pnew_cut,f_data),full_output=True,maxfev=15)
            a=res[0]
            print a
            #*** info about the fit: 
            infodict=res[2]
            mesg=res[3]
            flag=res[4]
#            print "\nPolynom of d-order: d=",d
#            print "optimization success (flag >0): flag=",flag
#            print "message: ",mesg
#            print "number of iteration",infodict['nfev']
#            print "initial guess (random):",init_values
#            print "result:",a
            #*** least squared error:
            e1=np.sum(infodict['fvec']**2)
#            print "least squared error:",e1
            de = abs((e1-e0)/e0)
#            print "error de:",de
            e0=e1
            #*** write to poly log file
            fid1.write('Polynom of '+str(d)+'-order\nInfo about the fit\n')
            fid1.write('optimization success (if flag >0): flag='+str(flag)+'\n')
            fid1.write('message:'+mesg+'\n')
            fid1.write('Initial guess:\na0(k)\n')
            for k in range(np.size(init_values)):
                fid1.write( '{:2.2e}\t'.format(init_values[k]))            
            fid1.write('\nOptimized results:\na(k)\n')
            for k in range(np.size(a)):
                fid1.write( '{:2.2e}\t'.format(a[k]))
            fid1.write('\nLeast squared error e1 = {:2.2e}\n'.format(e1))
            fid1.write('Error convergence : de = (e1-e0)/e0 = {:2.2e}\n'.format(de))
            
            
        ### Compute fitted polynom f_spline                
        P=np.poly1d(a)
#        print(P)
        f_spline=P(pnew_cut)
        error_spline = 100*infodict['fvec']**2

        fid.write("Chosen polynom of {:d}-order\nwith {:d} coefficients a(k):\n".format(d,d))
        for k in range(d):
            fid.write('{:2.2e}'.format(a[k])+'x^{:d}\t'.format(d-k))
        fid.write("\n")
        
        ### Derivative of the polynom P -> P_prime
        print '\nDerivative of the polynom'
        P_prime=np.polyder(P,m=1)
        d_prime=d-1
        f_spline_prime=P_prime(pnew_cut)
        fid.write( "polynom derivative coefficients"+"\n")
        for k in range(d_prime):
            fid.write('{:2.2e}'.format(P_prime.coeffs[k])+'x^{:d}\t'.format(d_prime-k))
        fid.write('\n')
        
        ### Log of the derivative P2
        ind=np.where(f_spline_prime>0)
        pnew_log=pnew_cut[ind]
        f_spline_prime=f_spline_prime[ind]
        ln_f_spline_prime = np.log(f_spline_prime)
        
        ### Linear regression of the log of the derivative
        print 'Linear regression of the log of the derivative'
        #*** initialize parameters
        init_value1=-3*np.random.rand(2)
        init_value1=np.power(10,init_value1)
        #*** minimization of least squared error:
        res1 = leastsq(error,init_value1,args=(pnew_log,ln_f_spline_prime),full_output=True,maxfev=15)
        a1=res1[0]
        #*** info about the fit: 
        infodict=res1[2]
        mesg=res1[3]
        flag=res1[4]
#        print "initial guess (random):",init_value1
#        print "result:",a1
#        print "optimization success (flag >0): flag=",flag
#        print "message: ",mesg
#        print "number of iteration",infodict['nfev']
        #*** compute linear regression:
        P1=np.poly1d(a1)
        ln_f_spline_prime_lin=P1(pnew_log)
        #*** write to logfile
        fid.write('Linear regression of the log of the derivative : ax+b\n')
        fid.write('Info about the fit\n')
        fid.write('optimization success (if flag >0): flag='+str(flag)+'\n')
        fid.write('message:'+mesg+'\n')
        fid.write('Results:\na\tb\n')
        fid.write( '{:2.2e}'.format(a1[0])+"\t"+'{:2.2e}'.format(a1[1])+"\n")
        
        return a1,d,pnew_log,f_spline,error_spline,ln_f_spline_prime,ln_f_spline_prime_lin