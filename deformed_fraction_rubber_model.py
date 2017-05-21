# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 20:49:17 2016

@author: Auriane
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Sep 08 10:59:16 2016

@author: Auriane
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FixedLocator,FormatStrFormatter,MultipleLocator
import itertools
from scipy.optimize import curve_fit,minimize,leastsq

plt.close('all')
#os.chdir("D:\\platzer\\Documents\\STAGE_AURIANE\\3_MACRO\\1_SAND_RUBBER_MIXTURE_MODEL\\")
num_trail=0
""" ###########################################################################
########################### Global variables ##################################
############################################################################"""

FRi = [0.05,0.10,0.20,0.30,0.50]
FRa = np.array(FRi)
pnew_length=1530
p_cut_HS=100.0
p_cut_LB=100.0
fnam_pure = ["NWH100R0C_new_dataexp_ep.txt","NWLB100R0C_new_dataexp_ep.txt"]
fnam_sand = ["LB","HS"]

p0_matrix=np.zeros((len(FRi),2))
F_matrix=np.zeros((len(FRi),2))
f_matrix=np.zeros((len(FRi),pnew_length))
erreur_m_list=[]
erreur_p_list=[]
pressure_list=[]

## PLOT TOOLS : 
#colormap = plt.cm.Greys
#couleur = [colormap(i) for i in np.linspace(0,0.8,len(FRi))]
savenam=['LB','HS']
markers=['o','^']
marker = itertools.cycle(('o', '^', 's', 'v', '*')) 
letter = itertools.cycle(('(a)','(b)','(c)'))#,'(d)','(e)'))
lines = itertools.cycle(("-","--","-.",":"))

colormap = plt.cm.Blues
couleur = [colormap(i) for i in np.linspace(0.5,1.0,len(FRi))]
nb_fig=5
plt.close()

""" ###########################################################################
################################ Functions ####################################
############################################################################"""

def error(a,x,data): # return difference between data and p(x), polynom of x with coefficients a(k)
    poly=np.poly1d(a)
    return poly(x)-data

def exp_etendu(x,p0,F,p_etoile,f_etoile):
    return (F-f_etoile)*(1-np.exp(-(x-p_etoile)/p0))+f_etoile
#expe_mod = Model(exp_etendu)

""" ###########################################################################
################################ Computation ##################################
############################################################################"""

for j0 in range(2): #[1]:#: #deux sables
    #*************************************************************************#
    ################################# LOG FILES ###############################
    #*** General log file
    fid = open(fnam_sand[j0]+'_logfile_'+str(num_trail)+'.txt','w')
    fid.write('Fitting p0 and F with logarithmic and polynomial derivative method\n')
    fid.write("Sand:"+fnam_sand[j0]+"\n\n")
    #*** Choosing order of polynom log file
    fid1 = open(fnam_sand[j0]+'_polyorder_logfile_'+str(num_trail)+'.txt','w')
    fid1.write("Sand:"+fnam_sand[j0]+"\n\n")
    #*************************************************************************#
    ################################ FIGURES ##################################
    # Fig 0: polynome et erreur
    fig=plt.figure(num=0+j0*nb_fig)
    gs = gridspec.GridSpec(1,3)
    ax01 = fig.add_subplot(gs[:2])
    ax02 = fig.add_subplot(gs[2])
    ax01.set_title(r"Polynomial fit of $f-f^{*}$ by $P(x)=\sum_{k=0}^{d} a_{k}x^{k}$")
    ax02.set_title("Least squared error\n"+r"between $f-f^{*}$ and $P$ ")
    ax02.set_ylabel('Error [%]')
    for axes in [ax01,ax02]:
        axes.set_xlabel("Pressure $p$ [kPa]")
        axes.set_xlim([50,550])
    # Fig 1:
    fig1=plt.figure(num=1+j0*nb_fig)
    ax1=fig1.add_subplot(111)
    ax1.set_title(r"Logarithm of the derivative of $f-f^{*}$ : $ln(P')$"+'\n'+r'Linear regression giving $p_0$ and $F$')
    ax1.set_ylabel(r"$ln(P')$ / linear regression ")
    # Fig 2:    
    fig2=plt.figure(num=2+j0*nb_fig)
    ax2=fig2.add_subplot(111)
    ax2.set_title(r"Exponential model : $f-f^{*} = (F-f^{*})(1-\exp(-\frac{p-p^*}{p_0}))$")
    
    for axes in [ax01,ax02,ax1,ax2]:
        axes.set_xlabel("Pression $p$ [kPa]",fontsize=14.0)
        axes.set_color_cycle(couleur)
        axes.grid()
#        axes.set_axis_bgcolor('lightgray')
    
    # Fig 3:    
    fig3=plt.figure(num=3+j0*nb_fig)
    gs3 = gridspec.GridSpec(2,1)
    ax30 = fig3.add_subplot(gs3[0])
    ax31 = fig3.add_subplot(gs3[1])
    for axes in [ax31,ax30]:
        axes.set_xlabel("$x_R$ [-]",fontsize=14.0)
        axes.set_xlim([0.0,0.55])
        axes.set_xticks([0.0,0.10,0.20,0.30,0.40,0.50])
        axes.set_xticklabels(axes.get_xticks(),size=13.0)
#        axes.grid()
    ax31.set_ylabel("$p_0$ [kPa]",fontsize=14.0)
    ax31.set_ylim([200,500])
    ax30.set_ylabel("$F$ [-]",fontsize=14.0)
    ax30.set_ylim([0.0,0.7])

    #Fig 4: 
    fig4,axarr=plt.subplots(3,1,num=4+j0*nb_fig,sharex=True)
#    fig4.set_size_inches(8.,12.,forward=True)
    axarr[-1].set_xlabel(r'Pressure $p$ [kPa]',fontsize=14.0)
    axarr[1].set_ylabel('Void ratios $e$ [-]',fontsize=14.0)
    for ax in axarr:
        ax.set_color_cycle(couleur)

    ### PURE SAND DATA ###
    
    datap=np.loadtxt(fnam_pure[j0])
    p0p = datap[0,0]
    p=datap[1:,0]
    e_void=datap[1:,1]
    
    erreur_p=np.zeros((len(FRi),pnew_length))
    erreur_m=np.zeros((len(FRi),pnew_length))
    pressure=np.zeros((len(FRi),pnew_length))
    iax4=-1
    for i0 in range(len(FRi)):#[3]:#
        #*********************************************************************#
        ####################### MIXTURE EXPERIMENTAL DATA #####################        
        FR=FRi[i0]
        ksi=FR/(1.0-FR)
        print 'Sand',fnam_sand[j0],'Rubber percentage',FR
        fid.write("\nrubber percentage: "+str(FR*100)+"%"+"\n")
        fid1.write("\nrubber percentage: "+str(FR*100)+"%"+"\n")
        fnam_mixture = ["NWH100R"+str(int(FR*100))+"C_new_dataexp_ep.txt","NWLB100R"+str(int(FR*100))+"C_new_dataexp_ep.txt"]
        datam=np.loadtxt(fnam_mixture[j0])
        p0m = datam[0,0]
        p_m=datam[1:,0]
        e_void_m=datam[1:,1]
#       print "length pm",len(p_m)
#       print "void 0",(e_void[0]-e_void_m[0])
        pmax=max(p[0],p_m[0])
        pmin=min(p_m[len(p_m)-1],p[len(p)-1])
        pnew,pnew_step=np.linspace(pmax,pmin,pnew_length,retstep=True)#min(len(p),len(p_m)))
#        print "mixture",len(p_m),"pure",len(p)
#        if FR==.10 and j0==0:
#            e_void_m=e_void_m-0.003#0.005#+(e_void[0]-e_void_m[0])   
        f=np.zeros(len(pnew))
        df=np.zeros(len(pnew))
        e_void_m_new=np.zeros(len(pnew))
        e_void_new=np.zeros(len(pnew))
        for k in range(len(pnew)):
            idx=np.argmin((p_m-pnew[k])**2)
            jdx=np.argmin((p-pnew[k])**2)
            e_void_m_new[k]=e_void_m[idx]
            e_void_new[k]=e_void[jdx]
            f[k]=-(e_void_m[idx]-e_void[jdx])/FR
        for j in range(len(pnew)-1):
            df[j]=(f[j+1]-f[j])/(pnew[j+1]-pnew[j])
        n_cut=np.argmin((pnew-[p_cut_HS,p_cut_LB][j0])**2)
        n_etoile=n_cut    
        p_etoile=pnew[n_etoile]
        f_etoile=f[n_etoile]
#        print "f_etoile",f_etoile
#        print "p_etoile",p_etoile
        
        #*********************************************************************#
        ########################## SPLINE FIT #################################
        pnew_cut=pnew[n_cut:]
        f_data=f[n_cut:]-f_etoile # curve to fit
        
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
                
        #*********************************************************************#
        ############## EXPONENTIAL MODEL WITH FITTED PARAMETERS ###############
        print 'Exponential model'
        if j0==0: #LB sand
            p0=-1/a1[0]
            F=f_etoile+p0*np.exp(a1[1]-(p_etoile/p0))
            p0_matrix[i0,j0]=p0
            F_matrix[i0,j0]=F
        else: #HS sand
            p0 = p0_matrix[i0,0]
            F = F_matrix[i0,0]
        f_exp=exp_etendu(x=pnew_cut,p0=p0,F=F,p_etoile=p_etoile,f_etoile=f_etoile)    
        #**** write to log file
        fid.write('\n'+r"Exponential model : $f-f^{*} = (F-f^{*})(1-\exp(-\frac{p-p^*}{p_0}))$"+'\n')
        fid.write( "p0=-1/a="+'{:2.2f}'.format(p0)+"\n")
        fid.write( "F="+'{:2.2f}'.format(F)+"\n")
        
        #*** Compute error
        print 'Computing relative error between exponential model and experimental data'
        epf_model=np.zeros(len(pnew))
        emf_model=np.zeros(len(pnew))
        for ii in range(len(pnew)):#range(n_cut,len(pnew)):#
            f_model=-exp_etendu(pnew[ii],p0,F,p_etoile,f_etoile)
			#*** pure from mixture:
            idx=np.argmin((p_m-pnew[ii])**2)
            em1=e_void_m[idx]
				#PATRICK
            epf_model[ii]=(em1-f_model*FR)
				#BOGDAN:
            # epf_model[ii]=(e1*(1+ksi)+f_model*ksi)/((1-f_model)*ksi+1)
			#*** mixture from pure:
            idx=np.argmin((p-pnew[ii])**2)
            e1=e_void[idx]
				#PATRICK:
            emf_model[ii]=(e1+FR*f_model)
				#BOGDAN:
            # emf_model[ii]=(e1*((1-f_model)*ksi+1)-f_model*ksi)/(1+ksi)
			#*** error
            erreur_m[i0,ii]=100*abs(emf_model[ii]-em1)/em1
            erreur_p[i0,ii]=100*abs(epf_model[ii]-e1)/e1
#            erreur_m[i2,ii]=abs(emf_model[ii]-em1)
#            erreur_p[i2,ii]=abs(epf_model[ii]-e1)
            pressure[i0,ii]=pnew[ii]
        erreur_m[i0,:n_cut-1]=0.
        erreur_p[i0,:n_cut-1]=0.0
        #*********************************************************************#
        ############################# PLOTS ###################################
        print 'Plotting...'
        ### Fig0:       
        ax01.plot(pnew_cut,f[n_cut:]-f_etoile,'ok',ms=1.5,mfc="w")        
        poly,=ax01.plot(pnew_cut,f_spline,'-',lw=1.5)
        txt=r'$d={:d}$'.format(d)        
        ax01.text(pnew_cut[-1],f_spline[-1],txt,va='top',ha='left',bbox=dict(fc='white',edgecolor=poly.get_color()))
        ax02.plot(pnew_cut,error_spline,'-',lw=1.5)
#        f_spline_model.plot_residuals(ax=ax02)
#        ### Fig1: logarithm derivation of f* + linear regression giving p0 and F    
        ax1.plot(pnew_log,ln_f_spline_prime,'--',lw=1.5)#,'o',ms=1.0)
        ax1.plot(pnew_log,ln_f_spline_prime_lin,'-k',lw=1.5)
#        ### Fig2:exponential model vs expe data for f-f*
        ax2.plot(pnew,f-f_etoile,'ok',mfc="w",ms=1.5,markevery=40)
        ax2.plot(pnew[n_cut:],f_exp-f_etoile)
        ### Fig4: exponential model vs expe data for void ratio
        if i0==0 or i0==2 or i0==4:
            iax4+=1
            ax4=axarr[iax4]
            PS,=ax4.plot(p,e_void,'o',mfc='white',mec='k',ms=9.0,markevery=40,alpha=0.7,label=r'$e(x_{R}=0,p)$')
            M,=ax4.plot(p_m,e_void_m,'^',mfc='w',mec='k',ms=9.0,markevery=40,alpha=0.7,label=r'$e(x_{R},p)$')
            
    #        n_cut=0
            pnew_cut=pnew[n_cut:]
            PS0,=ax4.plot(pnew_cut,epf_model[n_cut:],'--k',linewidth=1.5,label='$\\tilde{e}(x_{R},p)$')
            M0, =ax4.plot(pnew_cut,emf_model[n_cut:],'-k',linewidth=1.5,label='$\\tilde{e}(x_{R},p)$')
            x=[p_etoile,p_etoile]
            y=[ax4.get_ylim()[0],ax4.get_ylim()[1]]
            ax4.plot(x,y,':k',lw=1.0)
            if i0==0:
                ax4.text(x[0],y[1],r'$p^*$',ha='left',va='top',fontsize=14)
            ax4.text(0.99,0.98,letter.next(),ha='right',va='top',fontsize=18.0,transform=ax4.transAxes)
            ax4.text(0.01,0.1,r'$x_{R}='+str(int(FR*100))+'\%$',ha='left',va='center',fontsize=18.0,transform=ax4.transAxes)
    #        ax4.set_yticks(ax4.get_yticks()[1:-1:2])
    #        ax4.set_yticklabels(['{:.2f}'.format(ytick) for ytick in ax4.get_yticks()],size=16.0)
    #        ax4.set_xticklabels(['{:d}'.format(int(xtick)) for xtick in ax4.get_xticks()],size=16.0)
            majorLocator   = MultipleLocator([0.01,0.02,0.05,0.05,0.1][i0])
            majorFormatter = FormatStrFormatter('%.2f')
            ax4.yaxis.set_major_locator(majorLocator)
            ax4.yaxis.set_major_formatter(majorFormatter)
            ax4.set_yticklabels(ax4.get_yticks(),size=13.0)
            ax4.set_xticks([0,100,200,300,400,500,600])
            ax4.set_xticklabels(ax4.get_xticks(),fontsize=13.0)
#    axarr[0].legend(handles=[PS,PS0],loc="upper center",handletextpad=0.2,
#        ncol=2, mode="expand", borderaxespad=0.5,fontsize=13.0)
#    axarr[0].legend(handles=[PS,M,PS0],loc="upper right",handletextpad=0.2,borderaxespad=0.5,fontsize=11.0)
    fig4.tight_layout()
    ### Fig3: p0(xR) and F(xR)
    FRa==np.array(FRi)
    ax30.plot(FRa,F_matrix[:,j0],':sk')
    ax30.text(0.5,0.9,r'(a)',ha='center',va='center',fontsize=14.0,transform=ax30.transAxes)
    majorLocator   = FixedLocator([0.,0.2,0.4,0.6])
    majorFormatter = FormatStrFormatter('%.1f')
    minorLocator   = FixedLocator([0.1,0.3,0.5,0.7])
    ax30.yaxis.set_major_locator(majorLocator)
    ax30.yaxis.set_major_formatter(majorFormatter)
    ax30.yaxis.set_minor_locator(minorLocator)
    ax30.set_yticklabels(ax30.get_yticks(),size=13.0)
#    ax30.set_yticks([0,0.2,0.4,0.6])
#    ax30.set_yticklabels(['{:.1f}'.format(ytick) for ytick in ax30.get_yticks()],size=16.0)
#    ax30.set_yticks([0,0.2,0.4,0.6],minor=True)
    ax31.plot(FRa,p0_matrix[:,j0],':ok')
    ax31.text(0.5,0.9,r'(b)',ha='center',va='center',fontsize=14.0,transform=ax31.transAxes)
#    ax31.set_yticks([200,300,400,500],minor=True)
#    ax31.set_yticklabels(['{:d}'.format(int(ytick)) for ytick in ax31.get_yticks()],minor=False,size=16.0)
    majorLocator   = FixedLocator([200,300,400,500])
    majorFormatter = FormatStrFormatter('%d')
    minorLocator   = FixedLocator([250,350,450])
    ax31.yaxis.set_major_locator(majorLocator)
    ax31.yaxis.set_major_formatter(majorFormatter)
    ax31.yaxis.set_minor_locator(minorLocator)
    ax31.set_yticklabels(ax31.get_yticks(),size=13.0)
#    gs.tight_layout(fig)
    
    ax02.yaxis.set_label_position('right')
    ax02.yaxis.set_ticks_position('right')
    gs3.tight_layout(fig3)
    
    for figs in [fig3,fig4]:#fig,fig1,fig2,
#        tight_layout()
        figs.savefig(".\\figures\\"+fnam_sand[j0]+"_Fig"+str(figs.number-j0*nb_fig)+".png",format='png',bbox_inches='tight',dpi=300.0)
        figs.savefig(".\\figures\\"+fnam_sand[j0]+"_Fig"+str(figs.number-j0*nb_fig)+".eps",format='eps',bbox_inches='tight',dpi=300.0)

#        plt.close(figs)
    fid.write('Evolution of p0 with xR:\n'+str(p0_matrix)+'\n')
    fid.write('Evolution of F with xR:\n'+str(F_matrix)+'\n')
        
    fid.close()
    fid1.close()
    
#*****************************************************************************#
################### FIGURE 2 : ERROR MAP ######################################
    erreur_m_list.append(erreur_m)
    erreur_p_list.append(erreur_p)
    pressure_list.append(pressure)

#fig, axes = plt.subplots(nrows=2, ncols=2)
#for dat, ax in zip(data, axes.flat):
#    # The vmin and vmax arguments specify the color limits
#    im = ax.imshow(dat, vmin=0, vmax=2)
#
## Make an axis for the colorbar on the right side
#cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
#fig.colorbar(im, cax=cax)
#%%
y,x=np.meshgrid(np.ones(pnew_length),np.array(FRi))

fig2,ax2=plt.subplots(1,2,num=2*nb_fig+1,sharey=True)
cs1=ax2[0].contourf(x,pressure,erreur_m_list[0],cmap=plt.cm.YlOrBr,vmax=2.4,vmin=0.0)
cs2=ax2[1].contourf(x,pressure,erreur_m_list[1],cmap=plt.cm.YlOrBr,vmax=2.4,vmin=0.0)
cmin=min([cs1.get_clim()[0],cs2.get_clim()[0]])
print cmin
cmax=max([cs1.get_clim()[1],cs2.get_clim()[1]])
print cmax
cs1.set_clim(cmin,cmax)
cs2.set_clim(cmin,cmax)
cb=fig2.colorbar(cs2,ax=ax2.ravel().tolist(),orientation="horizontal")    

cb.set_clim(cmin,cmax)
#cb_ticks=np.linspace(cmin,cmax,10)
#cb.set_ticks(cb_ticks)
cb.set_label("error $\\mathcal{E}$ [%]",fontsize=14.0)
for ax in ax2:
    ax.set_xticks([0.0,0.1,0.2,0.3,0.4,0.5])
    ax.set_xticklabels(["{:d}".format(int(xtick*100)) for xtick in ax.get_xticks()],size=13.0)
    ax.set_yticks([0,100,200,300,400,500])
    ax.set_xlabel("Volume ratio of rubber $x_R$ [-]",ha='center',fontsize=14.0)
ax2[0].set_yticklabels(["{:d}".format(ytick) for ytick in ax2[0].get_yticks()],size=13.0)
ax2[0].set_title('(a): Hostun Sand',fontsize=14.0)
ax2[1].set_title('(b): Leighton Buzzard Sand',fontsize=14.0)
ax2[0].set_ylabel("Pressure $p$ [kPa]",fontsize=14.0)
fig2.savefig('.\\figures\\error_map.png',format='png',bbox_inches='tight',dpi=250.0)
fig2.savefig('.\\figures\\error_map.eps',format='eps',bbox_inches='tight',dpi=500.0) 