# -*- coding: utf-8 -*-
"""
Created on Thu Feb 04 15:08:15 2016

@author: platzer
"""

import numpy as np
import matplotlib.pyplot as plt
import itertools
plt.close('all')


FRi = [0.05,0.10,0.20,0.30,0.50,1.0]#[0.05]#
fnam_pure = ["NWH100R0C_new_dataexp_ep.txt","NWLB100R0C_new_dataexp_ep.txt"]
savenam=['HS','LB'] # order of plot: HS first (top (a)) and LB next (bottom (b))

#==============================================================================
# Plot settings
#==============================================================================
width,height = plt.figaspect(1)
fig2,axarr= plt.subplots(2,1,sharex=True,figsize=(1.3*width,1.5*height))
colormap = plt.cm.Greys
couleur = [colormap(i) for i in np.linspace(0.2,1.,len(FRi))]
marker1 = itertools.cycle(('^', 'o','v', 's','<','>'))
xtext=np.zeros(2)
ytext=np.zeros(2)

#==============================================================================
# Load and plot experimental results
#==============================================================================
for j in [0,1]:#range(1): #deux sables

    ### LOAD PURE SAND DATA ###
    datap=np.loadtxt(fnam_pure[j])
    p=datap[1:,0]
    e_void=datap[1:,1]
    xtext[j]=p[-1]
    ytext[j]=e_void[-1]

    ### Legend and labels ###
    ep_legend=[]
    labels=[]
    axarr[j].grid()

    ### PLOT PURE SAND DATA ###
    PS,=axarr[j].plot(p,e_void,ls='',marker='8',mfc='w',mec='k',ms=7.0,markevery=50,label=r'$x_{R}=0%$')
    ep_legend.append(PS)
    axarr[j].set_color_cycle(couleur)

    for i in range(len(FRi)):#range(1):
        FR=FRi[i]
        fnam_mixture = ["NWH100R"+str(int(FR*100))+"C_new_dataexp_ep.txt","NWLB100R"+str(int(FR*100))+"C_new_dataexp_ep.txt"]

        ### LOAD MIXTURE ###
        datam=np.loadtxt(fnam_mixture[j])
        p0m = datam[0,0]
        p_m=datam[1:,0]
        e_void_m=datam[1:,1]
        ### correction for HS sand at 10% (not made in the computation, can remove)
        if FR==.10 and j==0:
            print "true"
            print (e_void[0]-e_void_m[0])
            e_void_m=e_void_m-0.005#+(e_void[0]-e_void_m[0])

        ### PLOT MIXTURE ###
        M,=axarr[j].plot(p_m,e_void_m,ls='',mec='k',marker=marker1.next(),ms=7.0,markevery=40,alpha=1.0,label=r'$x_{R}= '+str(int(FR*100))+'\%$')
        ep_legend.append(M)
        labels.append(r'$x_{R}= '+str(int(FR*100))+'\%$')

ax1=axarr[0]
ax2=axarr[1]
#ax2.set_ylim(bottom=0.0)
## add (a) and (b) on the subplots:
ax1.text(0.05,0.1,r'(a)',ha='left',va='center',fontsize=18.0,transform=ax1.transAxes)
ax2.text(0.05,0.1,r'(b)',ha='left',va='center',fontsize=18.0,transform=ax2.transAxes)
## add an arrow showing decrease of rubber fraction:
ax1.text(0.91,0.5,r'$x_{R}$',ha='left',va='center',fontsize=14.0,transform=ax1.transAxes)
ax1.annotate('', xy=(0.9,0.9), xycoords='axes fraction',
                xytext=(0.9,0.1), textcoords='axes fraction',size=16.0,
#                bbox=dict(boxstyle="square",facecolor="w",edgecolor='k'),
                arrowprops=dict(arrowstyle="<-",lw=1.0,color='k',relpos=(0.0, 0.0))
                )
ax2.text(0.96,0.6,r'$x_{R}$',ha='left',va='center',fontsize=14.0,transform=ax2.transAxes)
ax2.annotate('', xy=(0.95,0.9), xycoords='axes fraction',
                xytext=(0.95,0.1), textcoords='axes fraction',size=16.0,
#                bbox=dict(boxstyle="square",facecolor="w",edgecolor='k'),
                arrowprops=dict(arrowstyle="<-",lw=1.0,color='k',relpos=(0.0, 0.0))
                )
## configure axis (lim and labels)
ax1.set_ylabel(r'void ratio $e$ [-]',fontsize=14.0)
ax2.set_ylabel(r'void ratio $e$ [-]',fontsize=14.0)
ax2.set_xlabel(r'pressure $p$ [kPa]',fontsize=14.0)
ax2.set_xlim(0,600)
ax1.set_yticks([0.4,0.5,0.6,0.7,0.8])
ax1.set_yticklabels(ax1.get_yticks(),size=12.0)
ax2.set_yticks([0.1,0.2,0.3,0.4,0.5,0.6])
ax2.set_yticklabels(ax2.get_yticks(),size=12.0)
#leg=ax2.legend(handles=ep_legend,loc=3,bbox_to_anchor=(0.01, 0.02, 0.98, .102),\
#           ncol=4, mode="expand", borderaxespad=0.,fontsize=12.0)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
fig2.tight_layout(h_pad = 0.05)
#fig2.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.05)
fig2.savefig('.\\figures\\plot_exp.png',format='png',dpi=500,bbox='tight')
fig2.savefig('.\\figures\\plot_exp.eps',format='eps',dpi=500,bbox='tight')
plt.show()