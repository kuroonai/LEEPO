#!/usr/bin/python
# -*- coding: utf-8 -*-
'''

@author:Naveen Kumar Vasudevan, 
        400107764,
        Doctoral Student, 
        The Xi Research Group, 
        Department of Chemical Engineering,
        McMaster University, 
        Hamilton, 
        Canada.
        
        naveenovan@gmail.com
        https://naveenovan.wixsite.com/kuroonai

This program was made based on the following research work: 

Xu, Xiaoming, Mansoor A. Khan, and Diane J. Burgess. "Predicting hydrophilic 
drug encapsulation inside unilamellar liposomes." 
International journal of pharmaceutics 423.2 (2012): 410-418.

DOI : https://doi.org/10.1016/j.ijpharm.2011.12.019



# Python code to calculate Encapuslation efficiecny of Liposomes

(L)iposome

    (E)ncapsulation
    
        (E)fficiency
        
            (P)redictor &
            
                (O)ptimizer
            
# Program Usage : python LEEPO-sphere.py mu sigma, sizerange ... 

# Example command : python LEEPO.py 190.0 38.0 1000 0.23 5.2 0.45 140 5 60.04

Enter the system parameters in the order of:

1. Mean                         (μ - nm {float})
2. Standard deviation           (σ - nm {float})
3. Size range                   (eg: for 1-1000 nm enter "1000")
4. Anticipated defect           (% {float})
5. Bilayer thickness           (d - nm)
6. Molecular area of lipid      (a - nm^2)
7. Lipid molar concentration    (c - mMol/L)
8. Total volume                 (V - ml)
9. Expected Encapsulation efficiency (%) [optional] - provide only if defect value optimization is required.

'''

import sys
import os
import matplotlib.pyplot as plt; plt.rcdefaults()
import math
from scipy import constants as cons
import time




#descrip='Enter the system parameters in the order of:\n\n1. Mean (mu - nm {float})\n2. Standard deviation (sigma - nm {float})\n3. Size range (eg: for 1-1000 nm enter "1000")\n4. Anticipated defect (% {float})\n5. membrane thickness (d - nm)\n6. Molecular area of lipid (a - nm^2)\n7. Lipid molar concentration (c - mMol/L)\n8. Total volume (V - ml)\n9. Expected Encapsulation efficiency (%) [optional]'
#parser = argparse.ArgumentParser(description=descrip, formatter_class=argparse.RawTextHelpFormatter)
#args = parser.parse_args()
timer=time.time()
if len(sys.argv) < 9 or len(sys.argv) > 10:
  print('\n\nSyntax error encountered:\n\nUsage:\n\nEnter the system parameters in the order of:\n\n1. Mean (μ - nm)\n2. Standard deviation (σ - nm)\n3. Size range (eg: for 1-1000 nm enter "1000")\n4. Anticipated defect (%)\n5. Bilayer thickness (d - nm)\n6. Molecular area of lipid (a - nm^2)\n7. Lipid molar concentration (c - mMol/L)\n8. Total volume (V - ml)\n9. Expected Encapsulation efficiency (%) [optional]\n\n')
  raise Exception("Syntax: *.py mu sd ...")

cwd = os.getcwd()

if os.name=='nt':folder='%s\\LEEPO-plots'%cwd
elif os.name=='posix':folder='%s/LEEPO-plots'%cwd


if os.path.isdir(folder):
    os.chdir(folder)
    
else:
    os.mkdir('LEEPO-plots')
    os.chdir(folder)

mu, sigma, sizerange, defect, d, molarea, c, V  = float(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),float(sys.argv[6]),float(sys.argv[7]),float(sys.argv[8])

if len(sys.argv) == 10:expEE = float(sys.argv[9])

print("\nInput values are as follows:\n")    
print("Mean (μ) - %.4f nm\n"%mu)
print("\tStandard deviation (σ) - %.4f nm\n" %sigma)
print("\t\tSize range - %d nm\n" %sizerange)
print("\t\t\tDefect  - %.4f %%\n" %defect)
print("\t\t\t\tLipid bilayer thickness (d) - %.4f nm\n" %d)
print("\t\t\t\t\tLipid molecular area (a) -  %.4f nm^2\n" %molarea)
print("\t\t\t\t\t\tLiposome molar concentration (c)- %.4f mMol/L\n" %c)
print("\t\t\t\t\t\t\tTotal volume (V) - %.4f ml\n" %V)
if len(sys.argv) == 10:print("\t\t\t\t\t\t\t\tExpected/Experimental EE- %.4f %%\n" %expEE)

EE=0
EE_back=0
iteration=0
max_defect="yes"

def plotEE(ev, ee):
    labels = 'Encapsulated', 'Free'
    sizes = [sum(ev), V-sum(ev)]
    colors = ['Green', 'Red']
    explode = (0.0, 0.1)   # explode 1st slice
 
    # Plot
    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
        autopct='%1.2f%%', shadow=True, startangle=100)
 
    plt.axis('equal')
    #plt.show()
    if abs(ee-expEE) >= abs(float(0.01)):
        plt.title("LEEPO-Predicted EE")
        plt.savefig("LEEPO-predicted.pdf", bbox_inches='tight')
        plt.close()
    elif abs(ee-expEE) <= abs(float(0.01)):
        plt.title("LEEPO-Optimized EE")
        plt.savefig("LEEPO-Optimised.pdf", bbox_inches='tight')
        plt.close()
        
    return()

def plotEE_old(ev, ee):
    labels = 'Encapsulated', 'Free'
    sizes = [sum(ev), V-sum(ev)]
    colors = ['Green', 'Red']
    explode = (0.0, 0.1)   # explode 1st slice
 
    # Plot
    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
        autopct='%1.2f%%', shadow=True, startangle=100)
 
    plt.axis('equal')
    #plt.show()
    plt.title("Old Model's - Predicted EE")
    plt.savefig("old model-predicted.pdf", bbox_inches='tight')
    plt.close()
            
    return()
    
def complotter(p1,p2,lst,lst1,colour,colour1,ylab,title,name):
    print(("\n\tPlotting %s")%title)
    plt.figure()
    plt.plot(p1,linestyle=lst, color=colour,lw=1.0, label='Predicted')
    plt.plot(p2,linestyle=lst1, color=colour1,lw=1.0, label='Optimized')
    plt.title(title)
    plt.xlabel('Vesticle size (nm)')
    plt.ylabel(ylab)
    plt.legend(loc='best')
    plt.savefig(name, bbox_inches='tight')
    plt.close('all')

def compareplot(mu, sigma, sizerange, defect,defect1, d, molarea, c, V):
    sigbymu=sigma/mu
    x=list(range(1,sizerange+1)) #particle side/dia in (nm)
    p=[0]*sizerange #probability of given size
    
    dx=[0]*sizerange
    v_in=[0]*sizerange #internal volume micro - L
    area_out=[0]*sizerange # outer surface area nm2
    area_in=[0]*sizerange # inner surface area nm2
    lnum_out=[0]*sizerange #lipid molecule number out
    lnum_in=[0]*sizerange #lipid molecule number in
    k=[0]*sizerange #lipid number per incident
    vesiclenumber=[0]*sizerange
    entrapvolume=[0]*sizerange
    
    dx1=[0]*sizerange
    v_in1=[0]*sizerange #internal volume micro - L
    area_out1=[0]*sizerange # outer surface area nm2
    area_in1=[0]*sizerange # inner surface area nm2
    lnum_out1=[0]*sizerange #lipid molecule number out
    lnum_in1=[0]*sizerange #lipid molecule number in
    k1=[0]*sizerange #lipid number per incident
    vesiclenumber1=[0]*sizerange
    entrapvolume1=[0]*sizerange
    
    print("\nplotting the predicted and optimised outputs for comparison\n")
    
    for i in range(0,sizerange):
        #probability with no unit
        p[i]=(1/math.sqrt(2 * math.pi * (sigbymu)**2 * x[i]**2) ) * math.exp((math.log(x[i])-math.log(mu))**2/(-2 * (sigbymu)**2))
        #defect in terms of the probability
      
    for i in range(0,sizerange):    
        dx[i] = defect*(p[i]/max(p))
        dx1[i] = defect1*(p[i]/max(p))
    for i in range(0,sizerange):
    
        #internal volume in (micro L)
        v_in[i]=( (4/3) * math.pi * (float(x[i])/2 - d)**3 * (1-dx[i]) + (4/3) * math.pi * (float(x[i])/2 - d/2)**3 * dx[i] ) * 10**-18
        v_in1[i]=( (4/3) * math.pi * (float(x[i])/2 - d)**3 * (1-dx1[i]) + (4/3) * math.pi * (float(x[i])/2 - d/2)**3 * dx1[i] ) * 10**-18
        #outer surface area
        area_out[i]= 4 * math.pi * (float(x[i])/2)**2
        area_out1[i]= 4 * math.pi * (float(x[i])/2)**2

        
        #inner surface area        
        if x[i] < d:
            area_in[i]=0
            area_in1[i]=0
               
        else:
            area_in[i]= 4 * math.pi * (float(x[i])/2 - d)**2 * (1-dx[i]) + 4 * math.pi * (float(x[i])/2 - d/2)**2 * dx[i]
            area_in1[i]= 4 * math.pi * (float(x[i])/2 - d)**2 * (1-dx1[i]) + 4 * math.pi * (float(x[i])/2 - d/2)**2 * dx1[i]

        lnum_out[i]=area_out[i]/molarea
        lnum_in[i]=area_in[i]/molarea
        
        lnum_out1[i]=area_out1[i]/molarea
        lnum_in1[i]=area_in1[i]/molarea
    
        k[i]=p[i]*(lnum_out[i]+lnum_in[i])
        k1[i]=p[i]*(lnum_out1[i]+lnum_in1[i])

    m=c*V*cons.N_A/10**6/sum(k)
    m1=c*V*cons.N_A/10**6/sum(k1)

    for i in range(0,sizerange):
        vesiclenumber[i]=m*p[i]
        entrapvolume[i]=vesiclenumber[i]*v_in[i]/1000
        vesiclenumber1[i]=m1*p[i]
        entrapvolume1[i]=vesiclenumber1[i]*v_in1[i]/1000
    
    
    complotter(dx,dx1,'-','-','r','g','Defect %','Defect Percentage','LEEPO-defect.pdf')
    complotter(v_in,v_in1,'-','-','r','g',r'Internal volume ($\mu$l)','Internal volume','LEEPO-V_internal.pdf')
    complotter(area_out,area_out1,'-','-','r','g','a_out (nm^2)','Outer Surface area','LEEPO-area_outer.pdf')
    complotter(area_in,area_in1,'-','-','r','g','a_in (nm^2)','Inner Surface area','LEEPO-area_inner.pdf')
    complotter(lnum_out,lnum_out1,'-','-','r','g','# units','Outer Lipid units','LEEPO-lin.pdf')
    complotter(lnum_in,lnum_in1,'-','-','r','g','# units','Inner Lipid units','LEEPO-lout.pdf')
    complotter(k,k1,'-','-','r','g','# units','Lipid units','LEEPO-lipid units.pdf')
    complotter(vesiclenumber,vesiclenumber1,'-','-','r','g','# units','Vesicle number','LEEPO-Vesiclenumber.pdf')
    complotter(entrapvolume,entrapvolume1,'-','-','r','g','Entraped volume (ml)','Entraped volume','LEEPO-entrapvol.pdf')
    print(("\nPlots completed and can be accessed at '%s'\n")%folder)
    
def max_defect (mu, sigma, sizerange, defect, d, molarea, c, V, EE):
    
    sigbymu=sigma/mu
    x=list(range(1,sizerange+1)) #particle side/dia in (nm)
    p=[0]*sizerange #probability of given size
    dx=[0]*sizerange
    v_in=[0]*sizerange #internal volume micro - L
    area_out=[0]*sizerange # outer surface area nm2
    area_in=[0]*sizerange # inner surface area nm2
    lnum_out=[0]*sizerange #lipid molecule number out
    lnum_in=[0]*sizerange
    
    k=[0]*sizerange 
    
    vesiclenumber=[0]*sizerange
    entrapvolume=[0]*sizerange
    
    EE_new=EE
    defect_old=defect
    iteration=0
    print('\nOptimal range determination step begins !!\n')
    while abs(EE_new-74.05) >= abs(float(0.01)):
        if EE_new < 74.05: defect=defect+(abs(EE_new-74.05)/10)
        elif EE_new > 74.05: defect=defect-(abs(EE_new-74.05)/10)
        #else: exit()
        for i in range(0,sizerange):
        #probability with no unit
            p[i]=(1/math.sqrt(2 * math.pi * (sigbymu)**2 * x[i]**2) ) * math.exp((math.log(x[i])-math.log(mu))**2/(-2 * (sigbymu)**2))
        #defect in terms of the probability
        for i in range(0,sizerange):    
            dx[i] = defect*(p[i]/max(p))
        for i in range(0,sizerange):
    
        #internal volume in (micro L)
            v_in[i]=( (4/3) * math.pi * (float(x[i])/2 - d)**3 * (1-dx[i]) + (4/3) * math.pi * (float(x[i])/2 - d/2)**3 * dx[i] ) * 10**-18
        #outer surface area
            area_out[i]= 4 * math.pi * (float(x[i])/2)**2
        
        #inner surface area        
            if x[i] < d:area_in[i]=0
    
            else:
                area_in[i]= 4 * math.pi * (float(x[i])/2 - d)**2 * (1-dx[i]) + 4 * math.pi * (float(x[i])/2 - d/2)**2 * dx[i]
    
            lnum_out[i]=area_out[i]/molarea
            lnum_in[i]=area_in[i]/molarea
    
            k[i]=p[i]*(lnum_out[i]+lnum_in[i])

        m=c*V*cons.N_A/10**6/sum(k)

        for i in range(0,sizerange):
            vesiclenumber[i]=m*p[i]
            entrapvolume[i]=vesiclenumber[i]*v_in[i]/1000

        EE_new=sum(entrapvolume)/V * 100
        print(("Iteration No: %d ##############################")%iteration)
        print(("\n\t\tThe encapsulation efficiency for the system at interation %d is %.4f %% , The defect is %.4f %%, Δ efficiency = %.4f %%\n")%(iteration,EE_new,abs(defect),abs(EE_new-expEE)))
        print(("\t\tThe entrapment volume is %.4f μl/μmol\n")%(sum(entrapvolume)/(c*V)*1000))

        iteration=iteration+1
        
        #if abs(EE_new-expEE) != abs(0.0001):exit 
    
    print(("\nThe optimal Encapsulation efficiency range we can expect before liposome degradation is %.4f - %.4f %%")%(EE,EE_new))
    print(("\n\nThe Optimal defect/void %% is between %.4f - %.4f %%\n")%(abs(defect_old),defect))
      
def encap_opt (mu, sigma, sizerange, defect, d, molarea, c, V, EE):
    
    sigbymu=sigma/mu
    x=list(range(1,sizerange+1)) #particle side/dia in (nm)
    p=[0]*sizerange #probability of given size
    dx=[0]*sizerange
    v_in=[0]*sizerange #internal volume micro - L
    area_out=[0]*sizerange # outer surface area nm2
    area_in=[0]*sizerange # inner surface area nm2
    lnum_out=[0]*sizerange #lipid molecule number out
    lnum_in=[0]*sizerange
    
    k=[0]*sizerange 
    
    vesiclenumber=[0]*sizerange
    entrapvolume=[0]*sizerange
    
    EE_new=EE
    iteration=0
    print('Optimization step begins !!\n')
    while abs(EE_new-expEE) >= abs(float(0.00001)):
        if EE_new < expEE: defect=defect+(abs(EE_new-expEE)/10)
        elif EE_new > expEE: defect=defect-(abs(EE_new-expEE)/10)
        #else: exit()
        for i in range(0,sizerange):
        #probability with no unit
            p[i]=(1/math.sqrt(2 * math.pi * (sigbymu)**2 * x[i]**2) ) * math.exp((math.log(x[i])-math.log(mu))**2/(-2 * (sigbymu)**2))
        #defect in terms of the probability
        for i in range(0,sizerange):    
            dx[i] = defect*(p[i]/max(p))
        for i in range(0,sizerange):
    
        #internal volume in (micro L)
            v_in[i]=( (4/3) * math.pi * (float(x[i])/2 - d)**3 * (1-dx[i]) + (4/3) * math.pi * (float(x[i])/2 - d/2)**3 * dx[i] ) * 10**-18
        #outer surface area
            area_out[i]= 4 * math.pi * (float(x[i])/2)**2
        
        #inner surface area        
            if x[i] < d:area_in[i]=0
    
            else:
                area_in[i]= 4 * math.pi * (float(x[i])/2 - d)**2 * (1-dx[i]) + 4 * math.pi * (float(x[i])/2 - d/2)**2 * dx[i]
    
            lnum_out[i]=area_out[i]/molarea
            lnum_in[i]=area_in[i]/molarea
    
            k[i]=p[i]*(lnum_out[i]+lnum_in[i])

        m=c*V*cons.N_A/10**6/sum(k)

        for i in range(0,sizerange):
            vesiclenumber[i]=m*p[i]
            entrapvolume[i]=vesiclenumber[i]*v_in[i]/1000

        EE_new=sum(entrapvolume)/V * 100
        print(("Iteration No: %d ##############################")%iteration)
        print(("\n\t\tThe encapsulation efficiency for the system at interation %d is %.4f %% , The defect is %.4f %%, Δ efficiency = %.4f %%\n")%(iteration,EE_new,defect,abs(EE_new-expEE)))
        print(("\t\tThe entrapment volume is %.4f μl/μmol\n")%(sum(entrapvolume)/(c*V)*1000))

        iteration=iteration+1
        
        #if abs(EE_new-expEE) != abs(0.0001):exit
    print('Optimization ends!!\n\n')
    print(('The optimized values are :\n\nExpected EE \t\t- %.4f %%\n\nEE predicted by model \t- %.4f %%\n\nEE Optimized by model \t- %.4f %%\n\nDefect Percentage \t- %.4f %%\n\nTotal Vesicles \t\t- %.4E\n\nTotal internal volume \t- %.4f ml\n\nEntrapment volume \t- %.4f μl/μmol\n\nNumber of optimization steps - %d steps\n\n')%(expEE,EE,EE_new,abs(defect),m,sum(entrapvolume),sum(entrapvolume)/(c*V)*1000,iteration-1))
    max_defect (mu, sigma, sizerange, defect, d, molarea, c, V, EE_new)
    
    print(('\n\nThe total execution time is  %s seconds\n\n')%(time.time()-timer))
    plotEE(entrapvolume, EE_new)
    
    compareplot(mu, sigma, sizerange, float(sys.argv[4]),defect, d, molarea, c, V)
    
def encap_old  (mu, sigma, sizerange, defect, d, molarea, c, V):
        
    sigbymu=sigma/mu

    x=list(range(1,sizerange+1)) #particle side/dia in (nm)

    p=[0]*sizerange #probability of given size
    dx=[0]*sizerange

    v_in=[0]*sizerange #internal volume micro - L

    area_out=[0]*sizerange # outer surface area nm2
    area_in=[0]*sizerange # inner surface area nm2
    
    lnum_out=[0]*sizerange #lipid molecule number out
    lnum_in=[0]*sizerange #lipid molecule number in
    
    k=[0]*sizerange #lipid number per incident 

    vesiclenumber=[0]*sizerange
    entrapvolume=[0]*sizerange
    
    

    for i in range(0,sizerange):
        #probability with no unit
        p[i]=(1/math.sqrt(2 * math.pi * (sigbymu)**2 * x[i]**2) ) * math.exp((math.log(x[i])-math.log(mu))**2/(-2 * (sigbymu)**2))
        #defect in terms of the probability
    for i in range(0,sizerange):    
        dx[i] = defect*(p[i]/max(p))
    for i in range(0,sizerange):
    
        #internal volume in (micro L)
        v_in[i]=( (4/3) * math.pi * (float(x[i])/2 - d)**3 * (1-dx[i]) + (4/3) * math.pi * (float(x[i])/2 - d/2)**3 * dx[i] ) * 10**-18
        #outer surface area
        area_out[i]= 4 * math.pi * (float(x[i])/2)**2
        
        #inner surface area        
        if x[i] < d:area_in[i]=0
    
        else:
            area_in[i]= 4 * math.pi * (float(x[i])/2 - d)**2 * (1-dx[i]) + 4 * math.pi * (float(x[i])/2 - d/2)**2 * dx[i]
    
        lnum_out[i]=area_out[i]/molarea
        lnum_in[i]=area_in[i]/molarea
    
        k[i]=p[i]*(lnum_out[i]+lnum_in[i])

    m=c*V*cons.N_A/10**6/sum(k)

    for i in range(0,sizerange):
        vesiclenumber[i]=m*p[i]
        entrapvolume[i]=vesiclenumber[i]*v_in[i]/1000

    EE=sum(entrapvolume)/V * 100
    
    print("\nW.R.T Existing model with no parameter for defect:\n\n")
    print(("\n\t\tThe encapsulation efficiency for the system is %.4f %%\n\n")%(EE))
    print(("\n\t\tThe entrapment volume is %.4f μl/μmol\n")%(sum(entrapvolume)/(c*V)*1000))
    print(('\n\t\tTotal Vesicles \t\t- %.4E\n\n\t\tTotal internal volume \t- %.4f ml\n\n\t\tEntrapment volume \t- %.4f μl/μmol\n\n')%(m,sum(entrapvolume),sum(entrapvolume)/(c*V)*1000))
    
      
    plotEE_old(entrapvolume, EE)


    
def encap (mu, sigma, sizerange, defect, d, molarea, c, V):
    
    sigbymu=sigma/mu

    x=list(range(1,sizerange+1)) #particle side/dia in (nm)

    p=[0]*sizerange #probability of given size
    dx=[0]*sizerange

    v_in=[0]*sizerange #internal volume micro - L

    area_out=[0]*sizerange # outer surface area nm2
    area_in=[0]*sizerange # inner surface area nm2
    
    lnum_out=[0]*sizerange #lipid molecule number out
    lnum_in=[0]*sizerange #lipid molecule number in
    
    k=[0]*sizerange #lipid number per incident 

    vesiclenumber=[0]*sizerange
    entrapvolume=[0]*sizerange
    
    

    for i in range(0,sizerange):
        #probability with no unit
        p[i]=(1/math.sqrt(2 * math.pi * (sigbymu)**2 * x[i]**2) ) * math.exp((math.log(x[i])-math.log(mu))**2/(-2 * (sigbymu)**2))
        #defect in terms of the probability
    for i in range(0,sizerange):    
        dx[i] = defect*(p[i]/max(p))
    for i in range(0,sizerange):
    
        #internal volume in (micro L)
        v_in[i]=( (4/3) * math.pi * (float(x[i])/2 - d)**3 * (1-dx[i]) + (4/3) * math.pi * (float(x[i])/2 - d/2)**3 * dx[i] ) * 10**-18
        #outer surface area
        area_out[i]= 4 * math.pi * (float(x[i])/2)**2
        
        #inner surface area        
        if x[i] < d:area_in[i]=0
    
        else:
            area_in[i]= 4 * math.pi * (float(x[i])/2 - d)**2 * (1-dx[i]) + 4 * math.pi * (float(x[i])/2 - d/2)**2 * dx[i]
    
        lnum_out[i]=area_out[i]/molarea
        lnum_in[i]=area_in[i]/molarea
    
        k[i]=p[i]*(lnum_out[i]+lnum_in[i])

    m=c*V*cons.N_A/10**6/sum(k)

    for i in range(0,sizerange):
        vesiclenumber[i]=m*p[i]
        entrapvolume[i]=vesiclenumber[i]*v_in[i]/1000

    EE=sum(entrapvolume)/V * 100
    
    print(("\n\t\tThe encapsulation efficiency for the system is %.4f %%\n\n")%(EE))
    print(("\n\t\tThe entrapment volume is %.4f μl/μmol\n")%(sum(entrapvolume)/(c*V)*1000))
    print(('\n\t\tDefect Percentage \t- %.4f%%\n\n\t\tTotal Vesicles \t\t- %.4E\n\n\t\tTotal internal volume \t- %.4f ml\n\n\t\tEntrapment volume \t- %.4f μl/μmol\n\n')%(defect,m,sum(entrapvolume),sum(entrapvolume)/(c*V)*1000))

    encap_old  (mu, sigma, sizerange, 0, d, molarea, c, V)
    
    if 'expEE' in globals() and abs(EE-expEE) >= abs(float(0.01)):
        encap_opt (mu, sigma, sizerange, defect, d, molarea, c, V, EE)
    
    else:return
    
        
            
    plotEE(entrapvolume, EE)
    
       
encap (mu, sigma, sizerange, defect, d, molarea, c, V)

              
exit()





























