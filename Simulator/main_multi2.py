__author__ = 'pjweggy'

from DWF_new import *
import numpy as np
import sys
import time

def GenSingleInfo(parentlist,prefix,ll,pho,Recomfrag,i,j):
    p1=parentlist[2*j]+1 ##data use index from 1-500
    p2=parentlist[2*j+1]+1  ##data use index from 1-500
    P1Seq=np.loadtxt(prefix+'_indv'+str(p1)+'_gen'+str(i)+'_Seq.dat') ##file name need to be unified
    P2Seq=np.loadtxt(prefix+'_indv'+str(p2)+'_gen'+str(i)+'_Seq.dat') ##
    if i == 0:
        if p1 <=nindv/2:
            P1AP=np.zeros((P1Seq.shape))
        else:
            P1AP=np.ones((P2Seq.shape))
        if p2 <=nindv/2:
            P2AP=np.zeros((P1Seq.shape))
        else:
            P2AP=np.ones((P2Seq.shape))
    else:
        P1AP=np.loadtxt(prefix+'_indv'+str(p1)+'_gen'+str(i)+'_AP.dat')
        P2AP=np.loadtxt(prefix+'_indv'+str(p2)+'_gen'+str(i)+'_AP.dat')
    DSeq,DAP,DRec = GenDescSingle(P1Seq,P2Seq,P1AP,P2AP,ll,pho,Recomfrag)
    np.savetxt(prefix+'_indv'+str(j+1)+'_gen'+str(i+1)+'_Seq.dat',DSeq,fmt='%d')
    np.savetxt(prefix+'_indv'+str(j+1)+'_gen'+str(i+1)+'_AP.dat',DAP,fmt='%d')
    np.savetxt(prefix+'_indv'+str(j+1)+'_gen'+str(i+1)+'_Rec.dat',DRec,fmt='%d')
    return 0


prefix="sim"
nindv=200

ll = 5000001

No = nindv #total population size: pop1 + pop2
Ne = 10000 #migration rate is 0.5
Ng = 5 #number of generations to evolve
r = 200000.0
pho = float(r)/(4*Ne*(ll-1)) ##recombination rate:f() of r,ll; different for ms and macs: ms use r=4Ne*(ll-1)*pho, macs use r=4Ne*pho
#phopara=pho*ll
Recomfrag=np.loadtxt(prefix+'_Recomrdc_reshape.dat')
if Recomfrag.ndim == 1:
    Recomfrag=np.vstack((Recomfrag,np.zeros(len(Recomfrag))))



indvlist=list()
for i in range(nindv):
    indvlist.append(i)

##start DWF process
chrs = 0
if chrs == 0:
    for i in range(Ng):
        parentlist = ChooseParent(No)
        #print parentlist
        if(i == 0):
            Parent = parentlist
        else:
            Parent = np.vstack((Parent,parentlist))
        for j in range(nindv):
            GenSingleInfo(parentlist,prefix,ll,pho,Recomfrag,i,j)
        print ('NOW begin generation {0}'.format(i+1))


    np.savetxt(prefix+'_Parent.dat',Parent,fmt='%d')
elif chrs == 1:
    for i in range(Ng):
        parentlist = np.loadtxt(prefix+'_Parent.dat')
        for j in range(nindv):
            GenSingleInfo(parentlist[i],prefix,ll,pho,Recomfrag,i,j)
        print ('NOW begin generation {0}'.format(i+1))
