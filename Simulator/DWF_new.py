#!/usr/bin/env python
__author__ = 'pjweggy'


import numpy as np
import random
#from scipy import *
import timeit
import sys
import pickle

start=timeit.timeit()

##dataset
##first 2Nf rows are females
##last 2Nm raws are males
##each line is a haplotype; nearby odd and even rows form one individual's two haplotypes

def ChooseParent(No):
    parentlist = list()
    for i in range(No):
        temp=random.sample(range(No),2)
        parentlist.append(temp[0])
        parentlist.append(temp[1])
    return parentlist

'''
def GenDescSeqonly(PSeq,PAP,parentlist,ll,pho,recompos):
    nr,nc = PSeq.shape
    DSeq = zeros((nr,nc))
    DAP = zeros((nr,nc))
    DRec =zeros((nr,nc))
    for i in range(len(parentlist)):
        pp = parentlist[i]
        DSeq[i,],DAP[i,],DRec[i,] = GenDescinfoSeqonly(PSeq[(2*pp):(2*pp+2),],PAP[(2*pp):(2*pp+2),],ll,pho,recompos)
    return DSeq,DAP,DRec

def GenDescinfoSeqonly(Pseq,Pseqap,ll,pho,recompos):
    ##recombination events
    r_c=random.binomial(ll,pho) ##recombination event counts
    r_pos=list()
    if r_c != 0:
        for i in range(r_c):
            r_pos.append(random.rand())
    r_pos=sorted(r_pos)

    nlen=Pseq.shape[1]
    k=random.randint(0,1)
    dseq=zeros((1,nlen))
    dap=list()
    r_point=0
    for i in range(nlen):
        if recompos[i]<r_pos(r_point):
            dap.append(k)
        else:
            k=1-k
            dap.append(k)
            r_point += 1
    for i in range(nlen):
        dseq[0,i]=Pseq[dap(i),i]
        dap
    return dseq,r_pos
'''

def GenDescSingle(P1Seq,P2Seq,P1AP,P2AP,ll,pho,Recomfrag):
    nr,nc=P1Seq.shape
    DSeq=np.zeros((nr,nc))
    DAP=np.zeros((nr,nc))
    DRec=np.zeros((nr,nc))
    for i in range(int(nr/2)):
        recompos=Recomfrag[i]
        
        DSeq[2*i,],DAP[2*i,],DRec[2*i,]=GenDescinfo(P1Seq[(2*i):(2*i+2),],P1AP[(2*i):(2*i+2),],ll,pho,recompos)
        DSeq[2*i+1,],DAP[2*i+1,],DRec[2*i+1,]=GenDescinfo(P2Seq[(2*i):(2*i+2),],P2AP[(2*i):(2*i+2),],ll,pho,recompos)
    return DSeq,DAP,DRec


def GenDesc(PSeq,PAP,parentlist,ll,pho,recompos):
    nr,nc = PSeq.shape
    DSeq = np.zeros((nr,nc))
    DAP = np.zeros((nr,nc))
    DRec = np.zeros((nr,nc))
    for i in range(len(parentlist)):
        pp = parentlist[i]
        DSeq[i,],DAP[i,],DRec[i,] = GenDescinfo(PSeq[(2*pp):(2*pp+2),],PAP[(2*pp):(2*pp+2),],ll,pho,recompos)
    return DSeq,DAP,DRec

##two sequences for one
def GenDescinfo(Pseq,Pseqap,ll,pho,recompos):
    random.seed()
    np.random.seed()
    nlen=Pseq.shape[1]
    ##1.decide which chr to inherit
    k=random.randint(0,1)
    ##2.decide recombination events
    r_c=np.random.binomial(ll,pho) ##recombination event counts
    r_pos=list()
    if r_c != 0:
        for i in range(r_c):
            pos=np.random.rand()*(recompos[-1]-recompos[0])+recompos[0]
            r_pos.append(pos)
    else:
        r_pos.append(0)
    r_pos=sorted(r_pos)
    r_pos.append(1)
    rec_event=list()
    r_point=0
    if r_pos[r_point] == 0:
        recc = 0
        rec_event = [k] * nlen
    else:
        recc=len(r_pos)-1
        for i in range(nlen):
            if recompos[i]<r_pos[r_point]:
                rec_event.append(k)
            else:
                k=1-k
                rec_event.append(k)
                r_point += 1

    #3.generate seq and ap according to recombination events
    dseq=np.zeros((1,nlen))
    dap=np.zeros((1,nlen))

    for i in range(nlen):
        dseq[0,i]=Pseq[rec_event[i],i]
        dap[0,i] = Pseqap[rec_event[i],i]
    #print dseq,dap
    print('there are {0} recombination events'.format(recc))
    return [dseq,dap,rec_event]




def GenDesctrace(Nf,Nm,pho,length):

##[parent1(female) ; parent2(male) ; main haplotype herited from parent1 ; main haplotype herited from parent2 ;
##recombination position from parent1 (0-length of haplotype) ;  recombination position from parent2]
    descinfo=np.zeros((1,6))
    descinfo[0,0]=random.randint(1,Nf)
    descinfo[0,1]=Nf+random.randint(1,Nm)
    descinfo[0,2]=random.randint(1,2)
    descinfo[0,3]=random.randint(1,2)

    temp=random.uniform(0,1)
    if temp <= pho*float(length):
        descinfo[0,4]=random.randint(1,length)
    else:
        descinfo[0,4]=0

    temp=random.uniform(0,1)
    if temp <= pho*float(length):
        descinfo[0,5]=random.randint(1,length)
    else:
        descinfo[0,5]=0
    return descinfo

def GenDescseq(Pseq,tracemat):
    l=Pseq.shape[1]
    dseq=np.zeros((2,l))

    if tracemat[4]==0:
        dseq[0]=Pseq[2*(tracemat[0]-1)+tracemat[2]-1]
        #print dseq[0]
        #print Pseq[2*(tracemat[0]-1)+tracemat[2]-1]
    else:
        dseq[0,0:tracemat[4]]=Pseq[2*(tracemat[0]-1)+tracemat[2]-1,0:(tracemat[4])]
        dseq[0,tracemat[4]:l]=Pseq[2*(tracemat[0]-1)+trans(tracemat[2])-1,tracemat[4]:l]
        #print dseq[0]
        #print Pseq[2*(tracemat[0]-1)+tracemat[2]-1,0:(tracemat[4])]
        #print Pseq[2*(tracemat[0]-1)+trans(tracemat[2])-1,tracemat[4]:l]
        #print Pseq[2*(tracemat[0]-1)+trans(tracemat[2])-1]

    if tracemat[5]==0:
        dseq[1]=Pseq[2*(tracemat[1]-1)+tracemat[3]-1]
    else:
        dseq[1,0:(tracemat[5])]=Pseq[2*(tracemat[1]-1)+tracemat[3]-1,0:(tracemat[5])]
        dseq[1,tracemat[5]:l]=Pseq[2*(tracemat[1]-1)+trans(tracemat[3])-1,tracemat[5]:l]
    return dseq



def trans(k):
    s=0
    if k==1:
        s=2
    elif (k==2):
        s=1
    return s

def Recevent_comb(rec_event1,rec_event2):
    nlen1=len(rec_event1)
    nlen2=len(rec_event2)
    if nlen1==nlen2:
        Rec_event=np.vstack((rec_event1,rec_event2))
    elif nlen1<nlen2:
        for i in range(nlen2-nlen1):
            rec_event1.append(0)
        Rec_event=np.vstack((rec_event1,rec_event2))
    elif nlen2<nlen1:
        for i in range(nlen1-nlen2):
            rec_event2.append(0)
        Rec_event=np.vstack((rec_event1,rec_event2))
    return Rec_event

###just one generation ahead
def traceback(Seq,ACPaint,Trace,Recomb,hap_ind):
    key='generation5'
    print('current generation haplotype is:')
    print(seq2str(Seq[key][hap_ind,]))
    print('current generation ancestor painting is:')
    print(seq2str(ACPaint[key][hap_ind,]))

    p_ind=int(Trace[key][hap_ind,])
    rec_event=Recomb[key][hap_ind,]

    key_p='generation4'
    print('parent generation haplotype is:')
    print(seq2str(Seq[key_p][2*(p_ind-1),]))
    print(seq2str(Seq[key_p][(2*(p_ind-1)+1),]))

    print('parent generation ancestor painting is:')
    print(seq2str(ACPaint[key_p][2*(p_ind-1),]))
    print(seq2str(ACPaint[key_p][(2*(p_ind-1)+1),]))

    if rec_event[-1]==0:
        print('recombination event: hap comes from left parent')
    else:
        print('recombination event: hap comes from right parent')
    a=rec_event[0:-1]
    #print a
    b=np.where(a==1)
    #print b[0]
    #print b[1]
    if len(b[0])==0:
        print('          no recombination happened!')
    else:
        #ss=seq2str(b[0])
        print('          recombination happens on positions:')
        print('          ',b[0])
    return 0

def seq2str(seq):
    s=''
    for x in seq:
        s=s+str(int(x))
    return s






