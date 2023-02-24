# Author: Yiming Zhang
import numpy as np
import sys
import argparse
from UseNumba import *
from PreProcessing import *
from UseTools import *
from HillClimbing import *


def get_parser():
    parser = argparse.ArgumentParser(description='PedFAM')
    parser.add_argument('-g', dest='generation', help='Number of Generations', type = int, required=True)
    parser.add_argument('-b', dest='block', help='Number of Blocks', type = int, required=True)
    parser.add_argument('-r', dest='rate', help='Recombination Rate', required=True)
    parser.add_argument('-c', dest='Chrom', help='Number of Chromosomes', type = int, required=True)
    parser.add_argument('-p', dest='panel', help='Number of Reference Panels', type = int, required=True)
    return parser

# Main function
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()    
    time_start_a = time.time()
    
    Num_generation = args.generation
    Num_Block = args.block
    recombination_r = float(args.rate)
    num_Chro = args.Chrom
    num_ref = args.panel
    
    # Preprocessing
    Prob_No_Switch_Happen_Chrom = np.zeros((num_Chro, Num_Block, Num_Block))
    Prob_Switch_Happen_Chrom = np.zeros((num_Chro, Num_Block, Num_Block))
    Prob_Ref_f_Chrom = np.zeros((num_Chro, num_ref, Num_Block, Num_Block))
    Prob_Ref_m_Chrom = np.zeros((num_Chro, num_ref, Num_Block, Num_Block))

    for i in range(num_Chro):
        Prob_No_Switch_Happen, Prob_Switch_Happen, Prob_Ref_f, Prob_Ref_m = PreProcess('example/AF_C'+str(i+1)+'.dat', 'example/test_1_C'+str(i+1)+'.dat', 'example/test_2_C'+str(i+1)+'.dat', 'example/ref_position_C'+str(i+1)+'.dat', Num_Block, recombination_r, num_ref)
        Prob_No_Switch_Happen_Chrom[i,:,:] = Prob_No_Switch_Happen
        Prob_Switch_Happen_Chrom[i,:,:] = Prob_Switch_Happen
        Prob_Ref_f_Chrom[i,:,:,:] = Prob_Ref_f
        Prob_Ref_m_Chrom[i,:,:,:] = Prob_Ref_m
    
    # Construct the Pedigree
    nodes = range((2**(Num_generation+1)) - 2, -1, -1)

    lists = []
    for i in nodes:
        lists.append(i)
    list = np.array(lists)
    root = creat(None, list,0, None)
    nodelists = list_preorder(root)

    mergeSort(nodelists)
    
    nodelist = typed.List()
    for i in nodelists:
        nodelist.append(i)
        
        
    # Start the Hill-Climbing Start Point 1    
    time_start1 = time.time()
    Identical_Settings = typed.List.empty_list(types.string)
    
    # Start point 1: Half 0 and half 1

    list_founder1 = np.zeros((2**(Num_generation-1)))
    list_founder2 = np.ones((2**(Num_generation-1)))

    list_founder = np.concatenate((list_founder1,list_founder2))
    founders_Opt1, Prob_Opt1, Identical_Settings = Infer(list_founder, Identical_Settings, Prob_Ref_f_Chrom, Prob_Ref_m_Chrom, Num_Block, Prob_No_Switch_Happen_Chrom, Prob_Switch_Happen_Chrom, nodelist, Num_generation, num_Chro)
    time_end1 = time.time()
    print('DP1 time cost:', time_end1-time_start1, 's')
    founder_OPT1 = founders_Opt1.astype('int32').tolist()
    print(founder_OPT1, Prob_Opt1)      
    
    
    # Start point 2: 0101...
    time_start1 = time.time()
    list_founder1 = np.zeros((2**(Num_generation)))
    for x in range(1, len(list_founder1), 2):
        list_founder1[x] = 1
    list_founder = list_founder1.copy()
    founders_Opt2, Prob_Opt2, Identical_Settings = Infer(list_founder, Identical_Settings, Prob_Ref_f_Chrom, Prob_Ref_m_Chrom, Num_Block, Prob_No_Switch_Happen_Chrom, Prob_Switch_Happen_Chrom, nodelist, Num_generation, num_Chro)
    time_end1 = time.time()
    print('DP2 time cost:', time_end1-time_start1, 's')
    founder_OPT2 = founders_Opt2.astype('int32').tolist()
    print(founder_OPT2, Prob_Opt2)
    if(Prob_Opt2 > Prob_Opt1):
        Prob_Opt = Prob_Opt2
        founders_Opt = founders_Opt2[:]
    else:
        Prob_Opt = Prob_Opt1
        founders_Opt = founders_Opt1[:]
    time_start1 = time.time()
    
    # Start point 3: All 0
    list_founder1 = np.zeros((2**(Num_generation)))
    list_founder = list_founder1[:]
    founders_Opt3, Prob_Opt3, Identical_Settings = Infer(list_founder, Identical_Settings, Prob_Ref_f_Chrom, Prob_Ref_m_Chrom, Num_Block, Prob_No_Switch_Happen_Chrom, Prob_Switch_Happen_Chrom, nodelist, Num_generation, num_Chro)
    time_end1 = time.time()
    print('DP3 time cost:', time_end1-time_start1, 's')
    founder_OPT3 = founders_Opt3.astype('int32').tolist()
    print(founder_OPT3, Prob_Opt3)
    if(Prob_Opt3 > Prob_Opt):
        Prob_Opt = Prob_Opt3
        founders_Opt = founders_Opt3[:]

    # Start point 4: All 1
    time_start1 = time.time()
    list_founder1 = np.ones((2**(Num_generation)))
    list_founder = list_founder1[:]
    founders_Opt4, Prob_Opt4, Identical_Settings = Infer(list_founder, Identical_Settings, Prob_Ref_f_Chrom, Prob_Ref_m_Chrom, Num_Block, Prob_No_Switch_Happen_Chrom, Prob_Switch_Happen_Chrom, nodelist, Num_generation, num_Chro)
    time_end1 = time.time()
    print('DP4 time cost:', time_end1-time_start1, 's')
    founder_OPT4 = founders_Opt4.astype('int32').tolist()
    print(founder_OPT4, Prob_Opt4)
    if(Prob_Opt4 > Prob_Opt):
        Prob_Opt = Prob_Opt4
        founders_Opt = founders_Opt4[:]

    founder_OPT = founders_Opt.astype('int32').tolist()
    print('final optimal:', founder_OPT, Prob_Opt)
    print('Total settings tested:', len(Identical_Settings))
    time_end_a = time.time()
    print('Total time cost:', time_end_a-time_start_a, 's')    