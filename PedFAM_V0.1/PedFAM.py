# Author: Yiming Zhang
import numpy as np
import random
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
    parser.add_argument('-s', dest='start', help='Number of Random Starting points', type = int, default=5)
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
    num_start = args.start
    
    # Preprocessing
    Prob_No_Switch_Happen_Chrom = np.zeros((num_Chro, Num_Block, Num_Block))
    Prob_Switch_Happen_Chrom = np.zeros((num_Chro, Num_Block, Num_Block))
    Prob_Ref_f_Chrom = np.zeros((num_Chro, num_ref, Num_Block, Num_Block))
    Prob_Ref_m_Chrom = np.zeros((num_Chro, num_ref, Num_Block, Num_Block))

    for i in range(num_Chro):
        Prob_No_Switch_Happen, Prob_Switch_Happen, Prob_Ref_f, Prob_Ref_m = PreProcess('example/AF_C'+str(i+1)+'.dat', 'example/Geno_C'+str(i+1)+'.dat', 'example/Position_C'+str(i+1)+'.dat', Num_Block, recombination_r, num_ref)
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
        
        
    # Start the Hill-Climbing   
    time_start1 = time.time()
    Identical_Settings = typed.List.empty_list(types.string)
    PROB_OPT = np.iinfo(np.int32).min
    FOUNDER_OPT = np.empty((2**(Num_generation)))

    for x in range(num_start):
        time_start1 = time.time()
        list_founder = np.asarray([random.randint(0,num_ref-1) for _ in range(2**(Num_generation))])
        founders_Opt, Prob_Opt, Identical_Settings = Infer(list_founder, Identical_Settings, Prob_Ref_f_Chrom, Prob_Ref_m_Chrom, Num_Block, Prob_No_Switch_Happen_Chrom, Prob_Switch_Happen_Chrom, nodelist, Num_generation, num_Chro, num_ref)
        time_end1 = time.time()
        print('Start', x+1, 'time cost:', time_end1-time_start1, 's')
        if founders_Opt is not None:
            founders_Opt = founders_Opt.astype('int32').tolist()
        print(founders_Opt, Prob_Opt)     
        if(Prob_Opt > PROB_OPT):
            PROB_OPT = Prob_Opt
            FOUNDER_OPT = founders_Opt.copy() 

    time_end_a = time.time()
    print('Total time cost:', time_end_a-time_start_a, 's')   
    print('Final:', FOUNDER_OPT, PROB_OPT)   