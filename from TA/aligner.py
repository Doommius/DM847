#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 17:04:39 2018

@author: alexander


run it like this:
    
    python3 aligner.py BLOSUM62.txt 0 globindequences.fasta



NOTE!! If have not imlemented a full version of the Needleman-Wunsch algorithm.
But I calculate the scoring matrices and the parts of the traceback matrices.
The alignments and the total blast scores are printed, when you run the program.


"""





import sys
import numpy as np


nm = ['HBB_HUMAN',
      'HBB_GORGO',
      'HBB_COLPO',
      'HBB_SEMEN',
      'HBB_MARFO',
      'HBB_CANLA',
      'HBB_URSMA',
      'HBB_LYNLY',]



def load_BLOSUM62(matrix_filename, gap_penalty):
    
    ## importing relevant rows of the BLOSUM62 matrix.
    with open(matrix_filename) as matrix_file:
        lst = []
        k = 0
        for m in matrix_file:
            if m.strip('\n').split()[0] == 'A' and k == 0:
                lst.append(m.strip('\n').split())
                k += 1
            else:
                lst.append(m.strip('\n').split()[1:])

    ## removing the description
    lst = lst[6:]

    # Converting the string numbers into integers
    for i,j in enumerate(lst):
        if i > 0:
            lst[i] = [int(m) for m in j]
            lst[i][-1] = int(gap_penalty)
        #print(j)
        
    # Converting the lists with integers into a matrix
    matrix = np.array(lst[1:])
    aa = lst[0]
    print(matrix)

    return aa, matrix



def load_FASTA(FASTA_filename):
    
    ## importing relevant rows of the BLOSUM62 matrix.
    with open(FASTA_filename) as fasta:
        lst = []
        names = []
        for m in fasta:
            if m.strip('\n')[0] == '>':
                lst.append([])
                names.append(m.strip('\n'))
            else:
                lst[-1] += m.strip('\n').split()

    for i in range(len(lst)):
        lst[i] = ''.join(lst[i])

    
    return lst, names




def fill_out_score_mtrx_equal_lens(sq0, sq1, mtrx, tr, aa, blosum62, n1, n2):
    
    for i in range(len(mtrx)-1):
        mtrx[i+1][i+1] =  mtrx[i][i]+blosum62[aa.index(sq0[i])][aa.index(sq1[i])]
        
    print(n1, sq0)
    print(n2, sq1)
    score = [int(m) for m in mtrx[np.diag_indices_from(mtrx)]]
    print(score[-1])
    
    
    print('\n\n')
    



def fill_out_score_mtrx(sq0, sq1, mtrx, tr, aa, blosum62, n1, n2, mis):
    
    print('\n\n')
    gap_i = -2
    gap_j = -2
    for i in range(len(mtrx)-1):
        for j in range(len(mtrx[0])-1):
            if mtrx[i+1][j+1] == 0:
                
                ## I am only doing anything if an insertion in made 
                if np.argmax([ mtrx[i][j]+blosum62[aa.index(sq0[j])][aa.index(sq1[i])], mtrx[i][j+1]+gap_j, mtrx[i+1][j]+gap_i]) == 2:
                    mtrx[i+1][j+1] = np.max([ mtrx[i][j]+blosum62[aa.index(sq0[j])][aa.index(sq1[i])], mtrx[i][j+1]+gap_j, mtrx[i+1][j]+gap_i])
                    tr[i+1][j+1] = 2
                    ## 2 = left
                    gab_i = 0
                    gab_j = 0
                if np.argmax([ mtrx[i][j]+blosum62[aa.index(sq0[j])][aa.index(sq1[i])], mtrx[i][j+1]+gap_j, mtrx[i+1][j]+gap_i]) == 0:
                    mtrx[i+1][j+1] = np.max([ mtrx[i][j]+blosum62[aa.index(sq0[j])][aa.index(sq1[i])], mtrx[i][j+1]+gap_j, mtrx[i+1][j]+gap_i])
                    tr[i+1][j+1] = 0
                    ## 0 = diag
    
    
    align = ''
    r = 0
    d = 0
    reg = 0
    reg += mis
    
    for i in range(len(sq0))[::-1]:
        if i == 0 and len(tr[0]) > len(tr) and mis > 0:
            align += '_'
            break
        if i-reg < 0:
            break
        if tr[i+d - reg][i+r] == 0:
            align += sq1[i-reg]
            
        if tr[i+d - reg][i+r] == 2 and mis > 0:
            align += sq1[i-reg]+'_'
            
            for b in tr[i+d - reg][:i+r][::-1]:
                if b == 2:
                    align += '_'
                    mis -= 1
                    r -= 1
                if b == 0:
                    break
   
    
    print(n1, sq0)
    #print(n2, sq1)
    align = align[::-1]
    print(n2, align)
    print(mtrx[-1][-1])


    print('\n\n')







if __name__ == '__main__':
    aa, blosum62 = load_BLOSUM62(sys.argv[1], sys.argv[2])
    seqs, names = load_FASTA(sys.argv[3])
    print(aa)
    #print(names)
    
    
    for i in range(1, len(seqs)):
        
        ## did like this because the human sequence is never the shortest
        mtrx = np.zeros(( len(seqs[i])+1, len(seqs[0])+1 ))
        tr = np.zeros(( len(seqs[i])+1, len(seqs[0])+1 ))
        
        ## initializing the penalty rows and columns
        mtrx[0][1:] += int(sys.argv[2])
        mtrx[:,0][1:] += int(sys.argv[2])
  
        
        if len(seqs[0]) == len(seqs[i]):
            fill_out_score_mtrx_equal_lens(seqs[0], seqs[i], mtrx, tr, aa, blosum62, nm[0], nm[i])
        
        if len(seqs[0]) > len(seqs[i]):
            fill_out_score_mtrx(seqs[0], seqs[i], mtrx, tr, aa, blosum62, nm[0], nm[i], len(seqs[0]) - len(seqs[i]))
        
        ## Does not matter in this case
        if len(seqs[0]) < len(seqs[i]):
            pass
            
            
    
    
    