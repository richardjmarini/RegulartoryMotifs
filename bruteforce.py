#!/usr/bin/env python

from time import time
from argparse import ArgumentParser
from itertools import product


def brute_force(sequences, k):


    iterations= 0

    # cycle through the list of all possible permutations
    # of possible offsets of length k across all the sequences

    found= []
    for offset in product( range(len(sequences[0])-k), repeat= len(sequences) ):  
        
        for i in range(1, len(offset)):

            # if the k-mer of length k from the first sequence
            # match the k-mer of length k from the next sequence
            # then we found a motif between the two sequences

            if sequences[0][offset[0]:offset[0]+k] == sequences[i][offset[i]:offset[i]+k]:  
                found.append(((offset, sequences[0][offset[0]:offset[0]+k]), iterations))

            iterations+= 1

        # if we did find a match between two k-mer's across any of the sequences
        # then return the offsets and the k-mer we found
        if len(found): 

            return found


if __name__ == '__main__':

    parser= ArgumentParser(description= "Regulartory Motifs")
    parser.add_argument("--input", default= "implanted.txt", help= "file of sequences")
    parser.add_argument("--k", type= int, default= 8, help= "k-mer length")
    parser.add_argument("--start", type= int, default= 0, help= "starting sequence")
    parser.add_argument("--end", type= int, default= 5, help= "ending sequence")

    args= parser.parse_args()

    print args

    fh= open(args.input, "r")
    sequences= fh.readlines()
    fh.close()

    # search time will grow expodentially so limit 
    # the sequences we're considering via start:end slice

    start_time= time()
    results= brute_force(sequences[args.start:args.end], args.k)
    end_time= time()

    if len(results):
        for result, iterations in results:
            (offset, pattern)= result
            print "Found motif %s at offset %s in %s iterations" % (pattern, offset, iterations)
    else:
        print "Could not find any motifs"

    print end_time - start_time, "seconds"
