#!/usr/bin/env python

from argparse import ArgumentParser
from itertools import product


def brute_force(sequences, k):
    
    for offset in product( range(len(sequences[0])-k+1), repeat= len(sequences) ):  

        notfound= False

        for i in range(1, len(offset)):

            if sequences[0][offset[0]:offset[0]+k] != sequences[i][offset[i]:offset[i]+k]:  
                notfound= True
                break 

        if not notfound:
            return (offset, sequences[0][offset[0]:offset[0]+k])


if __name__ == '__main__':

    parser= ArgumentParser(description= "Regulartory Motifs")
    parser.add_argument("--input", default= "implanted.txt", help= "file of sequences")
    parser.add_argument("--k", type= int, default= 8, help= "k-mer length")
    parser.add_argument("--start", type= int, default= 0, help= "starting sequence")
    parser.add_argument("--end", type= int, default= 3, help= "ending sequence")

    args= parser.parse_args()

    fh= open(args.input, "r")
    sequences= fh.readlines()
    fh.close()
    
    found= brute_force(sequences[args.start:args.end], args.k)
    if found:
        (offset, pattern)= found
        print "Found motif %s at offset %s" % (pattern, offset)
    else:
        print "Could not find any motifs"
