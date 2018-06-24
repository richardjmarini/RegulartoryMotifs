#!/usr/bin/env python

from argparse import ArgumentParser
from bruteforce import brute_force
from branchandbound import MotifSearch

if __name__ == '__main__':

    parser= ArgumentParser(description= "Regulartory Motifs")
    parser.add_argument("--input", default= "implanted.txt", help= "file of sequences")
    parser.add_argument("--k", type= int, default= 8, help= "k-mer length")
    parser.add_argument("--start", type= int, default= 0, help= "starting sequence")
    parser.add_argument("--end", type= int, default= 3, help= "ending sequence")
    parser.add_argument('--optimization', type= int, default= 0, help= "sets optimization level [0= bruteforce, 1= branchandbound")

    args= parser.parse_args()

    fh= open(args.input, "r")
    sequences= [line.replace('\n', '').upper() for line in fh.readlines()]
    fh.close()

    if args.optimization == 0:
    
        results= brute_force(sequences[args.start:args.end], args.k)
        if len(results):
            for result in results:
                (offset, pattern)= result
                print "Found motif %s at offset %s" % (pattern, offset)
        else:
            print "Could not find any motifs"

    elif args.optimization == 1:

        search= MotifSearch(args.k)

        # list the sequence, position and pattern of each motif we found
        for (sequence, position) in zip(range(args.start, args.end), search(sequences[args.start:args.end])):
             pattern= sequences[sequence][position:position + args.k]
             print "found motif", pattern, "at position", position, "in sequence", sequences[sequence]

    else:
        print stderr, "uknown optimization level"

