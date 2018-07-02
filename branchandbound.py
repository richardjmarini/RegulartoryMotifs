#!/usr/bin/env python

from time import time
from pprint import pprint
from argparse import ArgumentParser
from collections import Counter

class MotifSearch(object):
    """
    Searches for motifs using branch and bound method.
    """

    def __init__(self, k, debug= False):

        self.offsets= []
        self.k= k
        self.debug= debug

        self.iterations= 0
        self.max_score= 0

    def score(self, kmer, sequences):
        """
        sequences: 2 dimensional array of sequences
        """
  
        if self.debug:
            print "\nkmer starting positions=", kmer
            profile_matrix= []

        # compute the consensus score for each slice
        # of the motif matrix
        consensus_score= 0
        for i in range(self.k):
 
            # create and populate profile matrix 
            profile= dict(A= 0, C= 0, T= 0, G= 0)
            for x in range(len(kmer)):
                y= kmer[x]
                profile[sequences[x][y+i]]+= 1

                if self.debug:
                    print "pos= %s, x= %s, y=%s, y+i=%s, sequence= %s alignment= %s" % (i, x, y, y+i, sequences[x], sequences[x][y+i])

            consensus_score+= max(profile.values())
            if self.debug:
                profile_matrix.append(profile)
                print 

        if self.debug:
            print "profile matrix="
            pprint(profile_matrix)
            print "consensus_score=", consensus_score

        return consensus_score                 

    def next(self, sequences, offsets, max_score):

        self.iterations+= 1
    
        num_sequences= len(sequences)
        num_offsets= len(offsets)
        score= max_score

        # if we're not yet at a leaf node 
        if (num_offsets != num_sequences):
    
            # page: 111
            # optimisticScore <-- Score(s, i, DNA) + (t - i) * l
            opt_score= self.score(offsets, sequences) + (num_sequences - num_offsets) * self.k
   
            # page: 108 + bottom of page 107
            if (opt_score >= max_score):
                   
                # for each k-mer starting position
                for kmer in range(len(sequences[num_offsets]) - (self.k + 1)):
                    max_score= self.next(sequences, offsets + [kmer], max_score)
    
                score= max_score

        else:

            # if we finally reached a leaf node
            # and the current_score is better than the max_score
            # then save the offset and return the current score
            # otherwise, the max_score still stands

            current_score= self.score(offsets, sequences)
            if (current_score > max_score):
                self.offsets= offsets
                score= current_score

        return score
    
    def search(self, sequences):
    
        max_score= 0
 
        #: page 109
        # cycle through each k-mer starting position 
        for kmer in range(len(sequences[0]) - (self.k + 1)):
            max_score= self.next(sequences, [kmer], max_score)

        self.max_score= max_score
    
        return self.offsets

    def __call__(self, sequences):

        return self.search(sequences)


if __name__ == '__main__':

    parser= ArgumentParser(description= "Regulartory Motifs")
    parser.add_argument("--input", default= "implanted.txt", help= "file of sequences")
    parser.add_argument("--k", type= int, default= 8, help= "k-mer length")
    parser.add_argument("--start", type= int, default= 0, help= "starting sequence")
    parser.add_argument("--end", type= int, default= 7, help= "ending sequence")
    parser.add_argument('--debug', dest= 'debug', action= 'store_true', default= False, help= "show debug trace")


    args= parser.parse_args()
    print args
    debug= args.debug
    if debug:
        print "args", args.__dict__

    fh= open(args.input, "r")
    sequences= [line.replace('\n', '').upper() for line in fh.readlines()]
    fh.close()
    
    motif_search= MotifSearch(args.k, debug)

    start_time= time()
    # list the sequence, position and pattern of each motif we found
    for (sequence, position) in zip(range(args.start, args.end), motif_search(sequences[args.start:args.end])):
             
         pattern= sequences[sequence][position:position + args.k]

         print "found motif", pattern, "at position", position, "in sequence", sequences[sequence]
    end_time= time()

    print "max score", motif_search.max_score, "in", motif_search.iterations, "iterations", "in", end_time - start_time, "seconds"
    
