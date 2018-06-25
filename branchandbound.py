#!/usr/bin/env python

from argparse import ArgumentParser
from collections import Counter

class MotifSearch(object):
    """
    Searches for motifs using branch and bound method.
    """

    def __init__(self, k):

        self.offsets= []
        self.k= k

    def score(self, kmer, sequences):
        """
        sequences: 2 dimensional array of sequences
        """
        
        # compute the consensus score for each slice
        # of the motif matrix
        consensus_score= 0
        for i in range(self.k):
 
            # create and populate profile matrix 
            profile= dict(A= 0, C= 0, T= 0, G= 0)
            for x in range(len(kmer)):
                y= kmer[x]
                profile[sequences[x][y+i]]+= 1

            consensus_score+= max(profile.values())

        return consensus_score                 

    def next_vertex(self, sequences, offsets, max_score, skip):
    
        num_offsets= len(offsets)
        num_sequences= len(sequences)
        
        # if we're at a leaf node
        if (num_offsets == num_sequences):
    
            current_score= self.score(offsets, sequences)
    
            # if the current_score is max_score
            # then save the path we took to get here
            # and return the current_score
            # otherwise, the max_score still stands

            if (current_score > max_score):
                self.offsets= offsets
                return current_score
            else:
                return max_score
    
        else:
    
            # page: 111
            # optimisticScore <-- Score(s, i, DNA) + (t - i) * l
            opt_score= self.score(offsets, sequences) + (num_sequences - num_offsets) * self.k if num_offsets > 1 else self.k * num_sequences
   
            # page: 108 + bottom of page 107
            # if the optimistic score is worst than the max_score
            # then we can prune this branch and keep the max_score 
            # BYPASS 
            if (opt_score < max_score):
                skip= skip + 1
                return max_score

            else:
                   
                # for each starting positoin of the kmer 
                for kmer in range(len(sequences[num_offsets]) - (self.k + 1)):
                    max_score= self.next_vertex(sequences, offsets + [kmer], max_score, skip)
    
                return max_score
    
    def search(self, sequences):
    
        skip= 0
        max_score= 0
 
        #: page 109
        # cycle through each k-mer starting position 
        for kmer in range(len(sequences[0]) - (self.k + 1)):
            max_score= self.next_vertex(sequences, [kmer], max_score, skip)
    
        return self.offsets

    def __call__(self, sequences):

        return self.search(sequences)


if __name__ == '__main__':

    parser= ArgumentParser(description= "Regulartory Motifs")
    parser.add_argument("--input", default= "implanted.txt", help= "file of sequences")
    parser.add_argument("--k", type= int, default= 8, help= "k-mer length")
    parser.add_argument("--start", type= int, default= 0, help= "starting sequence")
    parser.add_argument("--end", type= int, default= 3, help= "ending sequence")

    args= parser.parse_args()

    fh= open(args.input, "r")
    sequences= [line.replace('\n', '').upper() for line in fh.readlines()]
    fh.close()
    
    search= MotifSearch(args.k)

    # list the sequence, position and pattern of each motif we found
    for (sequence, position) in zip(range(args.start, args.end), search(sequences[args.start:args.end])):

         pattern= sequences[sequence][position:position + args.k]

         print "found motif", pattern, "at position", position, "in sequence", sequences[sequence]

