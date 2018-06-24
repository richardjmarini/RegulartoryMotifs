#!/usr/bin/env python

from argparse import ArgumentParser
from collections import Counter

class MotifSearch(object):
    """
    Searches for motifs using branch and bound method.
    """

    def __init__(self, k):

        self.alignments= []
        self.k= k

        # scoring function: find consense string with max score
        self.score= lambda kmer, sequences: sum(
            [max(Counter(sequences[pos][val + i] for pos, val in enumerate(kmer)).values()) for i in range(self.k)]
        )

    def next_vertex(self, sequences, path, best, pruned):
    
        depth= len(path)
        l= len(sequences)
        
        if (depth == l):
    
            current= self.score(path, sequences)
    
            # if the current score is best
            # then save the path we took to get here
            # and return the current score
            # otherwise, the best score still stands

            if (current > best):
                self.alignments= path
                return current
            else:
                return best
    
        else:
    
            # optimisticScore <-- Score(s, i, DNA) + (t - i) * l
            opt_score= self.score(path, sequences) + (l - depth) * self.k if depth > 1 else self.k * l
    
            if (opt_score < best):
                pruned= pruned + 1
                return best

            else:

                for kmer in range(len(sequences[depth]) - (self.k + 1)):
                    next= tuple([i for i in path] + [kmer])
                    best= self.next_vertex(sequences, next, best, pruned)
    
                return best
    
    def branch_and_bound(self, sequences):
    
        pruned= 0
        best= 0
 
        # cycle through each k-mer starting position 
        for kmer in range(len(sequences[0]) - (self.k + 1)):
            best= self.next_vertex(sequences, [kmer], best, pruned)
    
        return self.alignments

    def __call__(self, sequences):

        return self.branch_and_bound(sequences)


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

