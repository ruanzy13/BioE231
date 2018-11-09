#!/usr/bin/python

"""

Needleman-Wunsch Aligner
Bioengineering 131/231, Fall 2018

Command-line script to read in the contents of a multi-FASTA containing two
sequences and a score matrix, perform a global alignment, and print the
resulting alignment matrix and optimal alignment to STDOUT.

"""

import os
import sys

class NWAligner:
    def __init__(self, score_matrix_fname):
        self.score_matrix, self.gap_penalty = self.load_score_matrix(score_matrix_fname)

    @staticmethod
    def load_score_matrix(fname):
        score_matrix = {}
        gap_penalty = None

        with open(fname) as fp:
            for line_num, line in enumerate(fp):
                # ignore comments in matrix file
                if line.startswith("#"):
                    continue
                line=line.rstrip()
                if line_num==0:
                    line=line.replace(" ", "")
                    name=[]
                    for each in line:
                        score_matrix[each]={}# this is a new dictionary with name as score_matrix['A'] or other character
                        name.append(each)
                    #print(name)
                elif line_num > 0:
                    #line=line.replace(" ", "")
                    value=[]
                    line = list(line)
                    i = 0
                    j = 0
                    while i < len(line):
                        if line[i] =="-":
                            line[i+1]="-"+line[i+1]
                            del line[i]
                        i+=1
                    while j < len(line):
                        if j < len(line)-1 and line[j]!=" " and line[j+1]!=" ":
                            line[j+1]=line[j]+line[j+1]
                            del line[j]
                        j+=1
                    for each in line:
                        if each !=" ":
                            value.append(int(each))
                    while i < len(line):
                        if line[i]=="-":
                            line[i+1]="-"+line[i+1]
                        if i< (len(line)-1):# there is a exception with value 11
                            if line[i+1]!=" " or line[i+1]!="-":
                                line[i+1]=line[i]+line[i+1]
                        elif line[i] != " ":
                            print(line[i])
                            #value.append(int(line[i]))
                        i+=1
                        print(value)
                    if len(value)!=1:
                        j = 0
                        while j < len(name):
                            score_matrix[name[line_num-1]][name[j]]=value[j]
                            j+=1
                    elif len(value)==1:
                        gap_penalty = value[0]

        return score_matrix,gap_penalty





    @staticmethod
    def load_FASTA(fname):
        """
        Input: (String) A path to a FASTA file containing exactly two sequences.
        Output: (List) A list containing two strings: one for each sequence.

        Example:

        >>> seqs = NWAligner.load_FASTA('example.fa')
        >>> seqs[0]
        'YAADSKATPGNPAFHQDEIFLARIAFIYQMWDGGQLKLIDYAPHHVMCEE'
        >>> seqs[1]
        'WVGQPNMKVQHWSNMKACCVKFITWTFIAPEKHACKWTETAYQADCDIIW'
        >>> len(seqs)
        2

        """

        seqs = []
        s=''
        with open(fname) as fp:
            for line in fp:
                line=line.rstrip()
                if line.startswith('>'):
                    if s!='':
                        seqs.append(s)
                        s=''
                else:
                    s=s+line
        seqs.append(s)#last entry  was not saved since we didn't find a new header
        if len(seqs)>2:
            print("There are more than two sequences in the file")
        return seqs
        ### TODO ###
        # Load FASTA file and return list of sequences.
        # Throw an error if there are more than two sequences in the file.







    def align(self, seq_x, seq_y, print_matrix = False):
        """
        Input: (Strings) Two sequences to be aligned (seq_x and seq_y).
               (Boolean) If print_matrix is True, print the dynamic programming
                         matrix before traceback.
        Output: (Tuple of strings) Two sequences, aligned.

        Example:

        >>> aligner = NWAligner('BLOSUM62')
        >>> seqs = aligner.load_FASTA('example.fa')
        >>> aligner.align(seqs[0], seqs[1])
        ('YAAD-SKATPGNPAF---HQDEIF--L-AR--IA-FIYQM-WDGGQLK-LIDYAPH-HVM-C---E-------E---',
         'W---VGQ--P-N--MKVQH----WSNMKA-CCV-KFI---TW------TFI--APEKH--ACKWTETAYQADCDIIW')

        """

        ###
        ### INITIALIZATION
        ###

        # create two empty matrices with sizes based on the input sequences.
        # one contains the dynamic programming matrix, the other contains
        # pointers we'll use during traceback
        matrix = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]#x is y-axis
        pointers = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]
        #score_matrix,gap_penalty = NWAligner.load_score_matrix(score_matrix_fname)########## This one needs more disscussion
        #Fill the top row of the matrix with scores for gaps
        i=0
        while i < (len(seq_y)+1):
            matrix[0][i]=0+i *self.gap_penalty
            i+=1
        #Fill the first column of the matrix with scores for gaps.
        j=0
        while j < (len(seq_x)+1):
            matrix[j][0]=0+j *self.gap_penalty
            j+=1
        #return(matrix)

        ### TODO ###
        # Fill the top row of the matrix with scores for gaps.
        # Fill the first column of the matrix with scores for gaps.








        ###
        ### RECURSION
        ###

        # fill the dynamic programming and pointer matrices
        for x in range(1, len(seq_x) + 1): #final x would be x=len(seq_x)
            for y in range(1, len(seq_y) + 1):
                match_score = self.score_matrix[seq_x[x - 1]][seq_y[y - 1]]#this is the score if align the two
                s_diag=matrix[x-1][y-1]+match_score
                s_lgap=matrix[x][y-1]+self.gap_penalty
                s_agap=matrix[x-1][y]+self.gap_penalty
                movement=[s_diag,s_lgap,s_agap]
                max_val=max(s_diag,s_lgap,s_agap)
                matrix[x][y]=max_val

                member=[]
                for nm,each in enumerate(movement,1):#1:s_diag 2:s_lgap 3:s_agap
                    if each == max_val:
                        member.append(nm)
                if len(member)==1:#to count how many of three is in max
                    if matrix[x][y]==s_lgap:
                        pointers[x][y]=-1# seq_x add a gap
                    elif matrix[x][y]==s_diag:
                        pointers[x][y]=2# 0 refers to exact match
                    elif matrix[x][y]==s_agap:
                        pointers[x][y]=1#seq_y add a gap
                elif len(member) == 2:
                    if 1 not in member:
                        if matrix[x-1][y] > matrix[x][y-1]:
                            pointers[x][y]=1#seq_y add a gap
                        else:
                            pointers[x][y]=-1# seq_x add a gap
                    elif 2 not in member:
                        if matrix[x-1][y-1] > matrix[x-1][y]:
                            pointers[x][y]=2# 0 refers to exact match
                        else:
                            pointers[x][y]=1#seq_y add a gap
                    else:
                        if matrix[x-1][y-1] > matrix[x][y-1]:
                            pointers[x][y]=2# 0 refers to exact match
                        else:
                            pointers[x][y]=-1# seq_x add a gap

                ### TODO ###
                # Take the maximum score of three possibilities:
                #   1) The element in the matrix diagonal from this one
                #      plus the score of an exact match
                #   2) The element to the left plus a gap penalty
                #   3) The element above plus a gap penalty
                # ... and set the current element (matrix[x][y]) equal to that
                #
                # Keep track of which of these choices you made by setting
                #   the same element (i.e., pointers[x][y]) to some value that
                #   has meaning to you.

        # print the dynamic programming matrix
        if print_matrix:
            for x in range(len(seq_x) + 1):
                print (" ".join(map(lambda i: str(int(i)), matrix[x])))

        ###
        ### TRACEBACK
        ###

        # starting from the bottom right corner, follow the pointers back
        x, y = len(seq_x), len(seq_y)

        # fill these lists with the aligned sequences
        align_x = []
        align_y = []

        while x > 0 or y > 0:
            move = pointers[x][y]
            if move==2:
                align_x.append(seq_x[x-1])
                align_y.append(seq_y[y-1])
                x-=1
                y-=1
            elif move==-1:
                align_x.append("-")
                align_y.append(seq_y[y-1])
                y-=1
            elif move==1:
                align_x.append(seq_x[x-1])
                align_y.append("-")
                x-=1



            ### TODO ###
            # Follow pointers back through the matrix to the origin.
            # Depending on which "move" you made at each element in the
            #   matrix, you'll either align seq_x to seq_y, seq_x to a gap, or
            #   seq_y to a gap.

        # flip the alignments, as they're reversed
        return ("".join(align_x[::-1]), "".join(align_y[::-1]))

    ###                                      ###
    ### NO NEED TO EDIT CODE BELOW THIS LINE ###
    ###                                      ###

if __name__ == '__main__':
    def usage():
        print ('usage: %s matrixfilename stringfilename')
        sys.exit(1)

    if len(sys.argv) != 3:
        usage()

    for fname in sys.argv[1:]:
        if not os.path.isfile(fname):
            print ('Can not open %s' % (fname,))
            usage()

    aligner = NWAligner(sys.argv[1])
    seqs = aligner.load_FASTA(sys.argv[2])
    result = aligner.align(seqs[0], seqs[1])

    print('>seq1\n%s\n>seq2\n%s' % (result[0], result[1]))
