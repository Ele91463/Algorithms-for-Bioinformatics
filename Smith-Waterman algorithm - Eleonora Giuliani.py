
"""Eleonora Giuliani 247161"""

"""
Implementation of the Smith-Waterman algorithm.
Instructions: run the code and enter the two sequences and the gap penality in the terminal

"""

"""Create a class to store metrics of the alignments"""
class Align:
    def __init__(self, align_seq1, align_seq2, bars, score, length, matches, mismatches, gaps, start_pos, end_pos):
        self.align_seq1 = align_seq1
        self.align_seq2 = align_seq2
        self.bars = bars
        self.score = score
        self.length = length
        self.matches = matches
        self.mismatches = mismatches
        self.gaps = gaps
        self.start_pos = start_pos
        self.end_pos = end_pos

    def __str__(self):
        return ("Alignment(score={}, length={}, matches={}, mismatches={}, gaps={})".
                format(self.score,
                       self.length,
                       self.matches,
                       self.mismatches,
                       self.gaps))


"""Define the substitution matrix"""
def init_subs_matrix(match_score, mismatch_score):
    subs_matrix = {
        'A': {'A': match_score, 'C': mismatch_score, 'G': mismatch_score, 'T': mismatch_score},    
        'C': {'A': mismatch_score, 'C': match_score, 'G': mismatch_score, 'T': mismatch_score},
        'G': {'A': mismatch_score, 'C': mismatch_score, 'G': match_score, 'T': mismatch_score},
        'T': {'A': mismatch_score, 'C': mismatch_score, 'G': mismatch_score, 'T': match_score}
    }
    return subs_matrix

"""Initialize the scoring matrix"""
def init_scoring_matrix(seq1, seq2):
    x_len = len(seq1) + 1
    y_len = len(seq2) + 1
    A = [[0 for x in range(y_len)] for x in range(x_len)]  
    return A

"""Initialize a traceback matrix"""
def init_traceback_matrix(seq1, seq2):
    x_len = len(seq1) + 1
    y_len = len(seq2) + 1
    B = [[0 for X in range(y_len)] for X in range(x_len)]
    return B


"""Fill the scoring matrix(record the choice)"""
def scoring(subs_matrix, A, B, seq1, seq2, g_pen):
   
    #initialize variables to keep track of the highest score and its position in the matrix
    score = 0
    pos = []

    x_len = len(seq1) + 1
    y_len = len(seq2) + 1
    
    #iterate through each element of the scoring matrix
    for i in range(1, x_len):
        for j in range(1, y_len):
            #calculate scores for the three possible directions: diagonal, up and right
            diag_score = A[i - 1][j - 1] + subs_matrix[seq1[i - 1]][seq2[j - 1]]  #match/mismatch
            up_score = A[i - 1][j] + g_pen    #gap delete
            horiz_score = A[i][j - 1] + g_pen  #gap insert

            #find the max_score 
            max_score = max(0, diag_score, up_score, horiz_score)

            #record the origin of the max_score
            if max_score == 0:
                B[i][j] = None
            elif max_score == diag_score:
                B[i][j] = "diagonal"
            elif max_score == up_score:
                B[i][j] = "up"
            elif max_score == horiz_score:
                B[i][j] = "horizontal"

            #update the max score in the cell 
            A[i][j] = max_score

            #update the score of the matrix and its position
            if max_score > score:
                score = max_score
                pos = [(i, j)]
            elif max_score == score:   
                pos.append((i,j))
                
    if score <= 0:
        return 0, [], B
    else:
        return score, pos, B
   

"""Starting at the highest score in the scoring matrix and ending at a matrix cell that
 has a score of 0, traceback based on the source of each score recursively to generate 
 the best local alignment"""

def traceback(seq1, seq2, A, B, init_pos, visited):
    align1 = []
    align2 = []
    bars = []
    i, j = init_pos
         
    while i > 0 and j > 0 and A[i][j] != 0:
        if (i, j) in visited:    #to keep track of the visited cell 
            break
        visited.add((i, j))
       
        if B[i][j] == "diagonal": 
            align1.append(seq1[i - 1])
            align2.append(seq2[j - 1])
            if seq1[i - 1] == seq2[j - 1]:
                bars.append("|")   #match
            else:
                bars.append(" ")   #mismatch
            i -= 1    #move to the previous cell in diagonal
            j -= 1
        elif B[i][j] == "up": 
            align1.append(seq1[i - 1])
            align2.append("-")    #gap
            bars.append(' ')
            i -= 1     #move to the cell above
        elif B[i][j] == "horizontal": 
            align1.append("-")    #gap
            align2.append(seq2[j - 1])
            bars.append(' ')
            j -= 1     #move to the left cell

    #it is necessary to reverse the alignment
    align1.reverse()
    align2.reverse()
    bars.reverse()

    #convert the list in string
    align_seq1 = ''.join(align1)   
    align_seq2 = ''.join(align2)
    bars = ''.join(bars)
    
    #compute the metrics
    alignment_length = len(align_seq1)
    matches = bars.count("|")
    gaps = align_seq1.count('-') + align_seq2.count('-')
    mismatches = alignment_length - matches - gaps

    return Align(align_seq1, align_seq2, bars, A[init_pos[0]][init_pos[1]], 
                 alignment_length, matches, mismatches, gaps, init_pos, (i, j))


"""perform traceback for each cell with the highest score"""
def tot_alignments2(sequence1, sequence2, A, B, pos):
    res = []
    visited = set()   #set to store visited cells during traceback to avoid sub-alignments
    x_len = len(A)    
    y_len = len(A[0])  
    for i in range(x_len - 1, 0, -1):       #start iteration from bottom-right cell
        for j in range(y_len - 1, 0, -1):
            init_pos = (i, j)
            if init_pos in pos and init_pos not in visited:
                res.append(traceback(sequence1, sequence2, A, B, init_pos, visited))
    return res


def tot_alignments(sequence1, sequence2, A, B):   #because the overlaps in trace-back are allowed I prefer to not use the set visited()
    res = []
    visited = set()
    x_len = len(A)
    y_len = len(A[0])
    for i in range(x_len - 1, 0, -1):       #start iteration from bottom-right cell
        for j in range(y_len - 1, 0, -1):
            init_pos = (i, j)
            if init_pos not in visited:
                res.append(traceback(sequence1, sequence2, A, B, (i, j), visited))
    return res


def is_matches(align_seq1, align_seq2, length_match):
    max_match = 0
    consecut = 0
    length_a = False

    for a, b in zip(align_seq1, align_seq2):
        if a == b:
            consecut += 1
            if consecut > max_match:
                max_match = consecut
            if consecut >= length_match:
                length_a = True
        else:
            consecut = 0
    
    return length_a, max_match

def filter_alignments(alignments, max_score):
    filtered = []
    threshold = 0.6 * max_score
    for alignment in alignments:
        if (alignment.score >= threshold and is_matches(alignment.align_seq1, alignment.align_seq2, 3)[0]):
            filtered.append(alignment)
    
    return sorted(filtered, key=lambda x: is_matches(x.align_seq1, x.align_seq2, 3)[1], reverse=True)

#request input 
seq1 = input("Enter the first sequence: ")
seq2 = input("Enter the second sequence: ")
match_score = int(input("Enter the match score: "))
mismatch_score = int(input("Enter the mismatch score: "))
g_pen = int(input("Enter the gap penalty: "))


#define the values of the substitution matrix
subs_matrix = init_subs_matrix(match_score, mismatch_score)
#initialize the scoring matrix
A = init_scoring_matrix(seq1, seq2)
#initialize the traceback matrix
B = init_traceback_matrix(seq1, seq2)
#perform scoring
score, pos, B = scoring(subs_matrix, A, B, seq1, seq2, g_pen)
#perform traceback for each highest score
res = tot_alignments(seq1, seq2, A, B, pos)

filt_sorted = filter_alignments(res, max_score)


#print the results
align_num = 1  
if filt_sorted == []:
    print("No possible alignments")
else:
    for x in filt_sorted:
        print("Alignment {}:".format(align_num))
        print("Sequence 1:", x.align_seq1)
        print("           ", x.bars)
        print("Sequence 2:", x.align_seq2)
        print("Score:", x.score)
        print("Length of the alignment:", x.length)
        print("N.matches:", x.matches)
        print("N.mismatches:", x.mismatches)
        print("N.gaps:", x.gaps)
        print(" ")
        align_num += 1 

    