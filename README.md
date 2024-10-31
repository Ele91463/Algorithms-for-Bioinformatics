# Algorithms-for-Bioinformatics

Implementation of the Smith-Waterman Algorithm:

1. Determine the substitution matrix and the gap penalty scheme
2. Construct a scoring matrix and initialize its first row and first column
   a) The size of the scoring matrix is (n+1)*(m+1)
   b) Note the 0-based indexing
3. Fill the scoring matrix using this schema
4. Starting at the highest score in the scoring matrix and ending at a matrix cell that has a score of 0, traceback based on the source of each score recursively to generate the best local alignment


Substitution matrix and gap penalty:
A substitution matrix assigns each pair of nucleotides (or aminoacids) a score.
Usually matches get positive scores, whereas mismatches get relatively lower scores.
A gap penalty function determines the score cost for opening or extending gaps.

Initialize the scoring matrix:
The dimensions of the scoring matrix are 1+length of each sequence
All the elements of the first row and of the first column are set to 0
The extra first row and first column make it possible to align one sequence to another at any position, and setting them to 0 makes the terminal gap free from penalty

Scoring:
Score each element of the matrix from left to right, top to bottom
The highest score among all the possible scores and 0 is used and the source of that score is recorded

Traceback:
Starting at the element with the highest score, traceback based on the source of each score recursively, until 0 is encountered
The segments that have the highest similarity score based on the given scoring system is generated in this process
To obtain the second best local alignment, apply the traceback process starting at the second highest score outside the trace of the best alignment

