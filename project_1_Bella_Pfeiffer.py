from typing import List
#codon table for future functions
CODON_TABLE = {
    'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S',
    'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y', 'UGU': 'C',
    'UGC': 'C', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAU': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T',
    'ACG': 'T', 'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGU': 'S',
    'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V',
    'GUG': 'V', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAU': 'D',
    'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
}


#task 1 ------
def HammingDistance(s: str, t: str) -> int:
    # verify input: check if s and t are the same length
    if len(s) != len(t):
        raise ValueError("Inputs are not the same length")
    
    # verify input: check if s and t contain only 'C', 'A', 'G', or 'T' characters
    valid_chars = {'C', 'A', 'G', 'T'}
    if not all(char in valid_chars for char in s) or not all(char in valid_chars for char in t):
        raise ValueError("Strings must contain only 'C', 'A', 'G', or 'T' characters")
    
    # initialize counter variable
    counter = 0
    # corrected the range and indexing
    for i in range(len(s)):
        if s[i] != t[i]:
            counter += 1
    return counter


#tester for HammingDistance
output1 = HammingDistance("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT")
print(output1)

#task 2 -----
def Translation(s: str) -> str:
    # verify input: check if s and t contain only 'C', 'A', 'G', or 'U' characters
    valid_chars = {'C', 'A', 'G', 'U'}
    if not all(char in valid_chars for char in s):
        raise ValueError("Strings must contain only 'C', 'A', 'G', or 'U' characters")
    
    # creating new string for storage
    protein_string = ""
    # each 3 = codon, so step = 3
    for i in range(0, len(s), 3):
        # define the codon using i
        codon = s[i:i+3]
        # find match in CODON_TABLE, add amino acid to new string
        amino_acid = CODON_TABLE.get(codon, '') 
        if amino_acid == 'Stop':  # stop translation when a stop codon is encountered
            break
        protein_string += amino_acid
    return protein_string

#tester for translation 
output2 = Translation("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA")
print(output2)

#task 3 ----
def FindingMotif(s: str, t: str) -> List[int]:
    #create list
    locations = []
   
    #check if t is longer than s
    if len(s) <= len(t):
        raise ValueError("Your first string cannot be shorter than your first string. Please try again.")
    # verify input: check if s and t contain only 'C', 'A', 'G', or 'T' characters
    valid_chars = {'C', 'A', 'G', 'T'}
    if not all(char in valid_chars for char in s) or not all(char in valid_chars for char in t):
        raise ValueError("Strings must contain only 'C', 'A', 'G', or 'T' characters")
    
    # find the position of the symbol you're searching for
    # go through s, looking for matches with t
    for i in range(len(s) - len(t) + 1):
        # check if the substring of `s` from `i` to `i + len(t)` matches `t`
        if s[i:i+len(t)] == t:
            # add the position (1-based index)
            locations.append(i + 1)
    # return the locations
    return locations

#tester for finding motif
output3 = FindingMotif("GATATATGCATATACTT", "ATAT")
print(output3)

#task 4 -------
def RNASplicing(s: str, introns: List[str]) -> str:
    # verify input: check if s and t contain only 'C', 'A', 'G', or 'T' characters
    valid_chars = {'C', 'A', 'G', 'T'}
    if not all(char in valid_chars for char in s):
        raise ValueError("Strings must contain only 'C', 'A', 'G', or 'T' characters")
    
     # remove introns from the DNA sequence
    for intron in introns:
        s = s.replace(intron, "")
    
    # transcribe DNA to RNA by replacing 'T' with 'U'
    rna = s.replace('T', 'U')
    
    # Translate RNA to protein
    protein = []
    for i in range(0, len(rna) - 2, 3):  # Step by 3 to read each codon
        codon = rna[i:i+3]
        amino_acid = CODON_TABLE.get(codon, '')
        if amino_acid == 'Stop':  # Stop codon
            break
        protein.append(amino_acid)
    return ''.join(protein)

# tester for rna splicing
introns = ["ATCGGTCGAA", "ATCGGTCGAGCGTGT"]
output4 = RNASplicing("ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG", introns)
print(output4)

#task 5 -------
def LongestCommonSubstring(k: List[str]) -> str:
    # verify input: check if s and t contain only 'C', 'A', 'G', or 'T' characters
    valid_chars = {'C', 'A', 'G', 'T'}
    for seq in k:
        if not all(char in valid_chars for char in seq):
            raise ValueError("Strings must contain only 'C', 'A', 'G', or 'T' characters")
    
    # find the shortest sequence for initial substring comparisons
    k.sort(key=len)
    base_seq = k[0]
    other_seqs = k[1:]
    
    # initialize variables for storing the longest common substring
    lcs = ""
    max_len = 0

    # i understand this is not the most efficient algorithm (O(n^2), this was the way it was demonstrated in BINF 527 so i used that model
    for i in range(len(base_seq)):
        for j in range(i + 1, len(base_seq) + 1):
            # store substring
            substr = base_seq[i:j]
            
            # check if substring is common in all other sequences
            if all(substr in seq for seq in other_seqs):
                # update ONLY if loop finds a longer common substring
                if len(substr) > max_len:
                    lcs = substr
                    max_len = len(substr)
    
    return lcs

#tester for longest common substring
seqs = ["GATTACA", "TAGACCA", "ATACA"]
output5 = LongestCommonSubstring(seqs)
print(output5)


#task 6 -------
def FindingSubsequence(s: str, t: str) -> List[int]:
    # verify input: check if s and t contain only 'C', 'A', 'G', or 'T' characters
    valid_chars = {'C', 'A', 'G', 'T'}
    if not all(char in valid_chars for char in s) or not all(char in valid_chars for char in t):
        raise ValueError("Strings must contain only 'C', 'A', 'G', or 'T' characters")
            
    # this solution returns a subset of the indices where the symbols appear; it only returns one solution as specified in the question
    indices = []
    t_index = 0
    
    # Traverse through `s` to find each character of `t`
    for s_index in range(len(s)):
        if t_index < len(t) and s[s_index] == t[t_index]:
            # Record the 1-based index
            indices.append(s_index + 1)
            # Move to the next character in `t`
            t_index += 1
        
        # Stop if we have found all characters of `t`
        if t_index == len(t):
            break
    
    return indices

#tester for finding subsequences 
output6 = FindingSubsequence("TATGCTAAGATC", "ACG")
print(output6)

#task 7 ------
def LongestCommonSubsequence(s: str, t: str) -> str:
    # verify input: check if s and t contain only 'C', 'A', 'G', or 'T' characters
    valid_chars = {'C', 'A', 'G', 'T'}
    if not all(char in valid_chars for char in s) or not all(char in valid_chars for char in t):
        raise ValueError("Strings must contain only 'C', 'A', 'G', or 'T' characters")

    #determine the lengths of s and t
    n, m = len(s), len(t) 
    #2D list with empty strings
    #will store lcsubsequence between s:i and t:j
    dp = [["" for _ in range(m + 1)] 
              for _ in range(n + 1)]
    #iterate through each character
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if s[i - 1] == t[j - 1]: #if characters at these positions match
                dp[i][j] = dp[i - 1][j - 1] + s[i - 1] #add char to subsequence ending @ that location
            else: #otherwise choose longer subsequence
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1], key=len)

    return dp[-1][-1] #bottom right cell will have lcsubsequence

#tester for longest common subsequence
output7 = LongestCommonSubsequence("AACTTGA", "ACACTGTGA")
print(output7)

#task 8 -------
def ShortestCommonSupersequence(s: str, t: str) -> str:
    # use previous function to use the longest common subsequence
    lcs = LongestCommonSubsequence(s, t)
    
    # work backwards from LCS to build the shortest common supersequence
    i = j = 0
    scs = []
    
    # traverese the LCS, interleaving characters from s and t
    for c in lcs:
        # add characters from s until we reach the next character in LCS
        while i < len(s) and s[i] != c:
            scs.append(s[i])
            i += 1
        # add characters from t until we reach the next character in LCS
        while j < len(t) and t[j] != c:
            scs.append(t[j])
            j += 1
        # add the LCS character itself and move both pointers
        scs.append(c)
        i += 1
        j += 1
    
    # append remaining characters in s and t
    scs.extend(s[i:])
    scs.extend(t[j:])
    
    return ''.join(scs)


#tester for shortest common supersequence
output8 = ShortestCommonSupersequence("ATCTGAT", "TGCATA")
print(output8)