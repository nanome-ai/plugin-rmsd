import numpy as np 
from nanome.util import Logs



# needleman wunsch algorithm
def global_align(res1,res2,gap_penalty = -1, mismatch_penalty = -1, match_reward = 2):
    
    seq1 = list(map(lambda res: res.type, res1.residues))
    seq2 = list(map(lambda res: res.type, res2.residues))
    res_list1 =list(res1.residues)
    res_list2 = list(res2.residues)
    Logs.debug("res_list1 is ",res_list1)
    Logs.debug("length of reslist1 is ",len(res_list1))
    #create the table
    m, n = len(seq1), len(seq2)
    score = np.zeros((m+1, n+1))      
    
    # file the first column and row
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i
    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j
    # fill the table wtih scores
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                match = score[i - 1][j - 1] + match_reward
            else:
                match = score[i - 1][j - 1] + mismatch_penalty
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert)
    Logs.debug("table filled")
    Logs.debug(score)
    # Traceback and compute the alignment 
    align1, align2 = '', ''
    # start from the bottom right cell
    i,j = m,n
    while i > 0 and j > 0: 
        Logs.debug(i," ",j)

        # Logs.debug("i and j are: ",i," ",j)
        score_current = score[i][j]
        score_diagonal = score[i-1][j-1]
        score_up = score[i][j-1]
        score_left = score[i-1][j]
        # two residuses match, only deselect when the selected atoms are not matched
        if score_current == score_diagonal + match_reward and seq1[i-1] == seq2[j-1]:
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        # two of the residues do not match, the deselect both
        elif score_current == score_diagonal + mismatch_penalty and seq1[i-1] != seq2[j-1]:
            for x in res_list1[i-1].atoms:
                x.selected = False
            for y in res_list2[j-1].atoms:
                y.selected = False
            i -= 1
            j -= 1
        # seq1 has an extra residue, deselect it
        elif score_current == score_left + gap_penalty:
            align1 += seq1[i-1]
            align2 += '---'
            for x in res_list1[i-1].atoms:
                x.selected = False
            i -= 1
        # seq2 has an extra residue, deselect it
        elif score_current == score_up + gap_penalty:
            align1 += '---'
            align2 += seq2[j-1]
            for x in res_list2[j-1].atoms:
                x.selected = False
            j -= 1
    Logs.debug("traceback done1")
    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq1[i-1]
        align2 += '---'
        for x in res_list1[i-1].atoms:
            x.selected = False
        i -= 1
    while j > 0:
        align1 += '---'
        align2 += seq2[j-1]
        for x in res_list2[j-1].atoms:
            x.selected = False
        j -= 1
    Logs.debug("traceback done2")
    Logs.debug("align1 is ",align1)
    Logs.debug("align2 is ",align2)
    # finalize(align1, align2)

    return res1,res2
def local_align(plugin):
    pass