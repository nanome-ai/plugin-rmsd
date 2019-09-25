import numpy as np 
from nanome.util import Logs
from math import ceil
from itertools import product

# needleman wunsch algorithm
# the param only_score was used for clustalW
def global_align(complex1,complex2,gap_penalty = -1, mismatch_penalty = 0, match_reward = 3, only_score = False):
    match_count = 0
    clustalW_score = 0
    selected_res1 = selected_res(complex1)
    selected_res2 = selected_res(complex2)

    # list of residues type of the complex
    seq1 = list(map(lambda res: res.type, selected_res1))
    seq2 = list(map(lambda res: res.type, selected_res2))

    # run the "smart occupancy selection method" on the residue lists of both complexes
    res_list1 =list(map(lambda a:select_occupancy(a),selected_res1))
    res_list2 =list(map(lambda a:select_occupancy(a),selected_res2))

    # create the table of global alignment
    m, n = len(seq1), len(seq2)
    shorter_len = min(m,n)
    score = np.zeros((m+1, n+1))      
    
    # file the first column and first row of the table
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i
    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j

    # fill the table wtih scores
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[ j-1]:
                match = score[i - 1][j - 1] + match_reward
            else:
                match = score[i - 1][j - 1] + mismatch_penalty
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert)

    # Traceback and compute the alignment 
    # aligns are the output sequences with gaps (both delete and insert)
    # finals are the output sequences that should be the same after the sequence alignment
    align1, align2 = '', ''
    final1, final2 = '', ''

    # start from the bottom right cell
    i,j = m,n

    # go left and up until touching the 1st row/column
    while i > 0 and j > 0: 
        score_current = score[i][j]
        score_diagonal = score[i-1][j-1]
        score_up = score[i][j-1]
        score_left = score[i-1][j]
        # two residuses match, only deselect when the selected atoms don't match (problem in the pdb file)
        if score_current == score_diagonal + match_reward and \
           seq1[i-1] == seq2[j-1] and seq1[i-1] != 'UNK' and seq2[j-1] != 'UNK':
            # align1 += seq1[i-1]
            # align2 += seq2[j-1]
            # final1 += seq1[i-1]
            # final2 += seq2[j-1]
            # clustalW_score += match_reward
            match1=list(map(lambda a:a.selected,res_list1[i-1].atoms))
            match2=list(map(lambda a:a.selected,res_list2[j-1].atoms))
            
            if match1 != match2 and not only_score:
                
                for x in res_list1[i-1].atoms:
                        x.selected = False

                for x in res_list2[j-1].atoms:
                        x.selected = False
            else:
                align1 += seq1[i-1]
                align2 += seq2[j-1]
                final1 += seq1[i-1]
                final2 += seq2[j-1]
                clustalW_score += match_reward
                match_count += 1
            i -= 1
            j -= 1

        # two of the residues do not match, deselect both
        elif score_current == score_diagonal + mismatch_penalty and \
             seq1[i-1] != seq2[j-1] or (seq1[i-1] == 'UNK' and seq2[j-1] == 'UNK'):
            if not only_score:
                for x in res_list1[i-1].atoms:
                    x.selected = False
                for y in res_list2[j-1].atoms:
                    y.selected = False
            clustalW_score += mismatch_penalty
            i -= 1
            j -= 1
            
        # seq1 has an extra residue, deselect it
        elif score_current == score_left + gap_penalty:
            align1 += seq1[i-1]
            align2 += '---'
            if not only_score:
                for x in res_list1[i-1].atoms:
                    x.selected = False
            clustalW_score += gap_penalty
            i -= 1

        # seq2 has an extra residue, deselect it
        elif score_current == score_up + gap_penalty:
            align1 += '---'
            align2 += seq2[j-1]
            if not only_score:
                for x in res_list2[j-1].atoms:
                    x.selected = False
            clustalW_score += gap_penalty
            j -= 1

    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq1[i-1]
        align2 += '---'
        if not only_score:
            for x in res_list1[i-1].atoms:
                x.selected = False
        clustalW_score += gap_penalty
        i -= 1
    while j > 0:
        align1 += '---'
        align2 += seq2[j-1]
        if not only_score:
            for x in res_list2[j-1].atoms:
                x.selected = False
        clustalW_score += gap_penalty
        j -= 1
    
    # return complex1,complex2
    # return clustalW_score
    if shorter_len != 0:
        rt = 1-(match_count/shorter_len)
    else:
        rt = 0
        Logs.debug("one of the complexes has no atom selected")
    return rt


def local_align(mobile,target):
    pass

# use the online server of t-coffee to get the multiple sequence alignment result.
# COMMENTED OUT BECAUSE VINCENT SAID WE CAN'T USE ANOTHER WEBSITE'S FUNCTION YET
# http://tcoffee.crg.cat/apps/tcoffee/do:regular
# input: a list of complexes
# output: the result of the t-coffee alignment, including the matching, mismatching, and the gaps
# def tcoffee(complexes):
#     letter_code = {"ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "ASX":"B", "CYS":"C",
#                    "GLU":"E", "GLN":"Q", "GLX":"Z", "GLY":"G", "HIS":"H", "ILE":"I",
#                    "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P", "SER":"S",
#                    "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V",}
    # # get the sequence of residues in all the complexes in the input
    # three_letter_list={}
    # for x in complexes:

    #     three_letter_list[x.name] = list(map(lambda res:res.type,x))
    # # deal with all the occupancy and un-full-selected residues things
    # # 1. get the one letter symbol of all the complexes, and format them
    # for x in three_letter_list:
    #     for j,y in enumerate(x):
    #         if y in letter_code:
    #             three_letter_list[x][j] = letter_code[y]
    # # paste the fasta sequences into the textbox
    # # click "submit" button on the webpage. (and follow it to the new webpage?)
    # # download the output file from the new webpage
    # # parse the result 
    # # deselect the non-common parts


# takes in a single residue
def select_occupancy(residue):
    occ_dict = {}
    for a in residue.atoms:
        if a._occupancy < 1:
            name = a.name
            if name in occ_dict:
                occ_dict[name][0].append(a)
                occ_dict[name][1].append(a._occupancy)
            else:
                occ_dict[name]=[[a],[a._occupancy]]

    for p in occ_dict:
        top_n = round(sum(occ_dict[p][1]))
        occ_dict[p][0].sort(key=lambda x: x._occupancy, reverse=True)
        occ_dict[p][0] = occ_dict[p][0][top_n:]
        for a in occ_dict[p][0]:
            a.selected = False

    return residue

# select the residues whose atoms are all selected.
def selected_res(complexes):
    residues = list(map(lambda a:a,complexes.residues))
    rt = []
    # if there's an unselected atom in the residue, don't include it in the list
    for residue in residues:
        selected_bool = True
        for atom in residue.atoms:
            if atom.selected == False:
                selected_bool = False
        if selected_bool:
            rt.append(residue)
    return rt

# fill the nd-matrix "score" with the scores recursively
# table is table_size, dim is the current dimension index
def fill_score(score, table, dim, loop_indices, reward_penalty,seq_list):
    if dim < 0:
        # compare all the sequences
        match_list = []
        for x in range(len(table)):
            match_list.append(seq_list[-x-1][loop_indices[-x-1]-1])

        match_index = list(loop_indices)
        match_index_diag = [x-1 for x in match_index]

        if all(elem == match_list[0] for elem in match_list):
            match = score[tuple(match_index_diag)] + reward_penalty[2]
        else:
            match = score[tuple(match_index_diag)] + reward_penalty[1]
       
        # question: if different number of dimensions don't match, are the gap penalties different?
        gap_points = get_gap_points(loop_indices)
        gap_results = [score[tuple(x)]+reward_penalty[0] for x in gap_points]
        score[tuple(loop_indices)] = max(match,max(gap_results))
        return
    
    else:
        for x in range(table[-dim-1]-1):
            loop_indices[-dim-1] += 1
            fill_score(score,table,dim-1,loop_indices, reward_penalty,seq_list) 
           
        loop_indices[-dim-1] = 0

# take in a point and return all the gap points
# used in multiple_global_align
# Ex. input：(i,j,k)
#     output:[(i,j,k-1), (i,j-1,k),(i,j-1,k-1),(i-1,j,k-1),(i-1,j-1,k)]
def get_gap_points(ijk):
    for x in ijk:
        if x == 0:
            # Logs.debug("point on the side or at the origin")
            return
    match_index = list(ijk)
    match_index_diag = [x-1 for x in match_index]
    comb_list = []
    for x in range(len(match_index)):
        comb_list.append([match_index[x],match_index_diag[x]])
    
    # use manhattan distance instead
    gap_points = []
    for x in range(len(ijk)):
        temp = match_index[:]
        temp[x] = match_index_diag[x]
        gap_points.append(temp)
    return gap_points

