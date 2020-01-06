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
    rest_types1 = list(map(lambda res: res.type, selected_res1))
    rest_types2 = list(map(lambda res: res.type, selected_res2))

    # run the "smart occupancy selection method" on the residue lists of both complexes
    res_list1 =list(map(lambda a:select_occupancy(a),selected_res1))
    res_list2 =list(map(lambda a:select_occupancy(a),selected_res2))

    # create the table of global alignment
    m, n = len(rest_types1), len(rest_types2)
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
            if rest_types1[i-1] == rest_types2[ j-1]:
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
           rest_types1[i-1] == rest_types2[j-1] and rest_types1[i-1] != 'UNK' and rest_types2[j-1] != 'UNK':
            # align1 += rest_types1[i-1]
            # align2 += rest_types2[j-1]
            # final1 += rest_types1[i-1]
            # final2 += rest_types2[j-1]
            # clustalW_score += match_reward
            match1=list(map(lambda a:a.selected,res_list1[i-1].atoms))
            match2=list(map(lambda a:a.selected,res_list2[j-1].atoms))
            
            if match1 != match2 and not only_score:
                
                for x in res_list1[i-1].atoms:
                        x.selected = False

                for x in res_list2[j-1].atoms:
                        x.selected = False
            else:
                align1 += rest_types1[i-1]
                align2 += rest_types2[j-1]
                final1 += rest_types1[i-1]
                final2 += rest_types2[j-1]
                clustalW_score += match_reward
                match_count += 1
            i -= 1
            j -= 1

        # two of the residues do not match, deselect both
        elif score_current == score_diagonal + mismatch_penalty and \
             rest_types1[i-1] != rest_types2[j-1] or (rest_types1[i-1] == 'UNK' and rest_types2[j-1] == 'UNK'):
            if not only_score:
                for x in res_list1[i-1].atoms:
                    x.selected = False
                for y in res_list2[j-1].atoms:
                    y.selected = False
            clustalW_score += mismatch_penalty
            i -= 1
            j -= 1
            
        # rest_types1 has an extra residue, deselect it
        elif score_current == score_left + gap_penalty:
            align1 += rest_types1[i-1]
            align2 += '---'
            if not only_score:
                for x in res_list1[i-1].atoms:
                    x.selected = False
            clustalW_score += gap_penalty
            i -= 1

        # rest_types2 has an extra residue, deselect it
        elif score_current == score_up + gap_penalty:
            align1 += '---'
            align2 += rest_types2[j-1]
            if not only_score:
                for x in res_list2[j-1].atoms:
                    x.selected = False
            clustalW_score += gap_penalty
            j -= 1

    # Finish tracing up to the top left cell
    while i > 0:
        align1 += rest_types1[i-1]
        align2 += '---'
        if not only_score:
            for x in res_list1[i-1].atoms:
                x.selected = False
        clustalW_score += gap_penalty
        i -= 1
    while j > 0:
        align1 += '---'
        align2 += rest_types2[j-1]
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


def local_align(complex1,complex2,gap_penalty = -2, mismatch_penalty = -1, match_reward = 3, only_score = False):
    match_count = 0
    clustalW_score = 0
    selected_res1 = selected_res(complex1)
    selected_res2 = selected_res(complex2)
    max_cell = [0,0]
    max_cell_value = 0
    stop_traceback = 0

    # list of residues type of the complex
    rest_types1 = list(map(lambda res: res.type, selected_res1))
    rest_types2 = list(map(lambda res: res.type, selected_res2))

    # run the "smart occupancy selection method" on the residue lists of both complexes
    res_list1 =list(map(lambda a:select_occupancy(a),selected_res1))
    res_list2 =list(map(lambda a:select_occupancy(a),selected_res2))

    # create the table of global alignment
    m, n = len(rest_types1), len(rest_types2)
    shorter_len = min(m,n)
    score = np.zeros((m+1, n+1))      
    
    # file the first column and first row of the table
    for i in range(0, m + 1):
        score[i][0] = 0
    for j in range(0, n + 1):
        score[0][j] = 0

    # fill the table wtih scores
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if rest_types1[i-1] == rest_types2[ j-1]:
                match = score[i - 1][j - 1] + match_reward
            else:
                match = score[i - 1][j - 1] + mismatch_penalty
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert,0)

            if score[i][j] > max_cell_value:
                max_cell_value = score[i][j]
                max_cell = [i,j]
   

    # Traceback and compute the alignment 
    # aligns are the output sequences with gaps (both delete and insert)
    # finals are the output sequences that should be the same after the sequence alignment
    align1, align2 = '', ''
    final1, final2 = '', ''

    i,j = m,n

    while i > max_cell[0]:
        # align1 += rest_types1[i-1]
        # align2 += '---'
        if not only_score:
            for x in res_list1[i-1].atoms:
                x.selected = False
        clustalW_score += gap_penalty
        i -= 1
    while j > max_cell[1]:
        # align1 += '---'
        # align2 += rest_types2[j-1]
        if not only_score:
            for x in res_list2[j-1].atoms:
                x.selected = False
        clustalW_score += gap_penalty
        j -= 1
    


    # start from the bottom right cell
    i,j = max_cell

    # go left and up until touching the 1st row/column
    while i > 0 and j > 0: 
        score_current = score[i][j]
        score_diagonal = score[i-1][j-1]
        score_up = score[i][j-1]
        score_left = score[i-1][j]
        # two residuses match, only deselect when the selected atoms don't match (problem in the pdb file)
        if score_current == 0:
            stop_traceback = 1
            break
        if score_current == score_diagonal + match_reward and \
           rest_types1[i-1] == rest_types2[j-1] and rest_types1[i-1] != 'UNK' and rest_types2[j-1] != 'UNK':
           
            match1=list(map(lambda a:a.selected,res_list1[i-1].atoms))
            match2=list(map(lambda a:a.selected,res_list2[j-1].atoms))
            
            if match1 != match2 and not only_score:
                
                for x in res_list1[i-1].atoms:
                        x.selected = False

                for x in res_list2[j-1].atoms:
                        x.selected = False
            else:

                align1 += rest_types1[i-1]
                align2 += rest_types2[j-1]
                final1 += rest_types1[i-1]
                final2 += rest_types2[j-1]
                clustalW_score += match_reward
                match_count += 1
            i -= 1
            j -= 1
           
        # two of the residues do not match, deselect both
        elif score_current == score_diagonal + mismatch_penalty and \
             rest_types1[i-1] != rest_types2[j-1] or (rest_types1[i-1] == 'UNK' and rest_types2[j-1] == 'UNK'):
            if not only_score:
                for x in res_list1[i-1].atoms:
                    x.selected = False
                for y in res_list2[j-1].atoms:
                    y.selected = False
            clustalW_score += mismatch_penalty
            i -= 1
            j -= 1
           
        # rest_types1 has an extra residue, deselect it
        elif score_current == score_left + gap_penalty:
            align1 += rest_types1[i-1]
            align2 += '---'
            if not only_score:
                for x in res_list1[i-1].atoms:
                    x.selected = False
            clustalW_score += gap_penalty
            i -= 1
            
        # rest_types2 has an extra residue, deselect it
        elif score_current == score_up + gap_penalty:
            align1 += '---'
            align2 += rest_types2[j-1]
            if not only_score:
                for x in res_list2[j-1].atoms:
                    x.selected = False
            clustalW_score += gap_penalty
            j -= 1
          
    # Finish tracing up to the top left cell
    while i > 0:
        align1 += rest_types1[i-1]
        align2 += '---'
        if not only_score:
            for x in res_list1[i-1].atoms:
                x.selected = False
        clustalW_score += gap_penalty
        i -= 1
    while j > 0:
        align1 += '---'
        align2 += rest_types2[j-1]
        if not only_score:
            for x in res_list2[j-1].atoms:
                x.selected = False
        clustalW_score += gap_penalty
        j -= 1
    Logs.debug("final1 is ",final1)
    Logs.debug("final2 is ",final2)
    
    # return complex1,complex2
    # return clustalW_score
    if shorter_len != 0:
        rt = 1-(match_count/shorter_len)
    else:
        rt = 0
        Logs.debug("one of the complexes has no atom selected")
    return rt

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



