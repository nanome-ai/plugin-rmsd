import numpy as np 
from nanome.util import Logs
from math import ceil
from itertools import product

# needleman wunsch algorithm
def global_align(complex1,complex2,gap_penalty = -1, mismatch_penalty = -1, match_reward = 3):

    selected_res1 = selected_res(complex1)
    selected_res2 = selected_res(complex2)

    # list of residues type of the complex
    seq1 = list(map(lambda res: res.type, selected_res1))
    seq2 = list(map(lambda res: res.type, selected_res2))

    # run the "smart occupancy selection method" on the residue lists of both complexes
    res_list1 =list(map(lambda a:select_occupancy(a),selected_res1))
    res_list2 =list(map(lambda a:select_occupancy(a),selected_res2))
    
    # Logs.debug("res_list1 is ",res_list1)
    Logs.debug("length of reslist1 is ",len(res_list1))
    Logs.debug("length of reslist2 is ",len(res_list2))

    # create the table of global alignment
    m, n = len(seq1), len(seq2)
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
    Logs.debug("table filled")
    # Logs.debug(score)
    # Traceback and compute the alignment 
    align1, align2 = '', ''
    final1, final2 = '', ''
    # start from the bottom right cell
    i,j = m,n
    path=[]
    while i > 0 and j > 0: 
        score_current = score[i][j]
        score_diagonal = score[i-1][j-1]
        score_up = score[i][j-1]
        score_left = score[i-1][j]
        # two residuses match, only deselect when the selected atoms are not matched
        if score_current == score_diagonal + match_reward and \
           seq1[i-1] == seq2[j-1] and seq1[i-1] != 'UNK' and seq2[j-1] != 'UNK':
            path.append("match "+seq1[i-1]+" "+str(score_current)+" ["+str(i)+", "+str(j)+"]")
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            final1 += seq1[i-1]
            final2 += seq2[j-1]
            # Logs.debug(seq1[i-1],' ',seq2[j-1])
            match1=list(map(lambda a:a.selected,res_list1[i-1].atoms))
            match2=list(map(lambda a:a.selected,res_list2[j-1].atoms))
            if match1 != match2:
                # Logs.debug("different atoms in the residues, please select backbone_only")
                # Logs.debug("The two residues are not same: ")
                # Logs.debug("type 1 is ",res_list1[i-1].type)
                # Logs.debug("type 2 is ",res_list2[j-1].type)
                # Logs.debug("atoms 1:   ",list(map(lambda a:a.name,res_list1[i-1].atoms)))
                # Logs.debug("occupancy: ",list(map(lambda a:a._occupancy,res_list1[i-1].atoms)))
                # Logs.debug("selected:  ",list(map(lambda a:a.selected,res_list1[i-1].atoms)))
                # Logs.debug("serial:    ",res_list1[i-1].serial)
                # Logs.debug("atoms 2:   ",list(map(lambda a:a.name,res_list2[j-1].atoms)))
                # Logs.debug("occupancy  ",list(map(lambda a:a._occupancy,res_list2[j-1].atoms)))
                # Logs.debug("selected:  ",list(map(lambda a:a.selected,res_list2[j-1].atoms)))
                # Logs.debug("serial:    ",res_list2[j-1].serial)
                
                for x in res_list1[i-1].atoms:
                        x.selected = False

                for x in res_list2[j-1].atoms:
                        x.selected = False

            i -= 1
            j -= 1

        # two of the residues do not match, the deselect both
        elif score_current == score_diagonal + mismatch_penalty and \
             seq1[i-1] != seq2[j-1] or (seq1[i-1] == 'UNK' and seq2[j-1] == 'UNK'):
            path.append("mismatch "+seq1[i-1]+" "+seq2[j-1]+" "+str(score_current)+" ["+str(i)+", "+str(j)+"]")
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
            path.append("gap align1 "+str(score_current)+" ["+str(i)+", "+str(j)+"]")
            for x in res_list1[i-1].atoms:
                x.selected = False
            i -= 1
        # seq2 has an extra residue, deselect it
        elif score_current == score_up + gap_penalty:
            align1 += '---'
            align2 += seq2[j-1]
            path.append("gap align2 "+str(score_current)+" ["+str(i)+", "+str(j)+"]")
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
    Logs.debug("len of align1 is ",len(align1))
    Logs.debug("len of align2 is ",len(align2))
    Logs.debug("final1 is ",final1)
    Logs.debug("final2 is ",final2)
    Logs.debug("len of final1 is ", len(final1))
    Logs.debug("len of final2 is ",len(final2))
    f= open("psa.txt","w")
    for x in path:
        f.write(x)
        f.write("\n")
    f.close()

    return complex1,complex2

def local_align(mobile,target):
    pass

# This function align multiple sequences. It uses dynamic programming method 
# and make an n-dimensional scoring matrix
def multi_global_align(complexes,gap_penalty = -1, mismatch_penalty = -1, match_reward = 3):
    Logs.debug("complexes are: ",complexes)
    
    # for every complex, only select the residues whose atoms are all selected
    selected_res_list = []
    for x in complexes:
        selected_res_list.append(selected_res(x))

    # a list of the lists of residues type of the complex
    seq_list = []
    for x in selected_res_list:
        seq_list.append(list(map(lambda res:res.type,x)))
    # run the "smart occupancy selection method" on the residue lists of both complexes
    res_lists = []
    for x in selected_res_list:
        res_lists.append(list(map(lambda a:select_occupancy(a),x)))
    
    for i,x in enumerate(res_lists):
        Logs.debug("length of reslist ",i," is ",len(x))

    # m, n = len(seq1), len(seq2)
    # score = np.zeros((m+1, n+1))      
    # create the table of global alignment
    len_list = []
    table_size = ()
    for x in seq_list:
        len_list.append(len(x))
        table_size = table_size + (len(x)+1,)
    score = np.zeros(table_size)
    # i is the dimension index, x is the dimension length
    for i,x in enumerate(table_size):
        # j is the index in one dimension
        for j in range(1,x):
            # create a tuple for index
            temp_index = np.zeros(len(table_size), dtype=int)
            temp_index[i] = j
            temp_index = tuple(temp_index)
            score[temp_index] = gap_penalty * j

    # list of indices in all dimensions
    score_index = np.zeros(len(table_size),dtype = int)
    reward_penalty = [gap_penalty, mismatch_penalty, match_reward]
    fill_score(score, table_size, len(table_size)-1, score_index  ,reward_penalty,seq_list)
    Logs.debug("scrore table size is: ",score.shape)
    # Traceback and compute the alignment  
    aligns = []
    finals = []
    for x in range(len(table_size)):
        aligns.append('')
        finals.append('')
        
    # start from the bottom right cell
    trace = [int(x-1) for x in table_size]
    path=[]
    while all(x > 0 for x in trace): 
        diag_trace = [x-1 for x in trace]
        score_current = score[tuple(trace)]
        score_diagonal = score[tuple(diag_trace)]
        gap_points = get_gap_points(trace) 
        # Logs.debug("score size", score.shape)
        # Logs.debug(trace)
        # Logs.debug([x[::-1] for x in gap_points])
        gap_scores = [score[tuple(x)]+reward_penalty[0] for x in gap_points]

        # two residuses match, only deselect when the selected atoms are not matched
        seqs_val = []
        for i,x in enumerate(diag_trace):
            seqs_val.append(seq_list[i][x])
        if score_current == score_diagonal + match_reward and \
            all(elem == seqs_val[0] for elem in seqs_val) and all(elem != 'UNK' for elem in seqs_val):
            path.append("match "+seqs_val[0]+" "+str(score_current)+" "+str(trace))
            for i, x in enumerate(aligns):
                aligns[i] += seqs_val[i]
                finals[i] += seqs_val[i]
            matches = []
            for i, x in enumerate(diag_trace):
                matches.append(list(map(lambda a:a.selected,res_lists[i][x].atoms)))
            if not all(elem == matches[0] for elem in matches):    
                # Logs.debug("The two residues are not same: ")
                # Logs.debug("type 1 is ",res_lists[0][trace[0]-1].type)
                # Logs.debug("type 2 is ",res_lists[1][trace[1]-1].type)
                # Logs.debug("atoms 1:   ",list(map(lambda a:a.name,res_lists[0][trace[0]-1].atoms)))
                # Logs.debug("occupancy: ",list(map(lambda a:a._occupancy,res_lists[0][trace[0]-1].atoms)))
                # Logs.debug("selected:  ",list(map(lambda a:a.selected,res_lists[0][trace[0]-1].atoms)))
                # Logs.debug("serial:    ",res_lists[0][trace[0]-1].serial)
                # Logs.debug("atoms 2:   ",list(map(lambda a:a.name,res_lists[1][trace[1]-1].atoms)))
                # Logs.debug("occupancy  ",list(map(lambda a:a._occupancy,res_lists[1][trace[1]-1].atoms)))
                # Logs.debug("selected:  ",list(map(lambda a:a.selected,res_lists[1][trace[1]-1].atoms)))
                # Logs.debug("serial:    ",res_lists[1][trace[1]-1].serial)
                for i,y in enumerate(res_lists):
                    for x in y[diag_trace[i]].atoms:
                        x.selected = False
            for i in range(len(trace)):
                trace[i] -= 1

            

        # two of the residues do not match, the deselect both
        elif score_current == score_diagonal + mismatch_penalty and \
             all(elem != seqs_val[0] for elem in seqs_val[1:]) or all(elem=='UNK' for elem in seqs_val):
            path.append("mismatch "+seqs_val[0]+" "+seqs_val[1]+" "+str(score_current)+" "+str(trace))
            for i,y in enumerate(res_lists):
                for x in y[diag_trace[i]].atoms:
                    x.selected = False
                
            for i in range(len(trace)):
                trace[i] -= 1
            
        # seq1 has an extra residue, deselect it
        elif score_current in gap_scores:
            # find the dimension that's been chosen
            # loop through all the dimensions and if it has gap, "---"
            temp=trace[:]
            temp2=0
            for i,x in enumerate(gap_scores):
                if score_current == x:
                    aligns[i] += seq_list[i][diag_trace[i]]
                    for x in res_lists[i][diag_trace[i]].atoms:
                        x.selected = False
                    trace[i] -= 1
                    temp2=i
                else:
                    aligns[i] += '---'
            path.append("gap align"+str(temp2+1)+" "+str(score_current)+" "+str(temp))

            # for i in range(len(trace)):
            #     trace[i] -= 1


    Logs.debug("traceback done1")
    # Finish tracing up to the top left cell
    for x in range(len(trace)):
        while trace[x] > 0:
            Logs.debug("trace[x] is ",trace[x])
            Logs.debug("len of seqlist[x] is ",len(seq_list[x]))
            aligns[x] += seq_list[x][trace[x]]
            for y in res_lists[x][trace[x]].atoms:
                y.selected = False
            for z in [y for y in range(len(trace)) if y != x] :
                aligns[z] += '---'
            
            trace[x] -= 1
    for x in range(len(aligns)):
        Logs.debug(aligns[x])
        Logs.debug("len of align",x+1," is ",len(aligns[x]))
    for x in range(len(finals)):
        Logs.debug(finals[x])
        Logs.debug("len of final",x+1," is ",len(finals[x]))
    f= open("msa.txt","w")
    for x in path:
        f.write(x)
        f.write("\n")
    f.close()
    return complexes


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
        # Logs.debug("name is   ",p)
        # Logs.debug("occupancy ",occ_dict[p][1])
        # Logs.debug("top_n is  ",top_n)
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
        # Logs.debug("loop indices is: ",loop_indices)
        for x in range(len(table)):
            # Logs.debug("loop indices is: ",loop_indices[x])
            # Logs.debug("seq list shape is: ",len(seq_list)," ",len(seq_list[0])," ",len(seq_list[1]))
            # Logs.debug("seq list is ",len(seq_list[-x-1]))
            # Logs.debug("loop idx is ",loop_indices[-x-1])
            # Logs.debug("-x-1 is ",-x-1)
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
# Ex. input：(i,j,k)
#     output:[(i,j,k-1), (i,j-1,k),(i,j-1,k-1),(i-1,j,k-1),(i-1,j-1,k)]
def get_gap_points(ijk):
    for x in ijk:
        if x == 0:
            Logs.debug("point on the side or at the origin")
            return
    match_index = list(ijk)
    match_index_diag = [x-1 for x in match_index]
    comb_list = []
    for x in range(len(match_index)):
        comb_list.append([match_index[x],match_index_diag[x]])
    # gap_points = list(product(*comb_list))
    # gap_points.remove(ijk)
    # gap_points.remove(tuple(match_index_diag))
    
    # use manhattan distance instead
    gap_points = []
    for x in range(len(ijk)):
        temp = match_index[:]
        temp[x] = match_index_diag[x]
        gap_points.append(temp)
    return gap_points

