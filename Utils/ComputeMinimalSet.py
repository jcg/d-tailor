'''
Created on Apr 9, 2012

@author: cambray

This script uses a monte carlo simulation to compute a set of sequences (with exactly one element for each combination) having minimal pair-wise distance between all of them.

INPUT: 
    - CSV file with all solutions
OUTPUT:
    - CSV file with final set of sequences selected
    - FASTA file with final set of sequences selected
    - CSV file containing the nucleotide distance matrix between all final sequences    

'''

import pickle
import csv
from random import randint, choice
from math import sqrt
from os import path
import sys

code = {'ata':'i', 'atc':'i', 'att':'i', 'atg':'m',
        'aca':'t', 'acc':'t', 'acg':'t', 'act':'t',
        'aac':'n', 'aat':'n', 'aaa':'k', 'aag':'k',
        'agc':'s', 'agt':'s', 'aga':'r', 'agg':'r',
        'cta':'l', 'ctc':'l', 'ctg':'l', 'ctt':'l',
        'cca':'p', 'ccc':'p', 'ccg':'p', 'cct':'p',
        'cac':'h', 'cat':'h', 'caa':'q', 'cag':'q',
        'cga':'r', 'cgc':'r', 'cgg':'r', 'cgt':'r',
        'gta':'v', 'gtc':'v', 'gtg':'v', 'gtt':'v',
        'gca':'a', 'gcc':'a', 'gcg':'a', 'gct':'a',
        'gac':'d', 'gat':'d', 'gaa':'e', 'gag':'e',
        'gga':'g', 'ggc':'g', 'ggg':'g', 'ggt':'g',
        'tca':'s', 'tcc':'s', 'tcg':'s', 'tct':'s',
        'ttc':'f', 'ttt':'f', 'tta':'l', 'ttg':'l',
        'tac':'y', 'tat':'y', 'taa':'*', 'tag':'*',
        'tgc':'c', 'tgt':'c', 'tga':'*', 'tgg':'w'}

def translate(seq):
    return "".join([code[seq[i:i+3]] for i in range(0,len(seq),3)])


#################
## Select the best set of sequences

def monte_carlo_min_dist(in_seq, pkl_path, pkl_import=False, random_period=5000, focus_period=500,convergence=1000):
    """
    Use monte-carlo scheme to find a set of sequence with minimal overall distance (nucleotide level)
    - in_seq : csv file as ouput by DB2CSV
    - pkl_path : path to read/write pkl
    - pkl_import : wether to import an existing pkl to start from
    - random_period : nbr of iteration when change are completely random
    - focus_period : nbr of iteration when changes are focuse on the 20% worst ids
    """
    seqs_nt = {}
    idnbr = {}
    h1 = csv.DictReader(open(in_seq))
    in_ids = []
    for data in h1:
        #print data
        if idnbr.has_key(data["des_solution_id"]):
            idnbr[data["des_solution_id"]] += 1
        else:
            idnbr[data["des_solution_id"]] = 1
        in_ids.append(data["des_solution_id"])
        data["des_solution_id"] = "%s_%i" % (data["des_solution_id"], idnbr[data["des_solution_id"]])
        seqs_nt[data["des_solution_id"]] = data["sequence"]

    #init set
    if pkl_import:
        current_set = pickle.load(open(pkl_import,'rb'))
        ids = current_set.keys()
        # check if new combinations
        if len(ids)<len(in_ids):
            for id_comb in in_ids:
                if not id_comb in ids:
                    current_set[id_comb]= randint(1, idnbr[id_comb])
    else:
        current_set = {}
        for combi in idnbr:
            # choose one sequence per combination
            jackpot = randint(1, idnbr[combi])
            current_set[combi] = jackpot
    #init distance
    distnt = {}
    idsumnt = {}
    total_dist_nt = 0
    ids = current_set.keys()
    allfullids = []
    for id_comb in ids:
        if idnbr[id_comb] > 1:
            allfullids.extend([(id_comb, nbr) for nbr in range(1, idnbr[id_comb]+1)])
    ids.sort()
    for i in range(len(ids)):
        fullid1= "%s_%s" % (ids[i], current_set[ids[i]])
        distnt[ids[i]] = {}
        for j in range(i+1, len(ids)):
            fullid2= "%s_%s" % (ids[j], current_set[ids[j]])
            distnt[ids[i]][ids[j]] = sum([1 for k in range(len(seqs_nt[fullid1])) if seqs_nt[fullid1][k]!=seqs_nt[fullid2][k]])        
            total_dist_nt += distnt[ids[i]][ids[j]]
            for id_comb in (ids[i],ids[j]):
                if idsumnt.has_key(id_comb):
                    idsumnt[id_comb]+=distnt[ids[i]][ids[j]]
                else:
                    idsumnt[id_comb]=distnt[ids[i]][ids[j]]
    totaltry    = 0
    currenttry  = 0
    toggle_pkl  = 0
    change_flag = True
    new_total_dist_nt = total_dist_nt
    new_idsumnt = idsumnt
    go_on = True
    while go_on:
        totaltry+=1
        currenttry+=1
        if change_flag:
            print "stop random"
            change_flag = False
            if not focus_period:
                random_flag = 1
                continue
            random_flag = -focus_period
            # identify max id_comb
            sumidnt = {}
            for id_comb in ids:
                if sumidnt.has_key(idsumnt[id_comb]):
                    sumidnt[idsumnt[id_comb]].append(id_comb)
                else:
                    sumidnt[idsumnt[id_comb]]=[id_comb]
            sums = sumidnt.keys()
            sums.sort(reverse=True)
            # get the % highest distance ids (those with several solutions available)
            maxnbr = 0.20 * len(ids)
            fullidsol = []
            for dist in sums:
                for id_comb in sumidnt[dist]:
                    fullidsol.extend([(id_comb, nbr) for nbr in range(1,idnbr[id_comb]+1) if (nbr != current_set[id_comb])])
                if len(fullidsol) >= maxnbr:
                    break
        #print fullidsol
        # choose a new sequence
        # equal probability to all solutions, irrespective of id_comb
        if random_flag < 0:
            combi, jackpot = fullidsol.pop(randint(0,len(fullidsol)-1))
            if not len(fullidsol):
                random_flag = 0
                currenttry  = 0
                print "start random"
            #print combi, jackpot
        elif random_flag > 0:
            combi, jackpot = choice(allfullids)
            while jackpot == current_set[combi]:
                combi, jackpot = choice(allfullids)
        else:
            combi, jackpot = choice(allfullids)
            while jackpot == current_set[combi]:
                combi, jackpot = choice(allfullids)
            print "start random"
            currenttry = 0
        # update status
        random_flag+=1
        new_id = "%s_%s" % (combi, jackpot)
        # recalculate distance
        new_distnt = {}
        for id_comb in ids:
            if id_comb != combi:
                fullid = "%s_%s" % (id_comb, current_set[id_comb])
                new_distnt[id_comb] = sum([1 for k in range(len(seqs_nt[fullid])) if seqs_nt[fullid][k]!=seqs_nt[new_id][k]])
                if  id_comb < combi:
                    new_total_dist_nt += new_distnt[id_comb] - distnt[id_comb][combi]
                    new_idsumnt[id_comb] += new_distnt[id_comb] - distnt[id_comb][combi]
                    new_idsumnt[combi] += new_distnt[id_comb] - distnt[id_comb][combi]
                else:
                    new_total_dist_nt += new_distnt[id_comb] - distnt[combi][id_comb]
                    new_idsumnt[id_comb] += new_distnt[id_comb] - distnt[combi][id_comb]
                    new_idsumnt[combi] += new_distnt[id_comb] - distnt[combi][id_comb]
        if (new_total_dist_nt < total_dist_nt):    
            print "%s (%s):\t%i\t%i" % (totaltry, currenttry, total_dist_nt, new_total_dist_nt-total_dist_nt)
            total_dist_nt = new_total_dist_nt
            idsumnt = new_idsumnt
            current_set[combi] = jackpot
            for id_comb in new_distnt.keys():
                if  id_comb < combi:
                    distnt[id_comb][combi] = new_distnt[id_comb]
                else:
                    distnt[combi][id_comb] = new_distnt[id_comb]
            output = open(pkl_path+str(toggle_pkl),'wb')
            toggle_pkl = abs(toggle_pkl-1)
            pickle.dump(current_set, output)
            output.close()
            currenttry = 0
        else:
            new_total_dist_nt = total_dist_nt
            new_idsumnt = idsumnt
        if random_flag > random_period:
            change_flag = True
        if convergence and currenttry>convergence:
            return
    return

def get_final_set_feats(seq_in, pkl, seq_out, get_distance=True, verbose=True):
    """
    filter initial dataset with selected data
    output minimal and updated dataset as csv
    if get_distance, also output the distance matrix between sequences (nt)
    statistics are for unique combinations of different sequences
    """
        
    seqs_nt = {}
    feats   = {}
    idnbr = {}
    current_set = {}
    h = csv.DictReader(open(seq_in))
    for data in h:
        if data.has_key(None):
            data.pop(None)
        if idnbr.has_key(data["des_solution_id"]):
            idnbr[data["des_solution_id"]] += 1
        else:
            idnbr[data["des_solution_id"]] = 1
        data["des_solution_id"] = "%s_%i" % (data["des_solution_id"], idnbr[data["des_solution_id"]])
        #data["sequence"] = data["sequence"][3:102]
        seqs_nt[data["des_solution_id"]] = data["sequence"]
        feats[data["des_solution_id"]] = data
        
    # load set
    if pkl:
        pkl_file = open(pkl,'rb')
        current_set = pickle.load(pkl_file)
        # modify data structure as needed ie if has not been extended
        if not type(current_set[current_set.keys()[0]]) is list:
            for id_comb in current_set:
                current_set[id_comb] = [current_set[id_comb]]
    else:
        current_set = {id_comb:range(1,idnbr[id_comb]+1) for id_comb in idnbr}
    #print current_set
    #
    ids = current_set.keys()
    ids.sort()
    total = 0.0
    if get_distance:
        hnt = open('%s_distance_matrix_nt.csv' % seq_out,'w')
        for id_comb in ids:
            for j in range(1,len(current_set[id_comb])+1):
                total+=1
                hnt.write(",%s_%s" % (id_comb, j))
        hnt.write("\n")
    #modify csv headers
    fields = ['final.id_comb', 'serie.id_comb']
    fields.extend(h.fieldnames)
    fields.insert(3,'repeat.id_comb')
    fields[4]='param.id_comb'
    hfasta = open('%s.fas' % seq_out,'w')
    hfeat = csv.DictWriter(open('%s_feats.csv' % seq_out,'w'), fieldnames=fields) 
    hfeat.writeheader()
    # calculate matrix
    alldistnt = []
    distnt = {}
    nbr   = 0.0
    total = total*(total-1)/2
    for i in range(len(ids)):
        for i2 in range(len(current_set[ids[i]])):
            finalid1 = "%s_%s" % (ids[i], i2+1)
            fullid1  = "%s_%s" % (ids[i], current_set[ids[i]][i2])
            #update infos for writting csv
            feats[fullid1].update({'final.id_comb':finalid1, 'repeat.id_comb':(i2+1),'param.id_comb':ids[i]})
            feats[fullid1].pop('des_solution_id')            
            hfeat.writerow(feats[fullid1])
            hfasta.write(">%s\n%s\n" % (finalid1, seqs_nt[fullid1].upper()))
            distnt[finalid1] = {}
            if get_distance:
                hnt.write("%s" % finalid1)
                for j in range(len(ids)):
                    for j2 in range(len(current_set[ids[j]])):
                        nbr+=1
                        finalid2 = "%s_%s" % (ids[j], j2+1)
                        fullid2  = "%s_%s" % (ids[j], current_set[ids[j]][j2])
                        if finalid2 == finalid1:
                            distnt[finalid1][finalid2] = 0
                        elif distnt.has_key(finalid2):
                            distnt[finalid1][finalid2] = distnt[finalid2][finalid1]
                        else:
                            distnt[finalid1][finalid2] = sum([1 for k in range(len(seqs_nt[fullid1])) if seqs_nt[fullid1][k]!=seqs_nt[fullid2][k]])
                            alldistnt.append(distnt[finalid1][finalid2])
                        hnt.write(",%s" % distnt[finalid1][finalid2])     
                hnt.write("\n")
    if get_distance:
        average_nt = sum(alldistnt)/total
        sd_nt      = sqrt(sum([(dist-average_nt)**2 for dist in alldistnt])/(total-1))
        hnt.write("average distance: %.2f nt +/- %.2f (sd)" % (average_nt, sd_nt))
        if verbose:
            print"\n\n"
            print"################### Summary ###################"
            print"number of combinations: %d" % (current_set.keys().__len__())
            print"average distance nt: %.2f +/- %.2f" % (average_nt, sd_nt)
            print"###############################################"
        return ((average_nt, sd_nt))
    return                

def dist_wrapper(path_seq, pkl_import=False):
    """
    streamline operations
    """
    
    path_seq = path_seq.replace('.csv','') 
        
    if pkl_import and path.exists(path_seq+".pkl0"):
        pkl_import = "%s.pkl0" % (path_seq)
    else:
        pkl_import=False
            
    monte_carlo_min_dist("%s.csv" % (path_seq),
                         "%s.pkl" % (path_seq),
                          pkl_import=pkl_import,
                          random_period=5000,
                          focus_period=500,
                          convergence=1000)

    get_final_set_feats("%s.csv" % (path_seq),
                        "%s.pkl0" % (path_seq),
                        "%s_min_set" % (path_seq),
                        get_distance=True,
                        verbose=True)
    return

if __name__ == "__main__":
    if len(sys.argv) == 2:
        db_file = sys.argv[1]
    else:
        db_file = "../testFiles/outputFiles/tfec_2.sqlite.generated_solutions.csv"

    dist_wrapper(path_seq=db_file,pkl_import=True)

    
