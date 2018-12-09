'''
Created on Oct 3, 2011

@author: jcg
'''

codon2aa_table =  { 'aaa' : 'lys',
                    'aag' : 'lys',
                    'tta' : 'leu',
                    'ttg' : 'leu',
                    'ctt' : 'leu',
                    'ctc' : 'leu',
                    'cta' : 'leu',
                    'ctg' : 'leu',
                    'gaa' : 'glu',
                    'gag' : 'glu',
                    'caa' : 'gln',
                    'cag' : 'gln',
                    'cat' : 'his',
                    'cac' : 'his',
                    'aat' : 'asn',
                    'aac' : 'asn',
                    'tct' : 'ser',
                    'tcc' : 'ser',
                    'tca' : 'ser',
                    'tcg' : 'ser',
                    'agt' : 'ser',
                    'agc' : 'ser',
                    'cgt' : 'arg',
                    'cgc' : 'arg',
                    'cga' : 'arg',
                    'cgg' : 'arg',
                    'aga' : 'arg',
                    'agg' : 'arg',
                    'tgg' : 'trp',
                    'gct' : 'ala',
                    'gcc' : 'ala',
                    'gca' : 'ala',
                    'gcg' : 'ala',
                    'tgt' : 'cys',
                    'tgc' : 'cys',
                    'ggt' : 'gly',
                    'ggc' : 'gly',
                    'gga' : 'gly',
                    'ggg' : 'gly',
                    'gat' : 'asp',
                    'gac' : 'asp',
                    'ttt' : 'phe',
                    'ttc' : 'phe',
                    'atg' : 'met',
                    'tat' : 'tyr',
                    'tac' : 'tyr',
                    'gtt' : 'val',
                    'gtc' : 'val',
                    'gta' : 'val',
                    'gtg' : 'val',
                    'act' : 'thr',
                    'acc' : 'thr',
                    'aca' : 'thr',
                    'acg' : 'thr',
                    'cct' : 'pro',
                    'ccc' : 'pro',
                    'cca' : 'pro',
                    'ccg' : 'pro',
                    'att' : 'ile',
                    'atc' : 'ile',
                    'ata' : 'ile',
                    'taa' : 'stop',
                    'tag' : 'stop',
                    'tga' : 'stop' }

from Data import *
from math import log, exp, sqrt
from random import randint, choice, random
import re
import string
from os import system, devnull, chdir, path
from subprocess import call, Popen, PIPE, check_output
import sys

######################
##
#  General Functions
##
######################

# validate a CDS region (i.e. begins with start codon, ends with stop codon, and does not have an in-frame stop codon in the middle of the sequence. 
def validateCDS(cds=""):
    #normalize CDS
    cds = cds.lower().replace('u','t')
    
    #cds length is multiple of 3
    if len(cds) % 3 != 0:
        return False
    
    #starts with a start codon (ATG, GTG, TTG)
    if cds[0:3] not in ('atg','gtg','ttg'):
        return False
    
    #stop with a stop codon
    if cds[-3:] not in ('taa','tag','tga'):
        return False
        
    #stop codon in the middle
    if len(set([cds[i:i+3] for i in range(3,len(cds)-3,3)]).intersection(['taa','tag','tga'])) != 0:
        return False
    
    return True
    
# translate codons to aa
def translateCDS(sequence):
    aa_seq = ""
    for i in range(0,len(sequence),3):
        aa_seq += codon2aa_table[sequence[i:(i+3)]]+" "
    
    return aa_seq
        

# return the indices of the array 'array' sorted
def argsort(array, reverse=False):
    return sorted(range(len(array)), key=array.__getitem__, reverse=reverse)

#return sequence complementary to seq
def complementary(seq):
    return str(seq).translate(string.maketrans("atcg", "tagc"))

#return random nucleotide (different than original)
def randomMutation(nucleotide):
    possible_mut = list(set(bases) - set(nucleotide))
    
    return choice(possible_mut)

#returns a random sequence with defined size
def randomSequence(size):
    sequence = ""
    
    for i in xrange(size):
        sequence += choice(bases)
    
    return sequence 
    
# count the number of different characters between str1 and str2 (hamming distance)
def diff(str1, str2):
    nbr=0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            nbr+=1
    return nbr    

def lin(x, y):
    """
    Summary
        Linear regression of y = ax + b
    Usage
        real, real, real, real, real = lin(list, list)
    Returns coefficients to the regression line "y=ax+b" from x[] and y[], 
            R^2 Value, sum of squares and mean error 
    """
    if len(x) != len(y):
        raise ValueError, 'unequal length'
    n = len(x)
    sx = sy = sxx = syy = sxy = 0.0
    for i, j in zip(x, y):
        sx = sx + i
        sy = sy + j
        sxx = sxx + i*i
        syy = syy + j*j
        sxy = sxy + i*j
    det = sxx * n - sx * sx
    a, b  = (sxy * n - sy * sx)/det, (sxx * sy - sx * sxy)/det
    mean_y = sy/n
    mean_error = residual = 0.0
    for i, j in zip(x, y):
        mean_error = mean_error + (j - mean_y)**2
        residual  = residual  + (j - a * i - b)**2
    r2  = 1 - residual/mean_error
    ssq = residual / (n-2)
    #var_a, var_b = ss * n / det, ss * Sxx / det
    return a, b, r2, ssq, mean_error

def appendLabelToDict(somedict,label):
    return dict(map(lambda (key, value): (label+str(key), value), somedict.items()))


def average(array):
    return sum(array)*1.0/len(array)

def stddev(array):
    avg = average(array)    
    return sqrt(average(map(lambda x: (x - avg)**2, array)))
    
#select randomly based on probability, prob_list = [('1',0.1),('2',0.5),('3',0.4)]
def pick_random_tuple(prob_list):
    r, s = random(), 0
    for num in prob_list:
        s += num[1]
        if s >= r:
            return num[0]
        
#select randomly based on probability, prob_list = [0.1,0.5,0.4]
def pick_random(prob_list):
    r, s, i = random(), 0, 0
    for num in prob_list:
        s += num
        if s >= r:
            return i
        i += 1        

    
######################
##
#  Scoring Functions
##
######################

def structureAnalysis(structure_file, propertyOfInterest="ss"):
    '''
        given a structure file it returns the position that are either single stranded "ss" or double stranded "ds"
    '''
    
    if path.exists("tmp/structures/"+structure_file+".ct"):
        if propertyOfInterest == "ds":
            output = Popen("perl 3rdParty/unafold/ss-count.pl tmp/structures/" + structure_file + ".ct | awk /[[:digit:]][[:blank:]]0/'{print $1}'", stdout=PIPE, shell=True).stdout.read()
            return [eval(k) for k in output.split()]
        elif propertyOfInterest == "ss":
            output = Popen("perl 3rdParty/unafold/ss-count.pl tmp/structures/" + structure_file + ".ct | awk /[[:digit:]][[:blank:]]1/'{print $1}'", stdout=PIPE, shell=True).stdout.read()
            return [eval(k) for k in output.split()]
        else:
            return []


def analyzeCodons(seq, positions = None, data_table = cai_table):
    '''
       given a sequence it returns a list with two elements: [ list_of_codons, list_of_codons_cai]
    '''    
    
    if positions == None:
        positions = range(0,len(seq),3)
    
    seq = seq.lower();    
    codons = []
    codons_cai = []
    for i in positions:
        codon = seq[i:i+3]        
        codons.append(codon)
        if data_table.has_key(codon):
            codons_cai.append(data_table[codon])
        else:
            codons_cai.append("NA")
    return [codons, codons_cai]


def get_alternate_codons(codon, data=tai_tuller, dist=0):
    """
        returns a alternate codon to codon
        data: dictionary with a map between codons and tAI
        dist: 0   --> only synonymous codon
              1-3 --> only codon with 1-3 nt difference from original
    """    
    if dist == 0:
        # return only syn codon
        return [(syn_cod, data[syn_cod]) for syn_cod in aa2codon_table[codon2aa_table[codon]] if syn_cod != codon]
    else:
        # return syn codon and codon 1 nt away
        return [(alt_cod, data[alt_cod]) for alt_cod in codons_list if (alt_cod != codon and diff(codon, alt_cod) <= dist)]

def get_tai(seq, data=tai_tuller):
    seq = seq.lower()
#    seq = seq.replace("t","u")
    return [(seq[i:i+3], data[seq[i:i+3]]) for i in range(0,len(seq),3)]
    
def analyze_tai(seq, window=21, data=tai_tuller, method="harmonic"):
    seq = seq.lower()
#   seq = seq.replace("t","u")
    scores   = [data[seq[i:i+3]] for i in range(0,len(seq),3)]
    smoothie = []
    if window > 1:
        if method == "geometric":
            for i in range(len(scores)-window+1):
                smoothie.append(1)
                for j in range(i,i+window):
                    smoothie[-1] *= scores[j]
                smoothie[-1] = smoothie[-1]**(1/float(window))
        if method == "harmonic":
            for i in range(len(scores)-window+1):
                smoothie.append(sum([1/v for v in scores[i:i+window]])/window) 
    return scores, smoothie


def analyze_hydropathy(seq):
    seq = seq.lower();
    score = 0
    len_sq = 0
    for i in range(0,len(seq),3):
        if hydropathy_index_table.has_key(seq[i:i+3]):
            score += (hydropathy_index_table[seq[i:i+3]])
            len_sq += 1
    score /= len_sq
    return score

def analyze_cai(seq):
    seq = seq.lower();
    score = 0
    len_sq = 0
    for i in range(0,len(seq),3):
        if cai_table.has_key(seq[i:i+3]):
            score += log(cai_table[seq[i:i+3]])
            len_sq += 1
    score /= len_sq
    return exp(score)

def analyze_bottleneck(sequence, window=20, data=tai_tuller, method="harmonic"):
    score, smooth = analyze_tai(sequence, window, data=data, method=method)
    
    return [score, smooth]

def analyze_bottleneck_pos(sequence, score, smooth):
    return smooth.index(max(smooth))+1

def analyze_bottleneck_abs_strength(sequence, score, smooth):
    return max(smooth)

def analyze_bottleneck_rel_strength(sequence, score, smooth):
    return (len(score)*max(smooth)/sum([1/v for v in score]))

def analyze_ntcontent(seq):
    
    seq = seq.replace('u', 't')
    nuc_freq = { 'NucleotideContentAT': 0, 
                 'NucleotideContentGC': 0, 
                 'NucleotideContentA' : 0, 
                 'NucleotideContentT' : 0, 
                 'NucleotideContentG' : 0, 
                 'NucleotideContentC' : 0}
    
    for i in range (len(seq)):
        nuc_freq['NucleotideContent'+seq[i].upper()] += 1
    
    nuc_freq['NucleotideContentA'] /= float(seq.__len__())
    nuc_freq['NucleotideContentT'] /= float(seq.__len__())
    nuc_freq['NucleotideContentG'] /= float(seq.__len__())
    nuc_freq['NucleotideContentC'] /= float(seq.__len__())    
    nuc_freq['NucleotideContentAT'] = (nuc_freq['NucleotideContentA']+nuc_freq['NucleotideContentT'])
    nuc_freq['NucleotideContentGC'] = (nuc_freq['NucleotideContentG']+nuc_freq['NucleotideContentC'])    
    
    return  nuc_freq

def analyze_pwm_score(seq, pwm):
    #print seq
    max_score = -99
    max_pos   = 0 
    
    pwm_length = len(pwm[pwm.keys()[0]])
    
    for pos in range(0,len(seq)-pwm_length+1):
        tmp_score = pwm_score(seq[pos:(pos+pwm_length)], pwm)
        
        if tmp_score > max_score:
            max_score = tmp_score
            max_pos = pos
            
    return { 'MotifScore' : max_score , 'MotifPosition' : max_pos }
        
def pwm_score(seq,pwm):
    score=0    
    for i in range(0,len(seq)):
        score += pwm[seq[i]][i]
    return  (score)

def analyze_riboden(seq, filename, init_rate=1):        
     
    chdir(project_dir)
    
    chdir(project_dir+"/3rdParty/TASEP/")
        
    s1_file = filename + ".seq"
    
    system("echo '" + str(seq).upper() + "' > " + s1_file)        
    fnull = open(devnull, 'w')   #this line is necessary to omit output generated by UNAFOLD             
    output = Popen("./executable_file "+s1_file+" "+str(init_rate)+" | tail -1 | awk '{print $4}'", stdout=PIPE, stderr=None, shell=True).stdout.read().strip()
                
    system("rm %s*" % s1_file)
    fnull.close()
    
    chdir(project_dir)
    
    if output == '':
        return None
    else:
        return float(output)
 

def analyze_terminator(seq):
    #use transtermHP
    #create input files
    file1 = 'test' + str(randint(0000,1000)) + '.fa';
    file2 = 'test' + str(randint(0000,1000)) + '-fake.coords';
    
    system("echo '>seq1\n" + str(seq) + "' > tmp/transterm_files/"+file1)
    system("echo 'fakegene1    1 2    seq1\nfakegene1    " + str(len(seq)-1) + " " + str(len(seq)) +"    seq1' > tmp/transterm_files/"+file2)
    #run transtermhp
    output = Popen("./3rdParty/transterm/transterm -p 3rdParty/transterm/expterm.dat tmp/transterm_files/"+file1+" tmp/transterm_files/"+file2+" 2> tmp/transterm_files/err.txt | awk /TERM/'{print $8}'", stdout=PIPE, stderr=None, shell=True).stdout.read()

    if output == '':
        return 'No'
    else:
        #print "Terminator output: " + output
        all_term_scores = [float(score) for score in output.split()]
        if max(all_term_scores) > 90:
            return 'Yes'
        elif max(all_term_scores) > 70:
            return 'Maybe'
        else:
            return 'No'               


def analyze_duplex_structure(seq1,seq2,filename):        
     
    chdir(project_dir)
    
    structure_pairs = {}
    data = {}
        
    r = randint(4,1000)
    
    s1_file = str(r) + "s1.seq"
    s2_file = str(r) + "s2.seq"

    system("echo '" + str(seq1) + "' > " + s1_file)
    system("echo '" + str(seq2) + "' > " + s2_file)        
    fnull = open(devnull, 'w')   #this line is necessary to omit output generated by UNAFOLD         
    call("./3rdParty/unafold/hybrid-min -n RNA -o "+filename+" "+ s1_file + " " + s2_file + " -t 37 -T 37", shell = True, stdout = fnull, stderr = fnull)     #code is necessary to omit output generated by UNAFOLD    
    
    if path.exists(project_dir+"/"+filename+".ct"):
        system("mv %s*.ct tmp/structures/" % filename)
        system("mv %s*.asc tmp/structures/" % filename)
        
        output_ds = Popen("cat tmp/structures/" + filename + ".ct | awk '{print $1 , $5}'", stdout=PIPE, shell=True).stdout.read()
        list_of_pairs = output_ds.split('\n')[1:-1]
        
        for pair_str in list_of_pairs:
            pair = pair_str.split()        
            structure_pairs[pair[0]] = pair[1]
                
    data['RNADuplexPairs'] = structure_pairs
    system("rm %s*" % filename)
    system("rm %s*.*" % r)
    fnull.close()
    
    return data 

def analyze_duplex_structure_rnafold(seq1,seq2,filename):
        
     
    chdir(project_dir)
    
    data = {}
        
    r = randint(4,1000)
    
    s12_file = str(r) + "s1.s2.seq"

    system("echo '>s1\n" + str(seq1) + "\n>s2\n"+str(seq2)+"\n' > " + s12_file)
        
    fnull = open(devnull, 'w')   #this line is necessary to omit output generated by UNAFOLD         
    call("./3rdParty/vienna/RNAduplex < "+s12_file+" > "+ filename + ".rnafold.out", shell = True, stdout = fnull, stderr = fnull)     #code is necessary to omit output generated by RNAfold    
    
    if path.exists(project_dir+"/"+filename+".rnafold.out"):
        system("mv %s*.rnafold.out tmp/structures/" % filename)
                        
    data['RNADuplexPairs'] = "NA"
#    system("rm %s*" % filename)
    system("rm %s*.*" % r)
    fnull.close()
    
    return data 

def analyze_duplex_mfe(filename,region = None):
    data = {}
         
    chdir(project_dir)
    
    if path.exists(project_dir+"/tmp/structures/"+filename+".ct"):
        try:
            res = check_output(["./3rdParty/unafold/ct-energy" , "tmp/structures/" + filename + ".ct"])
            
            if res != "":
                data['RNADuplexMFE'] = float(str(res).rstrip())
            else:
                data['RNADuplexMFE'] = 'NA'
                
        except NameError:
            data['RNADuplexMFE'] = 'NA'

        #system("rm tmp/structures/%s*" % filename)
    else:
        data['RNADuplexMFE'] = 'NA'
        
    return data

def analyze_duplex_mfe_rnafold(filename,region = None):
    data = {}
         
    chdir(project_dir)
    
    if path.exists(project_dir+"/tmp/structures/"+filename+".rnafold.out"):
        try:            
            output = check_output(["tail", "-n 1", "tmp/structures/" + filename + ".rnafold.out"]).rstrip()

            m = re.search('-?\d+\.\d+',output)            
            res = m.group(0)
            
            if res != "":
                data['RNADuplexRNAFoldMFE'] = float(res)
            else:
                data['RNADuplexRNAFoldMFE'] = 'NA'
                
        except NameError:
            data['RNADuplexRNAFoldMFE'] = 'NA'

        #system("rm tmp/structures/%s*" % filename)
    else:
        data['RNADuplexRNAFoldMFE'] = 'NA'
        
    return data

def analyze_duplex_ds(filename,seq1 = "", region = None):
    data = {}   
    
    chdir(project_dir)
    
    if region == None:
        if path.exists(project_dir+"/tmp/structures/"+filename+".ct"):            
            output_ds = Popen("perl 3rdParty/unafold/ss-count.pl tmp/structures/" + filename + ".ct | awk /[[:digit:]][[:blank:]]0/'{print $1}'", stdout=PIPE, shell=True).stdout.read()
            output = output_ds.split()
            data['RNADuplexDoubleStrandedBasesList_Mol1'] = [eval(k) for k in output if eval(k) <= len(seq1)]
            data['RNADuplexDoubleStrandedBases_Mol1'] = len(data['RNADuplexDoubleStrandedBasesList_Mol1'])
            data['RNADuplexDoubleStrandedBasesList_Mol2'] = [eval(k)-len(seq1) for k in output if eval(k) > len(seq1)]
            data['RNADuplexDoubleStrandedBases_Mol2'] = len(data['RNADuplexDoubleStrandedBasesList_Mol2'])
        else:
            data['RNADuplexDoubleStrandedBasesList_Mol1'] = 'NA'
            data['RNADuplexDoubleStrandedBases_Mol1'] = 'NA'
            data['RNADuplexDoubleStrandedBasesList_Mol2'] = 'NA'
            data['RNADuplexDoubleStrandedBases_Mol2'] = 'NA'

    #else:
    # TODO: get just DS bases from one region of the rna structure
        
    return data

def analyze_duplex_ss(filename,seq1 = "",region = None):
    data = {}       
      
    chdir(project_dir)
    
    if region == None:
        if path.exists(project_dir+"/tmp/structures/"+filename+".ct"):
            output_ss = Popen("perl 3rdParty/unafold/ss-count.pl tmp/structures/" + filename + ".ct | awk /[[:digit:]][[:blank:]]1/'{print $1}'", stdout=PIPE, shell=True).stdout.read()                    
            output = output_ss.split()
            data['RNADuplexSingleStrandedBasesList_Mol1'] = [eval(k) for k in output if eval(k) <= len(seq1)]
            data['RNADuplexSingleStrandedBases_Mol1'] = len(data['RNADuplexDoubleStrandedBasesList_Mol1'])
            data['RNADuplexSingleStrandedBasesList_Mol2'] = [eval(k)-len(seq1) for k in output if eval(k) > len(seq1)]
            data['RNADuplexSingleStrandedBases_Mol2'] = len(data['RNADuplexDoubleStrandedBasesList_Mol2'])
        else:
            data['RNADuplexSingleStrandedBasesList_Mol1'] = 'NA'
            data['RNADuplexSingleStrandedBases_Mol1'] = 'NA'
            data['RNADuplexSingleStrandedBasesList_Mol2'] = 'NA'
            data['RNADuplexSingleStrandedBases_Mol2'] = 'NA'    #else:
    # TODO: get just DS bases from one region of the rna structure
        
    return data


def analyze_structure_prob(seq,filename,window=50,region=[]):
    
    
    chdir(project_dir)
    
    structure_pairs = {}
    data = {}

    system("echo '>" + filename + "\n" + str(seq) + "' > " + filename + ".fa")        
    fnull = open(devnull, 'w')   #this line is necessary to omit output generated by UNAFOLD         
    call("3rdParty/vienna/RNAplfold -d2 -noLP -W "+str(window)+" -u 1 < " + filename + ".fa", shell = True, stdout = fnull, stderr = fnull)     #code is necessary to omit output generated by RNAplfold

    output_ss = Popen("cat " + filename + "_lunp | awk '{print $1 \"\t\" $2}'", stdout=PIPE, shell=True).stdout.read()
    l = output_ss.rstrip().split('\n')[2:]    
    
    for p in l:
        pair = p.split()
        if pair[1] != 'NA':
            structure_pairs[pair[0]]=float(pair[1])
        else:
            structure_pairs[pair[0]]='NA'
    
    if region != []:
        reg_avg = 0
        for pos in region:
            reg_avg += structure_pairs[str(pos)]
        reg_avg = reg_avg/len(region)
        data['StructureProb'] = reg_avg        
    else:
        data['StructureProb'] = sum(structure_pairs.values())/len(structure_pairs.keys())
        
    data['StructureProbList'] = structure_pairs 
        
    # remove tmp files
    system("rm %s*" % filename)
    #system("mv " + filename + "* tmp/unafold_files/")
    fnull.close()             
 
    return data

def analyze_ensemble(seq,filename,sample_size=100):
    
    chdir(project_dir)

    system("echo '>" + filename + "\n" + str(seq) + "' > " + filename + ".fa")                 

    output_ss = Popen("./3rdParty/vienna/RNAsubopt -d2 -noLP -s -p "+str(sample_size)+" < " + filename + ".fa | tail -n "+str(sample_size), stdout=PIPE, shell=True).stdout.read()
        
    l = output_ss.rstrip().split('\n')
    ens_st = []
    
    string_aux = ""
    
    for st in l:
        string_aux += str(seq) + "\n" + str(st) + "\n"
        
    system("echo '" + str(string_aux) + "' > " + filename + ".st")
    output_ss = Popen("./3rdParty/vienna/RNAeval -d2 < " + filename + ".st | perl -lne 'm/.* \((.*)\)$// print $1'", stdout=PIPE, shell=True).stdout.read().rstrip()
    ens_st = map(lambda x: float(x), output_ss.rstrip().split('\n')[2:])

    data = {}
    data['StructureEnsembleSample'] = ens_st
    data['StructureEnsembleSampleMean'] = average(ens_st)
    data['StructureEnsembleSampleSD'] = stddev(ens_st)
    
    system("rm %s*" % filename)
    
    # remove tmp files
    #system("rm %s*" % filename)
    #system("mv " + filename + "* tmp/unafold_files/")             
 
    return data

def analyze_structure(seq,filename,ensemble=False):
    
        
    chdir(project_dir)

    system("echo '" + str(seq) + "' > " + filename + ".seq")        
    fnull = open(devnull, 'w')   #this line is necessary to omit output generated by UNAFOLD
    if ensemble:         
        call("./3rdParty/unafold/UNAFold.pl -n RNA " + filename + ".seq", shell = True, stdout = fnull, stderr = fnull)     #code is necessary to omit output generated by UNAFOLD
    else:
        call("./3rdParty/unafold/hybrid-ss-min -n RNA " + filename + ".seq", shell = True, stdout = fnull, stderr = fnull)     #code is necessary to omit output generated by UNAFOLD
                
    if os.path.isfile(filename+".ct"):        
        system("mv %s*.ct tmp/structures/" % filename)
        # remove tmp files
    
    system("rm %s*" % filename)
    #system("mv " + filename + "* tmp/unafold_files/")
    fnull.close()             
 
    return 1

def analyze_structure_rnafold(seq,filename):
    
        
    chdir(project_dir)

    system("echo '" + str(seq) + "' > " + filename + ".seq")        
    fnull = open(devnull, 'w')   #this line is necessary to omit output generated by RNAfold
    call("./3rdParty/vienna/RNAfold -i " + filename +".seq --noPS -d2 > " + filename + ".rnafold.out", shell = True, stdout = fnull, stderr = fnull)     #code is necessary to omit output generated by RNAplfold
                
    if os.path.isfile(filename+".rnafold.out"):        
        system("mv %s*.rnafold.out tmp/structures/" % filename)
        # remove tmp files
    
    system("rm %s*" % filename)
    #system("mv " + filename + "* tmp/unafold_files/")
    fnull.close()             
 
    return 1

def analyze_structure_mfe(filename):
    data = {}
         
    chdir(project_dir)
    
    if path.exists(project_dir+"/tmp/structures/"+filename+".ct"):
        
        output = check_output(["./3rdParty/unafold/ct-energy" , "tmp/structures/" + filename + ".ct"]).rstrip()
        mfe_list = [float(a) for a in output.split('\n')]
        
        data['StructureMFE'] = mfe_list[0]
    else:
        data['StructureMFE'] = 0
        
    return data

def analyze_structure_mfe_rnafold(filename):
    data = {}
         
    chdir(project_dir)
    
    if path.exists(project_dir+"/tmp/structures/"+filename+".rnafold.out"):
        
        output = check_output(["tail", "-n 1" , "tmp/structures/" + filename + ".rnafold.out"]).rstrip()
        
        m = re.search('-?\d+\.\d+',output)
        data['StructureRNAFoldMFE'] = float(m.group(0))
        
    else:
        data['StructureRNAFoldMFE'] = "NA"
        
    return data


def analyze_structure_ds(filename,region = None):
    data = {}   
    
    chdir(project_dir)
    
    if region == None:
        if path.exists(project_dir+"/tmp/structures/"+filename+".ct"):
            output_ds = Popen("perl 3rdParty/unafold/ss-count.pl tmp/structures/" + filename + ".ct | awk /[[:digit:]][[:blank:]]0/'{print $1}'", stdout=PIPE, shell=True).stdout.read()            
            data['StructureDoubleStrandedList'] = [eval(k) for k in output_ds.split()]
            data['StructureDoubleStranded'] = len(data['StructureDoubleStrandedList'])
        else:
            data['StructureDoubleStrandedList'] = 'NA'
            data['StructureDoubleStranded'] = 'NA'
    #else:
    # TODO: get just DS bases from one region of the rna structure
        
    return data

def analyze_structure_ss(filename,region = None):
    data = {}       
      
    chdir(project_dir)
    
    if region == None:
        if path.exists(project_dir+"/tmp/structures/"+filename+".ct"):
            output_ss = Popen("perl 3rdParty/unafold/ss-count.pl tmp/structures/" + filename + ".ct | awk /[[:digit:]][[:blank:]]1/'{print $1}'", stdout=PIPE, shell=True).stdout.read()                    
            data['StructureSingleStrandedList'] = [eval(k) for k in output_ss.split()]
            data['StructureSingleStranded'] = len(data['StructureSingleStrandedList'])
        else:
            data['StructureSingleStrandedList'] = 'NA'
            data['StructureSingleStranded'] = 'NA'
    #else:
    # TODO: get just DS bases from one region of the rna structure
        
    return data

def analyze_structure_accessibility(filename):
    data = {}       
       
    chdir(project_dir)
    
    if path.exists(project_dir+"/tmp/structures/"+filename+".ct"):
        #output = Popen("perl 3rdParty/unafold/ss-count.pl -w tmp/structures/" + filename + ".ct | awk '{sum +=  $2; count += 1}  END{print sum/count}'", stdout=PIPE, shell=True).stdout.read()
        #data['StructureAccessibility'] = eval(output)                    
        output = Popen("perl 3rdParty/unafold/ss-count.pl -w tmp/structures/" + filename + ".ct | awk '{print $2}'", stdout=PIPE, shell=True).stdout.read()            
        data['StructureEnsembleAccessibilityBasesList'] = [eval(e) for e in output.split('\n') if e != '']
        data['StructureEnsembleAccessibility'] = sum(data['StructureEnsembleAccessibilityBasesList'])/len(data['StructureEnsembleAccessibilityBasesList'])
    else:
        data['StructureEnsembleAccessibility'] = 'NA'

        
    return data
                  

######################
##
#  Operators Functions
##
######################


def SimpleStructureOperator(sequence, structurefile, structure_range, mutable_region, cds_region, direction, keep_aa = True, ss_bases=None, ds_bases=None):
    '''
        Operator that given a sequence and structure, mutates the sequence to change structure
            seq: sequence in structure
            structure_file: a file containing the structure correspondent to sequence
            structure_range - start and end position to calculate structure - a tuple in the form (start, end)  
            mutable_region - a list with all bases that can be mutated
            cds_regions - a list of pairs with begin and end of CDSs - example: [(0,100), (150,300)]            
            direction: either increase ('+') or decrease ('-') single bases 
    ''' 
    
    if not mutable_region: #it's not possible to mutate
        return None
    
    if ss_bases == None:
        ss_bases = structureAnalysis(structurefile, "ss")
    
    if ds_bases == None:
        ds_bases = structureAnalysis(structurefile, "ds")
    
    if direction == '+':       
        #get double stranded bases    
        baseToMutate = [(b+structure_range[0]-1) for b in ds_bases if (b+structure_range[0]-1) in mutable_region]
    elif direction == '-':
        #get single stranded bases
        baseToMutate = [(b+structure_range[0]-1) for b in ss_bases if (b+structure_range[0]-1) in mutable_region]
    else:
        sys.stderr.write("Direction Unknown")

    mutated = False
    iteration = 0
    
    while not mutated and iteration <= 100:         
        #position to mutate
        index_to_mutate = baseToMutate.pop(randint(0,len(baseToMutate)-1)) if len(baseToMutate) != 0 else mutable_region.pop(randint(0,len(mutable_region)-1))
                    
        #mutate base
        if keep_aa == True and index_to_mutate >= cds_region[0] and index_to_mutate <= cds_region[1]:
            codon_p = (index_to_mutate-cds_region[0])/3
            initial_codon = sequence[(cds_region[0]+codon_p*3):(cds_region[0]+codon_p*3+3)]
            #codon position
            p1 = (cds_region[0]+codon_p*3)
            p2 = (cds_region[0]+codon_p*3+1)
            p3 = (cds_region[0]+codon_p*3+2)         
            alt_codons = ([c for c in aa2codon_table[codon2aa_table[initial_codon]] if c!= initial_codon])
            
            while len(alt_codons) != 0:
                rnd_alt_codon = alt_codons.pop(randint(0,len(alt_codons)-1))
                
                valid = True 
                if initial_codon[0] != rnd_alt_codon[0] and not(p1 in mutable_region):
                    valid = False
                if initial_codon[1] != rnd_alt_codon[1] and not(p2 in mutable_region):
                    valid = False
                if initial_codon[2] != rnd_alt_codon[2] and not(p3 in mutable_region):
                    valid = False    
                
                if valid == True:
                    mutated = True 
                    new_seq = sequence[:(cds_region[0]+codon_p*3)] + rnd_alt_codon + sequence[(cds_region[0]+codon_p*3+3):]
                    #print "------"
                    #print initial_codon
                    #print rnd_alt_codon
                    #print sequence
                    #print new_seq
            
        else:
            mutated = True
            new_seq = list(sequence)
            if direction == '+':
                comp = complementary(sequence[index_to_mutate])
            else:
                comp = randomMutation(sequence[index_to_mutate])
            new_seq[index_to_mutate] = comp
            #print sequence
            #print "".join(new_seq)     
            
        iteration+=1                
                
    return "".join(new_seq)

def SimpleBneckOperator(seq, direction="", distance=0):
    """
    Modify bottleneck in seq
    direction: + --> increase bottleneck
               - --> decrase bottleneck
    One codon is chosen at random among the corresponding tertile based on tai value
    This codon is randomly replaced by possible alternative codon in the appropriate tertile
    Distance control the codon similarity.
    distance: 0   --> only synonymous codon
              1-3 --> only codon with 1-3 nt difference from original
    """
    #print "dist: " , str(distance)
    if direction == "":
        return None    
          
    iteration = 0
    tai = get_tai(seq)
    alt_codons= []

    while len(alt_codons)==0 and iteration != 100:
        iteration += 1
        codonPosition2replace = randint(0,len(tai)-1)
        alt_codons = get_alternate_codons(tai[codonPosition2replace][0], dist=distance)
        
        #print "attempt to replace: " + str(codonPosition2replace) + "-" + str(tai[codonPosition2replace])
        
        if direction == "+":
            # randomly pick lower tai for the alternate solution given allowable distance
            alt_codons = [alt_cod for alt_cod in alt_codons if alt_cod[1] < tai[codonPosition2replace][1]]
        elif direction == "-":
            # randomly pick higher tai for the alternate solution given allowable distance
            alt_codons = [alt_cod for alt_cod in alt_codons if alt_cod[1] > tai[codonPosition2replace][1]]
            
        if len(alt_codons) != 0:
            randomAltCodon = randint(0,len(alt_codons)-1)
            #print "Bot operator: old codon -> " + str(tai[codonPosition2replace][0]) 
            #print "Bot operator: new codon -> " + str(alt_codons[randomAltCodon][0])

    
    if iteration == 100:
        return None
    
    return seq[:3*codonPosition2replace]+alt_codons[randomAltCodon][0]+seq[3*(codonPosition2replace+1):]

def SimpleNtContentOperator(seq, direction=0, nucleotide = [], mutable_region = None, cds_region = (0,9), keep_aa=True):
    if direction == 0:
        return seq
        
    mutated = False    
    seq = seq.lower() 
    
    #check direction to decide possible mutations
    if direction == '-':
        #select mutable position based on presence of nucleotide(s)        
        mutable_positions = [pos for pos in mutable_region if seq[pos] in set(nucleotide)]
    elif direction == '+':
        #select mutable position based on absence of nucleotide(s)            
        mutable_positions = [pos for pos in mutable_region if not (seq[pos] in set(nucleotide))]      

    while not mutated and mutable_positions.__len__() != 0 :
        #print mutable_positions
        rnd_pos = mutable_positions.pop(randint(0,mutable_positions.__len__()-1))
        
        if direction == '-':
            possible_mutations = list(set(bases) - set(nucleotide))
        elif direction == '+':
            possible_mutations = list(nucleotide)    
        
        while not mutated and possible_mutations.__len__() != 0:
                                     
            new_seq = seq[:rnd_pos]+possible_mutations.pop(randint(0,possible_mutations.__len__()-1))+seq[rnd_pos+1:]
            
            if keep_aa == True and rnd_pos >= cds_region[0] and rnd_pos <= cds_region[1]:
                #check if AA remains the same
                codon_p = (rnd_pos-cds_region[0])/3
                initial_codon = seq[(cds_region[0]+codon_p*3):(cds_region[0]+codon_p*3+3)]
                final_codon = new_seq[(cds_region[0]+codon_p*3):(cds_region[0]+codon_p*3+3)]
                
                #print "initial codon: " + str(initial_codon) + " AA: " + codon2aa_table[initial_codon]
                #print "final codon: " + str(final_codon) + " AA: " + codon2aa_table[final_codon]            
        
                if codon2aa_table[initial_codon] == codon2aa_table[final_codon]:
                    mutated = True
            else:
                mutated = True                
                    
    if mutated == False:
        return None
    
    return new_seq
          

def SimplePWMScoreOperator(seq, pwmnt, direction=0, mutable_region = None, keep_aa=False, max_iter=100):    
    if direction == 0:
        return seq
    elif direction == '+':
        direction = 1;
    elif direction == '-':
        direction = -1 
        
    iteration    = 0
    mutated = False
    new_seq = seq

    res = analyze_pwm_score(seq,pwmnt)
    motif_position = res['MotifPosition']

    motifsize = len(pwmnt[pwmnt.keys()[1]])

    while not mutated:
        iteration += 1
        if iteration == max_iter:
            return None 
        #draw position at random and check if better can be found
        pos = choice(list(set(mutable_region) & set(range(motif_position,motif_position+motifsize))))
        
        rel_pos = pos - motif_position
        current_value = pwmnt[new_seq[pos].replace('u','t')][rel_pos]
        base = randomMutation(seq[pos])
        if direction * (pwmnt[base.replace('u','t')][rel_pos] - current_value) > 0: #check whether the change goes in the right direction
            new_seq = seq[:pos]+base+seq[pos+1:]
            mutated = True
            if keep_aa:
                #get codon assuming the sequence is in frame
                cod_offset = (pos)%3
                cod        = seq[pos-cod_offset:(pos-cod_offset+3)]
                new_cod    = new_seq[pos-cod_offset:(pos-cod_offset+3)]
                if codon2aa_table[cod] != codon2aa_table[new_cod]:
                    mutated = False
                    new_seq = seq
                    
    return new_seq


def SimplePWMPositionOperator(seq, pwmnt, position=None, mutable_region = None, keep_aa=False, max_iter=100):    
    if position == None:
        return seq
        
    iteration = 0
    mutated = False
    new_seq = seq

    res = analyze_pwm_score(seq,pwmnt)
    motif_position = res['MotifPosition']

    motifsize = len(pwmnt[pwmnt.keys()[1]])

    while not mutated:
        iteration += 1
        if iteration == max_iter:
            return None 

        if choice([True,False]): #remove unwanted motif
            #print "remove current motif"
            direction = -1
            pos = choice(list(set(mutable_region) & set(range(motif_position,motif_position+motifsize))))
            tmp_motif_position=motif_position        
        else: #create desired motif
            #print "improve new motif"
            direction = 1
            pos = choice(list(set(mutable_region) & set(range(position,position+motifsize))))
            tmp_motif_position = position
            
        #print "tmp_motif_position: ",tmp_motif_position
        #print "pos: ", pos
        rel_pos = pos - tmp_motif_position
        current_value = pwmnt[new_seq[pos].replace('u','t')][rel_pos]
        base = randomMutation(seq[pos])
        if direction * (pwmnt[base.replace('u','t')][rel_pos] - current_value) > 0:
        # check whether the change goes in the right direction
            new_seq = seq[:pos]+base+seq[pos+1:]
            mutated = True
            if keep_aa:
                #get codon assuming the sequence is in frame
                cod_offset = (pos)%3
                cod        = seq[pos-cod_offset:(pos-cod_offset+3)]
                new_cod    = new_seq[pos-cod_offset:(pos-cod_offset+3)]
                if codon2aa_table[cod] != codon2aa_table[new_cod]:
                    mutated = False
                    new_seq = seq
                    
    #print "--------"
    #print analyze_pwm_score(seq,pwmnt) , seq
    #print analyze_pwm_score(new_seq,pwmnt) , new_seq 
    #print "--------"
                    
    return new_seq


def randomMutationOperator(sequence, keep_aa, mutable_region, cds_region, pos=None, n_mut = [1,2]):
    '''
        Operator that given a sequence, mutates the sequence randomly
            sequnce: sequence   
            mutable_region - a list with all bases that can be mutated
            cds_regions - a list of pairs with begin and end of CDSs - example: [(0,100), (150,300)]             
    '''    
    mutableCodonsPosition = [c for c in range(cds_region[0],cds_region[1],3) if set([c,c+1,c+2]).issubset(mutable_region)]
    mutableUTRPosition = list(set(mutable_region) - set(range(cds_region[0],cds_region[1])))      
    
    if mutableCodonsPosition == [] and mutableUTRPosition == []:
        sys.stderr.write("randomMutationOperator: No codons available for mutation\n")
        return None
    else:
        if keep_aa == True:
            if (mutableUTRPosition == []) or (mutableCodonsPosition != [] and choice([True,False])):
                return mutateCDS(sequence, keep_aa, mutableCodonsPosition, cds_region, pos, n_mut)
            else:
                return mutateAll(sequence, keep_aa, mutableUTRPosition, cds_region, pos, n_mut)
        else:
            return mutateAll(sequence, keep_aa, mutable_region, cds_region, pos, n_mut)
                
    
def mutateCDS(sequence, keep_aa, mutableCodonsPosition, cds_region, pos=None, n_mut = [1,2]):
    if keep_aa == True:                
        result = analyzeCodons(sequence,mutableCodonsPosition)
        
        n_mutations = choice(n_mut)
        
        codons = (result[0])        
        codons_ind = range(0,codons.__len__())
        
        mutated = False
        while codons_ind.__len__() != 0 and n_mutations > 0:
            rnd_ind       = codons_ind.pop(randint(0,codons_ind.__len__()-1))
            rnd_codon     = codons[rnd_ind]            
            alt_codons = [c for c in aa2codon_table[codon2aa_table[rnd_codon]] if c!= rnd_codon and codon2aa_table[c]!='stop']
            if alt_codons.__len__() != 0:                      
                mutated = True
                n_mutations -= 1
                new_codon = choice(alt_codons)
                real_codon_pos = mutableCodonsPosition[rnd_ind]
                codon_position = (real_codon_pos-cds_region[0])/3
                all_codons = analyzeCodons(sequence,range(cds_region[0],cds_region[1]+1,3))[0]
                all_codons[codon_position]=new_codon
                new_seq = sequence[:cds_region[0]] + ''.join(c for c in all_codons) + sequence[cds_region[1]+1:]
                sequence = new_seq
                
        if mutated == False:
            sys.stderr.write("RandomMutator: Not able to mutate sequence keeping AA\n")
            return None
        else:                  
            return new_seq

def mutateAll(sequence, keep_aa, mutable_region, cds_region, pos=None, n_mut = [1,2]):
    
    #####
    # Not necessary to keep AA
    #
    n_mutations = choice(n_mut)
    
    if mutable_region == []:
        return None
    
    while n_mutations!=0:
        n_mutations -= 1    
        if pos != None:
            intersect_mut = list(set(mutable_region) & set(pos))
            if intersect_mut != []:
                position_to_mutate = choice(intersect_mut)
            else:
                position_to_mutate = choice(mutable_region)
        else:
            position_to_mutate = choice(mutable_region)
        mutation = randomMutation(sequence[position_to_mutate])
                
        new_seq = sequence[:position_to_mutate] + mutation + sequence[position_to_mutate+1:]
        sequence = new_seq
    
    return new_seq               
    

def SimpleCAIOperator(sequence, cai_range, keep_aa, mutable_region, cds_regions, direction = '+'):
    '''
        Operator that given a sequence, mutates the sequence to change CAI
            sequnce: sequence 
            cai_range - start and end position to calculate cai - a tuple in the form (start, end)  
            mutable_region - a list with all bases that can be mutated
            cds_regions - a list of pairs with begin and end of CDSs - example: [(0,100), (150,300)]            
            direction: either increase ('+') or decrease ('-') CAI 
    '''    
    mutated = False
    mutableCodonsPosition = [c for c in range(cai_range[0],cai_range[1],3) if set([c,c+1,c+2]).issubset(mutable_region)]
                    
    if len(mutableCodonsPosition) == 0:
        sys.stderr.write("SimpleCAIOperator: No codons available for mutation\n")
        return None
    
    result = analyzeCodons(sequence,mutableCodonsPosition)

    codons = (result[0])    
    codons_cai = (result[1])    
    codons_ind = range(0,codons.__len__())
    
    while not mutated and codons_ind.__len__() != 0:
        
        rnd_ind       = codons_ind.pop(randint(0,codons_ind.__len__()-1))
        rnd_codon     = codons[rnd_ind]
        rnd_codon_cai = codons_cai[rnd_ind] 
        
        #select alternative codons
        if   keep_aa == True and direction == '+':
            alt_codons = [c for c in aa2codon_table[codon2aa_table[rnd_codon]] if c!= rnd_codon and cai_table[c] > rnd_codon_cai and codon2aa_table[c]!='stop']
        elif keep_aa == True and direction == '-':
            alt_codons = [c for c in aa2codon_table[codon2aa_table[rnd_codon]] if c!= rnd_codon and cai_table[c] < rnd_codon_cai and codon2aa_table[c]!='stop']
        elif keep_aa == False and direction == '+':
            alt_codons = list(k for k,v in cai_table.iteritems() if v>rnd_codon_cai and codon2aa_table[k]!='stop')
        elif keep_aa == False and direction == '-':
            alt_codons = list(k for k,v in cai_table.iteritems() if v<rnd_codon_cai and codon2aa_table[k]!='stop')
            
        if alt_codons.__len__() != 0:                      
            mutated = True
            new_codon = choice(alt_codons)
            #print "new: " + str(new_codon)                
        
    if mutated == False:
        sys.stderr.write("SimpleCAIOperator: Not able to mutate sequence\n")
        return None
    else:              
        #print "CAI operator: old_codon -> " + str(rnd_codon)
        #print "CAI operator: new_codon -> " + str(new_codon)    
        real_codon_pos = mutableCodonsPosition[rnd_ind]
        codon_position = (real_codon_pos-cai_range[0])/3
        all_codons = analyzeCodons(sequence,range(cai_range[0],cai_range[1]+1,3))[0]             
        all_codons[codon_position]=new_codon

                                     
        new_seq = sequence[:cai_range[0]] + ''.join(c for c in all_codons) + sequence[cai_range[1]+1:]
        return new_seq               
    
def hammingDistance(seq1,seq2):
    score_nt = 0

    for i in range(0,len(seq1)):
        if seq1[i] != seq2[i]:
            score_nt += 1
            
    return score_nt
    
def SimpleHydropathyIndexOperator(sequence, hi_range, keep_aa, mutable_region, cds_regions, direction = '+'):
    '''
        Operator that given a sequence, mutates the sequence to change CAI
            sequnce: sequence 
            cai_range - start and end position to calculate cai - a tuple in the form (start, end)  
            mutable_region - a list with all bases that can be mutated
            cds_regions - a list of pairs with begin and end of CDSs - example: [(0,100), (150,300)]            
            direction: either increase ('+') or decrease ('-') CAI 
    '''    
    mutated = False
    mutableCodonsPosition = [c for c in range(hi_range[0],hi_range[1],3) if set([c,c+1,c+2]).issubset(mutable_region)]
                    
    if len(mutableCodonsPosition) == 0:
        sys.stderr.write("SimpleHydropathyIndexOperator: No codons available for mutation\n")
        return None
    
    result = analyzeCodons(sequence,mutableCodonsPosition,data_table = hydropathy_index_table)
    codons = (result[0])    
    codons_hi = (result[1])    
    codons_ind = range(0,codons.__len__())
    
    while not mutated and codons_ind.__len__() != 0:
        
        rnd_ind       = codons_ind.pop(randint(0,codons_ind.__len__()-1))
        rnd_codon     = codons[rnd_ind]
        rnd_codon_hi  = codons_hi[rnd_ind] 
        
        #select alternative codons
        if   keep_aa == True and direction == '+':
            alt_codons = [c for c in codons_list if c!= rnd_codon and hydropathy_index_table[c] > rnd_codon_hi and codon2aa_table[c]!='stop' and hammingDistance(c, rnd_codon)==1]
        elif keep_aa == True and direction == '-':
            alt_codons = [c for c in codons_list if c!= rnd_codon and hydropathy_index_table[c] < rnd_codon_hi and codon2aa_table[c]!='stop' and hammingDistance(c, rnd_codon)==1]
        elif keep_aa == False and direction == '+':
            alt_codons = list(k for k,v in cai_table.iteritems() if v>rnd_codon_hi and codon2aa_table[k]!='stop')
        elif keep_aa == False and direction == '-':
            alt_codons = list(k for k,v in cai_table.iteritems() if v<rnd_codon_hi and codon2aa_table[k]!='stop')
            
        if alt_codons.__len__() != 0:                      
            mutated = True
            new_codon = choice(alt_codons)
            #print "new: " + str(new_codon)
        
    if mutated == False:
        sys.stderr.write("SimpleHydropathyIndexOperator: Not able to mutate sequence\n")
        return None
    else:          
        all_codons = analyzeCodons(sequence,range(hi_range[0],hi_range[1]+1,3))[0]             
        all_codons[rnd_ind]=new_codon 
                                     
        new_seq = sequence[:hi_range[0]] + ''.join(c for c in all_codons) + sequence[hi_range[1]+1:]
        return new_seq

######################
##
#  Validation Functions
##
######################

def look_for_RBS(seq):
    max_score = -99
    max_pos   = -1

    #print seq
        
    #look for putative start sites
    for i in range(0,len(seq)-2,1):
        codon = seq[i:(i+3)] 
        if codon in ('atg','ttg','gtg'): #there is a start codon
            start = max([0,i - 17]) #not farther than 17 nucleotides
            stop  = max([0,i - 4])  #at least 3 nucleotides away
            
            aux_seq = seq[start:stop]                        
            res = analyze_pwm_score(aux_seq,pwm_rbs)
            score = res['MotifScore']
            position = res['MotifPosition']
            
            #print aux_seq , "\t" , score , "\t" , position
            
            #return start codon
            position=i
            #print str(i) + "\t" + aux_seq + "\t" + str([score, start + position])
            
            if score >= max_score:
                max_score = score 
                #max_pos   = start + position
                #return start codon
                max_pos   = position
                
    return (max_score, max_pos)

def look_for_promoters(seq):
    #look for putative promoters                    
    res = analyze_pwm_score(seq,pwm_pro15)
    score1 = res['MotifScore']
    position1 = res['MotifPosition']
    res = analyze_pwm_score(seq,pwm_pro16)
    score2 = res['MotifScore']
    position2 = res['MotifPosition']
    res = analyze_pwm_score(seq,pwm_pro17)
    score3 = res['MotifScore']
    position3 = res['MotifPosition']
    res = analyze_pwm_score(seq,pwm_pro18)
    score4 = res['MotifScore']
    position4 = res['MotifPosition']
    res = analyze_pwm_score(seq,pwm_pro19)
    score5 = res['MotifScore']
    position5 = res['MotifPosition']
    
    all_scores = [(score1, position1, 15),(score2, position2, 16),(score3, position3, 17),(score4, position4, 18),(score5, position5, 19)]         
                                               
    return max(all_scores, key=lambda x: x[0])

def look_for_terminators(seq):
    #use transtermHP
    #create input files
    file1 = 'test' + str(randint(0000,1000)) + '.fa';
    file2 = 'test' + str(randint(0000,1000)) + '-fake.coords';
    
    system("echo '>seq1\n" + str(seq) + "' > tmp/transterm_files/"+file1)
    system("echo 'fakegene1    1 2    seq1\nfakegene1    " + str(len(seq)-1) + " " + str(len(seq)) +"    seq1' > tmp/transterm_files/"+file2)
    #run transtermhp
    output = Popen("./3rdParty/transterm/transterm -p 3rdParty/transterm/expterm.dat tmp/transterm_files/"+file1+" tmp/transterm_files/"+file2+" 2> tmp/transterm_files/err.txt | awk /TERM/'{print $8}'", stdout=PIPE, stderr=None, shell=True).stdout.read()

    if output == '':
        return 0
    else:
        #print "Terminator output: " + output
        all_term_scores = [float(score) for score in output.split()] 
        return max(all_term_scores)   

#### WHAT's THIS for?

def get_coli_cds():
    """
    Return of dictionary of E. coli CDSs indexed by name
    """
    cds = {}
    bname = {}
    eck = {}
    h = open('/Users/cambray/Dropbox/WorkStuff/Sequences/Genomes/k12_all_cds.csv')
    h.readline()
    for l in h:
        l = l.strip()
        gene, blattner, eckid, seq = l.split(",")
        cds[gene] = seq
        bname[gene] = blattner
        eck[gene] = eckid
    h.close()
    return [cds, bname, eck]

def clean_cds(cds):
    to_drop = []
    for gene in cds[0]:
        if len(cds[0][gene])%3 != 0:
            to_drop.append(gene)
        elif cds[0][gene][1:3] != "tg" or cds[0][gene][0:3] == "ctg":
            to_drop.append(gene)
        else:
            for i in range(0,len(cds[0][gene])-3,3):
                if cds[0][gene][i:i+3] in scod:
                    to_drop.append(gene)
                    break
    for gene in to_drop:
        gene, cds[0].pop(gene)
    return cds

def get_core_cds(cds, core_path="/Users/cambray/Dropbox/WorkStuff/Sequences/Genomes/core_genome_coli_29_MicroScope_11_10_26.csv"):
    core = []
    h = open(core_path)
    h.readline()
    for l in h:
        l = l.strip()
        splitline = l.split(",")
        core.append(splitline[5].replace("'",""))
    h.close()
    cds[0] = dict([(k, v) for (k, v) in cds[0].items() if k in core])
    return cds        

def regression_tai(cds, out_path="", win=21, length=90):
    if out_path:
        h = open(out_path, "w")
    h.write("gene,"+",".join(["s%i" % i for i in range(2, length/3+2)])+","+",".join(["w%i" % i for i in range(2, length/3-win+3)])+",slope,intercept,r2, sumsq, merror\n")
    for gene in cds:
        seq = cds[gene][3:length+3].lower()
        if not len(seq)<length:
            score, smoothie = analyze_tai(seq, window=win)
            reg    = lin(range(2,len(smoothie)+2), smoothie)
            if out_path:
                h.write("%s,%s,%s,%s\n" % (gene,
                                           ",".join([str(n) for n in score]),
                                           ",".join([str(n) for n in smoothie]),
                                           ",".join([str(n) for n in reg])))
    if out_path:
        h.close()

def bneck_tai(cds, out_path="", win=21, data=tai_tuller):
    bname = cds[1]
    cds   = cds[0] 
    if out_path:
        h = open(out_path, "w")
        h.write("gene,bname,abs_pos,abs_strength,rel_pos,rel_strength,gene_length,gene_strength\n")
    for gene in cds:
        if len(cds[gene])/3 > win:
            scores, smooths = analyze_tai(cds[gene][:-3], window=win, method="harmonic", data=tai_tuller)
            bneck = max(smooths) 
            pos   = smooths.index(bneck)
            gene_tai = sum([1/v for v in scores])/len(scores)
            if out_path:
                h.write("%s,%s,%i,%.3f,%.3f,%.3F,%i,%.3f\n" % (gene, bname[gene], pos ,bneck, float(pos)/(len(cds[gene])/3-win), bneck/gene_tai, len(scores), gene_tai))
    if out_path:
        h.close()
    return
                
                
    
def construct_pssm(cds, length=90, out_path="", prob=None):
    """
    Construct Position Specific Scoring Matrices with log-likelihood values
    length: size of analyzed region from start, in bp  (sequences that are not this size are discarded)
    prob : a dict of bases with a priori expected probabilities  
    """
    cds = cds[0]
    if not prob:
        prob = {"a":0.25, "t":0.25, "g":0.25, "c":0.25}
    m = {"a":[0]*length, "t":[0]*length, "g":[0]*length, "c":[0]*length}
    tot_gene = 0.0
    for gene in cds:
        if len(cds[gene]) >= length:
            tot_gene += 1
            for i in range(length):
                m[cds[gene][i]][i] += 1
    for k in m:
        m[k] = [log((v/tot_gene)/prob[k]) for v in m[k]]
    if out_path:
        h = open(out_path, "w")
        h.write(","+",".join([str(i) for i in range(1,length+1)])+"\n")
        for b in ["a", "t", "g", "c"]:
            h.write(b+","+",".join(["%.2f" % v for v in m[b]])+"\n")
        h.close()
    return m


if __name__ == '__main__':
    pass