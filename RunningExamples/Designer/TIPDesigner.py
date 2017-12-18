'''
Created on Nov 12, 2012

@author: jcg
'''

import Functions
import sys
from SequenceDesigner import SequenceDesigner
from Features import CAI,Bottleneck,Structure,NucleotideContent,HydropathyIndex
from DesignOfExperiments.Design import FullFactorial,RandomSampling,Optimization, CustomDesign
from Data import codon2aa_table, project_dir
from random import randint

class TIPDesigner(SequenceDesigner):
    
    def __init__(self, name, seed, design, dbfile, createDB=True):
        SequenceDesigner.__init__(self, name, seed, design, dbfile, createDB)
        self.designErrors = {}
        
    def configureSolution(self, solution):
        '''
        Solution configuration
        '''
        if solution.sequence == None:
            return 0

        #Populate solution with desired features
        solution.mutable_region=range(92,188)
        solution.cds_region=(89,1022)
        #We want to keep AA
        keep_aa = True
    
        #CAI
        cai_obj = CAI.CAI(solution=solution,label="cds",args= {  'cai_range' : (89,187), 'mutable_region': range(92,188), 'cds_region':(89,1022), 'keep_aa' : keep_aa } )
        #Bottleneck
        bot_obj = Bottleneck.Bottleneck(solution=solution, label="cds", args = {'bottleneck_range' : (89,1022), 'mutable_region' : range(92,188), 'keep_aa' : keep_aa})
        bot_pos = Bottleneck.BottleneckPosition(bot_obj)
        bot_rls = Bottleneck.BottleneckRelativeStrength(bot_obj)
        bot_obj.add_subfeature(bot_pos)
        bot_obj.add_subfeature(bot_rls)
        #NucleotideContent (AT)
        nc_obj = NucleotideContent.NucleotideContent(solution=solution,label="cds",args= { 'ntcontent_range' : (92,109), 'mutable_region' : range(92,188),'cds_region' : (89,1022), 'keep_aa' : keep_aa } )
        at_obj = NucleotideContent.NucleotideContentAT(nc_obj)
        nc_obj.add_subfeature(at_obj)
        #Hydropathy Index
        hi_obj = HydropathyIndex.HydropathyIndex(solution=solution,label="cds",args= {  'hi_range' : (116,148), 'mutable_region': range(92,188), 'cds_region':(89,1022), 'keep_aa' : keep_aa } )    
    
        #Structures
        #[-30,30]
        st1_obj = Structure.Structure(solution=solution,label="utrCds",args= { 'structure_range' : (59,118), 'mutable_region' : range(92,188),'cds_region' : (89,1022), 'keep_aa' : keep_aa } )
        st_mfe = Structure.StructureMFE(st1_obj)
        st_ds  = Structure.StructureDoubleStranded(st1_obj)
        st_ss  = Structure.StructureSingleStranded(st1_obj)
        #st_acc = Structure.StructureEnsembleAccessibility(st1_obj)
        st1_obj.add_subfeature(st_mfe)
        st1_obj.add_subfeature(st_ds)
        st1_obj.add_subfeature(st_ss)
        #st1_obj.add_subfeature(st_acc)    
        #[01,60]
        st2_obj = Structure.Structure(solution=solution,label="fivepCds",args= { 'structure_range' : (89,148), 'mutable_region' : range(92,188),'cds_region' : (89,1022), 'keep_aa' : keep_aa } )
        st_mfe = Structure.StructureMFE(st2_obj)
        st_ds  = Structure.StructureDoubleStranded(st2_obj)
        st_ss  = Structure.StructureSingleStranded(st2_obj)
        #st_acc = Structure.StructureEnsembleAccessibility(st2_obj)
        st2_obj.add_subfeature(st_mfe)
        st2_obj.add_subfeature(st_ds)
        st2_obj.add_subfeature(st_ss)
        #st2_obj.add_subfeature(st_acc)
        #[31,90]
        st3_obj = Structure.Structure(solution=solution,label="threepCds",args= { 'structure_range' : (119,178), 'mutable_region' : range(92,188),'cds_region' : (89,1022), 'keep_aa' : keep_aa } )
        st_mfe = Structure.StructureMFE(st3_obj)
        st_ds  = Structure.StructureDoubleStranded(st3_obj)
        st_ss  = Structure.StructureSingleStranded(st3_obj)
        #st_acc = Structure.StructureEnsembleAccessibility(st3_obj)
        st3_obj.add_subfeature(st_mfe)
        st3_obj.add_subfeature(st_ds)
        st3_obj.add_subfeature(st_ss)
        #st3_obj.add_subfeature(st_acc)
        
        solution.add_feature(cai_obj)
        solution.add_feature(bot_obj)
        solution.add_feature(nc_obj)
        solution.add_feature(hi_obj)
        solution.add_feature(st1_obj) 
        solution.add_feature(st2_obj)
        solution.add_feature(st3_obj)
        
    
    def validateSolution(self, solution):
        '''
        Solution validation tests
        '''
        return True
        
        if solution.sequence == None:
            return 0
        
        valid = True
        e_type=None
      
        variable_region = solution.sequence[solution.mutable_region[0]:(solution.mutable_region[-1]+1)]
                
        #No Promoters
        promoter_putative_region = solution.sequence[(solution.mutable_region[0]-30):(solution.mutable_region[-1]+31)]
        (score, position, spacer) = Functions.look_for_promoters(promoter_putative_region)
        if score >= 0.3990166: #0.95 percentile for Promoter PWM scores
            valid = False
            position = (solution.mutable_region[0]-30) + position
            e_type = "prom"
            sys.stderr.write("SolutionValidator: High Promoter score\n")
                
        #No RBS
        rbs_putative_region = solution.sequence[(solution.mutable_region[0]-12):(solution.mutable_region[-1]+11)]
        (score, position) = Functions.look_for_RBS(rbs_putative_region)
        if score >= 0.820515: #0.95 percentile for RBS PWM scores
            valid = False
            position = (solution.mutable_region[0]-12) + position
            e_type = "rbs"      
            sys.stderr.write("SolutionValidator: High RBS score: "+str(score)+" - "+str(position)+"\n")    
        
        #No Terminator
        term_putative_region = solution.sequence[(solution.mutable_region[0]-30):(solution.mutable_region[-1]+31)]
        score = Functions.look_for_terminators(term_putative_region)
        if score >= 90: #90% confidence from transtermHP
            e_type = "term"
            valid = False
            sys.stderr.write("SolutionValidator: High Terminator score\n")    
        
        #No stop codon
        for i in range(0,len(variable_region),3):
            if codon2aa_table[variable_region[i:i+3]]=='stop':
                sys.stderr.write("SolutionValidator: Stop codon found\n")
                e_type = "stopcodon"
                valid = False
                position = solution.mutable_region[0]+i
                break            
                    
        #No restriction enzymes
        rsite_putative_region = solution.sequence[(solution.mutable_region[0]-6):(solution.mutable_region[-1])]
        if 'ggtacc' in rsite_putative_region:
            e_type = "rsite"
            sys.stderr.write("SolutionValidator: Restriction enzyme found\n")
            position = (solution.mutable_region[0]-6)+rsite_putative_region.index('ggtacc')
            valid = False
        
        if 'ggatcc' in rsite_putative_region:
            e_type = "rsite"
            sys.stderr.write("SolutionValidator: Restriction enzyme found\n")
            position = (solution.mutable_region[0]-6)+rsite_putative_region.index('ggatcc')
            valid = False                
            
        #No bases run of more than 5nt
        repeat_putative_region = solution.sequence[(solution.mutable_region[0]-6):(solution.mutable_region[-1]+7)]
        for b_run in ('aaaaaa','tttttt','gggggg','cccccc'):
            if b_run in repeat_putative_region:
                sys.stderr.write("SolutionValidator: Bases run for more than 5 found\n")            
                valid = False
                position = (solution.mutable_region[0]-6) + repeat_putative_region.index(b_run)
                e_type = "repeats"


        solution.valid = valid
                
        if valid == True:
            return (valid,None,None)
        else:
            return (valid,position,e_type)

    def switchConfiguration(self,solution, keep_aa=True):
        if solution != None and solution.sequence != None:
            for feature in solution.features.values():
                feature.keep_aa = keep_aa
                for sub_feature in feature.subfeatures.values():
                    sub_feature.keep_aa = keep_aa


    def additionalConfigurationPreMutation(self, solution):
        '''
        This method is executed before the mutation iteration happens and can be used to set additional mutational properties
        
        Here we are implementing a stochastic process that decides if AA within CDS should be kept or not
        '''
        pass
        #change configuration to not keep_aa
        if randint(1, 100) <= 1: # 1% of the times, don't keep AA
            #Change features configuration
            self.switchConfiguration(solution, False)
        
    def additionalConfigurationPostMutation(self, solution):
        '''
        This method is executed after the mutation iteration happens and can be used to set additional mutational properties
        '''
        self.switchConfiguration(solution, True)

    
if __name__ == '__main__':
    #Exemplary Seed sequence from which mutants should be derived (find all seed in ./testFiles/auxiliaryData/TIP_seed_sequences.csv)
    seed = 'agggcccaagttcacttaaaaaggagatcaacaatgaaagcaatttaggtactgaaacatcttaatcatgcacataaggaggtaccataatgtcagtcggcgtttgttacccgaaatgctgctggccaggaacgcacgctacacaatactacgacctcaggactaaaccgcagattatcgttggaacggcaggtggatccggtggatcgggcgggtcatactaccatcaccaccaccaccatttggagtcggaaaatttatacttccaatccggttccgccggatctgctgccggttccggtgaatttagcaaaggagaagaacttttcaccggagtagtcccgattctggttgaattagatggtgatgttaatgggcacaaattttctgtccgtggagagggtgaaggtgatgctacaaacggaaaactcacccttaaattcatttgcacaacgggaaaactaccagtaccgtggccaacactggtcactactctgacctatggtgttcaatgcttttcccgttatccggatcacatgaaacggcatgactttttcaagagtgccatgcccgaaggttatgtacaggaacgcactatatctttcaaagatgacgggacctacaagacgcgtgctgaagtcaagtttgaaggtgatacccttgttaatcgtatcgagttaaagggtattgattttaaagaagatggaaacattcttggacacaaactcgagtacaactttaactcacacaatgtatacatcacggcagacaaacaaaagaatggaatcaaagctaacttcaaaattcgccacaacgttgaagatggttccgttcaactagctgagcattatcaacaaaatactccaattggagatggacctgtccttttaccagacaaccattacctgtcgacacaatctgtcctttcgaaagatcccaacgaaaagcgtgaccacatggtccttcttgagtttgtaactgctgctgggattacacatggcatggatgagctctacaaataa'

    #Design Methodology and thresholds
    design_param = {
                   "cdsCAI":                        { 'type' : 'REAL'    , 'thresholds' : { '1': (0,0.28), '2': (0.28,0.44), '3': (0.44,1)} },
                   "utrCdsStructureMFE":            { 'type' : 'REAL'    , 'thresholds' : { '1': (-32,-11.5), '2': (-11.5,-7.5), '3': (-7.5,0)} },
                   "fivepCdsStructureMFE":          { 'type' : 'REAL'    , 'thresholds' : { '1': (-32,-14), '2': (-14,-9), '3': (-9,0)} },
                   "threepCdsStructureMFE":         { 'type' : 'REAL'    , 'thresholds' : { '1': (-32,-16.5), '2': (-16.5,-11.5), '3': (-11.5,0)} },                   
                   "cdsBottleneckPosition":         { 'type' : 'INTEGER' , 'thresholds' : { '1': (1,33), '2': (34,290)} },
                   "cdsBottleneckRelativeStrength": { 'type' : 'REAL'    , 'thresholds' : { '1': (1.0,1.48), '2': (1.48,2.1)} },
                   "cdsNucleotideContentAT":        { 'type' : 'REAL'    , 'thresholds' : { '1': (0,0.58), '2': (0.58,1) } },                   
                   "cdsHydropathyIndex":            { 'type' : 'REAL'    , 'thresholds' : { '1': (-5,-1), '2': (-1,1), '3': (1,5)} },
                   }
    
    
    
    #In the TIP project cdsBottleneckPositionLevel = 2 and cdsBottleneckRelativeStrength = 2 is not possible so these should be removed from listDesigns
    featureList = ["cdsCAI","utrCdsStructureMFE","fivepCdsStructureMFE","threepCdsStructureMFE","cdsBottleneckPosition","cdsBottleneckRelativeStrength","cdsNucleotideContentAT","cdsHydropathyIndex"]
    design = FullFactorial(featureList,design_param)
    botPosIndx = featureList.index('cdsBottleneckPosition')
    botStrIndx = featureList.index('cdsBottleneckRelativeStrength')
    newListDesigns = [combination for combination in design.listDesigns if not (combination.split(".")[botPosIndx]=='2' and combination.split(".")[botStrIndx]=='2')]    
    design.listDesigns = newListDesigns 
    ######                
    
    tip_designer = TIPDesigner("tip", seed, design, project_dir+"/testFiles/outputFiles/tip_db", createDB=True)
    tip_designer.run(selection="directional")