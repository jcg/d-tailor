'''
Created on Nov 1, 2012

@author: jcg
'''

from SequenceAnalyzer import SequenceAnalyzer
from Features import CAI,StructureRNAFold,RNADuplexRNAFold
from Functions import validateCDS
from Data import project_dir

class TranslationFeaturesEcoliAnalyzer(SequenceAnalyzer):
    
    '''
    Initializes class that analyzes sequence features 
    '''
    
    def __init__(self, input_file, input_type):        
        SequenceAnalyzer.__init__(self,input_file,input_type)         
            
    def configureSolution(self, solution):
        solution.valid = validateCDS(solution.sequence[49:])
        
        if solution.valid:        
            #CAI
            cai_obj = CAI.CAI(solution=solution,label="cds",args= {  'cai_range' : (49,len(solution.sequence)) } )
            
            #Look for RBS
            dup_obj1 = RNADuplexRNAFold.RNADuplexRNAFoldRibosome(solution1=solution, label="sd16s", args = { 'rnaMolecule1region' : (25,48) })
            dup_mfe = RNADuplexRNAFold.RNADuplexRNAFoldMFE(dup_obj1)
            dup_obj1.add_subfeature(dup_mfe)
            
            #MFE [-30,30]
            st1_obj = StructureRNAFold.StructureRNAFold(solution=solution,label="utr",args= { 'structure_range' : (49-30,49+30) } )
            st_mfe = StructureRNAFold.StructureRNAFoldMFE(st1_obj)
            st1_obj.add_subfeature(st_mfe)                    
            
            solution.add_feature(cai_obj)
            solution.add_feature(dup_obj1)
            solution.add_feature(st1_obj)
            
    
    def outputStart(self):
        print "gene_name,sd_hyb_energy,mfe_structure,cai"
    
    def output(self, solution):
        if solution.valid:
            print ",".join([solution.solid,str(solution.scores['sd16sRNADuplexRNAFoldMFE']),str(solution.scores['utrStructureRNAFoldMFE']),str(solution.scores['cdsCAI'])])

if __name__ == '__main__':
    
    seqAnalyzerTest = TranslationFeaturesEcoliAnalyzer(project_dir+"/testFiles/genomes/partial_ecoli_genome.csv","CSV")
    #seqAnalyzerTest = TranslationFeaturesEcoliAnalyzer(project_dir+"/testFiles/genomes/ecoli_genome.csv","CSV")
    seqAnalyzerTest.run()
    
