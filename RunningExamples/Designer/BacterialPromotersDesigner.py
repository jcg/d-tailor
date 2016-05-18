'''
Created on Dec 22, 2012

@author: jcg
'''

import Data
from SequenceDesigner import SequenceDesigner
from Features import NucleotideContent,Motif
from DesignOfExperiments.Design import FullFactorial,Optimization,CustomDesign


class BacterialPromotersDesigner(SequenceDesigner):
    
    def __init__(self, name, seed, design, dbfile, createDB=True):
        SequenceDesigner.__init__(self, name, seed, design, dbfile, createDB)
        self.max_sol_counter = 100000
        self.max_iterations  = 100
        
    def configureSolution(self, solution):
        '''
        Solution configuration
        '''
                
        if solution.sequence == None:
            return 0
        
        #Populate solution with desired features
        solution.mutable_region=range(0,len(solution.sequence)) # whole region
        solution.cds_region = None
        solution.keep_aa = False
        
        # The entire promoter + 5' utr will have 75 nucleotides (and TSS will be at nucleotide 50)
        #   UP element     -35          spacer        -10                 5' utr
        # UUUUUUUUUUUUUUU MMMMMM SSSSSSSSSSSSSSSSSSS DDDDDD AAAA OOOOOOOOOOOOOOOOOOOOOOOOO
        
        #UP element (-50,-36)
        up_obj = NucleotideContent.NucleotideContent(solution=solution,label="up",args= {  'ntcontent_range' : (0,14), 'mutable_region' : range(0,15) } )
        upat_obj = NucleotideContent.NucleotideContentAT(up_obj)
        #-35 motif (-35,-30)
        m35_obj = Motif.Motif(solution=solution,label="m35",args= {  'motif_range' : (15,20), 'pwm' : Data.pwm_35, 'mutable_region' : range(15,21) } )
        m35score_obj = Motif.MotifScore(m35_obj)
        #-10 motif (-10,-5)
        m10_obj = Motif.Motif(solution=solution,label="m10",args= {  'motif_range' : (40,45), 'pwm' : Data.pwm_10 , 'mutable_region' : range(40,46) } )
        m10score_obj = Motif.MotifScore(m10_obj)
        #lacI operator (-6, +25)
        mlacI_obj = Motif.Motif(solution=solution,label="mlacI",args= {  'motif_range' : (0,74), 'pwm' : Data.pwm_lacI , 'mutable_region' : range(0,74) } )
        mlacIscore_obj = Motif.MotifScore(mlacI_obj)
        mlacIpos_obj = Motif.MotifPosition(mlacI_obj)
            
        solution.add_feature(upat_obj)
        solution.add_feature(m35score_obj)
        solution.add_feature(m10score_obj)                    
        solution.add_feature(mlacIscore_obj)
        solution.add_feature(mlacIpos_obj)

if __name__ == '__main__':
    #Seed sequence from which mutants will be derived
    seed='cggcttaactcgagagctgacctgcttctacctagccttacagtggtaacgacccaatctgcgtagcgcaacgca'
    
    #Design Methodology and thresholds
    design_param = {  "upNucleotideContentAT" : { 'type' : 'REAL' , 
                                                  'thresholds' : { '1': (0,0.25), '2': (0.25,0.75), '3': (0.75,1)} },
                      "m35MotifScore"         : { 'type' : 'REAL' , 
                                                  'thresholds' : { '1': (-12.0,-6.81), '2': (-6.81,0.63), '3': (0.63,11.0)} },
                      "m10MotifScore"         : { 'type' : 'REAL' , 
                                                  'thresholds' : { '1': (-12.0,-8.19), '2': (-8.19,0.32), '3': (0.32,11.0)} },
                      "mlacIMotifScore"       : { 'type' : 'REAL' , 
                                                  'thresholds' : { '1': (-2,2), '2': (2,6), '3': (6,12)} },
                      "mlacIMotifPosition"    : { 'type' : 'INTEGER' , 
                                                  'thresholds' : { '1': 10, '2': 19, '3': 44} } }        
            
    design = FullFactorial(["upNucleotideContentAT","m35MotifScore","m10MotifScore","mlacIMotifScore","mlacIMotifPosition"],design_param)
    #design = Optimization(["upNucleotideContentAT","m35MotifScore","m10MotifScore","mlacIMotifScore","mlacIMotifPosition"],design_param,"1.3.3.1.1")
    #design = CustomDesign(["upNucleotideContentAT","m35MotifScore","m10MotifScore","mlacIMotifScore","mlacIMotifPosition"],design_param,["1.1.1.1.1","2.2.2.2.2","3.3.3.3.3"])
    
    tfec_designer = BacterialPromotersDesigner("bpec", seed, design, Data.project_dir+"/testFiles/outputFiles/bpec", createDB=True)
    tfec_designer.run(selection="directional")
    
    
    
    
    
      
