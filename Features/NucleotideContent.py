'''
Created on Nov 16, 2011

@author: jcg
'''

from Features.Feature import Feature
import Functions
from uuid import uuid4

class NucleotideContent(Feature):
    """
    Nucleotide Content Feature
        solution - solution where nucleotide content should be computed
        label - some label to append to the name
        hi_range - start and end position to calculate nucleotide content - a tuple in the form (start, end)  
        mutable_region - a list with all bases that can be mutated
        cds_region - a pair with begin and end of CDSs - example: (0,100)
        keep_aa - boolean option indicating if in the design mode amino acids should be kept        
    """
    def __init__(self, nucleotideContentObject = None, solution=None, label="", args = { 'ntcontent_range' : (0,9), 
                                                                                         'mutable_region' : None, 
                                                                                         'cds_region' : None, 
                                                                                         'keep_aa' : True }):
        if nucleotideContentObject == None: #create new instance
            #General properties of feature
            Feature.__init__(self, solution=solution, label=label)
            #Specifics of this Feature
            self.ntcontent_range    = args['ntcontent_range']
            self.sequence           = solution.sequence[self.ntcontent_range[0]:self.ntcontent_range[1]+1]
            self.mutable_region     = args['mutable_region'] if args.has_key('mutable_region') else solution.mutable_region
            self.cds_region         = args['cds_region']    if args.has_key('cds_region') else solution.cds_region
            self.keep_aa            = args['keep_aa']        if args.has_key('keep_aa') else solution.keep_aa
            self.set_scores()
            self.set_level()
        else:
            Feature.__init__(self, nucleotideContentObject)
            self.ntcontent_range    = nucleotideContentObject.ntcontent_range
            self.sequence           = nucleotideContentObject.sequence
            self.mutable_region     = nucleotideContentObject.mutable_region
            self.cds_region         = nucleotideContentObject.cds_region
            self.keep_aa            = nucleotideContentObject.keep_aa
            self.scores             = nucleotideContentObject.scores
    
    def set_scores(self, scoring_function = Functions.analyze_ntcontent):
        self.scores = Functions.appendLabelToDict(scoring_function(self.sequence), self.label)
            
    def mutate(self, operator=Functions.SimpleNtContentOperator):
        if not self.targetInstructions:
            return None
        new_seq = operator(self.solution.sequence, self.targetInstructions['direction'], self.nucleotides, self.mutable_region, self.cds_region, keep_aa=self.keep_aa)
        if not new_seq:
            return None                
        return Solution.Solution(sol_id=str(uuid4().int), sequence=new_seq, cds_region = self.cds_region, mutable_region = list(self.mutable_region), parent=self.solution, design=self.solution.designMethod)

class NucleotideContentAT(NucleotideContent):
    """
    Check AT content
    """
    def __init__(self,  nucleotideContentObject):
        NucleotideContent.__init__(self,nucleotideContentObject)
        self.nucleotides = ['a','t']
        self.set_level()
        
class NucleotideContentGC(NucleotideContent):
    """
    Check GC content
    """
    def __init__(self,  nucleotideContentObject):
        NucleotideContent.__init__(self,nucleotideContentObject)
        self.nucleotides = ['g','c']
        self.set_level()
        
class NucleotideContentA(NucleotideContent):
    """
    Check A content
    """
    def __init__(self,  nucleotideContentObject):
        NucleotideContent.__init__(self,nucleotideContentObject)                
        self.nucleotides = ['a']
        self.set_level()
        
class NucleotideContentT(NucleotideContent):
    """
    Check T content
    """
    def __init__(self,  nucleotideContentObject):
        NucleotideContent.__init__(self,nucleotideContentObject)
        self.nucleotides = ['t']
        self.set_level()
        
class NucleotideContentG(NucleotideContent):
    """
    Check G content
    """
    def __init__(self,  nucleotideContentObject):
        NucleotideContent.__init__(self,nucleotideContentObject)
        self.nucleotides = ['g']
        self.set_level()
        
class NucleotideContentC(NucleotideContent):
    """
    Check C content
    """
    def __init__(self,  nucleotideContentObject):
        NucleotideContent.__init__(self,nucleotideContentObject)
        self.nucleotides = ['c']
        self.set_level()
    
import Solution   