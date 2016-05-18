'''
Created on Feb 27, 2012

@author: jcg
'''

from Features.Feature import Feature
import Functions
from uuid import uuid4

class HydropathyIndex(Feature):
    """
    HydropathyIndex Feature
        solution - solution where hydropathy index should be computed
        label - some label to append to the name
        hi_range - start and end position to calculate hydropathy index - a tuple in the form (start, end)  
        mutable_region - a list with all bases that can be mutated
        cds_region - a pair with begin and end of CDSs - example: (0,100)
        keep_aa - boolean option indicating if in the design mode amino acids should be kept        
    """
    def __init__(self, hiObject = None, solution = None, label="", args = { 'hi_range' : (0,59), 
                                                                             'mutable_region' : None, 
                                                                             'cds_region' : None, 
                                                                             'keep_aa' : True }):
        if hiObject == None: #create new instance
            #General properties of feature
            Feature.__init__(self, solution=solution, label=label)
            #Specifics of this Feature
            self.hi_range          = args['hi_range']
            self.sequence          = solution.sequence[self.hi_range[0]:(self.hi_range[1]+1)]
            self.mutable_region    = args['mutable_region'] if args.has_key('mutable_region') else solution.mutable_region
            self.cds_region        = args['cds_region']    if args.has_key('cds_region') else solution.cds_region
            self.keep_aa           = args['keep_aa']        if args.has_key('keep_aa') else solution.keep_aa
            self.set_scores()
            self.set_level()
        else: #copy instance
            Feature.__init__(self, hiObject)
            self.hi_range           = hiObject.hi_range
            self.sequence           = hiObject.sequence
            self.mutable_region     = hiObject.mutable_region
            self.cds_region         = hiObject.cds_region
            self.keep_aa            = hiObject.keep_aa
            self.codons_hi          = hiObject.codons_hi
            self.scores             = hiObject.scores

    
    def set_scores(self, scoring_function=Functions.analyze_hydropathy):     
        self.scores[self.label+"HydropathyIndex"] = scoring_function(self.sequence)
                                                                     
    def mutate(self, operator=Functions.SimpleHydropathyIndexOperator):
        if not self.targetInstructions:
            return None
        new_seq = operator(self.solution.sequence, self.hi_range, self.keep_aa, self.mutable_region, self.cds_region, self.targetInstructions['direction'])
        if not new_seq:
            return None             
        return Solution.Solution(sol_id=str(uuid4().int), sequence=new_seq, cds_region = self.cds_region, mutable_region = self.mutable_region, parent=self.solution, design=self.solution.designMethod)
        
import Solution