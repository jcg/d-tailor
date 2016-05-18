'''
Created on Dec 15, 2012

@author: jcg
'''

from Features.Feature import Feature
import Functions
from uuid import uuid4


class StructureProb(Feature):
    """
    StructureProb Feature
        solution - solution where structure accessibility probability should be computed
        label - some label to append to the name of structure file
        structure_range - start and end position to calculate structure - a tuple in the form (start, end)  
        mutable_region - a list with all bases that can be mutated
        mutable_region - a list with all bases that can be mutated
        cds_region - a pair with begin and end of CDSs - example: (0,100)
        keep_aa - boolean option indicating if in the design mode amino acids should be kept
    """    
    def __init__(self, structureObject = None, solution = None, label="", args = { 'structure_range' : (0,59),
                                                                                   'acc_region' : [],
                                                                                   'window' : 50, 
                                                                                   'mutable_region' : None, 
                                                                                   'cds_region' : None, 
                                                                                   'keep_aa' : True }):
        
        #General properties of feature
        Feature.__init__(self, solution=solution, label=label)
        #Specifics of this Feature
        self.structurefile      = solution.solid + label
        self.structure_range    = args['structure_range']
        self.acc_region         = args['acc_region']        if args.has_key('acc_region') else []
        self.window             = args['window']            if args.has_key('window') else 50
        self.sequence           = solution.sequence[self.structure_range[0]:(self.structure_range[1]+1)]
        self.mutable_region     = args['mutable_region']    if args.has_key('mutable_region') else solution.mutable_region
        self.cds_region         = args['cds_region']        if args.has_key('cds_region') else solution.cds_region
        self.keep_aa            = args['keep_aa']           if args.has_key('keep_aa') else solution.keep_aa
        self.set_scores()
        self.set_level()                    
        
                            
    def set_scores(self, scoring_function=Functions.analyze_structure_prob): 
        #compute structure
        self.scores.update(Functions.appendLabelToDict(scoring_function(self.sequence,self.structurefile,self.window,self.acc_region), self.label))
                                            
    def mutate(self, operator=Functions.SimpleStructureOperator):
        if not self.targetInstructions:
            return None    
        
        if self.scores.has_key(self.label+'StructureProbList'):
            ss_bases = []
            ds_bases = []
            for base in self.scores[self.label+'StructureProbList'].keys():
                if self.scores[self.label+'StructureProbList'][base] < 0.5:
                    ds_bases.append(int(base))
                else:
                    ss_bases.append(int(base))
        else: 
            ss_bases = ds_bases = None        
          
        new_seq = operator(self.solution.sequence, self.solution.solid+self.label, self.structure_range, list(self.mutable_region), self.cds_region, self.targetInstructions['direction'], ss_bases=ss_bases, ds_bases=ds_bases)
        if not new_seq:
            return None   
        
        return Solution.Solution(sol_id=str(uuid4().int), sequence=new_seq, cds_region = self.cds_region, mutable_region = self.solution.mutable_region, parent=self.solution, design=self.solution.designMethod)
    
import Solution
