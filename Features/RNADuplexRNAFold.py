'''
Created on Jan 24, 2012

@author: jcg
'''

from Features.Feature import Feature
import Functions

class RNADuplexRNAFold(Feature):
    """
    RNADuplex Feature
        solution1 - molecule1 to hybridize
        solution2 - molecule2 to hybridize
        label - some label to append to the name
        rnaMolecule1/2region - start and end position to hybridize molecules - a tuple in the form (start, end)  
        NOTE: Design mode not implemented for this feature        
    """
    def __init__(self, duplexObject=None, solution1=None, solution2=None, label="", args = { 'rnaMolecule1region' : None,                                                                                            
                                                                                             'rnaMolecule2region' : None }):
        if duplexObject == None:           
            #General properties of feature
            Feature.__init__(self, solution=solution1, label=label)
            #Specifics of this Feature
            self.solution            = solution1
            self.solution2           = solution2
            self.duplexFile          = str(solution1.solid) + '-' + str(solution2.solid) + label
            self.rnaMolecule1region  = args['rnaMolecule1region']
            self.rnaMolecule1seq     = self.solution.sequence[self.rnaMolecule1region[0]:(self.rnaMolecule1region[1]+1)]
            self.mutable_region      = args['mutable_region1']   if args.has_key('mutable_region1') else None
            self.cds_regions         = args['cds_regions1']      if args.has_key('cds_regions1') else None
            self.keep_aa             = args['keep_aa1']          if args.has_key('keep_aa1') else False
            self.rnaMolecule2region  = args['rnaMolecule2region']
            self.rnaMolecule2seq     = self.solution2.sequence[self.rnaMolecule2region[0]:(self.rnaMolecule2region[1]+1)]
            self.mutable_region2     = args['mutable_region2']   if args.has_key('mutable_region2') else None
            self.cds_regions2        = args['cds_regions2']      if args.has_key('cds_regions2') else None
            self.keep_aa2            = args['keep_aa2']          if args.has_key('keep_aa2') else False
            self.set_scores()
            self.set_level()
        else:
            #General properties of feature
            Feature.__init__(self, duplexObject)
            #Specifics of this Feature
            self.solution           = duplexObject.solution
            self.solution2           = duplexObject.solution2
            self.duplexFile          = duplexObject.duplexFile
            self.rnaMolecule1region  = duplexObject.rnaMolecule1region
            self.rnaMolecule1seq     = duplexObject.rnaMolecule1seq 
            self.mutable_region     = duplexObject.mutable_region
            self.cds_regions        = duplexObject.cds_regions
            self.keep_aa            = duplexObject.keep_aa
            self.rnaMolecule2region  = duplexObject.rnaMolecule2region
            self.rnaMolecule2seq     = duplexObject.rnaMolecule2seq
            self.mutable_region2     = duplexObject.mutable_region2
            self.cds_regions2        = duplexObject.cds_regions2
            self.keep_aa2            = duplexObject.keep_aa2
            self.scores              = duplexObject.scores
            
                           
    
    def set_scores(self, scoring_function = Functions.analyze_duplex_structure_rnafold):
        self.scores.update(Functions.appendLabelToDict(scoring_function(self.rnaMolecule1seq,self.rnaMolecule2seq,self.duplexFile), self.label))        
        pass
    
    def mutate(self):        
        return Feature.randomMutation(self, mutable_region=self.mutable_region)


class RNADuplexRNAFoldRibosome(RNADuplexRNAFold):
    """
    RNADuplexRibosome Feature
        solution1 - molecule1 being hybridized (given by user)
        solution2 - molecule2 being hybridized (16S rRNA)
        label - some label to append to the name
        rnaMolecule1/2region - start and end position to hybridize molecules - a tuple in the form (start, end)  
        mutable_region - a list with all bases that can be mutated
        cds_region - a pair with begin and end of CDSs - example: (0,100)
        keep_aa - boolean option indicating if in the design mode amino acids should be kept        
    """
    def __init__(self, solution1=None, label="", args = { 'rnaMolecule1region' : None,
                                                          'mutable_region1' : None, 
                                                          'cds_region1' : None,
                                                          'keep_aa1' : True } ):
        ribRNA = Solution.Solution(sol_id = '16S', sequence='acctcctta')
        args.update({ 'rnaMolecule2region' : (0,8) })
        RNADuplexRNAFold.__init__(self,solution1=solution1, solution2=ribRNA,label=label, args=args)
        self.solution = solution1
        self.keep_aa          = args['keep_aa'] if args.has_key('keep_aa') else solution1.keep_aa
        self.mutable_region   = args['mutable_region'] if args.has_key('mutable_region') else solution1.mutable_region
        self.cds_regions      = args['cds_regions'] if args.has_key('cds_regions') else solution1.mutable_region
            
    
class RNADuplexRNAFoldMFE(RNADuplexRNAFold):
    """
    Manipulate the duplex MFE
    """
    def __init__(self, duplexObject, label = ""):
        RNADuplexRNAFold.__init__(self,duplexObject=duplexObject)
        self.label = self.label + label
        self.set_scores()
        self.set_level()      
        
    def set_scores(self, scoring_function=Functions.analyze_duplex_mfe_rnafold):    
        self.scores.update(Functions.appendLabelToDict(scoring_function(self.duplexFile), self.label))
        

import Solution