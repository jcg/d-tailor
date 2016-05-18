'''
Created on Nov 16, 2011

@author: jcg
'''

from Features.Feature import Feature
import Functions
from uuid import uuid4

class Motif(Feature):
    """
    Motif Feature
        solution - solution where motif score should be searched
        label - some label to append to the name
        pwm - position weights matrix (PWM) of the motif to be analyzed
        motif_range - start and end position where we should look for the motif - a tuple in the form (start, end)  
        mutable_region - a list with all bases that can be mutated
        cds_region - a pair with begin and end of CDSs - example: (0,100)
        keep_aa - boolean option indicating if in the design mode amino acids should be kept        
    """
    def __init__(self, motifObject=None, solution=None, label="", args = { 'pwm' : None,
                                                         'motif_range' : (0,9), 
                                                         'mutable_region' : None, 
                                                         'cds_region' : None, 
                                                         'keep_aa' : True }):
        
        if motifObject == None: #create new instance        
            #General properties of feature
            Feature.__init__(self, solution=solution, label=label)
            #Specifics of this Feature
            self.pwm                = args['pwm']
            self.motif_range        = args['motif_range']
            self.sequence           = solution.sequence[self.motif_range[0]:self.motif_range[1]+1]
            self.mutable_region     = [i-self.motif_range[0] for i in set(range(self.motif_range[0],self.motif_range[1]+1)) & set(args['mutable_region'])] if args.has_key('mutable_region') else solution.mutable_region
            self.cds_region         = args['cds_region']    if args.has_key('cds_region') else solution.cds_region
            self.keep_aa            = args['keep_aa']        if args.has_key('keep_aa') else solution.keep_aa        
            self.set_scores()
            self.set_level()        
        else: #copy instance
            Feature.__init__(self, solution=motifObject.solution, label=motifObject.label)
            self.pwm                = motifObject.pwm
            self.motif_range        = motifObject.motif_range 
            self.sequence           = motifObject.sequence 
            self.mutable_region     = motifObject.mutable_region 
            self.cds_region         = motifObject.cds_region 
            self.keep_aa            = motifObject.keep_aa         
            self.scores             = motifObject.scores
    
    def set_scores(self, scoring_function = Functions.analyze_pwm_score):
        self.scores = Functions.appendLabelToDict(scoring_function(self.sequence,self.pwm), self.label)


class MotifScore(Motif):
    """
    Manipulate the motif score
    """
    def __init__(self, motifObject, label = ""):
        Motif.__init__(self,motifObject)
        self.label = self.label + label
        self.set_level()
    
    def mutate(self, operator=Functions.SimplePWMScoreOperator):
        if not self.targetInstructions:
            return None
        
        new_seq = list(self.solution.sequence)
        mutated_seq = operator(self.sequence, self.pwm, self.targetInstructions['direction'], self.mutable_region, keep_aa=self.keep_aa)
        
        if mutated_seq == None:
            return None
        else:
            new_seq[self.motif_range[0]:self.motif_range[1]+1] = list(mutated_seq)
        new_seq = "".join(new_seq) 
                
        return Solution.Solution(sol_id=str(uuid4().int), sequence=new_seq, cds_region = self.cds_region, mutable_region = self.mutable_region, parent=self.solution, design=self.solution.designMethod)
    
    def defineTarget(self,desiredSolution):
        '''
        Function that determines if a target wasn't hit, and if not updates targetDirections 
        '''
        if desiredSolution == None:
            return True
    
        targetLevel = desiredSolution[self.label+self.__class__.__name__+"Level"]
        currentLevel = self.level
        
        # Check if there's an associated target position
        if desiredSolution.has_key(self.label+"MotifPositionLevel") and self.solution.levels[self.label+"MotifPositionLevel"] != desiredSolution[self.label+self.__class__.__name__+"Level"]:
            #first we need to set the right position level
            return False
        else:
            if currentLevel == targetLevel:
                return False
            elif currentLevel > targetLevel:
                self.targetInstructions['direction'] = '-' #decrease
            else:
                self.targetInstructions['direction'] = '+' #increase
                
            return True
    
class MotifPosition(Motif):
    """
    Manipulate the motif score
    """
    def __init__(self, motifObject, label = ""):
        Motif.__init__(self,motifObject)
        self.label = self.label + label
        self.set_level()
    
    def defineTarget(self,desiredSolution):
        '''
        Function that determines if a target wasn't hit, and if not updates targetDirections 
        '''
        if desiredSolution == None:
            return True
        
        targetLevel = desiredSolution[self.label+self.__class__.__name__+"Level"]
        currentLevel = self.level            
        
        if currentLevel != targetLevel:
            desiredMotifPosition = self.solution.designMethod.thresholds[self.label+self.__class__.__name__][targetLevel]
            self.targetInstructions['position'] = desiredMotifPosition
            return True
        else:
            return False

                    
    def mutate(self, operator=Functions.SimplePWMPositionOperator):
        if not self.targetInstructions:
            return None
        
        new_seq = list(self.solution.sequence)
        mutated_seq = operator(self.sequence, self.pwm, self.targetInstructions['position'], self.mutable_region, keep_aa=self.keep_aa)
        if mutated_seq == None:
            return None
        else:
            new_seq[self.motif_range[0]:self.motif_range[1]+1] = list(mutated_seq)
        new_seq = "".join(new_seq) 
                
        return Solution.Solution(sol_id=str(uuid4().int), sequence=new_seq, cds_region = self.cds_region, mutable_region = self.mutable_region, parent=self.solution, design=self.solution.designMethod)
    
import Solution   