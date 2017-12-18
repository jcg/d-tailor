'''
Created on Nov 16, 2011

@author: jcg
'''

from Features.Feature import Feature
import Functions
import sys
from uuid import uuid4
from random import choice

class Bottleneck(Feature):
    """
    Bottleneck Feature as described in doi:10.1186/gb-2011-12-2-r12
        solution - solution where bottleneck should be computed
        label - some label to append to the name
        bottleneck_range - start and end position to calculate bottleneck - a tuple in the form (start, end)  
        mutable_region - a list with all bases that can be mutated
        keep_aa - boolean option indicating if in the design mode amino acids should be kept
    """
    def __init__(self, bottleneckObject = None, solution= None, label="", args = { 'bottleneck_range' : (0,59), 
                                                                                   'mutable_region' : None , 
                                                                                   'keep_aa' : True}):
        if bottleneckObject == None: #create new instance
            #General properties of feature
            Feature.__init__(self, solution=solution, label=label)
            #Specifics of this Feature
            self.bottleneck_range = args['bottleneck_range']
            self.sequence         = solution.sequence[self.bottleneck_range[0]:self.bottleneck_range[1]+1]        
            self.mutable_region   = args['mutable_region'] if args.has_key('mutable_region') else solution.mutable_region
            self.cds_region       = args['cds_region']    if args.has_key('cds_region') else solution.cds_region
            self.keep_aa          = args['keep_aa']        if args.has_key('keep_aa') else solution.keep_aa             
            self.set_scores()
            self.set_level() 
            self.calculate_mutable_segments()      
        else: #copy instance
            Feature.__init__(self, bottleneckObject)            
            self.bottleneck_range = bottleneckObject.bottleneck_range
            self.sequence         = bottleneckObject.sequence        
            self.mutable_region   = bottleneckObject.mutable_region
            self.keep_aa          = bottleneckObject.keep_aa
            self.segmentMutation  = bottleneckObject.segmentMutation
            self.bot_scores       = bottleneckObject.bot_scores
            self.bot_smooth       = bottleneckObject.bot_smooth
            self.scores           = bottleneckObject.scores
                          
    
    def set_scores(self, scoring_function=Functions.analyze_bottleneck):
        self.scores     = {}
        #compute bottleneck
        res = scoring_function(self.sequence)
        self.bot_scores = res[0]
        self.bot_smooth = res[1]
        
    def calculate_mutable_segments(self):
        self.segmentMutation = {}   
        
        if self.mutable_region == None:
            return 
        
        for level in self.solution.designMethod.thresholds[self.label+"BottleneckPosition"].keys():
            segment = self.solution.designMethod.thresholds[self.label+"BottleneckPosition"][level]
            real_start = self.bottleneck_range[0] + (segment[0]-1)*3
            real_stop  = self.bottleneck_range[0] + (segment[1]-1)*3
            
            #check codons available for mutation
            mutableCodonsPosition = [c for c in range(real_start,real_stop+1,3) if set([c,c+1,c+2]).issubset(self.mutable_region)]
                
            self.segmentMutation[level] = mutableCodonsPosition                                                                         
    
    def mutate(self, operator=Functions.SimpleBneckOperator):
        if not self.targetInstructions:
            return None
                
        #check codons available for mutation
        mutableCodonsPosition = self.segmentMutation[self.targetInstructions["segment"]]
                
        if len(mutableCodonsPosition) == 0:
            sys.stderr.write("Bottleneck: No codons available for mutation - " + str(self.__class__.__name__) + "\n")
            return None
        
        mut_codons = [self.solution.sequence[cod_init:cod_init+3] for cod_init in mutableCodonsPosition]
        
        subseq = "".join(mut_codons) 
        mutated = operator(subseq, self.targetInstructions["direction"],distance = (0 if self.keep_aa==True else 0) )
        
        if mutated == None or len(mutated) != len(subseq) or mutated == subseq:
            return None

        #replace original sequence with mutated sequence
        new_seq = list(self.solution.sequence)
        # reconstruct final sequence after mutations
        j=0
        for i in mutableCodonsPosition:
            new_seq[i:i+3] = mutated[j:j+3]
            j += 3                    
        
        new_seq = "".join(new_seq)
        
        return Solution.Solution(sol_id=str(uuid4().int), sequence=new_seq, cds_region = self.solution.cds_region, mutable_region = self.mutable_region, parent=self.solution, design=self.solution.designMethod)    
    
class BottleneckPosition(Bottleneck):
    """
    Manipulate the position of bottleneck
    """
    def __init__(self,  bottleneckObject):
        Bottleneck.__init__(self, bottleneckObject)
        self.parentFeature = bottleneckObject
        self.set_scores()
        self.set_level()

    def set_scores(self, scoring_function=Functions.analyze_bottleneck_pos):    
        self.scores[self.label+"BottleneckPosition"] = scoring_function(self.sequence, self.bot_scores, self.bot_smooth)

    def defineTarget(self,desiredSolution):   
        self.targetInstructions = {}
        target_level = desiredSolution[self.label+'BottleneckPositionLevel']
        if target_level != self.level: #check if there is a level to achieve 
            
            levelsToMutate = [level for level in self.segmentMutation.keys() if self.segmentMutation[level] != []]           
            rndLevel = choice(levelsToMutate)
  
            if rndLevel == target_level:
                self.targetInstructions['segment'] = rndLevel
                self.targetInstructions['direction'] = '+'
            else: 
                self.targetInstructions['segment'] = rndLevel
                self.targetInstructions['direction'] = '-'
                
            return True

        return False
    

class BottleneckRelativeStrength(Bottleneck):
    """
    Manipulate the strength of bottleneck
    """
    def __init__(self, bottleneckObject):
        Bottleneck.__init__(self, bottleneckObject)
        self.parentFeature = bottleneckObject
        self.set_scores()
        self.set_level()
        
    def set_scores(self, scoring_function=Functions.analyze_bottleneck_rel_strength):    
        self.scores[self.label+"BottleneckRelativeStrength"] = scoring_function(self.sequence, self.bot_scores, self.bot_smooth)

    def defineTarget(self,desiredSolution):
        
        if desiredSolution == None:
            return True
        
        self.targetInstructions = {}
        #check if there is a target
        target_level_for_strength = desiredSolution[self.label+'BottleneckRelativeStrengthLevel']
        bneck_position_level      = self.parentFeature.subfeatures[self.label+'BottleneckPosition'].level       
        target_level_for_position = desiredSolution[self.label+'BottleneckPositionLevel']
        if  target_level_for_strength != self.level and bneck_position_level == target_level_for_position:
            
            other_segments = [segment for segment in self.segmentMutation.keys() if self.segmentMutation[segment] != [] and segment != target_level_for_position]           
            
            if  int(target_level_for_strength)-int(self.level) > 0:
                
                if len(other_segments) == 0 or choice([True,False,False]): #increase bottleneck strength in target segment
                    self.targetInstructions['segment'] = target_level_for_position
                    self.targetInstructions['direction'] ='+' #increase
                else: #decrease bottleneck strength in segments other than target                    
                    self.targetInstructions['segment'] = choice(other_segments)
                    self.targetInstructions['direction'] ='-' #decrease                                                
            else:
                if len(other_segments) == 0 or choice([True,False,False]): #decrease bottleneck strength in target segment
                    self.targetInstructions['segment'] = target_level_for_position
                    self.targetInstructions['direction'] ='-' #decrease
                else: #increase bottleneck strength in segments other than target                    
                    self.targetInstructions['segment'] = choice(other_segments)
                    self.targetInstructions['direction'] ='+' #increase
                
            return True

        return False                  

import Solution