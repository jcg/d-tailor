'''
Created on Jan 30, 2012

@author: jcg
'''

from Functions import pwm_score, randomSequence
import Data

def pwmScoresSample(pwm="",sampleSize=1000):
    """
        Generates 'sampleSize' random samples and score them using PWM in pwmFile 
    """
    scoresList = []
    motifsize = len(pwm[pwm.keys()[0]])
    
    for _ in xrange(sampleSize):
        random_seq = randomSequence(motifsize)
        scoresList.append(pwm_score(random_seq, pwm))
        
    return scoresList 
        
    
    
if __name__ == "__main__":    
    pwm_lacI = pwmScoresSample(Data.pwm_lacI,50000)
    print pwm_lacI , "\n", max(pwm_lacI) ,  min(pwm_lacI), sum((pwm_lacI))/50000 
    
