'''
Created on Jan 30, 2012

@author: jcg
'''

import math
from os import path

def pfm2pwm(pfmFile=""):
    """
        Converts a PFM in text format and vertical orientation to a dictionary in PWM format
    """
    pwm = {}
    
    project_dir = path.dirname(__file__)
    
    f = open(project_dir+"/"+pfmFile, 'r')
    keys = f.readline().split()
    
    for i in range(0,len(keys)):
        keys[i] = keys[i].lower()

    
    #initialize hash map
    for i in range(0,len(keys)):
        pwm[keys[i]] = []
    
    while True:
        values = f.readline().split()
        for i in range(0,len(values)):
            values[i] = eval(values[i])
            
        n = sum(values)        
        
        if not values:
            break
        else:            
            for i in range(0,len(keys)):
                base_weight = math.log((values[i]+math.sqrt(n)*0.25) / (n+math.sqrt(n)) / 0.25, 2)
                pwm[keys[i]].append(base_weight)
                
    return pwm
    
    
if __name__ == "__main__":
    #pwm_sd = pfm2pwm("../testFiles/pwm/e.coli/sd/ecoli_sd.pwm")
    #print pwm_sd
    
    pwm_lacI = pfm2pwm("../testFiles/pwm/e.coli/operators/ecoli_lacI.pwm")
    print pwm_lacI
    
    #pwm_35 = pfm2pwm("../testFiles/pwm/e.coli/promoters/ecoli_s70_35.pwm")
    #print pwm_35
    
    #pwm_10 = pfm2pwm("../testFiles/pwm/e.coli/promoters/ecoli_s70_10.pwm")
    #print pwm_10
    
    
