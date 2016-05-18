'''
Created on Nov 29, 2011

@author: jcg
'''

import sqlite3
import sys
#import time

def DBStatistics(db_file):
    """
        Prints the following DB statistics and info: Name    SolutionsGenerated    DesignsFound    SolutionsGenerated/DesignsFound
    """
    #### Create connection to DB
    con = sqlite3.connect(db_file)
    con.isolation_level = None
    con.row_factory = sqlite3.Row
    cur = con.cursor()

    #get statistics
    cur.execute("select (select count(1)  from desired_solution where status = 'DONE') as 'Designs Found', (select count(1)  from generated_solution) as 'Solutions Generated'")
    
    #print header
    print "DB name\tNum. of generated solutions\tNum. of combinations\tRatio"
    result = cur.fetchone()
    if result['Solutions Generated'] != 0:
            print db_file , "\t" , result['Solutions Generated'] , "\t" , result['Designs Found'] , "\t" , result['Designs Found']/float(result['Solutions Generated'])
        
    con.close()
    
    
if __name__ == '__main__':
    pass
    if len(sys.argv) == 2:
        db_file = sys.argv[1]
    else:
        db_file = "../testFiles/outputFiles/bpec.sqlite"
    
    #One-time
    DBStatistics(db_file)
    
    #Polling every 5 seconds
    #max_it = 5
    #while max_it != 0:
    #    DBStatistics(db_file)
    #    max_it -= 1
    #    time.sleep(5)