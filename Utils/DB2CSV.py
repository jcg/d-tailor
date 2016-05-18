'''
Created on Feb 29, 2012

@author: jcg
'''

import sqlite3
import sys

def DB2CSV(db_file):
    print "Generating CSV files for " , db_file , "...",

    #### Create connection to DB
    con = sqlite3.connect(db_file)
    con.isolation_level = None
    con.row_factory = sqlite3.Row
    cur = con.cursor()

    #### Write generated solutions
    output_sol = open(db_file+".generated_solutions.csv",'w')
    
    #header
    cur.execute("PRAGMA table_info(generated_solution)")
    output_sol.write(",".join([(v['name']) for v in cur.fetchall()])+"\n") 
     
    #rows    
    cur.execute("SELECT * FROM generated_solution")
    
    result = cur.fetchall()
    
    for r in result:
        output_sol.write(",".join([str(v) for v in r])+"\n")
    
    output_sol.close()
    
    #### Write design map    
    output_design = open(db_file+".design_list.csv",'w')
    
    cur.execute("PRAGMA table_info(desired_solution)")
    output_design.write(",".join([(v['name']) for v in cur.fetchall()])+"\n") 
     
    #rows    
    cur.execute("SELECT * FROM desired_solution")
    
    result = cur.fetchall()
    
    for r in result:
        output_design.write(",".join([str(v) for v in r])+"\n")
    
    output_design.close()
    
    con.close()
    
    print "Done!"
    
    
if __name__ == '__main__':
    pass

    if len(sys.argv) == 2:
        db_file = sys.argv[1]
    else:
        db_file = "../testFiles/outputFiles/bpec.sqlite"
        
    DB2CSV(db_file)