'''
Created on Nov 29, 2011

@author: jcg
'''

import sqlite3
import os
import sys

def DB2FASTA(db_file,split=False):
    print "Generating FASTA file(s) for " , db_file , "...",
    
    if not os.path.exists(db_file):
        print "File does not exist"
        return 0

    #### Create connection to DB
    con = sqlite3.connect(db_file)
    con.isolation_level = None
    con.row_factory = sqlite3.Row
    cur = con.cursor()

    #### Write FASTA file
    if not split:
        output_sol = open(db_file+".generated_solutions.fa",'w')
         
    #rows    
    cur.execute("SELECT * FROM generated_solution")
    
    result = cur.fetchall()
    
    for row in result:
        
        if split:
            output_sol = open(db_file+"."+row['generated_solution_id']+".fa",'w')
        
        seq = row['sequence'].upper()

        output_sol.write(">"+row['generated_solution_id']+" | "+row['des_solution_id']+"\n")
        output_sol.write(seq+"\n")
        
        if split:
            output_sol.close()
    
    if not split:
        output_sol.close()
    
    con.close()
    
    print "Done!"
    
if __name__ == '__main__':
    pass
    
    if len(sys.argv) == 2:
        db_file = sys.argv[1]
    else:
        db_file = "../testFiles/outputFiles/bpec.sqlite"

    DB2FASTA(db_file,split=False)
