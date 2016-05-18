'''
Created on Apr 4, 2012

@author: jcg
'''

from sqlite3 import connect,Row
import sys


def DBKinetics(db_file,resolution=50):
    """
        Prints a time series of desired solutions found as a function of number of generated solutions
    """
    #### Create connection to DB
    con = connect(db_file)
    con.isolation_level = None
    con.row_factory = Row
    cur = con.cursor()
    
    cur.execute("select count(1) from generated_solution")
    total_sol = cur.fetchone()[0]

    #print "Generated Solutions \t Desired Solutions Found"

    it = range(0,total_sol,total_sol/resolution)
    #it = range(0,total_sol,100)
    it.append(total_sol)

    #get statistics
    print "Num. of generated solutions\tNum. of combinations"
    for i in it:
        cur.execute("select count(DISTINCT code) as c from (select des_solution_id as code from generated_solution LIMIT " + str(i) + ") where code <> ''")
        #cur.execute("select min(rowid) as gen_sol from generated_solution where des_solution_id<>'' group by des_solution_id order by rowid")
        
        result = cur.fetchone()
    
        print i,"\t",result['c']
        #print result['c']
    
    con.close()
    
    
def DBKinetics2(db_file):
    """
        Prints a time series of generated solutions per target combination
    """
    #### Create connection to DB
    con = connect(db_file)
    con.isolation_level = None
    con.row_factory = Row
    cur = con.cursor()
    
    cur.execute("select min(rowid) as gen_sol from generated_solution where des_solution_id<>'' group by des_solution_id order by gen_sol")
        
    result = cur.fetchall()
    
    print "Num. of combinations\tNum. of generated solutions"
    for i in range(0,len(result)):
        print i+1,"\t",result[i]['gen_sol']
        #print result[i]['gen_sol']
    
    con.close()
    
if __name__ == '__main__':
    pass
    if len(sys.argv) == 2:
        db_file = sys.argv[1]
    else:
        db_file = "../testFiles/outputFiles/bpec.sqlite"
    
    DBKinetics(db_file) #number of combinations found per generated solution
    #DBKinetics2(db_file) #number of generated solutions per combination found 
        