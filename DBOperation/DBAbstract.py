'''
Created on Nov 1, 2012

@author: jcg
'''

class DBAbstract(object):
    '''
    Abstract class for DNA Sequence Designer Database
    '''

    def DBInit(self):
        '''
        Initialize database
        returns: Nothing
        '''
        pass
        
    def DBGetSolution(self, solution_id):
        '''
        Get solution given solution_id
        returns: a dictionary with a solution with all attributes
        '''
        pass
        
    def DBGetDesiredSolution(self):
        '''
        Get a desired solution that wasn't found yet
        returns: a dictionary with a desired solution or None
        '''
        pass
    
    def DBChangeStatusDesiredSolution(self,desired_solution_id,status):
        '''
        Change the status of a desired solution
        '''
        pass
        
    
    def DBGetClosestSolution(self, desiredSolution):
        '''
        Get a solution that is closer to the desired solution
        returns: a dictionary with a solution with all attributes
        '''
        pass
    
    def DBCheckDesign(self, desired_solution_id):
        '''
        Get the status of a solution to design
        returns: a boolean with the result of status == 'DONE'        
        '''
        pass
    
    def DBInsertSolution(self, solution, desired_solution_id):
        '''
        Insert solution into database
        returns: Nothing
        '''
        pass
    
    def DBCloseConnection(self):
        '''
        Closes connection to DB
        '''
        pass
        