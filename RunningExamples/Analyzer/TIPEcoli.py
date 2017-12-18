'''
Created on Jan 30, 2012

@author: jcg
'''

import Solution 
from Features import CAI,Bottleneck,Structure,NucleotideContent,HydropathyIndex
from Data import project_dir

#read CSV file with genes
wt_csv_file = project_dir+"/testFiles/auxiliaryData/taniguchi_pa_dataset.csv"
#wt_csv_file = 'random_60nt.csv'

f = open(wt_csv_file, 'r')
f.readline() #read header
#HEADER
print 'name,cdsCAI,cdsBottleneckPosition,cdsBottleneckRelativeStrength,cdsNucleotideContentAT,cdsHydropathyIndex,utrCdsStructureMFE,fivepCdsStructureMFE,threepCdsStructureMFE,PA,mrna,TE'



while True:
    values = f.readline().replace('\"','').replace('\n','').split(',')
 
    if values == ['']:
        break
    
    gene_name = values[0]
    gene_seq  = values[1].lower()
    gene_pa = values[2]
    gene_mrna = values[3]
    gene_te = values[4]
    
    if len(gene_seq[49:]) % 3 == 0: #and len(gene_seq[49:]) >= 100:    
        gene_solution = Solution.Solution(sol_id=gene_name, sequence=gene_seq) 
        
        #NucleotideContent
        nc_obj = NucleotideContent.NucleotideContent(solution=gene_solution,label="cds",args= { 'ntcontent_range' : (52,69) } )
        at_obj = NucleotideContent.NucleotideContentAT(nc_obj)
        nc_obj.add_subfeature(at_obj)

        #CAI
        cai_obj = CAI.CAI(solution=gene_solution,label="cds",args= {  'cai_range' : (52,len(gene_seq)+1) } )
        
        #Bottleneck
        bot_obj = Bottleneck.Bottleneck(solution=gene_solution, label="cds", args = {'bottleneck_range' : (49,len(gene_solution.sequence)) })
        bot_pos = Bottleneck.BottleneckPosition(bot_obj)
        bot_rls = Bottleneck.BottleneckRelativeStrength(bot_obj)
        bot_obj.add_subfeature(bot_pos)
        bot_obj.add_subfeature(bot_rls)
        
        #Hydropathy Index
        try:
            hi_obj = HydropathyIndex.HydropathyIndex(solution=gene_solution,label="cds",args= {  'hi_range' : (76,108) } )            
        except:
            gene_solution.scores['cdsHydropathyIndex']="NA"        
                        
        #Structures
        #[-30,30]
        st1_obj = Structure.Structure(solution=gene_solution,label="utrCds",args= { 'structure_range' : (19,78), 'mutable_region' : range(92,188),'cds_region' : (49,len(gene_solution.sequence)), 'keep_aa' : False } )
        st_mfe = Structure.StructureMFE(st1_obj)
        st_ds  = Structure.StructureDoubleStranded(st1_obj)
        st_ss  = Structure.StructureSingleStranded(st1_obj)
        #st_acc = Structure.StructureEnsembleAccessibility(st1_obj)
        st1_obj.add_subfeature(st_mfe)
        st1_obj.add_subfeature(st_ds)
        st1_obj.add_subfeature(st_ss)
        #st1_obj.add_subfeature(st_acc)    
        #[01,60]
        st2_obj = Structure.Structure(solution=gene_solution,label="fivepCds",args= { 'structure_range' : (49,108), 'mutable_region' : range(92,188),'cds_region' : (49,len(gene_solution.sequence)), 'keep_aa' : False } )
        st_mfe = Structure.StructureMFE(st2_obj)
        st_ds  = Structure.StructureDoubleStranded(st2_obj)
        st_ss  = Structure.StructureSingleStranded(st2_obj)
        #st_acc = Structure.StructureEnsembleAccessibility(st2_obj)
        st2_obj.add_subfeature(st_mfe)
        st2_obj.add_subfeature(st_ds)
        st2_obj.add_subfeature(st_ss)
        #st2_obj.add_subfeature(st_acc)
        #[31,90]
        st3_obj = Structure.Structure(solution=gene_solution,label="threepCds",args= { 'structure_range' : (79,138), 'mutable_region' : range(92,188),'cds_region' : (49,len(gene_solution.sequence)), 'keep_aa' : False } )
        st_mfe = Structure.StructureMFE(st3_obj)
        st_ds  = Structure.StructureDoubleStranded(st3_obj)
        st_ss  = Structure.StructureSingleStranded(st3_obj)
        #st_acc = Structure.StructureEnsembleAccessibility(st3_obj)
        st3_obj.add_subfeature(st_mfe)
        st3_obj.add_subfeature(st_ds)
        st3_obj.add_subfeature(st_ss)
        #st3_obj.add_subfeature(st_acc)
        
        gene_solution.add_feature(cai_obj)
        gene_solution.add_feature(bot_obj)
        gene_solution.add_feature(nc_obj)
        gene_solution.add_feature(hi_obj)
        gene_solution.add_feature(st1_obj) 
        gene_solution.add_feature(st2_obj)
        gene_solution.add_feature(st3_obj)
        
        print gene_solution.solid + ',' + \
              str(gene_solution.scores['cdsCAI']) + "," + \
              str(gene_solution.scores['cdsBottleneckPosition']) + "," + \
              str(gene_solution.scores['cdsBottleneckRelativeStrength']) + "," + \
              str(gene_solution.scores['cdsNucleotideContentAT']) + "," + \
              str(gene_solution.scores['cdsHydropathyIndex']) + "," + \
              str(gene_solution.scores['utrCdsStructureMFE']) + "," + \
              str(gene_solution.scores['fivepCdsStructureMFE']) + "," + \
              str(gene_solution.scores['threepCdsStructureMFE']) + "," + str(gene_pa) + "," + str(gene_mrna) + "," + str(gene_te)                                           

    else:
        print gene_name + ',NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA'