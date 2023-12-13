import os
import sys
import openpyxl as opx
from openpyxl import Workbook
import pandas as pd




        
def _update_excel(pi, requested_by, department, gene, edit, edit_size, injection_core, cage_number, ngs_date, success_num, submitted_num, notes):
    #open excel file
    animal_model_xl_dir = "Z:\ResearchHome\Groups\millergrp\home\common\Python\_pete\_AGES ( Automated Graphical Emailer System)\complete_animal_models.xlsx"
    animal_model_xl ="complete_animal_models.xlsx" 
    
    '''
    Investigator	Department	Gene name	Project type	insert size	Start date	End date	Days to completion	Initials of who completed project	Notes	CAGE Project #	

    '''
    
    workbook = opx.load_workbook(animal_model_xl_dir)
    
    mice_sheet = workbook.active
    
    investigators = str(pi +"," + requested_by)
    
    new_row = [investigators, department, gene, edit, edit_size, '', ngs_date, '', 'PMH', 'N/A', cage_number, notes]
    
    max_row = mice_sheet.max_row
   
    
    
    mice_sheet.append(new_row)
    workbook.save(animal_model_xl)
    workbook.close()
    
    df = pd.read_excel(animal_model_xl)
    print(f"max-row: {max_row}")
    print(df.tail())
    
    
    #go to 'mice' worksheet
    
    #find last row
    
    #update pi, dept, gene, edit, size, start?, end?, initials of who completed, cage#, 
    
    #move to 2nd sheet
    
    #update ngs date, cage#, submitted_num, num_success, injection, core           
    
    
    
    return
    



table_rows = [('2665904', 'Kanneganti, Thirumala-Devi', 'Baskaran, Yogi', 'Baskaran, Yogi', 'CAGE84', 'RIPK3-mBFP2_KI', '18', 'Well of Plate', 'Tail Snip/Toe Snip', 'CBT', 'Yes', '1', '4', 'KO', '5', 'TCU', '111423', ['CAGE84_hDGCR6L_F_R_short'],'testing_notes')]

        #table_rows = self.data



for row in table_rows:
    proj_data = list(row)
    
    print(proj_data)
    #unpack proj_data
    (srm_number,
    pi,
    requested_by,
    entered_by,
    cage_number,
    gene,
    sample_num,
    sample_format,
    sample_type,
    department,
    success,
    success_num,
    submitted_num,
    edit,
    edit_size,
    injection_core,
    ngs_date,
    cage_programs,
    notes) = proj_data
    
    #update excel sheet
_update_excel(pi, requested_by, department, gene, edit, edit_size, injection_core, cage_number, ngs_date, success_num, submitted_num, notes)