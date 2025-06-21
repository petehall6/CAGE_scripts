import pandas as pd

#Finds copy number using Depmap downloadable csv by converting to hd5 for quicker loading and parsing.
#returns "Number of copies: {copy number}" boardered by *
#dataset found @ https://depmap.org/portal/download/custom/
#Copy Number (Absolute)
#PMH 1/2023
#TODO multiple cell line entry?  Save results to core projects folder?

from cage_config import txtcolors



KI_list = ['Codon Change','Correction','Gene Tagging','Point Mutation', 'Targeted Deletion', 'Targeted Integration', 'Other']
top_btm_boarder = "*"*33
line_boarder = "*"+" "*31 + "*"
top_btm_boarder_flag = "*"*79
line_boarder_flag = "*"+" "*77 + "*"

def getCopyNumberInLine(cell_line, gene, proj_objective):
    
    print("Checking copy number information.\n")
    cell_line =str(cell_line).strip()
    gene = str(gene).strip().replace("-","")
    proj_objective = str(proj_objective).strip()
    
    
    df = pd.read_hdf("copy_number_data.h5")
    
    
    #strips hyphens for easier parsing.  Trying to remove hyphen in the hdf5 will give an unique index error
    #stripping here does not impact performance    
    df.columns= df.columns.str.replace("-","")

    #filters take a function so lambda was used to create anyonmous functions.
    cell_match = list(filter(lambda x: cell_line == x, df['cell_line']))
    gene_match = list(filter(lambda x: gene == x, df.columns))

    if len(cell_match) > 0:
        try:
            print("\n\nCell line match: ",cell_match[0])
            print("Gene match: ",gene_match[0])                      

        
            #This parses through df
            copyNumber = str(df.loc[df['cell_line'] == cell_match[0]][gene_match[0]])
            #copy number returns large string with cell and table information.  slice though to get just copy number        
            copyNumberSlice = copyNumber[5:8]
            
            #flag knock in projects at 2 or more genes
            if int(copyNumberSlice) >= 2 and proj_objective in KI_list:
                print(f"{txtcolors.WARNING}")
                print(top_btm_boarder_flag)
                print(line_boarder_flag)
                print(f"* BE ADVISED. THERE ARE{copyNumberSlice} COPIES PRESENT.  IS THIS TOO MANY FOR A KNOCK IN?  *")
                print(line_boarder_flag)
                print(top_btm_boarder_flag)
                input("PRESS ENTER TO CONTINUE.\n")
                print(f"{txtcolors.END}")
                
                
            #no flag
            else:
                print(f"{txtcolors.GOOD}")                       
                print(top_btm_boarder)
                print(line_boarder)
                print(f"*\tNumber of copies: {copyNumberSlice}"+" "*3+"*")
                print(line_boarder)
                print(top_btm_boarder)
                #input("Press any key to continue.\n")
                print(f"{txtcolors.END}")
                return
        
        except:
            print(f"{txtcolors.WARNING}")
            print("Problem finding cell line: {} or gene: {}".format(cell_line, gene))
            print("Probably an alias issue.")
            print(f"{txtcolors.END}")
            
    else:
        print(f"{txtcolors.FINE}")
        print("Copy number information not found.")
        print(f"{txtcolors.END}")
        return

def getCopyNumberTools():
    
    #print("Connecting to dataframe.  This will take a few moments.")
    data_store = "copy_number_data.h5"
    df = pd.read_hdf(data_store)
    #strips hyphens for easier parsing.  Trying to remove hyphen in the hdf5 will give an unique index error
    #stripping here does not impact performance    
    df.columns= df.columns.str.replace("-","")
    
    cell_line = input("Enter cell line: ").upper()
    gene = input("Enter gene: ").upper()

    #filters take a function so lambda was used to create anyonmous functions.
    cell_matches = list(filter(lambda x: cell_line in x, df['cell_line']))
    gene_matches = list(filter(lambda x: gene in x, df.columns))

    #if cell and gene matches found, list approximate matches
    if len(cell_matches) > 0 and len(gene_matches) > 0:
        while True:
            try:
                print("\nCell line matches:\n")
                for (i,j) in enumerate(cell_matches, start=1):
                    print(i,j)
                    
                cell_choice = int(input("\nEnter cell choice number: "))
                if cell_choice > len(cell_matches):
                    print("Invalid choice. Please re-select cell line.")
                    continue          
                
                print("\nGene matches: \n")
                for (i,j) in enumerate(gene_matches, start=1):
                    print(i,j)
                    
                gene_choice = int(input("\nEnter gene choice number: "))
                if gene_choice > len(gene_matches):
                    print("\nInvalid choice. Please re-select gene.\n")
                    continue
            except ValueError:
                    print("Invalid choice. Please re-select cell line and gene.")
                    continue
            
            cell_selection = cell_matches[cell_choice-1]
            gene_selection = gene_matches[gene_choice-1]

            print("\n\nCell choice: ",cell_selection)
            print("Gene choice: ",gene_selection)

            confirm = input("Press 'y' to confirm: ").upper()
                        

            if confirm == 'Y' or 'YES':
            #This finds copy number

                try:
                    copyNumber = str(df.loc[df['cell_line'] == cell_selection][gene_selection])
                    #copy number returns large string with cell and table information.  slice though to get just copy number        
                    copyNumberSlice = copyNumber[5:8]
                    #make it pretty
                    print(top_btm_boarder)
                    print(line_boarder)
                    print(f"*\tNumber of copies: {copyNumberSlice}"+" "*3+"*")
                    print(line_boarder)
                    print(top_btm_boarder)
                    return
        
                except KeyError:
                    print("Problem finding cell line: {} or gene: {}".format(cell_line, gene))
            #incorrect match choosen by user.  reset and repick
            elif confirm != 'Y':
                print("Resetting cell line and gene.  Please re-pick")
                del cell_line, gene
                getCopyNumberTools()
                    
    else:
        if not cell_matches:
            print("Cell line not found in database.")
        else:
            print("Gene was not found in database.")
        return None




