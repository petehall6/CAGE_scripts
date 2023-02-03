import pandas as pd


print("Connecting to dataframe.  This will take a few moments.")
data_store = "tiny.h5"
df = pd.read_hdf(data_store)
#strips hyphens for easier parsing.  Trying to remove hyphen in the hdf5 will give an unique index error
#stripping here does not impact performance    
df.columns= df.columns.str.replace("-","")

cell_lines = list((input("Enter cell lines separated by (,): ").strip(" ")).split(","))   

gene = input("Enter gene: ").upper()


for lines in cell_lines:
    #filters take a function so lambda was used to create anyonmous functions.
    cell_matches = list(filter(lambda x: lines in x, df['cell_line']))

    #if cell and gene matches found, list approximate matches
    if len(cell_matches) > 0:
        while True:
            try:
                print("Cell line matches:\n")
                for (i,j) in enumerate(cell_matches, start=1):
                    print(i,j)
                    
                cell_choice = int(input("\nEnter cell choice number: "))
                if cell_choice > len(cell_matches):
                    print("Invalid choice. Please re-select cell line.")
                    continue      
                
            except ValueError:
                    print("Invalid choice. Please re-select cell line and gene.")
                    continue
            
            cell_selection = cell_matches[cell_choice-1]          
            print("Cell choice: ",cell_selection)         

            confirm = input("Press 'y' to confirm: ").upper()

            if confirm == 'Y' or 'YES':
            #This finds copy number

                try:
                    copyNumber = str(df.loc[df['cell_line'] == cell_selection])
                    #copy number returns large string with cell and table information.  slice though to get just copy number        
                    #copyNumberSlice = copyNumber[5:8]
                    #make it pretty
                    top_btm_boarder = "*"*33
                    line_boarder = "*"+" "*31 + "*"
                    print(top_btm_boarder)
                    print(line_boarder)
                    print(f"*\tNumber of copies: a number goes here"+" "*3+"*")
                    print(line_boarder)
                    print(top_btm_boarder)
                    break
    
                except KeyError:
                    print("Problem finding cell line")
            elif confirm != "Y":
                print("bye")
            """   #incorrect match choosen by user.  reset and repick
            elif confirm != 'Y':
                print("Resetting cell line and gene.  Please re-pick")
                del cell_line, gene
                getCopyNumber(df) """
                    
    else:
        if not cell_matches:
            print("Cell line not found in database.")
        else:
            print("Gene was not found in database.")

#TODOGENE LOOP GOES HERE
