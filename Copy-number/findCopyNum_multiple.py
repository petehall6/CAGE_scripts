import pandas as pd

#Finds copy number using Depmap downloadable csv by converting to hd5 for quicker loading and parsing.
#returns "Number of copies: {copy number}" boardered by *
#dataset found @ https://depmap.org/portal/download/custom/
#Copy Number (Absolute)

#TODO multiple cell line entry?  Save results to core projects folder?


def initDataFrame():
    print("Connecting to dataframe.  This will take a few moments.")
    data_store = "copy_number_data.h5"
    df = pd.read_hdf(data_store)
    #strips hyphens for easier parsing.  Trying to remove hyphen in the hdf5 will give an unique index error
    #stripping here does not impact performance    
    df.columns= df.columns.str.replace("-","")
    
    return df

def getCopyNumber(df):
    
    cell_line = input("Enter cell line: ").upper()
    gene = input("Enter gene: ").upper()

    #filters take a function so lambda was used to create anyonmous functions.
    cell_matches = list(filter(lambda x: cell_line in x, df['cell_line']))
    gene_matches = list(filter(lambda x: gene in x, df.columns))

    #if cell and gene matches found, list approximate matches
    if len(cell_matches) > 0 and len(gene_matches) > 0:
        while True:
            try:
                print("Cell line matches:\n")
                for (i,j) in enumerate(cell_matches, start=1):
                    print(i,j)
                    
                cell_choice = int(input("\nEnter cell choice number: "))
                if cell_choice > len(cell_matches):
                    print("Invalid choice. Please re-select cell line.")
                    continue          

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

            print("Cell choice: ",cell_selection)
            print("Gene choice: ",gene_selection)

            confirm = input("Press 'y' to confirm: ").upper()
                        

            if confirm == 'Y' or 'YES':
            #This finds copy number

                try:
                    copyNumber = str(df.loc[df['cell_line'] == cell_selection][gene_selection])
                    #copy number returns large string with cell and table information.  slice though to get just copy number        
                    copyNumberSlice = copyNumber[5:8]
                    #make it pretty
                    top_btm_boarder = "*"*33
                    line_boarder = "*"+" "*31 + "*"
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
                getCopyNumber(df)
                    
    else:
        if not cell_matches:
            print("Cell line not found in database.")
        else:
            print("Gene was not found in database.")
        return None

def main():
    
    df = initDataFrame()

    cn = getCopyNumber(df)


if __name__ == "__main__":
    main()












