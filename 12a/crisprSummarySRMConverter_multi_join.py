import os
import pandas as pd

from xlutils.copy import copy
import pathlib
import xlrd
import xlwt
import shutil

# NOTE: If running standalone, check main() for User variables for os.paths
DEFAULT_TEMPLATE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "Converted-Picked-Template.xls"
)


def CRISPRSummaryToSeries(targetSummaryPath):
    """
        Converts standard CRISPR_summary_w_primers_NGS.txt file to
        a Pandas Series containing all of the data that would be included in
        the Converted-Picked.xls file uploaded to SRM2.

        arguments:
            targetSummaryPath: Expects os.path string to the TXT file. Does
                not have to be titled "CRISPR_summary_w_primers_NGS.txt"
    """
    df = pd.read_csv(targetSummaryPath, delimiter="\t")
    chosenDF = df[df.Picked == True]  # Slice to only include 'picked' rows

    # Generate a sliced DF from TXT file containing data relevant to final file
    df_cols_for_SRMTemplate = [0, 1]  # Init list with gRNA, sequence
    df_cols_for_SRMTemplate.extend(range(12, 16))  # fwd_oligo_name ... rev_oligo
    df_cols_for_SRMTemplate.extend(range(19, 23))  # fwd_primer_name ... rev_primer
    df_cols_for_SRMTemplate = [
        "Name",
        "gRNA",
        "fwd_oligo_name",
        "fwd_oligo",
        "rev_oligo_name",
        "rev_oligo",
        "Fwd_Primer_Name",
        "Fwd_Primer",
        "Rev_Primer_Name",
        "Rev_Primer",
    ]
    srmDF = chosenDF.loc[
        :, df_cols_for_SRMTemplate
    ].reset_index()  # New DF with only chosen info relevant to SRM
    srmRow = [
        pd.to_datetime("today").strftime("%m/%d/%Y")
    ]  # Init list w/ created date in string format
    srmRow.append(srmDF["Name"][0].split(".")[0])  # Fetch project number from gRNA name
    srmRow.append(len(srmDF.index))  # Add # of gRNAs
    srmRow = pd.Series(srmRow)  # Convert list to pd Series

    # Sort srmDF alphabetically and reset the index
    srmDF = srmDF.sort_values("Name").reset_index(drop=True)

    # For loop below fills remaining col with data from srmDF to create
    #   row that will be filled in for resulting XLS file
    for index, row in srmDF.iterrows():  # index is unused placeholder. We want 'row'
        srmRow = srmRow.append(
            row.reset_index(drop=True)[1:], ignore_index=True
        )  # Add each gRNA row to the end of the template row
    srmRow.fillna(value="", inplace=True)

    # print('Input TXT File "' + str(targetSummaryPath) + '" converted to Pandas Series.')
    return srmRow


def CRISPRSummarySeriesToXLS(seriesToUse, outputPath, templatePath=DEFAULT_TEMPLATE):
    """
        Writes a Pandas series (seriesToUse) to row 4 (index 3) in the provided
        template XLS file (templatePath), then saves the resulting XLS file
        to the provided outputPath.

        arguments:
            seriesToUse: Pandas series containing one element per cell to be written
                to row 4 (index 3) in template file.
            outputPath: os.path string leading to desired output, including filename
            templatePath: os.path string leading to a template XLS file
    """

    def writeCellAndCopyOldFormatting(newSheet, col, val):
        oldCellRef = newSheet._Worksheet__rows.get(3)._Row__cells.get(
            col
        )  # NOTE: All in row 3. Rework if generalizing
        newSheet.write(3, col, val)  # NOTE: Taking from srmRow. Rework if generalizing
        newCellRef = newSheet._Worksheet__rows.get(3)._Row__cells.get(col)
        if val != "" and newCellRef is not None:
            try:
                newCellRef.xf_idx = oldCellRef.xf_idx
            except:
                None
        # print('wrote cell #' + str(col) + ' in row #' + str(3))

    # print('Series to XLS start')

    with xlrd.open_workbook(filename=templatePath, formatting_info=True) as rb:
        wb = copy(rb)
        sheet = wb.get_sheet(0)
        numCols = rb.sheet_by_index(0).ncols

        for c in range(len(seriesToUse)):
            writeCellAndCopyOldFormatting(sheet, c, seriesToUse[c])

        for c in range(len(seriesToUse), numCols):
            writeCellAndCopyOldFormatting(sheet, c, "")

        wb.save(outputPath)
    # print("Pandas Series converted to XLS File.")
    print("\n  Output Filename: " + str(outputPath))
    print("Joining Completed")
    return wb


def TXT_To_SRM_XLS(inputFile, output):
    # print("Converting TXT to an XLS for SRM2...")
    # print("Template XLS File: " + str(template))
    finalRow = CRISPRSummaryToSeries(inputFile)
    return CRISPRSummarySeriesToXLS(finalRow, output)
    # Function over, output is "wb" from CRISPRSummarySeriesToXLS



def multiTXTtoXLS(
    file_list,
    outputFilename="Converted-Picked-Joined-Py.xls",
    templateToUse=None,
):
    """
        Takes multiple CRISPR_Summary.txt file paths via list as input, joins to form SRM XLS
            file in the same folder as the first supplied CRISPR_Summary.txt file_list[0].
    """
    print(file_list)
    
    print(f"\n\nMulti text list[0]: {file_list[0]}")
    
    
    
    #get names of projects
    #get gNRA names
    #get sequences
    #combine into df
    
    project_names = []
    #iniate seies with first project
    #joined_series = pd.Series
    
    
    joined_series = CRISPRSummaryToSeries(file_list[0])
    project_names.append(joined_series.iloc[1])
    num_of_gRNA = joined_series.iloc[2]
    
    
    
    for file in file_list[1:]:
        series = CRISPRSummaryToSeries(file)
        joined_series = joined_series.append(series.iloc[3:], ignore_index=True)
        project_names.append(series.iloc[1])
        num_of_gRNA = num_of_gRNA + series.iloc[2]
    
    
    joined_series.iloc[1] = ' and '.join(project_names)
    joined_series.iloc[2] = num_of_gRNA
    
    print(joined_series.to_string())
    
    
    outputPath = os.path.join(os.path.dirname(file_list[0]), outputFilename)
    if templateToUse is not None:
        newXLS = CRISPRSummarySeriesToXLS(joined_series, outputPath, templatePath=templateToUse)
    else:
        newXLS = CRISPRSummarySeriesToXLS(joined_series, outputPath)
    return newXLS 




def dualTXTtoXLS(
    summary_text1,
    summary_text2,
    outputFilename="Converted-Picked-Joined-Py.xls",
    templateToUse=None,
):
    """
        Takes two CRISPR_Summary.txt file paths as input, joins to form SRM XLS
            file in the same folder as the first supplied CRISPR_Summary.txt.
    """
    firstSeries = CRISPRSummaryToSeries(summary_text1)
    secondSeries = CRISPRSummaryToSeries(summary_text2)
    joinedSeries = firstSeries.append(secondSeries.iloc[3:], ignore_index=True)
    joinedSeries.iloc[1] = firstSeries.iloc[1] + " and " + secondSeries.iloc[1]
    joinedSeries.iloc[2] = firstSeries.iloc[2] + secondSeries.iloc[2]
    outputPath = os.path.join(os.path.dirname(summary_text1), outputFilename)
    if templateToUse is not None:
        newXLS = CRISPRSummarySeriesToXLS(
            joinedSeries, outputPath, templatePath=templateToUse
        )
    else:
        newXLS = CRISPRSummarySeriesToXLS(joinedSeries, outputPath)
    return newXLS

def SRMJoiner():
    filename = "CRISPR_summary_w_primers_NGS.txt"
    first_folder = ""
    second_folder = ""
    try:
        # Code block below can be used to specify other filename.
        # Commented out to default to default filename on 11/16/2020 - JS

        # use_default_filename_choice = ""
        # while use_default_filename_choice not in ["Y", "YES", "N", "NO"]:
        #     print("Use default filename 'CRISPR_summary_w_primers_NGS.txt'? (y/n)")
        #     use_default_filename_choice = input("--> ").strip().upper()
        # if use_default_filename_choice in ["N", "NO"]:
        #     print("Please enter the intended filename:")
        #     filename = input("--> ")
        # elif use_default_filename_choice in ["Y", "YES"]:
        #     filename = "CRISPR_summary_w_primers_NGS.txt"

        print("  Please enter the full path to the folder of the first file.")
        print("  * NOTE: This will serve as the output folder! *")
        first_folder = pathlib.PureWindowsPath(
            input("    --> ")
            .strip()
            .replace(
                r"Z:\ResearchHome\Groups\millergrp\home",
                r"\research_jude\rgs01_jude\groups\millergrp\home",
            )
            .replace(
                r"Z:\ResearchHome\Groups\millergrp\projects",
                r"\research_jude\rgs01_jude\groups\millergrp\projects",
            )
        ).as_posix()
        print("\n  Please enter the full path to the folder of the second file.")
        second_folder = pathlib.PureWindowsPath(
            input("    --> ")
            .strip()
            .replace(
                r"Z:\ResearchHome\Groups\millergrp\home",
                r"\research_jude\rgs01_jude\groups\millergrp\home",
            )
            .replace(
                r"Z:\ResearchHome\Groups\millergrp\projects",
                r"\research_jude\rgs01_jude\groups\millergrp\projects",
            )
        ).as_posix()

        firstFile = pathlib.PurePosixPath(first_folder).joinpath(filename)
        secondFile = pathlib.PurePosixPath(second_folder).joinpath(filename)
        dualTXTtoXLS(firstFile, secondFile)
        print("\n")
        return None
    except EOFError:
        print("\n")
        return None
    except FileNotFoundError:
        print("\nCould not find the files. Please try again.")
        print("You entered:")
        print(f"  1: {firstFile}")
        print(f"  2: {secondFile}")
        print()
        return None


def repick(current_project_path):
    try:
        filenames = ["CRISPR_summary.txt", "CRISPR_summary_w_primers_NGS.txt"]
        folder = current_project_path

        columns_to_show = [
            "Name",
            "gRNA",
            "long_0",
            "long_1",
            "long_2",
            "short_0",
            "BP",
            "Distance_from_BP",
            "Picked",
        ]

        example_input_file = os.path.join(folder, "CRISPR_summary.txt")
        example_df = pd.read_csv(example_input_file, sep="\t")
        example_df.index += 1

        print(example_df[columns_to_show])
        print(
            "\nPlease enter the index of ALL guides you would like picked, separated by a space."
        )
        print("Note: This will erase current picks.")
        guides = input("  --> ").strip()
        guides = guides.split(" ")
        # print(guides)

        formatted_guides = []
        for choice in guides:
            formatted_guides.append(int(choice) - 1)
        # print(formatted_guides)

        print("\nYou chose the following guides:")
        for choice in formatted_guides:
            print(f"  * {example_df.iloc[choice]['Name']}")

        print("\nConfirm picks? (y/n)")
        confirmation = input("  --> ").strip().lower()
        if confirmation in ["y", "yes"]:
            for filename in filenames:
                print(f"Repicking for {os.path.basename(filename)}...")
                input_fp = os.path.join(folder, filename)
                output_fp = os.path.join(folder, filename)
                backup_fp = os.path.join(folder, f"{filename}.bak")
                input_df = pd.read_csv(input_fp, sep="\t")
                input_df.index += 1
                shutil.copy(input_fp, backup_fp)
                new_df = input_df.copy(deep=True).assign(Picked="False")
                for choice in formatted_guides:
                    new_df.at[choice + 1, "Picked"] = True
                print(new_df[columns_to_show])
                new_df.to_csv(output_fp, sep="\t", index=False)
                print("\nMade .bak of original and saved updated summary.")
                print("\n")

            print(
                "Done repicking!\n  ** It is highly recommended to run primer command for new picks. **"
            )
        return None
    except EOFError:
        print("\n")
        return None
    except FileNotFoundError:
        print("\nCould not find the file. Please try again.")
        print("You entered:")
        print(f"  {current_project_path}")
        print()
        return None
    except ValueError:
        print(
            "You got a value error. I don't know what to do with that, I am sorry. Please try again."
        )
        return None


def singleConverter():
    try:
        print("  Please enter the full path to the project folder.")
        print("  * NOTE: This will serve as the output folder! *\n")
        filename = "CRISPR_summary_w_primers_NGS.txt"
        folder = pathlib.PureWindowsPath(
            input("    --> ")
            .strip()
            .replace(r"Z:\ResearchHome\Groups", r"\research_jude\rgs01_jude\groups")
            .replace(
                r"Z:\ResearchHome\Groups\millergrp\projects",
                r"\research_jude\rgs01_jude\groups\millergrp\projects",
            )
            # .replace("CORE PROJECTS", '"CORE PROJECTS"')
            # .split("\\")
            # + '"'
        ).as_posix()
        input_fp = pathlib.PurePosixPath(folder).joinpath(filename)
        output_fp = pathlib.PurePosixPath(folder).joinpath(
            "Converted-Picked-Manual-Py.xls"
        )
        TXT_To_SRM_XLS(input_fp, output_fp)
        return None
    except EOFError:
        print("\n")
        return None
    except FileNotFoundError:
        print("\nCould not find the file. Please try again.")
        print("You entered:")
        print(f"  {input_fp}")
        print()
        return None


def main():
    # USER VARIABLES
    # CRISPR_summary_w_primers_NGS.txt file path = TXTInput
    # targetPath = os.path.dirname(__file__)
    targetPath = os.path.abspath(
        r"Z:\ResearchHome\Groups\millergrp\home\common\CORE PROJECTS\hCBT-Miller-CAGE733"
    )
    targetFilename = "CRISPR_summary_w_primers_NGS" + ".txt"  # Filename
    TXTInput = os.path.join(targetPath, targetFilename)

    # Converted-Picked.xls template for headers/formatting
    templatePath = os.path.dirname(__file__)
    templateFilename = "Converted-Picked-Template" + ".xls"
    template = os.path.join(templatePath, templateFilename)

    # Output XLS file
    outputPath = targetPath
    outputFilename = "Converted-Picked-Py" + ".xls"
    output = os.path.join(outputPath, outputFilename)

    # Run actual function
    TXT_To_SRM_XLS(TXTInput, output, template)
    return None


if __name__ == "__main__":
    main()
