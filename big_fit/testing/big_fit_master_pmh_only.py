import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os
import shutil
import glob
import warnings
import xlwings as xw
import re

try:
    from natsort import natsort_keygen
    
except:
    os.system('pip install natsort')
    from natsort import natsort_keygen

inputCSV = 'input.csv'
ngs_dir = r'Z:\ResearchHome\Groups\millergrp\home\common\NGS'
script_dir = os.getcwd()
plt.rcParams['savefig.directory'] = (os.path.join(os.environ['USERPROFILE'], 'Desktop'))
holding_dir = os.path.join(os.getcwd(),'holding')
compiled_excel = 'compiled_data.xlsx'

#turns off openpyxl data valiadtion warning
warnings.simplefilter(action='ignore', category=UserWarning)

project_title=''

graph_loop = True

def find_csv():
    #clear holding_dir
    def clear_holding_dir():
        os.chdir(holding_dir)
        
        for file in os.scandir(holding_dir):
            os.remove(file)
        os.chdir(script_dir)
    
    clear_holding_dir()
    #get input from csv
    input_df = pd.read_csv(inputCSV,dtype=object)

    input_df = input_df.astype(str)
    
    #puts columns into list to containerize the data eaiser
    input_projects = input_df.values.tolist()
    
    print('\nScanning for csv\n')

    for target_proj in input_projects:
        
        cage_proj, ngs_date = target_proj
        cage_proj = str(cage_proj).strip()
        
        #appends 0 in from of months Jan->Sept
        if len(ngs_date) < 6:
            ngs_date = str('0' + str(ngs_date))
        
        
        joined_dir = os.path.join(ngs_dir,ngs_date,'joined').replace('\\\\','\\')
        os.chdir(joined_dir)

        #search joined folder for all_indel
        for file in glob.glob(f'{joined_dir}/*/{cage_proj}.csv',recursive=True):
            print(file)
            if cage_proj in file:
                shutil.copy(file, holding_dir)
    
    
    print('\nCSV files copied over.\n')

def compile_to_excel():
    print('\nCompiling data into Excel.\n')
    
    def _format_excel(compiled_df):
                    
            os.chdir(script_dir)
        
            comparison_df = pd.read_csv(inputCSV,dtype=object,usecols=['Project Name'])
            comparison_groups = comparison_df.values.tolist()
            
            #load workbook
            template_excel = 'template.xlsx'
            app = xw.App(visible=False)
            wb = xw.Book(template_excel)
            ws1 = wb.sheets['Sheet1']
            ws2 = wb.sheets['dont_touch']
            #update workbook
            ws1.range('B2').options(index=False,header=False).value = compiled_df

            #close and save as updated_form
            wb.save(compiled_excel)
            wb.close()
            app.quit()
            
    tmp_df_list=[]
    
    #instainate big df will add Sample Name column after comipling all csvs
    compiled_columns = ['CAGE#','Well#','Gene','Guide','0bp','In-frame', 'Out-of-Frame','Comparison_Group','Replicate','Graph_Group']
    compiled_df = pd.DataFrame(columns=[compiled_columns])
    
    os.chdir(holding_dir)
    
    #will return dir entry so can access name via attribute
    csv_list = os.scandir()
    
    #parse data from each csv
    for csv in csv_list:
        
        tmp_columns = ['Name','Sample','0bp','In-frame','Out-of-frame']
        
        #copy only the columns needed and append to list to concat into compiled_df later
        tmp_df = pd.read_csv(csv)
        
        tmp_df = tmp_df[tmp_columns]
        tmp_df.rename(columns={'Name':'Well#','Sample':'Guide'},inplace=True)
        
        #drop rows with no sample/guide info
        tmp_df = tmp_df.dropna(subset=['Guide'])
        
        cage_num = csv.name.split('_')[0]
        gene = csv.name.split('_')[1]
        tmp_df.insert(0,'CAGE#',cage_num,True)
        tmp_df.insert(2,'Gene',gene,True)
        
        tmp_df_list.append(tmp_df)
    
    #compile all csv df's and reset index, save as excel.  Saving excel happens in _format_excel(), template xlsx is overwritten
    compiled_df = pd.concat(tmp_df_list,ignore_index=True)
    compiled_df.index = compiled_df.index +1

    _format_excel(compiled_df)
    pass

def get_scores():
  
    replicates_found = False
    
    print('\nChoose time points and graphing groups in the compiled_data.xlsx.\n')
    #TODO
    input('\n\nPress Enter to continue\n\n')
    
    #Read in compiled excel, drop rows without a comparison selected and reset index for readability
    columns = ['CAGE#', 'Gene', 'Guide','Out-of-frame', 'Comparison_Group','Replicate','Graph_Group']
    score_df = pd.read_excel(compiled_excel, usecols=columns)
    score_df = score_df[score_df['Comparison_Group'].notna()]
    
    score_df.reset_index(drop=True,inplace=True)
    score_df.index = score_df.index + 1
    
    #need to loop through each graph_group and get the fitness scores
    #get unique graph groups
    graph_groups = score_df['Graph_Group'].unique().tolist()
    
    excel_df = pd.DataFrame()
    
    for group in graph_groups:
        #fitness_scores = pd.DataFrame()        
        group_df = score_df[score_df['Graph_Group'] == group].copy()
        try:
            group_df.loc[:,'Guide'] = group_df['Guide'].apply(lambda x: re.search(r'g\d*', x).group())
        except:
            print("************************Please double check guide names have been added to the guide column in the compiled_data.xlsx*************************")
        #Stats versions
        if not group_df['Replicate'].isnull().all():
            print("Running Stats")
            replicates_found = True
            
            #generate fitness scores
            #break out each initial and final time point into seperate df's then merge to have 1 flat combined df.
            #break out each replicate into seperate df and concat later
            init_df = group_df[['CAGE#','Gene','Guide','Out-of-frame','Comparison_Group','Replicate']][group_df['Comparison_Group'].str.contains('initial')]
            init_df.rename(columns={'Out-of-frame':'init_oof'},inplace=True) #need to have unique column for merge
            
            final_df = group_df[['CAGE#','Gene','Guide','Out-of-frame','Comparison_Group','Replicate']][group_df['Comparison_Group'].str.contains('final')]
            final_df.rename(columns={'Out-of-frame':'final_oof'},inplace=True)
            
            rep1_df = group_df[['CAGE#','Gene','Guide','Replicate','Out-of-frame','Comparison_Group']][group_df['Replicate'].str.contains('Rep1')]
            rep2_df = group_df[['CAGE#','Gene','Guide','Replicate','Out-of-frame','Comparison_Group']][group_df['Replicate'].str.contains('Rep2')]
            rep3_df = group_df[['CAGE#','Gene','Guide','Replicate','Out-of-frame','Comparison_Group']][group_df['Replicate'].str.contains('Rep3')]

            df_list = [rep1_df, rep2_df, rep3_df]

            #add init/final oof columns
            #rep df layout: CAGE# GENE Replicate Init_oof Final_oof
            for df in df_list:
                df.insert(4,'init_oof',True)
                df.insert(5,'final_oof',True)
                df['init_oof'] = df['Out-of-frame'][df['Comparison_Group'].str.contains('initial')]
                df['final_oof'] = df['Out-of-frame'][df['Comparison_Group'].str.contains('final')]
                df.drop(columns=['Out-of-frame','Comparison_Group'],inplace=True)
                df = df.groupby('CAGE#').first().reset_index()

            #combine rep dfs and clean up formatting
            combo_df = pd.concat(df_list).reset_index(drop=True)

            combo_df['Replicate'] = combo_df['Replicate'].apply(lambda x: re.search(r'Rep\d*', x).group())
            combo_df = combo_df.groupby(['CAGE#','Gene','Guide','Replicate']).agg({'init_oof':'first', 'final_oof':'last'}).reset_index()

            combo_df['fitness_score'] = (combo_df['final_oof'] / combo_df['init_oof']).round(2)
            combo_df.sort_values(by=['CAGE#','Gene','Guide','Replicate'],inplace=True)
            

            #calculate average and standard dev
            #by applying agg to the fitness_score column I get two separate columns and a single index.  if fitness_score is a agg() parameter I would get a multi_indexed df
            graphing_results_df = combo_df.groupby(['CAGE#', 'Gene', 'Guide'])['fitness_score'].agg(['mean', 'std']).reset_index() #*will give mean and std
            graphing_results_df.sort_values(by=['mean'],inplace=True)
            graphing_results_df.rename(columns={'mean':'fitness_score'},inplace=True)
            graphing_results_df.insert(3,'Comparison_Group',True)
            cat_cols = ['CAGE#','Guide']
            graphing_results_df['Comparison_Group'] = graphing_results_df[cat_cols].apply(lambda x: '_'.join(x.values.astype(str)),axis=1)
            graphing_results_df.drop_duplicates(inplace=True)

            avg_df= combo_df.groupby('Guide')['fitness_score'].mean().reset_index()
            avg_df.rename(columns={'fitness_score':'avg_fit_score'},inplace=True)
            
            stdev_df= combo_df.groupby('Guide')['fitness_score'].std().reset_index()
            stdev_df.rename(columns={'fitness_score':'stdev'},inplace=True)
            
            scores_df = combo_df.merge(avg_df, on='Guide').merge(stdev_df, on='Guide')
            scores_df.reset_index(inplace=True, drop=True)
            scores_df.index = scores_df.index + 1
            
            excel_columns = ['CAGE#','Gene','Guide','Replicate','init_oof','final_oof','fitness_score','avg_fit_score','stdev']
            scores_df = scores_df[excel_columns]

        #Non stats version
        else:
            print("Non stats version")
            replicates_found = False
            try:
                group_df['Guide'].apply(lambda x: re.search(r'g\d*', x).group())

            except:
                print("************************Please double check guide names have been added to the guide column in the compiled_data.xlsx*************************")
            #generate fitness scores
            #break out each initial and final time point into seperate df's then merge to have 1 flat combined df.#
            init_df = group_df[['CAGE#','Gene','Guide','Out-of-frame']][group_df['Comparison_Group'].str.contains('initial')]
            init_df['Comparison_Group'] = group_df['Comparison_Group'][group_df['Comparison_Group'].str.contains('initial')].str.replace('.initial','')
            #remove any days in listed in the comparison group so the initial and final df's can be merged
            init_df['Comparison_Group'] = init_df['Comparison_Group'].str.replace(r'\.(d\d*)\.?', lambda x: '.',regex=True, case=False)
            
            init_df.rename(columns={'Out-of-frame':'init_oof'},inplace=True) #need to have unique column for merge


            final_df = group_df[['CAGE#','Gene','Guide','Out-of-frame']][group_df['Comparison_Group'].str.contains('final')]
            final_df['Comparison_Group'] = group_df['Comparison_Group'][group_df['Comparison_Group'].str.contains('final')].str.replace('.final','')
            final_df['Comparison_Group'] = final_df['Comparison_Group'].str.replace(r'\.(d\d*)\.?', lambda x:'.', regex=True, case=False)
            
            final_df.rename(columns={'Out-of-frame':'final_oof'},inplace=True)


            graphing_results_df = init_df.merge(final_df,on=['CAGE#','Gene','Guide','Comparison_Group'])
            
            graphing_results_df['Guide'] = graphing_results_df['Guide'].str.strip()
            graphing_results_df['Comparison_Group'] = graphing_results_df['Comparison_Group'].str.strip().str.rstrip('.')

            graphing_results_df.insert(5,'fitness_score',(graphing_results_df['final_oof'] / graphing_results_df['init_oof']).round(2))
            graphing_results_df.drop(columns=['init_oof','final_oof'],inplace=True)
            
            graphing_results_df.index = graphing_results_df.index + 1
            
            
            print(f"graphing_results_df\n{graphing_results_df}\n\nPress Enter to continue\n\n")
            
            scores_df = graphing_results_df.copy()
            #scores_df.drop(columns=['Comparison_Group'],inplace=True)

        graph_scores(graphing_results_df,replicates_found)
        
        excel_df = pd.concat([excel_df,scores_df])

    print(excel_df )
    excel_df.to_excel('fitness_scores.xlsx',index=False)
        
def graph_scores(graphing_results_df, replicates_found):
        
    def _rotate_labels(plot):
        labels = []
        
        for label in plot.get_xticklabels():            
            text = label.get_text()
            guide = text.split('.')
            
            try:            
                wrapped_name = guide[0]+'\n'+guide[1] + '.' + guide[2] + '\n' + guide[3]
            except:
                wrapped_name = guide[0]+'\n'+guide[1] + '.' + guide[2]
            labels.append(wrapped_name)            
            fa_plot.xaxis.set_tick_params(which='major', pad=2)
        plot.xaxis.set_tick_params(which='both', pad=0)
        plot.set_xticklabels(labels, rotation=45, size=8.0, ha='right', rotation_mode='anchor')
    
    def _find_max_y(fa_score_df,replicates_found):
        
        if replicates_found == True:
            y_max = fa_score_df['fitness_score'].max() + fa_score_df['std'].max()
        else:
            y_max = fa_score_df['fitness_score'].max()
            
        if y_max > 1:
        
            y_buffer = 1.05
            
            y_upper_bound = y_buffer * y_max
            
        else:
            y_upper_bound = 1
        
        return y_upper_bound
    
    
    #TODO rework with regex so order inside of comparison group doesnt matter ie cell.guide.date vs guide.date.cell
    #CAGE#.GENE.guide cell type
    
    
    #if cell type is detecte
    if graphing_results_df['Comparison_Group'].str.split('.').str.len().max() == 3:
        
        graphing_results_df['guide_num'] = graphing_results_df['Comparison_Group'].str.extract(r'\.(g\d*)\.?').astype(str)
        print("got guide")
        print(graphing_results_df)
        tmp_df = pd.DataFrame(columns=['cell_type'])
        
        #removes CAGE# and g# no matter where it is in the comparison_group string
        tmp_df['cell_type'] = graphing_results_df['Comparison_Group'].str.replace(r'CAGE\d*\.', lambda x:'', regex=True, case=False)
        tmp_df['cell_type'] = tmp_df['cell_type'].str.replace(r'g\d', lambda x:'', regex=True, case=False)
        graphing_results_df['cell_type'] = tmp_df['cell_type'].str.replace(".","")
        
                    
        graphing_results_df['Guide'] = graphing_results_df['Comparison_Group'].str.split('.').str[0] + '.' + graphing_results_df['Gene'] + "." + graphing_results_df['guide_num'] +"."+ graphing_results_df['cell_type']
        
    else:
        
        graphing_results_df['guide_num'] = graphing_results_df['Comparison_Group'].str.extract(r'\.(g\d*)\.?').astype(str)

        graphing_results_df['Guide'] = graphing_results_df['Comparison_Group'].str.split('.').str[0] + '.' + graphing_results_df['Gene'] + '.' + graphing_results_df['guide_num']


    sort_choice =input("Choose graph sorting method (gene, guide, cage#, fitness, cell type): ").lower().strip().replace(" ","")
            
    while sort_choice not in ['gene','guide','cage#','fitness','celltype']:
        sort_choice = input("Invalid choice, please choose gene, guide, cage# or cell type: ")
    
    if sort_choice == 'gene':
        graphing_results_df.sort_values(by=['Gene'],inplace=True,key=natsort_keygen())
        
    elif sort_choice == 'guide':
        graphing_results_df.sort_values(by=['Guide'],inplace=True,key=natsort_keygen())
        
    elif sort_choice == 'cage#':
        graphing_results_df.sort_values(by=['CAGE#'],inplace=True,key=natsort_keygen())
        
    elif sort_choice == 'fitness':
        graphing_results_df.sort_values(by=['fitness_score'],inplace=True)
    elif sort_choice == 'celltype':
        graphing_results_df.sort_values(by=['Guide'],key = lambda x: x.str.split(" ").str[1], inplace=True)

    
        print("\n*******Graphing the following results*************\n")
        print(graphing_results_df)
        
    try:
        fa_score_df = graphing_results_df[['Guide','fitness_score','std']]
    except:
        fa_score_df = graphing_results_df[['Guide','fitness_score']]

    
    y_limit= _find_max_y(fa_score_df,replicates_found)
    
    bar_num = len(fa_score_df['Guide'].to_list())
    
    #*Stats Version
    if replicates_found == True:
        
        err_list = fa_score_df['std'].tolist()
    
        #tries to standardize the bar width.  Bar width is relative to the plot area so the more bars the smaller the width.
        #Fewer bars look weirdly wide.  Set width to 0.3 if 4 or fewer bars looks better, more than five increase
        #width so they don't look too thin
        if bar_num < 4:
            bar_width = 0.3
        else:
            bar_width = 0.5
        fa_plot = fa_score_df.plot.bar(
                        x='Guide',
                        y='fitness_score',
                        yerr='std',
                        capsize=10,
                        rot=0,
                        legend=False,
                        ylim=(0, y_limit), 
                        color='#008ccf',
                        align = 'center',
                        label = 'yes',
                        #sets bars in front of grid lines, effect doesn't work unless at least 3
                        zorder = 3, 
                        width = bar_width
                    )    
        
        #writes labels on top of bars
        
        for patch in enumerate(fa_plot.patches):
            
            #patch returns a tuple (bar#,rectangle)
            #hacky way to make this work without re-writting a bunch of stuff

            #sets rectangle to bar
            bar = patch[1]
            #gets appropriate error index from err_list
            err_buffer = float(err_list[patch[0]])
            
            
            fa_plot.text(
                #align middle
                bar.get_x() + bar.get_width() / 2,
                #set label slightly above bar
                bar.get_height() + 0.03 + err_buffer,
                #label string, and keep 2 decimals
                '{:.2f}'.format(bar.get_height()),
                #horizontal alignement
                ha='center',
                color='black',
                weight='bold',
                size=12
            )
            
    #*Non-stats version
    elif replicates_found == False:
                
        #tries to standardize the bar width.  Bar width is relative to the plot area so the more bars the smaller the width.
        #Fewer bars look weirdly wide.  Set width to 0.3 if 4 or fewer bars looks better, more than five increase
        #width so they don't look too thin
        
        
        if bar_num < 4:
            bar_width = 0.3
        else:
            bar_width = 0.5
        
        fa_plot = fa_score_df.plot.bar(
                        x='Guide',
                        y='fitness_score',
                        rot=0,
                        legend=False,
                        ylim=(0, y_limit), 
                        color='#008ccf',
                        align = 'center',
                        label = 'yes',
                        #sets bars in front of grid lines, effect doesn't work unless at least 3
                        zorder = 3, 
                        width = bar_width
                    )
        
        #writes labels on top of bars
        for bar in fa_plot.patches:
            fa_plot.text(
                #align middle
                bar.get_x() + bar.get_width() / 2,
                #set label slightly above bar
                bar.get_height() + 0.02,
                #label string, and keep 2 decimals
                '{:.2f}'.format(bar.get_height()),
                #horizontal alignement
                ha='center',
                color='black',
                weight='bold',
                size=12
            )

    #turns of the top and right border of the graphing area
    fa_plot.spines['top'].set_visible(False)
    fa_plot.spines['right'].set_visible(False)
    
    #draw grid and set behind the bars
    fa_plot.grid(axis='y', which='both', zorder=3)
    fa_plot.yaxis.set_minor_locator(AutoMinorLocator(2))
    #fa_plot.minorticks_on()
    #_wrap_labels(fa_plot,width=5)
    plt.xticks(fontsize=12)
    if bar_num > 3:
        _rotate_labels(fa_plot)


    #more y axis minor ticks
    plt.xlabel('Guide',fontsize=15, weight='bold',labelpad=10)
    plt.ylabel('Fitness Score',fontsize=15, weight='bold',labelpad=10)
    plt.title(f'Fitness Scores {project_title}', fontsize=20, weight='bold', pad=10)

    plt.show()
    
def main(graph_loop):

    find_csv() #comment this line out if you have already pulled your data and entered it into the compiled_data.xlsx
    compile_to_excel() #comment this line out if you have already pulled your data and entered it into the compiled_data.xlsx
    
    while graph_loop == True:
        get_scores()
        graph_loop = True if input('\n\n\nRun graphing step again? (y/n): ').lower().strip() == 'y' else False
    else:
        quit()
if __name__ == "__main__":
    main(graph_loop)
