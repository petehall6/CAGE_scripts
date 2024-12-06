import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os
import shutil
import glob
import warnings
import xlwings as xw
import re

inputCSV = 'input.csv'
ngs_dir = r'Z:\ResearchHome\Groups\millergrp\home\common\NGS'
script_dir = os.getcwd()
plt.rcParams['savefig.directory'] = (os.path.join(os.environ['USERPROFILE'], 'Desktop'))
holding_dir = os.path.join(os.getcwd(),'holding')
compiled_excel = 'compiled_data.xlsx'

#turns off openpyxl data valiadtion warning
warnings.simplefilter(action='ignore', category=UserWarning)

project_title=''


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
        
        cage_proj, ngs_date, comparison_group = target_proj
        cage_proj = str(cage_proj).strip()
        #appends 0 in from of months Jan->Sept
        if len(ngs_date) < 6:
            ngs_date = str('0' + str(ngs_date))
        
        
        joined_dir = os.path.join(ngs_dir,ngs_date,'joined').replace('\\\\','\\')
        os.chdir(joined_dir)

        #search joined folder for all_indel
        for file in glob.glob(f'{joined_dir}/{cage_proj}*/*all_indel*.csv',recursive=True):
            #print(file)
            if cage_proj in file:
                shutil.copy(file, holding_dir)
    
    print('\nCSV files copied over.\n')

def compile_to_excel():
    print('\nCompiling data into Excel.\n')
    
    def _format_excel(compiled_df):
                    
            os.chdir(script_dir)
        
            #TODO only read the input file once and stop being lazy
            comparison_df = pd.read_csv(inputCSV,dtype=object,usecols=['Project Name','Comparison Groups'])
            comparison_groups = comparison_df.values.tolist()
            
            #load workbook
            template_excel = 'template.xlsx'
            app = xw.App(visible=False)
            wb = xw.Book(template_excel)
            ws1 = wb.sheets['Sheet1']
            ws2 = wb.sheets['dont_touch']
            #update workbook
            ws1.range('B2').options(index=False,header=False).value = compiled_df
            
            #*Transfer comparison groups to excel
            #get unique CAGE#s
            
            time_point_list=[]
            rep_list=[]
            
            #cage_nums = sorted(set(compiled_df['CAGE#'].tolist()))
            #TODO comparisons to time point list
            
            for comparison in comparison_groups:
                
                grp = comparison[0] + '_' + comparison[1]
                
                time_point_list.append(f'{grp}_init_tp')
                time_point_list.append(f'{grp}_final_tp')
                
                for rep in range(1,4):
                    rep_list.append(f'{grp}_rep{rep}')

                
            time_point_df = pd.DataFrame(time_point_list)
            rep_df = pd.DataFrame(rep_list)
            
            ws2.range('A2').options(index=False,header=False).value = time_point_df
            ws2.range('C2').options(index=False,header=False).value = rep_df
            
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
        
        #print(csv.name)
        
        
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
    #print(compiled_df.head(20))
    
    _format_excel(compiled_df)
    pass

def get_scores():
    
    replicates_found = False
    
    print('\nChoose comparison groups in the compiled_data.xlsx.\n')
    #TODO uncomment
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
    
    for group in graph_groups:
        fitness_scores = pd.DataFrame()        
        group_df = score_df[score_df['Graph_Group'] == group].copy()
        group_df.loc[:,'Guide'] = group_df['Guide'].apply(lambda x: re.search(r'g\d*', x).group())
        
        #Stats versions
        if not group_df['Replicate'].isnull().all():
            print("Running Stats")
            replicates_found = True
            
            
            #generate fitness scores
            #break out each initial and final time point into seperate df's then merge to have 1 flat combined df.
            #break out each replicate into seperate df and concat later
            
            
            init_df = group_df[['CAGE#','Gene','Guide','Out-of-frame','Comparison_Group','Replicate']][group_df['Comparison_Group'].str.contains('init')]
            init_df.rename(columns={'Out-of-frame':'init_oof'},inplace=True) #need to have unique column for merge
            
            final_df = group_df[['CAGE#','Gene','Guide','Out-of-frame','Comparison_Group','Replicate']][group_df['Comparison_Group'].str.contains('final')]
            final_df.rename(columns={'Out-of-frame':'final_oof'},inplace=True)
            
            rep1_df = group_df[['CAGE#','Gene','Guide','Replicate','Out-of-frame','Comparison_Group']][group_df['Replicate'].str.contains('rep1')]
            rep2_df = group_df[['CAGE#','Gene','Guide','Replicate','Out-of-frame','Comparison_Group']][group_df['Replicate'].str.contains('rep2')]
            rep3_df = group_df[['CAGE#','Gene','Guide','Replicate','Out-of-frame','Comparison_Group']][group_df['Replicate'].str.contains('rep3')]

            df_list = [rep1_df, rep2_df, rep3_df]
        
            #add init/final oof columns
            #rep df layout: CAGE# GENE Replicate Init_oof Final_oof
            
            for df in df_list:
                df.insert(4,'init_oof',True)
                df.insert(5,'final_oof',True)
                df['init_oof'] = df['Out-of-frame'][df['Comparison_Group'].str.contains('init')]
                df['final_oof'] = df['Out-of-frame'][df['Comparison_Group'].str.contains('final')]
                df.drop(columns=['Out-of-frame','Comparison_Group'],inplace=True)

            #have to implicitly do groupby.  df_list isnt updating with new dfs for some reason
            rep1_df = rep1_df.groupby('CAGE#').first().reset_index()
            rep2_df = rep2_df.groupby('CAGE#').first().reset_index()
            rep3_df = rep3_df.groupby('CAGE#').first().reset_index()
            
            rep_df_list = [rep1_df, rep2_df, rep3_df]

            #combine rep dfs and clean up formatting
            combo_df = pd.concat(rep_df_list).reset_index(drop=True)
            
            combo_df['fitness_score'] = (combo_df['final_oof'] / combo_df['init_oof']).round(2)
            combo_df.sort_values(by=['CAGE#','Gene','Guide','Replicate'],inplace=True)
            
            #calculate average and standard dev
            #by applying agg to the fitness_score column I get two separate columns and a single index.  if fitness_score is a agg() parameter I would get a multi_indexed df
            graphing_results_df = combo_df.groupby(['CAGE#', 'Gene', 'Guide'])['fitness_score'].agg(['mean', 'std']).reset_index() #*will give mean and std
            graphing_results_df.sort_values(by=['mean'],inplace=True)
            graphing_results_df.rename(columns={'mean':'fitness_score'},inplace=True)
            
            avg_df= combo_df.groupby('Guide')['fitness_score'].mean().reset_index()
            avg_df.rename(columns={'fitness_score':'avg_fit_score'},inplace=True)
            
            stdev_df= combo_df.groupby('Guide')['fitness_score'].std().reset_index()
            stdev_df.rename(columns={'fitness_score':'stdev'},inplace=True)
            
            excel_df = combo_df.merge(avg_df, on='Guide').merge(stdev_df, on='Guide')
            excel_df.reset_index(inplace=True, drop=True)
            excel_df.index = excel_df.index + 1
            
            excel_columns = ['CAGE#','Gene','Guide','Replicate','init_oof','final_oof','fitness_score','avg_fit_score','stdev']
            excel_df = excel_df[excel_columns]
            
            #fitness_scores.append(excel_df)
        
        #Non stats version
        else:
            print("Non stats version")
            replicates_found = False
            
            #generate fitness scores
            #break out each initial and final time point into seperate df's then merge to have 1 flat combined df.
            init_df = group_df[['CAGE#','Gene','Guide','Out-of-frame']][group_df['Comparison_Group'].str.contains('init')]
            init_df.rename(columns={'Out-of-frame':'init_oof'},inplace=True) #need to have unique column for merge

            final_df = group_df[['CAGE#','Gene','Guide','Out-of-frame']][group_df['Comparison_Group'].str.contains('final')]
            final_df.rename(columns={'Out-of-frame':'final_oof'},inplace=True)


            graphing_results_df = init_df.merge(final_df,on=['CAGE#','Gene','Guide'])
            graphing_results_df.insert(4,'fitness_score',(graphing_results_df['final_oof'] / graphing_results_df['init_oof']).round(2))
            graphing_results_df.drop(columns=['init_oof','final_oof'],inplace=True)
            
            graphing_results_df.index = graphing_results_df.index + 1
            
            graphing_results_df.sort_values(by=['fitness_score'],inplace=True)
            
            graphing_results_df.to_excel('fitness_scores.xlsx')
        
            comparison_df= score_df['Comparison_Group']
            
            fitness_scores = pd.concat([fitness_scores,graphing_results_df,comparison_df],axis=1)
        
        fitness_scores.to_excel('fitness_scores.xlsx',index=False)
        graph_scores(graphing_results_df,replicates_found)

        

def graph_scores(graphing_results_df,replicates_found):
        
    def _rotate_labels(plot):
        labels = []
        
        for label in plot.get_xticklabels():            
            text = label.get_text()
            Guide = text.split('.')
            wrapped_name = Guide[0]+'\n'+Guide[1]+"."+Guide[2]
            labels.append(wrapped_name)            
            fa_plot.xaxis.set_tick_params(which='major', pad=2)
        plot.xaxis.set_tick_params(which='both', pad=0)
        plot.set_xticklabels(labels, rotation=45, size=8.0, ha='right', rotation_mode='anchor')
    
    def _find_max_y(fa_score_df):
        
        y_max = fa_score_df['fitness_score'].max()
        
        if y_max > 1:
        
            y_buffer = 1.05
            
            y_upper_bound = y_buffer * y_max
            
        else:
            y_upper_bound = 1
        
        return y_upper_bound
    print("\n*******Graphing the following results*************\n")
    print(f"{graphing_results_df}\n")
    
    graphing_results_df['Guide'] = graphing_results_df[['CAGE#','Gene','Guide']].agg('.'.join, axis=1)

    try:
        fa_score_df = graphing_results_df[['Guide','fitness_score','std']]
    except:
        fa_score_df = graphing_results_df[['Guide','fitness_score']]
    
    y_limit= _find_max_y(fa_score_df)
    
    bar_num = len(fa_score_df['Guide'].to_list())
    
    
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
    _rotate_labels(fa_plot)

    #more y axis minor ticks
    plt.xlabel('Guide',fontsize=15, weight='bold',labelpad=10)
    plt.ylabel('Fitness Score',fontsize=15, weight='bold',labelpad=10)
    plt.title(f'Fitness Scores {project_title}', fontsize=20, weight='bold', pad=10)

    plt.show()
    
def main():

    #find_csv() #comment this line out if you have already pulled your data and entered it into the compiled_data.xlsx
    #compile_to_excel() #comment this line out if you have already pulled your data and entered it into the compiled_data.xlsx
    get_scores()

if __name__ == "__main__":
    main()

print('Graphing Complete.')