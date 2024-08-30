import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os
import shutil
import glob
import warnings
import xlwings as xw
import textwrap

inputCSV = 'input.csv'
ngs_dir = r'Z:\ResearchHome\Groups\millergrp\home\common\NGS'
script_dir = os.getcwd()
plt.rcParams['savefig.directory'] = (os.path.join(os.environ['USERPROFILE'], 'Desktop'))
holding_dir = os.path.join(os.getcwd(),'holding')
compiled_excel = 'compiled_data_stats.xlsx'

#turns off openpyxl data valiadtion warning
warnings.simplefilter(action='ignore', category=UserWarning)

project_title=''

def main():
    #find_csv()
    #compile_to_excel()
    get_scores()

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
        for file in glob.glob(f'{joined_dir}/{cage_proj}*/*all_indel*.csv',recursive=True):
            #print(file)
            if cage_proj in file:
                shutil.copy(file, holding_dir)
    
    print('\nCSV files copied over.\n')

def compile_to_excel():
    print('\nCompiling data into Excel.\n')
    
    def _format_excel(compiled_df):
                    
            os.chdir(script_dir)
            
            #compiled_df.reset_index(drop=True, inplace=True)
            
            #*Transfer compiled df to excel
            
            #load workbook
            template_excel = 'template.xlsx'
            app = xw.App(visible=False)
            wb = xw.Book(template_excel)
            display_sheet = wb.sheets['Sheet1']
            functions_sheet = wb.sheets['dont_touch']
            time_plot_sheet = wb.sheets['time_plot']
            #update workbook
            display_sheet.range('B2').options(index=False,header=False).value = compiled_df
            time_plot_sheet.range('B2').options(index=False,header=False).value = compiled_df
            #*Transfer comparison groups to excel
            #get unique CAGE#s
            
            time_point_list=[]
            rep_list=[]
            
            cage_nums = sorted(set(compiled_df['CAGE#'].tolist()))
            
            for num in cage_nums:
                time_point_list.append(f'{num}_init_tp')
                time_point_list.append(f'{num}_final_tp')
                
                for rep in range(1,4):
                    rep_list.append(f'{num}_rep{rep}')

                
            time_point_df = pd.DataFrame(time_point_list)
            rep_df = pd.DataFrame(rep_list)
            
            functions_sheet.range('A2').options(index=False,header=False).value = time_point_df
            functions_sheet.range('C2').options(index=False,header=False).value = rep_df
            
            #close and save as updated_form
            wb.save(compiled_excel)
            wb.close()
            app.quit()
            
    
    tmp_df_list=[]
    
    #instainate big df will add Sample Name column after comipling all csvs
    compiled_columns = ['CAGE#','Well#','Gene','Guide_Name','0bp','In-frame', 'Out-of-Frame','Comparison_Group','Replicate']
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
        tmp_df.rename(columns={'Name':'Well#','Sample':'Guide_Name'},inplace=True)
        
        cage_num = csv.name.split('_')[0]
        gene = csv.name.split('_')[1]
        tmp_df.insert(1,'CAGE#',cage_num,True)
        tmp_df.insert(2,'Gene',gene,True)
        
        tmp_df_list.append(tmp_df)
    
    #compile all csv df's and reset index, save as excel.  Saving excel happens in _format_excel(), template xlsx is overwritten
    compiled_df = pd.concat(tmp_df_list,ignore_index=True)
    compiled_df.index = compiled_df.index +1
    #print(compiled_df.head(20))
    
    _format_excel(compiled_df)

def get_scores():
    
    def _fitness_scores(score_df):
                #Stats versions
        if pd.isna(score_df.loc[1,'Replicate']) == False:
            print("Running Stats")
            replicates_found = True
            
            
            #generate fitness scores
            #break out each initial and final time point into seperate df's then merge to have 1 flat combined df.
            #break out each replicate into seperate df and concat later
            
            
            init_df = score_df[['CAGE#','Gene','Guide_Name','Out-of-frame','Comparison_Group','Replicate']][score_df['Comparison_Group'].str.contains('init')]
            init_df.rename(columns={'Out-of-frame':'init_oof'},inplace=True) #need to have unique column for merge
            
            final_df = score_df[['CAGE#','Gene','Guide_Name','Out-of-frame','Comparison_Group','Replicate']][score_df['Comparison_Group'].str.contains('final')]
            final_df.rename(columns={'Out-of-frame':'final_oof'},inplace=True)
            
            rep1_df = score_df[['CAGE#','Gene','Guide_Name','Replicate','Out-of-frame','Comparison_Group']][score_df['Replicate'].str.contains('rep1')]
            rep2_df = score_df[['CAGE#','Gene','Guide_Name','Replicate','Out-of-frame','Comparison_Group']][score_df['Replicate'].str.contains('rep2')]
            rep3_df = score_df[['CAGE#','Gene','Guide_Name','Replicate','Out-of-frame','Comparison_Group']][score_df['Replicate'].str.contains('rep3')]
            
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

            
            combo_df = pd.concat(rep_df_list).reset_index()
            combo_df = combo_df.drop(columns=['index'])
            
            combo_df.insert(5,'fitness_score',(combo_df['final_oof'] / combo_df['init_oof']).round(2))
        
            #calculate average and standard dev
            #by applying agg to the fitness_score column I get two separate columns and a single index.  if fitness_score is a agg() parameter I would get a multi_indexed df
            graphing_results_df = combo_df.groupby(['CAGE#','Gene','Guide_Name'])['fitness_score'].agg(['mean','std']).reset_index() #*will give mean and std...good enough to send to plt
            graphing_results_df.sort_values(by=['mean'],inplace=True)
            graphing_results_df.rename(columns={'mean':'fitness_score'},inplace=True)
            
            
            avg_df= combo_df.groupby('Guide_Name')['fitness_score'].mean().reset_index()
            avg_df.rename(columns={'fitness_score':'avg_fit_score'},inplace=True)
            #avg_df.drop(columns=['Guide_Name'],inplace=True)
            
            stdev_df= combo_df.groupby('Guide_Name')['fitness_score'].std().reset_index()
            stdev_df.rename(columns={'fitness_score':'stdev'},inplace=True)
            #stdev_df.drop(columns=['Guide_Name'],inplace=True)
            
            excel_df = pd.concat([combo_df,avg_df,stdev_df],axis=1,ignore_index=False)
            #results.reset_index(inplace=True, drop=True)
            
            #print(excel_df.reset_index())
            
            excel_df.sort_values(by=['CAGE#'],inplace=True)

            print(excel_df.head(20))    
            excel_df.to_excel('fitness_scores_stats.xlsx',index=False)
        
        #Non stats version
        else:
            print("Non stats version")
            replicates_found = False
            
            #generate fitness scores
            #break out each initial and final time point into seperate df's then merge to have 1 flat combined df.
            init_df = score_df[['CAGE#','Gene','Guide_Name','Out-of-frame']][score_df['Comparison_Group'].str.contains('init')]
            init_df.rename(columns={'Out-of-frame':'init_oof'},inplace=True) #need to have unique column for merge

            final_df = score_df[['CAGE#','Gene','Guide_Name','Out-of-frame']][score_df['Comparison_Group'].str.contains('final')]
            final_df.rename(columns={'Out-of-frame':'final_oof'},inplace=True)

            graphing_results_df = init_df.merge(final_df,on=['CAGE#','Gene','Guide_Name'])
            
            graphing_results_df.insert(4,'fitness_score',(graphing_results_df['final_oof'] / graphing_results_df['init_oof']).round(2))
            graphing_results_df.drop(columns=['init_oof','final_oof'],inplace=True)
            
            graphing_results_df.index = graphing_results_df.index + 1
            
            graphing_results_df.sort_values(by=['fitness_score'],inplace=True)
            
            graphing_results_df.to_excel('fitness_scores.xlsx')

        #fitness_graph_scores(graphing_results_df,replicates_found)
    
    def _fold_change(fold_change_df):
        
        def __format_df(input_df):
            
            print(input_df)
            

            formated_df = input_df.pivot_table(index='Guide_Name',
                                               columns='Day',
                                               values='Out-of-frame',
                                               aggfunc='first')
            
            
            formated_df.reset_index(inplace=True)
            formated_df.rename_axis(None, axis=1, inplace=True)
            format_columns = ['Guide_Name','Day 3','Day 7','Day 14','Day 21']
            formated_df = formated_df[format_columns]
            #input(formated_df)
            
            normalized_df = formated_df
            
            normalized_df['Day 21'] = formated_df['Day 21'] / formated_df['Day 3']
            normalized_df['Day 14'] = formated_df['Day 14'] / formated_df['Day 3']
            normalized_df['Day 7'] = formated_df['Day 7'] / formated_df['Day 3']
            normalized_df['Day 3'] = formated_df['Day 3'] / formated_df['Day 3']
            
            
            
        
            normalized_df['Average'] = normalized_df[['Day 3','Day 7','Day 14','Day 21']].mean(axis=1)
            #formated_df.columns = columns
            #formated_df.reset_index(inplace=True)
            
            #formated_df.plot.line(x='Guide_Name',y=['Day 3', 'Day 7','Day 14','Day 21'])
            
            
            input(normalized_df)
            
            #input(formated_df)            
            
            plt.show()
            
            
            return None
        
        #print(fold_change_df)
        
        oof_fold_tmp_df = fold_change_df[['Guide_Name','Day','Out-of-frame']]
        
        #print(oof_fold_tmp_df)
        
        #in_fold_tmp_df = fold_change_df['Guide_Name','Day','In-frame']
        
        #total_fold_tmp_df = fold_change_df['Guide_Name','Day','Total-Indel']
        
        __format_df(oof_fold_tmp_df)
        
        
        
        
        
        
        #fold_change_graph(oof_fold_df, in_fold_df, total_fold_df)
        
        return None
    
    replicates_found = False
    
    print('\nChoose comparison groups in the compiled_data.xlsx.\n')
    #TODO uncomment
    #input('\n\nPress Enter to continue\n\n')
    
    #Read in compiled excel, drop rows without a comparison selected and reset index for readability
    score_columns = ['CAGE#', 'Gene', 'Guide_Name','Out-of-frame', 'Comparison_Group','Replicate']
    score_df = pd.read_excel(compiled_excel, usecols=score_columns, sheet_name='Sheet1')
    score_df = score_df.dropna(axis=0, how='all')

    score_df.reset_index(drop=True,inplace=True)
    score_df.index = score_df.index + 1
    
    fold_change_columns = ['Gene','Guide_Name','Day','In-frame','Out-of-frame','Total-Indel']
    fold_change_df = pd.read_excel(compiled_excel, usecols=fold_change_columns, sheet_name='time_plot')
    fold_change_df.dropna(inplace=True)
    
    #print(fold_change_df)
    
    #_fitness_scores(score_df)
    
    _fold_change(fold_change_df)

    
    
def fold_change_graph(oof_fold_df,in_fold_df,total_fold_df):
    
    return None
    
def fitness_graph_scores(graphing_results_df,replicates_found):
        
    def _rotate_labels(plot):
        labels = []
        
        for label in plot.get_xticklabels():            
            text = label.get_text()
            guide_name = text.split('.')
            wrapped_name = guide_name[0]+'\n'+guide_name[1]+"."+guide_name[2]
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
    
    print(graphing_results_df.head(20))
    
    graphing_results_df['Guide_Name'] = graphing_results_df[['CAGE#','Gene','Guide_Name']].agg('.'.join, axis=1)

    try:
        fa_score_df = graphing_results_df[['Guide_Name','fitness_score','std']]
    except:
        fa_score_df = graphing_results_df[['Guide_Name','fitness_score']]
    
    y_limit= _find_max_y(fa_score_df)
    
    bar_num = len(fa_score_df['Guide_Name'].to_list())
    
    
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
                        x='Guide_Name',
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
                        x='Guide_Name',
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
    #plt.autoscale()
    
    plt.show()
    

if __name__ == "__main__":
    main()

print('Graphing Complete.')