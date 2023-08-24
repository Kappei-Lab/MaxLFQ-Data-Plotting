#!/usr/bin/env python
# coding: utf-8


#Importing libraries 

import pandas as pd
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats
from scipy.stats import skew
import statsmodels.stats.multitest as multitest
from matplotlib.transforms import Bbox
import math
import os
import ast
import argparse
import json
import time
import requests
import sys
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from itertools import chain
from matplotlib.offsetbox import AnchoredText
from pylab import *
import urllib.parse
import urllib.request
import io

plt.rcParams['axes.unicode_minus']=False
plt.rcParams['pdf.fonttype']='42'
plt.rcParams['font.family'] = 'Arial'


def read_protein_pgroups(proteingroups_path,directory):
    

    """ Read the protein group file and create a directory for results
    Args:
    proteingroups_path: location of the protein groups file (just the file name)
    directory: current directory
    Returns:
    pgroups: dataframe containining the proteins groups data read from the file
    """


    global new_directory # directory name for analysis data
    global file_metadata

    analysis_indices = [int(x.split("_")[1]) for x in os.listdir(".") if x.startswith('Analysis_') and os.path.isdir(x)] # Find indices of Analysis folders
    
    if len(analysis_indices)>0:
        folder_no = max(analysis_indices)+1
    else:
        folder_no = 1

    new_directory=directory+"Analysis_"+str(folder_no) # create the directory name as "Analysis_1"

    os.mkdir(new_directory) # Now create the directory

    file_metadata = new_directory+os.sep+"Analysis_metadata.txt"

    pgroups=pd.read_table(proteingroups_path,sep="\t",dtype={'Protein names':str,'Gene names':str,'Only identified by site': str,'Reverse':str,'Potential contaminant':str})
    print("Initial dimensions of the table are ",pgroups.shape)
    print("\nColumn names are below:- \n\n",list(pgroups.columns))

    f = open(file_metadata, "a")
    f.write("\n..........Beginning of Analysis........\n")
    f.write("\nInitial dimensions of the table are : {},{}".format(pgroups.shape[0],pgroups.shape[1]))
    f.close()

    return pgroups

def remove_unreliable_proteins(pgroups,file_metadata):
    """
    This function removes the rows containing '+' in any of the columns: 'Only identified by site',
    'Reverse','Potential contaminant'
    Args:
    pgroups: dataframe to be processed
    file_metadata: name of the metadata file
    Returns:
    pgroup2: cleaned dataframe

    """
    pgroups2=pgroups.copy()

    if 'Only identified by site' in pgroups2.columns.values: # Check whether the column exist
        if len(pgroups2['Only identified by site'])!=(pgroups2['Only identified by site'].isnull().sum()):
            pgroups2=pgroups2[~(pgroups2['Only identified by site']=="+")]
    if 'Reverse' in pgroups2.columns.values: # Check whether the column exist
        if len(pgroups2['Reverse'])!=(pgroups2['Reverse'].isnull().sum()):
                pgroups2=pgroups2[~(pgroups2['Reverse']=="+")]
    if 'Potential contaminant' in pgroups2.columns.values: # Check whether the column exist
        if len(pgroups2['Potential contaminant'])!=(pgroups2['Potential contaminant'].isnull().sum()):
            pgroups2=pgroups2[~(pgroups2['Potential contaminant']=="+")]

    print("After cleaning dimensions of the table are ",pgroups2.shape)

    f = open(file_metadata, "a")
    f.write("\nAfter cleaning dimensions of the table are : {},{}".format(pgroups2.shape[0],pgroups2.shape[1]))
    f.close()

    return pgroups2


def connect(host='https://www.uniprot.org/uploadlists/'):
    """This function checks if a connection can be established with the host(webpage given)
    Args: webpage address
    Returns: True if able to connect otherwise false
    """
    try:
        urllib.request.urlopen(host) #Python 3.x
        return True
    except:
        return False

def get_name(Ids):
    """This function gives the genes names based on their protein Id input from Uniprot
    Args: Protein Ids
    Returns: Dataframe (table) with two columns "From" and "To" with protein Ids in columns "From" and their gene names on column "To"
    Note: Protein Ids that do not have a gene name, are not in the table that is returned
    """
    url = 'https://www.uniprot.org/uploadlists/'
    print(Ids)
    query= " ".join(Ids.values)
    params = {
    'from': 'ACC+ID',
    'to': 'GENENAME',
    'format': 'tab',
    'query': query
    }
    
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    rawData = pd.read_table(io.StringIO(response.decode('utf-8')))
    print(response.decode('utf-8'))
    
    
    return rawData

def clean_names(pgroups,plot_dict):
    """
    This function replaces the null values in the 'Protein names' column by the first item
    in Majority protein IDs' separated by ';'
    Then splits the  'Gene names' separated by ';' and select the first item
    """
    pnames = pgroups['Gene names'].copy()
    missing_bool = pgroups['Gene names'].isnull()
    first_id = pgroups['Majority protein IDs'].apply(lambda x:x.split(';')[0])
    
    if(plot_dict['Gene_names_from_Uniprot']):
        connected=connect(host='https://www.uniprot.org/uploadlists/')
        if(not(connected)):
            raise ConnectionError("Not Connected to the Internet")
        missing_gene_names=get_name(first_id[missing_bool])
        found=[]
        for i in range(len(first_id)):
            found.append(first_id.iloc[i] in missing_gene_names.From.values)
        mgn_dict={missing_gene_names["From"].loc[i]:missing_gene_names["To"].loc[i] for i in range(len(missing_gene_names))}
        first_id[found]=first_id[found].replace(mgn_dict)
        
    pnames[missing_bool] = first_id[missing_bool]
    
    pgroups['Gene names'] = pnames
    
    all_names=pgroups['Gene names'].str.split(';', expand=True)
    all_custom_proteins = list(chain.from_iterable(plot_dict['custom_groups']))
    custom_names_bool=all_names.isin(all_custom_proteins)
    in_custom=custom_names_bool.any(axis=1)

    for i in all_names:
        for j in all_names[in_custom].index:
            if(custom_names_bool[in_custom].loc[j,i]):
                all_names.loc[j,0]= all_names[custom_names_bool].loc[j,i]
    
    pgroups['Gene names'] = all_names[0]

    return pgroups

def select_lfq_cols(pgroups2,file_metadata):

	"""
	This function selects the lfq columns for control and treatment and removes the rows with
	large number of missing values
	Args:
	pgroups2: Input dataframe
	Returns:
	pgroup3: filtered dataframe

	"""


	print("LFQ intensity columns are \n")
	colnames=pgroups2.columns.values
	lfq_cols = [col for col in colnames if "LFQ" in col]
	print_df=pd.DataFrame({'LFQ':lfq_cols,'no':range(1,len(lfq_cols)+1)})
	print(print_df)

	is_t_correct='n'
	t_counter=1
	while is_t_correct!='y':

		treatment_ind = [int(x) for x in input("\nEnter the order of treatment columns (integers separated by space) ").split()]
		treatment_ind=[i-1 for i in treatment_ind]

		global treatment_colnames
		treatment_colnames=pd.Series(lfq_cols)[treatment_ind]
		print("\nSelected treatment columns are: \n")
		print(treatment_colnames)
		is_t_correct=input("\nIs this correct, enter y/n: ")
		t_counter+=1
		if t_counter>3:
			print("Maximum attempts exceeded. Please start from the beginning")
			break
			
	if is_t_correct=='y':

		is_c_correct='n'
		c_counter=1
		while is_c_correct!='y':

			control_ind = [int(x) for x in input("\nEnter the order of control columns (integers separated by space) ").split()]
			control_ind=[i-1 for i in control_ind]

			global control_colnames
			control_colnames=pd.Series(lfq_cols)[control_ind]
			print("\nSelected control columns are: \n")
			print(control_colnames)
			is_c_correct=input("\nIs this correct, enter y/n: ")
			c_counter+=1
			if c_counter>3:
				print("Maximum attempts exceeded. Please start from the beginning")
				break
		
	if is_t_correct=='y' and is_c_correct=='y':

		# Extract treatment and control LFQ values to separate dataframes
		treatment_LFQ=pgroups2[treatment_colnames]
		control_LFQ=pgroups2[control_colnames]

		nonzeros_control=(control_LFQ!=0).sum(axis=1)
		values_control=pd.Series(nonzeros_control.value_counts())
		values_control=values_control.sort_index()
		print("Number of non-zero values in 0(none)",*range(1,len(control_colnames)+1)," out of the",len(control_colnames),"control columns")
		print(values_control)

		nonzeros_treatment=(treatment_LFQ!=0).sum(axis=1)
		values_treatment=pd.Series(nonzeros_treatment.value_counts())
		values_treatment=values_treatment.sort_index()
		print("Number of non-zero values in 0(none)",*range(1,len(control_colnames)+1),"out of the",len(treatment_colnames),"treatment columns")
		print(values_treatment)


		# Delete the rows with zero values above the cut-off

		n=int(input("Enter the minimum number of non-zero values required either in control or treatment for each protein: "))
		pgroups3=pgroups2[~((nonzeros_control<n) & (nonzeros_treatment<n))]
		print("Rows has been removed ! \nDimensions of the dataframe are ",pgroups3.shape)
		
		f = open(file_metadata, "a")
		f.write("\nSelected treatment columns are :")
		for col in treatment_colnames:
			f.write(" {}".format(col))
		
		f.write("\nSelected control columns are :")
		for col in control_colnames:
			f.write(" {}".format(col))
		f.write("\nNumber of non-zero values in 0(none),1 ... n out of the n control columns are \n{}".format(values_control))
		f.write("\nNumber of non-zero values in 0(none),1 ... n out of the n treatment columns are \n{}".format(values_treatment))
		f.write("\nMinimum number of non-zero values selected is {}".format(n))
		
		f.write("\nAfter filtering dimensions of the table are : {},{}".format(pgroups3.shape[0],pgroups3.shape[1]))
				
		f.close()
		return pgroups3
		
	else:
		return None

def log_impute_col(column,plot_dict,nth_iter,file_metadata):

	"""
	This function log2 transform a column and replaces null values using lowest 5th percentile
	Args:
	column: An lfq intensity column
	impute_perc: lowest percentage of values to be used for imputation calculation
	std_shift: number of stds by which mean of the normal distribution to be shifted
	nth_iter: nth time this imputation is being performed
	file_metadata: name of the metadata_file
	Returns:
	log2_column: imputed and log2 transformed column

	"""
		
	col_len=len(column) # length the column to be imputed
	zeroes_len=sum(column==0) # number of zero values in the column
	

	if nth_iter==1:
		print(column.name,col_len,zeroes_len)
		f = open(file_metadata,"a")
		f.write("\n{} : {}".format(column.name,zeroes_len))
		f.close()

	if (zeroes_len!=0):
		column=column.replace(0,np.nan) # Change zeroes to nan before taking the log
	log2_column=np.log2(column)
	non_none_column=log2_column.dropna() # create a new column containing the non-null values

	percentile_nth=np.percentile(non_none_column,plot_dict['impute_perc']) # value of nth percentile
	selection_nth=non_none_column[non_none_column<=percentile_nth] # values below nth percentile
	selection_mean=np.mean(selection_nth) # mean of the values below nth percentile
	selection_std=np.std(selection_nth) # std of the values below nth percentile
	if plot_dict['impute_seed'] is not None:
		np.random.seed(nth_iter*100+plot_dict['impute_seed'])
	replacement_sample=np.random.normal((selection_mean+plot_dict['std_shift']*selection_std),selection_std,zeroes_len)


	# # save the distribution of values before and after imputation
	fig=plt.figure(figsize=(18,6))
	ax1=fig.add_subplot(1,3,1)
	ax1.hist(non_none_column)
	ax1.axvline(x=percentile_nth,zorder=2,linestyle='--',color='grey')
	ax1.set_title("Distribution before imputation")
	ax2=fig.add_subplot(1,3,2)
	ax2.hist(replacement_sample)
	ax2.set_title("Distribution of imputed values")
	# shuffle the above sample
	np.random.shuffle(replacement_sample)
	# Go through the nan values replace them with new values
	log2_column[log2_column.isnull()]=replacement_sample
	ax3=fig.add_subplot(1,3,3)
	ax3.hist(log2_column)
	ax3.set_title("Distribution after imputation")


	imputation_dir=new_directory+os.sep+"imputation_figures"

	if os.path.isdir(imputation_dir): # Check whether the directory already exist
		pass
	else:
		os.mkdir(imputation_dir)


	if nth_iter==1:
		imputed_figure_name_1=imputation_dir+os.sep+column.name+'.png'
		fig.savefig(imputed_figure_name_1,transparent=True,bbox_inches="tight",dpi=150)
		
	plt.close(fig)

	return log2_column

def multiple_impute(pgroups3,plot_dict,nth_iter,file_metadata):

	"""
	This function log2 transform all the lfq columns of a dataframe
	Args:
	pgroups3: A filtered dataframe from the previous step
	impute_perc: lowest percentage of values to be used for imputation calculation
	std_shift: number of stds by which mean of the normal distribution to be shifted
	nth_iter: nth this imputation is being performed
	file_metadata: name of the metadata_file
	Returns:
	pgroups: Transformed Dataframe

	"""
	print('\nImputation {} is in progress'.format(nth_iter))

	if nth_iter==1:
		print("\nColumn length and the number of zeroes in the column are\n ")
		
		f = open(file_metadata,"a")
		f.write("\n\nNumber of imputed values in the selected columns are")
		f.close()

	pgroups=pd.DataFrame()
	pgroups["Gene names"]=pgroups3['Gene names']

	for i in range(len(control_colnames)):
		
		
		name="control_log2_"+control_colnames.iloc[i]+'_'+str(nth_iter)
		column=pgroups3[control_colnames.iloc[i]]
		pgroups[name]=log_impute_col(column,plot_dict,nth_iter,file_metadata)
		
	for i in range(len(treatment_colnames)):
		name="treatment_log2_"+treatment_colnames.iloc[i]+'_'+str(nth_iter)
		
		column=pgroups3[treatment_colnames.iloc[i]]
		pgroups[name]=log_impute_col(column,plot_dict,nth_iter,file_metadata)

	return pgroups

def t_test(df,nth_iter):
	"""
	This function perform a ttest for each row of a dataframe between imputed treatment and control columns
	Args:
	df: A Dataframe containing imputed log2 transfored imputed values
	nth_iter: nth iteration of a number of imputations
	Returns:
	t_test_output: A dataframe containing Gene names, difference in means, and pvalues

	"""
	t_test_output=pd.DataFrame()
	t_test_output['Gene names']=df['Gene names']

	new_control_colnames=['control_log2_'+name for name in control_colnames+"_"+str(nth_iter)]
	new_treatment_colnames=['treatment_log2_'+name for name in treatment_colnames+"_"+str(nth_iter)]

	result =stats.ttest_ind(df[new_treatment_colnames],df[new_control_colnames],axis=1)
	mean_diffs = df[new_treatment_colnames].mean(axis=1) -  df[new_control_colnames].mean(axis=1)

	t_test_output["diff_"+str(nth_iter)]=mean_diffs
	t_test_output["pval_"+str(nth_iter)]=result[1]
	t_test_output["padj_"+str(nth_iter)]=multitest.multipletests(result[1],method="fdr_bh")[1]

	return t_test_output

def mean_of_imputed(pgroups_imputed):

	"""
	This function calculates mean of the imputed values for each replicates
	Args:
	pgroups_imputed: A dataframe containing multiple imputed values for each replicates
	"""
		
	colnames=pgroups_imputed.columns.values
	pgroups_imputed_sel=pgroups_imputed[[s for s in colnames if 'log2' in s]].copy()

	for col in control_colnames:

		cols_sel=[s for s in colnames if ('control_log2' in s) and (col in s)]
		new_col_name='mean_control_log2_'+col
		new_col=pgroups_imputed_sel[cols_sel].mean(axis=1)
		pgroups_imputed_sel.loc[:,new_col_name]=new_col.values

	for col in treatment_colnames:

		cols_sel=[s for s in colnames if ('treatment_log2' in s) and (col in s)]
		new_col_name='mean_treatment_log2_'+col
		new_col=pgroups_imputed_sel[cols_sel].mean(axis=1)
		pgroups_imputed_sel.loc[:,new_col_name]=new_col.values

	col_sel=[col for col in pgroups_imputed_sel.columns.values if 'mean' in col]
		
		
	return pgroups_imputed_sel[col_sel]

def iterate_impute_t_test(pgroups3,plot_dict,file_metadata):
	"""
	This function perform imputation and t-tests, n_iter number of times
	Args:
	pgroups3: filtered dataframe
	n_iter: number of times to repeat the imputation and t-tests
	Returns
	pgroups4: A dataframe with imputed values, ttest results
	"""
	pgroups4=pgroups3.copy()
	for i in range(1,plot_dict['n_iter']+1):
		pgroups3_imputed=multiple_impute(pgroups3,plot_dict,i,file_metadata)
		t_test_output=t_test(pgroups3_imputed,i)
		pgroups4=pd.concat([pgroups4,pgroups3_imputed.drop(['Gene names'],axis=1)],axis=1,sort=False)
		pgroups4=pd.concat([pgroups4,t_test_output.drop(['Gene names'],axis=1)],axis=1,sort=False)

	return pgroups4

def compile_convert(pgroups4):

	"""
	This function calculate the averages and standard deviations of diff, pval & padj
	Args:
	pgroups4: Dataframe containing diff, pval and padj
	Returns
	final_df: A dataframe containing the averages and standard deviations
	"""
    
	match_1 = [s for s in pgroups4.columns if "mean_control_log2_" in s]
	control_df=pgroups4[match_1]

	match_2 = [s for s in pgroups4.columns if "mean_treatment_log2_" in s]
	treatment_df=pgroups4[match_2]

	diff_cols=[s for s in pgroups4.columns if "diff" in s]
	diff_df=pgroups4[diff_cols]
	pval_cols=[s for s in pgroups4.columns if "pval" in s]
	pval_df=-np.log10(pgroups4[pval_cols])
	padj_cols=[s for s in pgroups4.columns if "padj" in s]
	padj_df=-np.log10(pgroups4[padj_cols])


	final_df=pd.DataFrame()
	final_df['Gene names']=pgroups4['Gene names']
    
	final_df['Control_mean_log2']=control_df.mean(axis=1)
	final_df['Treatment_mean_log2']=treatment_df.mean(axis=1)
    
	final_df['mean_log2_fold_change']=diff_df.mean(axis=1)
	final_df['std_log2_fold_change']=diff_df.std(axis=1)

	final_df['mean_-log10_pval']=pval_df.mean(axis=1)
	final_df['std_-1og10_pval']=pval_df.std(axis=1)

	final_df['mean_-log10_padj']=padj_df.mean(axis=1)
	final_df['std_-1og10_padj']=padj_df.std(axis=1)

	return final_df

def create_points(final_df,plot_dict):
	"""
	This function create a dataframe for points to be plotted
	Args:
	final_df: final dataframe containing results of statistical tests
	plot_dict: plotting parameters
	output:
	points_df: dataframe containing information for plotting different types of points
	"""

	n_custom_groups=len(plot_dict['custom_groups'])
	number_total_groups=3+n_custom_groups # total number of different types of points
	points_df=pd.DataFrame()
	points_df['up_bool']=((final_df['mean_log2_fold_change']>plot_dict['foldchange_cutoff'])& (final_df['mean_-log10_pval']>-np.log10(plot_dict['pval_cutoff'])))
	points_df['down_bool']=((final_df['mean_log2_fold_change']<-plot_dict['foldchange_cutoff'])& (final_df['mean_-log10_pval']>-np.log10(plot_dict['pval_cutoff'])))

	for i in range(n_custom_groups):
		points_df['custom_bool_'+str(i+1)]=final_df['Gene names'].isin(plot_dict['custom_groups'][i])

	# Remove up or down regulated based on the user input
	if (not plot_dict['upreg']):
		points_df['up_bool']=False
	if (not plot_dict['downreg']):
		points_df['down_bool']=False
		
		
	# remove custom group proteins from up & down regulated one

	points_df['sel_1']=points_df['up_bool'].copy()
	points_df['sel_2']=points_df['down_bool'].copy()
	for i in range(n_custom_groups):
		points_df['sel_1']=points_df['sel_1'] & (~points_df['custom_bool_'+str(i+1)])
		points_df['sel_2']=points_df['sel_2'] & (~points_df['custom_bool_'+str(i+1)])

	points_df['sel_3']=(~points_df['sel_1']) & (~points_df['sel_2'])

	for i in range(n_custom_groups):
		points_df['sel_3']=points_df['sel_3'] & (~points_df['custom_bool_'+str(i+1)])
		
	return points_df

def float_string(value):
	"""
	convert a floating point into integer-string
	"""
	if pd.isnull(value):
		return value
	else:
		return str(int(value))

def create_labels_legends(final_df,points_df,plot_dict):
	"""
	This function create a dataframe with labels and legends for the points
	Args:
	final_df: results of statistical tests
	num_legend_df_final: dataframe containing labels and legends
	"""
	num_legend_df=final_df.copy()

	up_num=num_legend_df['mean_log2_fold_change'][points_df['up_bool']].rank(ascending=False)
	num_legend_df['up_num']=up_num

	down_num=num_legend_df['mean_log2_fold_change'][points_df['down_bool']].rank(ascending=True)
	num_legend_df['down_num']=down_num
    
	start_num_up = 1
	start_num_down = 1
	if pd.notnull(num_legend_df['up_num'].max()):
		start_num_up=num_legend_df['up_num'].max()+1
	if pd.notnull(num_legend_df['down_num'].max()):
		start_num_down=num_legend_df['down_num'].max()+1


	for j in range(len(plot_dict['custom_groups'])):
		col_name='custom_bool_'+str(j+1)
		custom_df = pd.DataFrame(num_legend_df['mean_log2_fold_change'][points_df[col_name] & (~(points_df['up_bool']|points_df['down_bool']))])
		custom_sort_bool= custom_df['mean_log2_fold_change'] > 0
		custom_up_ord=custom_df[custom_sort_bool].sort_values(by=['mean_log2_fold_change'],ascending = False)
		custom_down_ord=custom_df[~custom_sort_bool].sort_values(by=['mean_log2_fold_change'],ascending = True)
		custom_up_ord['order']=np.arange(start_num_up,start_num_up+custom_up_ord.shape[0])
		custom_down_ord['order']=np.arange(start_num_down,start_num_down+custom_down_ord.shape[0])
		custom_df_ordered = pd.concat([custom_up_ord,custom_down_ord])
		if pd.notnull(custom_up_ord['order'].max()):
			start_num_up=custom_up_ord['order'].max()+1
		if pd.notnull(custom_down_ord['order'].max()):
			start_num_down=custom_down_ord['order'].max()+1

        
		custom_name='custom_'+str(j+1)+"_num"
		num_legend_df[custom_name]=custom_df_ordered['order']
		#start_num = start_num + custom_df_ordered.shape[0]

	columns=['up_num','down_num']+[('custom_'+str(i+1)+'_num') for i in range(len(plot_dict['custom_groups']))]
	num_legend_df.sort_values(by=columns,inplace=True)
	for col in columns:
		num_legend_df[col]=num_legend_df[col].apply(float_string)
		
	num_legend_df['color']=plot_dict['color_up']
	num_legend_df['color']=num_legend_df['color'].where(~points_df['sel_2'],plot_dict['color_down'])

	for i in range(len(plot_dict['custom_groups'])):
		col='custom_bool_'+str(i+1)
		num_legend_df['color']=num_legend_df['color'].where(~points_df[col],plot_dict['custom_colors'][i])


	num_legend_df_final=num_legend_df[(num_legend_df[columns].isnull().sum(axis=1)!=len(columns))].copy()

	def get_num(row):
		return row[~row.isnull()].values[0]

	num_legend_df_final['number']= num_legend_df_final[columns].apply(get_num,axis=1)
	num_legend_df_final['legend']=num_legend_df_final['number']+'. '+num_legend_df_final['Gene names']

	if(any(num_legend_df_final['mean_log2_fold_change']>0) and any(num_legend_df_final['mean_log2_fold_change']<0)):
		num_legend_df_final["row"] = np.nan
		num_legend_df_final["col"] = np.nan
    
		num_legend_df_final.loc[(num_legend_df_final['mean_log2_fold_change']>0),'row']=num_legend_df_final.loc[(num_legend_df_final['mean_log2_fold_change']>0),'number'].apply(lambda x:(int(x)-1)%15).values
		num_legend_df_final.loc[(num_legend_df_final['mean_log2_fold_change']>0),'col']=num_legend_df_final.loc[(num_legend_df_final['mean_log2_fold_change']>0),'number'].apply(lambda x:(int(x)-1)//15).values
    
		num_legend_df_final.loc[(num_legend_df_final['mean_log2_fold_change']<0),'row']=(num_legend_df_final.loc[(num_legend_df_final['mean_log2_fold_change']<0),'number'].apply(lambda x:int(x)-1)%15)+17
		num_legend_df_final.loc[(num_legend_df_final['mean_log2_fold_change']<0),'col']=num_legend_df_final.loc[(num_legend_df_final['mean_log2_fold_change']<0),'number'].apply(lambda x:int(x)-1)//15
	elif(any(num_legend_df_final['mean_log2_fold_change']>0) or any(num_legend_df_final['mean_log2_fold_change']<0)):
		row_numbers = pd.Series(range(num_legend_df_final.shape[0]))
		num_legend_df_final['row']=row_numbers.apply(lambda x:x%32).values
		num_legend_df_final['col']=row_numbers.apply(lambda x:x//32).values

	return num_legend_df_final

def save_fig(final_df,plot_dict,file_name):

	"""
	This function create and save the volcano plot
	Args:
	final_df: A dataframe containing mean values and standard deviations of diff, pval & padj
	plot_dict: A dictionary containing information required for plotting

	"""

	fig=plt.figure(figsize=(10,10))
	ax=fig.add_subplot(1,1,1)
    
	up_genes=((final_df['mean_log2_fold_change']>plot_dict['foldchange_cutoff'])& (final_df['mean_-log10_pval']>-np.log10(plot_dict['pval_cutoff'])))
	down_genes=((final_df['mean_log2_fold_change']<-plot_dict['foldchange_cutoff'])& (final_df['mean_-log10_pval']>-np.log10(plot_dict['pval_cutoff'])))

	## create the points_df dataframe

	points_df=create_points(final_df,plot_dict)
    
	## Create the labels and legends dataframe
	if(any(points_df["up_bool"]) or any(points_df["down_bool"]) or any(points_df["custom_bool_1"])):
		num_legend_df_final=create_labels_legends(final_df,points_df,plot_dict)
    

	## Plot the points with the right color and size
	if(plot_dict['Show_LFQ_intensity']):
		c_values_t=final_df['Treatment_mean_log2'][points_df['sel_1']]
		c_values_c=final_df['Control_mean_log2'][points_df['sel_2']]
		if(plot_dict['upreg'] and (~plot_dict['downreg'])):
			vmax,vmin=np.percentile(final_df['Treatment_mean_log2'],95),final_df['Treatment_mean_log2'].min()
		elif (plot_dict['downreg'] and (~plot_dict['upreg'])):
			vmax,vmin=np.percentile(final_df['Control_mean_log2'],95),final_df['Control_mean_log2'].min()
		else:
			vmax,vmin=max(np.percentile(final_df['Treatment_mean_log2'],95),np.percentile(final_df['Control_mean_log2'],95)),min(final_df['Treatment_mean_log2'].min(),final_df['Control_mean_log2'].min())
            
		#Changing color gradient of points based on intensity value
		a2=ax.scatter(final_df['mean_log2_fold_change'][points_df['sel_2']],final_df['mean_-log10_pval'][points_df['sel_2']],edgecolors="black",
                        c=c_values_c,cmap=plt.set_cmap(plot_dict['colormap_down']),s=250,alpha=1,zorder=2,vmin=vmin, vmax=vmax)

		a1=ax.scatter(final_df['mean_log2_fold_change'][points_df['sel_1']],final_df['mean_-log10_pval'][points_df['sel_1']],edgecolors="black",
                        c=c_values_t,cmap=plt.set_cmap(plot_dict['colormap_up']),s=250,alpha=1,zorder=2,vmin=vmin, vmax=vmax)
		if(plot_dict["upreg"]==False):
			a1=ax.scatter(final_df['mean_log2_fold_change'][up_genes],final_df['mean_-log10_pval'][up_genes],edgecolors="black",
				c=final_df['Treatment_mean_log2'][up_genes],cmap=plt.set_cmap(plot_dict['colormap_up']),s=50,alpha=1,zorder=2,vmin=vmin, vmax=vmax)
		if(plot_dict["downreg"]==False):
			a2=ax.scatter(final_df['mean_log2_fold_change'][down_genes],final_df['mean_-log10_pval'][down_genes],edgecolors="black",
                    c=final_df['Control_mean_log2'][down_genes],cmap=plt.set_cmap(plot_dict['colormap_down']),s=50,alpha=1,zorder=2,vmin=vmin, vmax=vmax)

		#Plotting Color bars for gradients
		if(any(up_genes)):
			cax1= inset_axes(ax,
                       width="30%",  # width = 5% of parent_bbox width
                       height="2%",  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(0.68, 0.03, 1, 1),
                       bbox_transform=ax.transAxes
                       )
			plt.colorbar(a1,cax=cax1,orientation='horizontal')
            
		if(any(down_genes)):
			cax2= inset_axes(ax,
                    width="30%",
                    height="2%",
                    loc='lower left',
                    bbox_to_anchor=(0.01, 0.03, 1, 1),
                    bbox_transform=ax.transAxes
                    )
        
			plt.colorbar(a2,cax=cax2, orientation='horizontal')

		for i in range(len(plot_dict['custom_groups'])):
			a3=ax.scatter(final_df['mean_log2_fold_change'][points_df['custom_bool_'+str(i+1)] & (final_df['mean_log2_fold_change']>0)],final_df['mean_-log10_pval'][points_df['custom_bool_'+str(i+1)] & (final_df['mean_log2_fold_change']>0)],edgecolors="black",c=final_df['Treatment_mean_log2'][points_df['custom_bool_'+str(i+1)] & (final_df['mean_log2_fold_change']>0)],cmap=plt.set_cmap(plot_dict['custom_colormap'][i]),s=250,alpha=1,zorder=3,vmin=vmin, vmax=vmax)
			a4=ax.scatter(final_df['mean_log2_fold_change'][points_df['custom_bool_'+str(i+1)] & (final_df['mean_log2_fold_change']<0)],final_df['mean_-log10_pval'][points_df['custom_bool_'+str(i+1)] & (final_df['mean_log2_fold_change']<0)],edgecolors="black",c=final_df['Control_mean_log2'][points_df['custom_bool_'+str(i+1)] & (final_df['mean_log2_fold_change']<0)],cmap=plt.set_cmap(plot_dict['custom_colormap'][i]),s=250,alpha=1,zorder=3,vmin=vmin, vmax=vmax)           
        
	else:
		ax.scatter(final_df['mean_log2_fold_change'][points_df['sel_1']],final_df['mean_-log10_pval'][points_df['sel_1']],edgecolors="black",
				color=plot_dict['color_up'],s=250,alpha=1,zorder=2)
		ax.scatter(final_df['mean_log2_fold_change'][points_df['sel_2']],final_df['mean_-log10_pval'][points_df['sel_2']],edgecolors="black",
			color=plot_dict['color_down'],s=250,alpha=1,zorder=2)
		for i in range(len(plot_dict['custom_groups'])):
			ax.scatter(final_df['mean_log2_fold_change'][points_df['custom_bool_'+str(i+1)]],final_df['mean_-log10_pval'][points_df['custom_bool_'+str(i+1)]],edgecolors="black",
				color=plot_dict['custom_colors'][i],s=250,alpha=1,zorder=3)
		if(plot_dict["upreg"]==False):
			a1=ax.scatter(final_df['mean_log2_fold_change'][up_genes],final_df['mean_-log10_pval'][up_genes],edgecolors="black",
				color=plot_dict['color_up'],s=50,alpha=1,zorder=2)
		if(plot_dict["downreg"]==False):
			a2=ax.scatter(final_df['mean_log2_fold_change'][down_genes],final_df['mean_-log10_pval'][down_genes],edgecolors="black",
				color=plot_dict['color_down'],s=50,alpha=1,zorder=2)
    
	ax.scatter(final_df['mean_log2_fold_change'][points_df['sel_3']],final_df['mean_-log10_pval'][points_df['sel_3']],edgecolors="black",color=plot_dict['base_color'],s=50,alpha=0.3,zorder=1)

	## plot the error bars for the significant points

	#error_bar_bool=(~num_legend_df_final['up_num'].isnull())|(~num_legend_df_final['down_num'].isnull())
	if(any(points_df["up_bool"]) or any(points_df["down_bool"]) or any(points_df["custom_bool_1"])):
		ax.errorbar(x=num_legend_df_final['mean_log2_fold_change'][points_df['sel_1']],
					y=num_legend_df_final['mean_-log10_pval'][points_df['sel_1']],
					yerr=num_legend_df_final['std_-1og10_pval'][points_df['sel_1']],
					xerr=num_legend_df_final['std_log2_fold_change'][points_df['sel_1']],
					ls='none',zorder=0,ecolor=plot_dict['color_up'],elinewidth=.5,capsize=2)

		ax.errorbar(x=num_legend_df_final['mean_log2_fold_change'][points_df['sel_2']],
					y=num_legend_df_final['mean_-log10_pval'][points_df['sel_2']],
					yerr=num_legend_df_final['std_-1og10_pval'][points_df['sel_2']],
					xerr=num_legend_df_final['std_log2_fold_change'][points_df['sel_2']],
					ls='none',zorder=0,ecolor=plot_dict['color_down'],elinewidth=.5,capsize=2)

		for i in range(len(plot_dict['custom_groups'])):
			if plot_dict['custom_errorbar'][i]==True:

				bool_column = 'color'
				custom_bool = (num_legend_df_final[bool_column]==plot_dict['custom_colors'][i])

				ax.errorbar(x=num_legend_df_final['mean_log2_fold_change'][custom_bool],
					y=num_legend_df_final['mean_-log10_pval'][custom_bool],
					yerr=num_legend_df_final['std_-1og10_pval'][custom_bool],
					xerr=num_legend_df_final['std_log2_fold_change'][custom_bool],
					ls='none',zorder=0,ecolor=plot_dict['custom_colors'][i],elinewidth=.5,capsize=2)

	## Number the points

		for i in range(0,num_legend_df_final.shape[0]):
			ax.text(num_legend_df_final['mean_log2_fold_change'].iloc[i],
					num_legend_df_final['mean_-log10_pval'].iloc[i],num_legend_df_final['number'].iloc[i],ha='center',va='center_baseline',
					fontsize=8,family=plot_dict['plot_font'],color=plot_dict['number_color'])

	# Plot the legends #
		if (any(num_legend_df_final['mean_log2_fold_change']>0) and any(num_legend_df_final['mean_log2_fold_change']<0)):
			ax.text(1.01,0.99,plot_dict['legend_upreg'],fontsize=plot_dict['legend_size'],color=plot_dict['color_up'],transform=ax.transAxes,verticalalignment='top',
				family=plot_dict['plot_font'],weight="bold")
			ax.text(1.01,0.96-(.03*16),plot_dict['legend_downreg'],fontsize=plot_dict['legend_size'],color=plot_dict['color_down'],transform=ax.transAxes,verticalalignment='top',
				family=plot_dict['plot_font'],weight="bold")
		elif any(num_legend_df_final['mean_log2_fold_change']>0):
			ax.text(1.01,0.99,plot_dict['legend_upreg'],fontsize=plot_dict['legend_size'],color=plot_dict['color_up'],transform=ax.transAxes,verticalalignment='top',
				family=plot_dict['plot_font'],weight="bold")
		elif any(num_legend_df_final['mean_log2_fold_change']<0):
			ax.text(1.01,0.99,plot_dict['legend_downreg'],fontsize=plot_dict['legend_size'],color=plot_dict['color_down'],transform=ax.transAxes,verticalalignment='top',
				family=plot_dict['plot_font'],weight="bold")
        
		for i in range(num_legend_df_final.shape[0]):
			ax.text(1.01+(.20*num_legend_df_final['col'].iloc[i]),0.96-(.03*num_legend_df_final['row'].iloc[i]),num_legend_df_final['legend'].iloc[i],
					fontsize=plot_dict['legend_size'],color=num_legend_df_final['color'].iloc[i],transform=ax.transAxes,verticalalignment='top',
					family=plot_dict['plot_font'])

	ax.axvline(linewidth=4,color="silver",zorder=0,alpha=.5)

	ax.axhline(linewidth=.5,alpha=0.5,color="black",zorder=0,linestyle='--',y=-np.log10(plot_dict['pval_cutoff']),dashes=(10,10))
	ax.axvline(linewidth=.5,alpha=0.5,color="black",zorder=0,linestyle='--',dashes=(10,10),x=plot_dict['foldchange_cutoff'])
	ax.axvline(linewidth=.5,alpha=0.5,color="black",zorder=0,linestyle='--',dashes=(10,10),x=-plot_dict['foldchange_cutoff'])

	ax.tick_params(axis='both',labelsize=15)

	if (plot_dict['x_lower_limit'] is not None) & (plot_dict['x_upper_limit'] is not None):
		ax.set_xlim(plot_dict['x_lower_limit'],plot_dict['x_upper_limit'])
	else:
		max_fc=final_df['mean_log2_fold_change'].max()
		min_fc=final_df['mean_log2_fold_change'].min()
		x_lim=max(abs(max_fc),abs(min_fc))
		if(x_lim%1>0.8):
			x_lim=math.ceil(x_lim)+0.5
		else:
			x_lim=math.ceil(x_lim)
		ax.set_xlim(-x_lim,x_lim)

	if (plot_dict['y_lower_limit'] is not None) & (plot_dict['y_upper_limit'] is not None):
		ax.set_ylim(plot_dict['y_lower_limit'],plot_dict['y_upper_limit'])

	if(plot_dict['upper_right_text']!=''):
			at1 = AnchoredText(plot_dict['upper_right_text'],loc='upper right', prop=dict(size=12,fontweight="bold",color=plot_dict['color_up'],family=plot_dict['plot_font']))
			at1.patch.set_boxstyle("round,pad=0.2,rounding_size=0.2")
			at1.patch.set_fc("none")
			ax.add_artist(at1)
	if(plot_dict['upper_left_text']!=''):
			at2 = AnchoredText(plot_dict['upper_left_text'],loc='upper left', prop=dict(size=12,fontweight="bold",color=plot_dict['color_down'],family=plot_dict['plot_font']))
			at2.patch.set_boxstyle("round,pad=0.2,rounding_size=0.2")
			at2.patch.set_fc("none")
			ax.add_artist(at2)
			
	
	ax.set_xlabel(plot_dict['xlabel'],color='black',fontsize=plot_dict['xylabel_size'],family=plot_dict['plot_font'])
	ax.set_ylabel(plot_dict['ylabel'],color='black',fontsize=plot_dict['xylabel_size'],family=plot_dict['plot_font'])
	ax.set_title(plot_dict['title'],color='black',fontsize=plot_dict['title_size'],family=plot_dict['plot_font'])

	if(any(points_df["up_bool"]) or any(points_df["down_bool"]) or any(points_df["custom_bool_1"])):
		if (any(num_legend_df_final['mean_log2_fold_change']>0) and any(num_legend_df_final['mean_log2_fold_change']<0)):
			pw=0.2*max(max(num_legend_df_final['col'])+1,max((len(plot_dict['legend_upreg'])/14),(len(plot_dict['legend_downreg'])/14)))
			ax.add_patch(plt.Rectangle((1.005,1.002),pw,-0.49,facecolor='lightgrey', alpha=0.2,
						clip_box=True,clip_on=False,linewidth = 0,transform=ax.transAxes))
			ax.add_patch(plt.Rectangle((1.005,0.488),pw,-0.49,facecolor='lightgrey', alpha=0.2,
						clip_box=True,clip_on=False,linewidth = 0,transform=ax.transAxes))
		elif(any(num_legend_df_final['mean_log2_fold_change']>0) or any(num_legend_df_final['mean_log2_fold_change']<0)):
			pw=0.2*max(max(num_legend_df_final['col'])+1,max((len(plot_dict['legend_upreg'])/14),(len(plot_dict['legend_downreg'])/14)))
			ax.add_patch(plt.Rectangle((1.005,1.002),pw,-1.005,facecolor='lightgrey', alpha=0.2,
				clip_box=True,clip_on=False,linewidth = 0,transform=ax.transAxes))

	no=1
	figure_name_1=new_directory+os.sep+file_name+'_LFQ_volcano_plot_'+str(no)+'.pdf'
	figure_name_2=new_directory+os.sep+file_name+'_LFQ_volcano_plot_'+str(no)+'.png'


	while os.path.isfile(figure_name_1) or os.path.isfile(figure_name_1) :
		no+=1
		figure_name_1=new_directory+os.sep+file_name+'_LFQ_volcano_plot_'+str(no)+'.pdf'
		figure_name_2=new_directory+os.sep+file_name+'_LFQ_volcano_plot_'+str(no)+'.png'

	fig.savefig(figure_name_1,transparent=True,bbox_inches="tight",dpi=300)
	fig.savefig(figure_name_2,transparent=True,bbox_inches="tight",dpi=300)




def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("protein_groups_file", help="Location of the protein group file",
               type=str)
	parser.add_argument("lfq_plot_info", help="Location of the plot_info",
               type=str)
	args = parser.parse_args()

	directory=os.getcwd()+os.sep

	print('Current directory is \n',directory)

	#print(args.lfq_plot_info)
	#print(args.protein_groups_file)
    
	with open(args.lfq_plot_info) as f:
		plot_dict = json.load(f)
    
	ls=directory.split(os.sep)
	file_name=ls[len(ls)-2]

	pgroups=read_protein_pgroups(args.protein_groups_file,directory)

	writer=None
	writer=pd.ExcelWriter(new_directory+os.sep+file_name+"_Processed_ProteinGroups.xlsx")
	pgroups.to_excel(writer,'Original')


	pgroups2 = remove_unreliable_proteins(pgroups,file_metadata)
        
	pgroups2 = clean_names(pgroups2,plot_dict)

	pgroups2.to_excel(writer,'Filtered_1')

	pgroups3=select_lfq_cols(pgroups2,file_metadata)
	pgroups3.to_excel(writer,'Filtered_2')

	pgroups4=iterate_impute_t_test(pgroups3,plot_dict,file_metadata)
	pgroups4a=mean_of_imputed(pgroups4)
	pgroups4=pd.concat([pgroups4,pgroups4a],axis=1)

	final_df=compile_convert(pgroups4)

	pgroups5=pd.concat([pgroups4,final_df.drop(['Gene names'],axis=1)],axis=1)

	pgroups6=pgroups5.sort_values(by=['mean_log2_fold_change'],ascending=False)
	pgroups6.to_excel(writer,'ttest_log_imputed')

	colnames=pgroups6.columns.values
	selected_columns=['Protein IDs','Gene names','Protein names','Number of proteins','Peptides','Razor + unique peptides',
        'Unique peptides','Sequence coverage [%]', 'Unique + razor sequence coverage [%]',
        'Unique sequence coverage [%]', 'Mol. weight [kDa]',
        'Sequence length']
	match_1 = [s for s in colnames if "mean_control_log2_" in s]
	selected_columns.extend(match_1)

	match_2 = [s for s in colnames if "mean_treatment_log2_" in s]
	selected_columns.extend(match_2)

	selected_columns.extend(['mean_log2_fold_change', 'std_log2_fold_change', 'mean_-log10_pval', 'std_-1og10_pval'])

	pgroups_final=pgroups6[selected_columns]
	print("Saving final data into excel sheet ...")
	pgroups_final.to_excel(writer,'Final')
	writer.save()
	print("Done!")


	f=open(file_metadata,"a")
	f.write("\n\nPlotting parameters used are below:\n")

	for key,val in plot_dict.items():
		f.write("\n{}:{}".format(key,val))

	f.close()

	print("Plotting data ... ")
	save_fig(final_df,plot_dict,file_name)
	print("Done!")

	f=open(file_metadata,"a")
	t = time.localtime()
	current_time = time.strftime("%Y-%b-%d %H:%M:%S", t)

	f.write("\n\nAnalysis finished at {}".format(current_time))
	f.write("\n\nOutput files were saved to:\n{}".format(os.path.abspath(new_directory)))
    
	f.write("\n\n.........End of Analysis........")
	f.close()

if __name__ == "__main__":
    main()

