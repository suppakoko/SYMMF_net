# coding: utf-8
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import numpy as np
import seaborn as sns
from sklearn import metrics
from scipy.stats import ks_2samp
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap
from sklearn import preprocessing
import os
import random
import copy
from scipy import stats
import sys
import argparse

my_parser = argparse.ArgumentParser(description='Optional app description')
my_parser.add_argument('-np', help="cut-off P-value for SYMMF-net correlation", default='1e-13')
my_parser.add_argument('-iter', help="iteration count value for AUC calculation", default='1000')
my_parser.add_argument('-ecdf', help="cut-off ECDF value for selection of top SYMMFs", default='0.2')
my_parser.add_argument('-p', help="cut-off P-value for AUC calculation it should be bigger than 0.001", type=float)

# Execute parse_args()
args = my_parser.parse_args()

if args.p == None:
	print_string = '''
	***************************************
	Thank you for using SYMMF-net analysis!
	***************************************
	How to use: After finishing NMF analysis using NMF_analysis.R, 
	copy the SYMMF_net_analysis.py to the same folder and 
	excute the following command with the cut-off P-value.

	ex) python SYMMF_net_analysis.py -p [cut-off P-value]


	if you want to see more options
	python SYMMF_net_analysis.py -h or --help


	'''
	print(print_string)
	sys.exit()

## filter settings
val_it = int(args.iter)
ecdf_cut_v = float(args.ecdf)
nwp_val = float(args.p)

#microbe variable declare
dif_auc_val = 0
network_cut_off = 0

#cut off pval of correlations
whole_network_pval_cutoff_l = [float(args.np)]
microbe_contribution_cuttoff = 0
minmum_auc = 0
ecdf_cut_str = str(ecdf_cut_v).replace('.','')

option_values = '''
######################
### option values  ###
######################
Iteration of AUC calculation : %d
ECDF cut-off : %f
SYMMF-net edge cut-off : %f
AUC calculation cut-off : %s 

''' %(val_it, ecdf_cut_v, nwp_val, str(whole_network_pval_cutoff_l[0]))
print(option_values)

#lists for coefmatrix
data_l = []
name_l = []
mic_l = []
 
#microbiome count matrix read and put into the matrix
mic_count_mat_raw = open('micro_count_mat.txt','r')
mcm_r = mic_count_mat_raw.read()
mcm_ls = mcm_r.strip().replace('g__','').split('\n') #mcm_r split by line
mcm_col_l = mcm_ls[0].split(' ') # counts number of samples
groupl = list(set(mcm_col_l)) # how many sample groups
col_num = len(groupl) # how many kinds of samples?
mcm_mat = np.zeros((len(mcm_ls)-1,len(mcm_col_l))) # make blank matrix species count x number of samples
mcm_row_l = [] #list for species list gathering

# purity comparison table
purity_table_read = open('Clustering_purity_value.txt','r')
purity_table_data = purity_table_read.read()
purity_lines = purity_table_data.strip().replace(' ','\t').replace('\n','\t').split('\t')
purity_float = []
for x in purity_lines:
    purity_float.append(float(x))
purity_array = np.array(purity_float[:-3])
purity_table = purity_array.reshape(23,3)
df = pd.DataFrame(data=purity_table, columns=['K-label','NMF','Raw'] )
th = np.linspace(0, 2*np.pi, 128)
mpl.style.use('default')
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_title('Purity comparison', color='C0')
ax.set_xlabel('Number of meta-feature (Km)')
ax.set_ylabel('Clustering Purity (%)')
ax.plot(df['K-label'], df['NMF'], 'C1', label='NMF')
ax.plot(df['K-label'], df['Raw'], 'C2', label='Raw')
ax.legend()
fig.savefig('purification_comparison_plot.png', dpi=300)
print('Drawing Purification comparision plot is done!')

# data insert to empty dataframe
for ln in range(0,len(mcm_ls)-1):
    mcm_l = mcm_ls[ln+1].split(' ') #split one count data from mcm_ls
    mcm_row_l.append(mcm_l[0]) # species list gathering
    #filling mcm_mat with count data 
    for val_id in range(0, len(mcm_l)-1):
        mcm_mat[ln,val_id] = mcm_l[val_id+1]

#mcm_df is species count dataframe without redundance
mcm_df = pd.DataFrame(mcm_mat, index=mcm_row_l, columns=mcm_col_l) 
mcm_df = mcm_df.groupby(mcm_df.index, axis=0).mean()
kscol_num = col_num * 3
mcm_count_ks = np.zeros((len(mcm_row_l), kscol_num)) #matrix number of species X number of samples
mcm_count_g_l = []
mcm_count_col_l = []

mcm_ks_df = pd.DataFrame()
if os.path.isfile("microbe_AUC_calculation.xlsx"):
    pre_mcm_ks_df = pd.read_excel("microbe_AUC_calculation.xlsx")
    header = list(pre_mcm_ks_df.columns)[1:]
    gl = []
    for x in range(0,len(header),3):
        gnam = header[x]
        gl.append(gnam)
        gl.append(gnam)
        gl.append(gnam)
    labl = list(pre_mcm_ks_df.loc[0,][1:])
    microbe_l = list(pre_mcm_ks_df.iloc[2:,0])
    mcm_ks_df = pre_mcm_ks_df.iloc[2:,1:]
    mcm_ks_df.columns = [gl, labl]
    mcm_ks_df.index = microbe_l
else:
    #data extract and calculate K-S from microbial count matrix
    for sp_o in range(len(mcm_row_l)): #iteration of species counts
        sp_data = mcm_df.loc[mcm_row_l[sp_o],] #counts data load by species name
        n_col = 0
        for g_idx1 in range(col_num): #col_num is sample group counts
            gn1 = groupl[g_idx1]
            gn2 = 'others'
            gn1_data = list(sp_data[gn1])
            gn2_data = list(sp_data[sp_data.index != gn1])

            #generate ratio and choose dominant
            gn1_data_avg = sum(gn1_data)/len(gn1_data)
            gn2_data_avg = sum(gn2_data)/len(gn2_data)

            #some avg value has a zero ratio is 3 other
            if gn1_data_avg * gn2_data_avg != 0:
                gn12_ratio = gn1_data_avg/gn2_data_avg
            else:
                gn12_ratio = 3
            dominant = 0
            if gn12_ratio == 3:
                if gn1_data_avg == 0:
                    if gn2_data_avg == 0:
                        dominant = 'nan'
                    else:
                        dominant = 0
                else:
                    dominant = 1
            elif gn12_ratio > 1:
                dominant = 1
            else:
                dominant = 0
            data_sele = sp_data
            #in dslab is datalabeling list has replaced column names selected and others
            dslab = []
            dsi = list(data_sele.index)
            for dt in dsi:
                if dt != gn1:
                    dslab.append('others')
                else:
                    dslab.append(gn1)
            dsval = list(data_sele)
            new_dslab = copy.deepcopy(dslab)
            fpr1, tpr1, th1 = metrics.roc_curve(dslab, dsval, pos_label = gn1)
            fpr2, tpr2, th2 = metrics.roc_curve(dslab, dsval, pos_label = 'others')
            roc_result = [metrics.auc(fpr1, tpr1),metrics.auc(fpr2, tpr2)]
            roc_result.sort(reverse=True)
            namel = []
            namel.append(gn1)
            namel.append('others')
            namel.sort()
            #AUC permutation
            p_count = 0
            list_AUC_p_val = []
            while p_count != val_it:
                random.shuffle(new_dslab)
                fpr, tpr, th = metrics.roc_curve(new_dslab, dsval, pos_label = gn1)
                auc_value = metrics.auc(fpr, tpr)
                if auc_value < 0.5:
                    auc_value = 1-auc_value
                list_AUC_p_val.append(auc_value)
                p_count += 1
            list_AUC_arr = np.array(list_AUC_p_val, dtype=float)
            AUC_p_value = 1-((stats.percentileofscore(list_AUC_arr, roc_result[0]))/100)
            gns = gn1
            for val_idx in range(3):
                if val_idx == 0:
                    mcm_count_ks[sp_o, n_col] = float(AUC_p_value)
                    if len(mcm_count_g_l) < kscol_num: #kscol_num is number of groups * 2
                        mcm_count_col_l.append('p-value')
                        mcm_count_g_l.append(gns)
                elif val_idx ==1:
                    mcm_count_ks[sp_o, n_col] = float(roc_result[0])
                    if len(mcm_count_g_l) < kscol_num:
                        mcm_count_col_l.append('AUC')
                        mcm_count_g_l.append(gns)
                else:
                    mcm_count_ks[sp_o, n_col] = dominant
                    if len(mcm_count_g_l) < kscol_num:
                        mcm_count_col_l.append('Dominant')
                        mcm_count_g_l.append(gns)
                n_col += 1
    mcm_ks_df = pd.DataFrame(mcm_count_ks, index=mcm_row_l, columns=[mcm_count_g_l,mcm_count_col_l])
    writer3 = ExcelWriter('microbe_AUC_calculation.xlsx') #mcm_count_df writing to Excel file
    mcm_ks_df.to_excel(writer3, 'microbe_AUC', index=True)
    writer3.save()
print('Each microbe AUC value generation is done!')

#data_l have blank list of number of samples
temp_l = []
for num in range(col_num):
    temp_l = []
    data_l.append(temp_l)
    del temp_l

k1n_l = [] #kln_l is K-K cluster name list
rc_data_l = [] #raw coefficient mat
rcol_l = []
col_l = []
out_str = ''
for numb in range(2,26):
    f_name = 'group_coefmat_%d.txt' %(numb) #coefmatrix read
    mf_nam = 'extracted_microb_feat_by_%d.txt' %(numb) #micrbe list
    r_cf_f_nam = 'raw_coefmat_%d.txt' %(numb) #raw coefmatrix read
    coefmat_op = open(f_name,'r')
    micr_l_raw = open(mf_nam, 'r')
    r_cf_f_raw = open(r_cf_f_nam, 'r')

    #microbioal list extract and make a DataFrame
    mi_ls = micr_l_raw.read() #microbe list readfile 
    mi_l = mi_ls.replace('g__','').strip().split('\n') #microbe count list
    temp_l = []
    for k_l in mi_l:
        temp_l.append(','.join(k_l.split(' ')[1:]).strip(','))
    mic_l.append(temp_l)
    del temp_l

    #file read and parsed into DataFrame input format
    coefmat_r = coefmat_op.read()
    rcoefmat_r = r_cf_f_raw.read()
    coefmat_ls = coefmat_r.replace('g__','').strip().split('\n')
    rcoefmat_ls = rcoefmat_r.replace('g__','').strip().split('\n')
    rcol_l = rcoefmat_ls[0].split(' ')
    col_l = coefmat_ls[0].split(' ')
    k_name = 'K-%d' %(numb)
    for dl in rcoefmat_ls[1:]:
        sd = dl.strip().split(' ')
        kk_cl_name = '%s_%s' %(k_name, sd[0])
        k1n_l.append(kk_cl_name)
        rc_data_l.append(sd[1:])
    for dataline in coefmat_ls[1:]:
        data_in = dataline.strip().split(' ')
        for num in range(col_num):
            data_l[num].append(float(data_in[num+1]))
                               
#color palete setting for condition group color coding
cmap = ListedColormap(sns.color_palette("bright",col_num))
cmap_l = cmap.colors.as_hex()
group_index_l = []

for nml in range(col_num):
    temp_l = []
    for ixn in range(len(rcol_l)):
        if rcol_l[ixn] == groupl[nml]:
            temp_l.append(ixn)
    group_index_l.append(temp_l)
    del temp_l
    
kscg_l = []
ksccol_l = []
rceof = np.array(rc_data_l, dtype=float)
rcoef_df = pd.DataFrame(rceof, index=k1n_l, columns=rcol_l)
raw_coefmat = np.zeros((len(k1n_l), kscol_num))

r_df_coefmat = pd.DataFrame()
if os.path.isfile("MMF_AUC_calculation.xlsx"):
    pre_r_df_coefmat = pd.read_excel("MMF_AUC_calculation.xlsx")
    header = list(pre_r_df_coefmat.columns)[1:]
    gl = []
    for x in range(0,len(header),3):
        gnam = header[x]
        gl.append(gnam)
        gl.append(gnam)
        gl.append(gnam)
    labl = list(pre_r_df_coefmat.loc[0,][1:])
    mmf_l = list(pre_r_df_coefmat.iloc[2:,0])
    r_df_coefmat = pre_r_df_coefmat.iloc[2:,1:]
    r_df_coefmat.columns = [gl, labl]
    r_df_coefmat.index = mmf_l
else:
    #data extract and calculate K-S
    for kk_order in range(len(k1n_l)):
        kk_data = rcoef_df.iloc[kk_order,]
        n_col = 0
        for g_idx1 in range(col_num):
            gn1 = groupl[g_idx1]
            gn2 = 'others'
            gn1_data = list(kk_data[gn1])
            gn2_data = list(kk_data[kk_data.index != gn1])
            gn1_data_avg = sum(gn1_data)/len(gn1_data)
            gn2_data_avg = sum(gn2_data)/len(gn2_data)
            #make a ratio
            if gn1_data_avg * gn2_data_avg != 0:
                gn12_ratio = gn1_data_avg/gn2_data_avg
            else:
                gn12_ratio = 3
            dominant = 0
            if gn12_ratio == 3:
                if gn1_data_avg == 0:
                    if gn2_data_avg == 0:
                        dominant = 'nan'
                    else:
                        dominant = 0
                else:
                    dominant = 1
            elif gn12_ratio > 1:
                dominant = 1
            else:
                dominant = 0
            data_sele = kk_data
            dslab = []
            for dts in list(data_sele.index):
                if dts != gn1:
                    dslab.append('others')
                else:
                    dslab.append(gn1)
            dsval = list(data_sele)
            namel = []
            namel.append(gn1)
            namel.append('others')
            namel.sort()
            gns = gn1
            fpr1, tpr1, th1 = metrics.roc_curve(dslab, dsval, pos_label = gn1)
            fpr2, tpr2, th2 = metrics.roc_curve(dslab, dsval, pos_label = gn2)
            roc_result = [metrics.auc(fpr1, tpr1),metrics.auc(fpr2, tpr2)]
            roc_result.sort(reverse=True)
            new_dslab = copy.deepcopy(dslab)
            roc_result.sort(reverse=True)

            #AUC permutation
            p_count = 0
            list_AUC_p_val = []
            while p_count != val_it:
                random.shuffle(new_dslab)
                fpr, tpr, th = metrics.roc_curve(new_dslab, dsval, pos_label = gn1)
                auc_value = metrics.auc(fpr, tpr)
                if auc_value < 0.5:
                    auc_value = 1-auc_value
                list_AUC_p_val.append(auc_value)
                p_count += 1
            list_AUC_arr = np.array(list_AUC_p_val, dtype=float)
            AUC_p_value = 1-((stats.percentileofscore(list_AUC_arr, roc_result[0]))/100)

            for val_idx in range(3):
                if val_idx == 0:
                    raw_coefmat[kk_order, n_col] = float(AUC_p_value)
                    if len(kscg_l) < kscol_num:
                        ksccol_l.append('p-value')
                        kscg_l.append(gns)
                elif val_idx == 1:
                    raw_coefmat[kk_order, n_col] = float(roc_result[0])
                    if len(kscg_l) < kscol_num:
                        ksccol_l.append('AUC')
                        kscg_l.append(gns)
                else:
                    raw_coefmat[kk_order, n_col] = dominant
                    if len(kscg_l) < kscol_num:
                        ksccol_l.append('Dominant')
                        kscg_l.append(gns)
                n_col += 1

    print('Calculation of Coefficient map AUC value done!')
    #Coefficient map AUC caluculation result writing to Excel file
    r_df_coefmat = pd.DataFrame(raw_coefmat, index=k1n_l, columns=[kscg_l, ksccol_l])
    writer_rcoef = ExcelWriter('MMF_AUC_calculation.xlsx')
    r_df_coefmat.to_excel(writer_rcoef, 'MMF_AUC_cal', index=True)
    writer_rcoef.save()

df_input = dict(zip(col_l, data_l))
# generate DataFrame
df_coefmat = pd.DataFrame(df_input, index=k1n_l)
# make max value series
df_max = df_coefmat.max(axis=1)
# normalize all data with max values
norm_df_coefmat = df_coefmat.divide(df_max, axis=0)

# gathering class specific fold changed K
# list for group stat reshare by secies
sp_nm_l = []
sp_nm_grcl_l = []
sp_nm_stat_l = []
def mlg(index_l, lb):
    rl =''''''
    for kk_n in index_l:
        kk_info = kk_n.split('_')
        bk = int(str(kk_info[0]).replace('K-',''))
        sk = int(str(kk_info[1]).replace('K',''))
        group_n = '%s_%s' %(bk, sk)
        mil = mic_l[bk-2][sk-1].strip(',')
        mic_c = mil.split(',')
        for spcies in mic_c:
            if spcies not in sp_nm_l:
                sp_nm_l.append(spcies)
                temp_l = []
                temp_l.append(group_n)
                sp_nm_grcl_l.append(temp_l)
                del temp_l
                temp_l = []
                temp_l.append(lb)
                sp_nm_stat_l.append(temp_l)
                del temp_l
            else:
                mc_idx = sp_nm_l.index(spcies)
                if group_n not in sp_nm_grcl_l[mc_idx]:
                    sp_nm_grcl_l[mc_idx].append(group_n)
                if lb not in sp_nm_stat_l[mc_idx]:
                    sp_nm_stat_l[mc_idx].append(lb)
        wstr = '%s\t%s\n' %(lb, mil)
        rl += wstr
    return rl

pheno_l = []
mic_nl = []
mic_cl = []
final_result = ''''''
sp_nm_l = []
sp_nm_grcl_l = []
sp_nm_stat_l = []
cluster_df = pd.DataFrame()
col_g_name_l = []
cluster_df_col_name = []

for x in range(col_num):
    lox = norm_df_coefmat.iloc[:,x].name
    loy = 'others'
    #rati_colap_cols is mean value of non x labeled columns
    rati_colap_cols = norm_df_coefmat.iloc[:,norm_df_coefmat.columns != lox].mean(axis='columns')
    coef_rati = norm_df_coefmat.iloc[:,x].divide(rati_colap_cols, axis=0)
    col_ratio = '%s_ratio' %(lox)
    col_stat = '%s_stat' %(lox)
    cluster_df[col_ratio] = coef_rati
    high_label = 'H-%s' %(lox)
    low_label = 'L-%s' %(lox)
    true_false = coef_rati > 1
    r1 = mlg(list(norm_df_coefmat.iloc[:,x][coef_rati>1].index), high_label)
    r2 = mlg(list(norm_df_coefmat.iloc[:,x][coef_rati<=1].index), low_label)
    
    #make dominant and hl list
    dominant_list = []
    hl_list = []
    for x in true_false:
        if x == True:
            hl_list.append(high_label)
            dominant_list.append(1)
        else:
            hl_list.append(low_label)
            dominant_list.append(0)
    hl_series = pd.Series(hl_list)
    hl_series.index = true_false.index
    cluster_df[col_stat]= hl_series
    dom_series = pd.Series(dominant_list)
    dom_series.index = true_false.index
    dom_stat = '%s_Dominant' %(lox)
    cluster_df[dom_stat]= dom_series
    cluster_df_col_name.append('ratio')
    col_g_name_l.append(lox)
    cluster_df_col_name.append('status')
    col_g_name_l.append(lox)
    cluster_df_col_name.append('Dominant')
    col_g_name_l.append(lox)
    ned = r_df_coefmat.loc[:,lox].iloc[:,0:2]
    cluster_df = pd.concat([cluster_df, ned], axis=1)
    cluster_df_col_name.append('p-value')
    col_g_name_l.append(lox)
    cluster_df_col_name.append('AUC')
    col_g_name_l.append(lox)
    #final result inclde both reduced and induced species list
    final_result += r1 + r2
 
#only group name in columns
ocngl = list(set(col_g_name_l))
col_label = ['ratio','status','Dominant','p-value','AUC']
print('Generation of ratio between each groups is done!')

#species list by KK cluster
spl_kk_cl = []
for num in range(len(k1n_l)):
    kk_na = k1n_l[num].split('_')
    pk = int(kk_na[0].replace('K-',''))-2
    sk = int(kk_na[1].replace('K',''))-1
    spl_kk_cl.append(mic_l[pk][sk])
mil_se = pd.Series(spl_kk_cl, index=k1n_l)
cluster_df = pd.concat([cluster_df,mil_se], axis=1)
col_names = list(cluster_df.columns)

#cluster_df_col_name.append('Cluster_include')
col_g_name_l.append('common')
cluster_df_col_name.append('Species_list')
cluster_df.columns = [col_g_name_l,cluster_df_col_name]
#cluster_df writing to Excel file
#writer = ExcelWriter('NMF_result_comparison_table.xlsx')
#cluster_df.to_excel(writer, 'NMF_result_comparison', index=True)
#writer.save()

#making pheno count list
phmic_l = final_result.strip().split('\n')
for phmic in phmic_l:
    inf = phmic.split('\t')
    stat = inf[0]
    micl = inf[1].strip(',').split(',')
    if stat not in pheno_l:
        pheno_l.append(stat)
        temp_l = []
        temp_l2 = []
        mic_nl.append(temp_l)
        mic_cl.append(temp_l2)
        st_idx = pheno_l.index(stat)
        for mic_n in micl:
            if mic_n not in mic_nl[st_idx]:
                mic_nl[st_idx].append(mic_n)
                mic_cl[st_idx].append(1)
            else:
                mic_idx = mic_nl[st_idx].index(mic_n)
                mic_cl[st_idx][mic_idx] += 1
    else:
        st_idx = pheno_l.index(stat)
        for mic_n in micl:
            if mic_n not in mic_nl[st_idx]:
                mic_nl[st_idx].append(mic_n)
                mic_cl[st_idx].append(1)
            else:
                mic_idx = mic_nl[st_idx].index(mic_n)
                mic_cl[st_idx][mic_idx] += 1

#writing misummary file
new_mi_df = pd.DataFrame()
origin_mi_index = []
for stat_name in ocngl:
    temp_df = mcm_ks_df[stat_name]
    auc_v_list_temp = str(list(temp_df['AUC']))
    #write AUC value each mi_auc for box plot data generation
    auc_v_list_temp = auc_v_list_temp[1:-1].split(', ')
    write_temp = '\t'.join(auc_v_list_temp)
    write_str_temp = 'mi_%s\t%s\n' %(stat_name, write_temp)
    temp_write = open('mi_%s.data' %(stat_name), 'w')
    temp_write.write(write_str_temp)
    temp_write.close()
    temp_df.columns = ['p-value','AUC', 'Dominant']
    temp_df = temp_df.reset_index()
    temp_df['index'] = temp_df['index'] + ':' + stat_name
    new_mi_df = pd.concat([new_mi_df, temp_df], axis = 0)
new_mi_df = new_mi_df.set_index('index')

#final summary table
new_df = pd.DataFrame()
for stat_name in ocngl:
    # write_temp collecting k-cluster AUC value
    write_temp = str(list(cluster_df[stat_name]['AUC']))
    write_temp = write_temp[1:-1].split(', ')
    temp_auc_list = '\t'.join(write_temp)
    write_auc_string = 'kk_%s\t%s\n' %(stat_name, temp_auc_list)
    nmf_auc_write = open('kk_%s.data' %(stat_name),'w')
    nmf_auc_write.write(write_auc_string)
    nmf_auc_write.close()
    temp_df = cluster_df[[stat_name,'common']]
    temp_df.columns = ['ratio','status','Dominant','p-value','AUC','Species_list']
    temp_df = temp_df.reset_index()
    temp_df['index'] = temp_df['index'] + ':' + stat_name
    new_df = pd.concat([new_df, temp_df], axis = 0)
#writer_new_df_unfiltered = ExcelWriter('nmf_summary_table_others_unfiltered.xlsx')
#new_df.to_excel(writer_new_df_unfiltered, 'final_summary_table', index=True)
#writer_new_df_unfiltered.save()
    
new_df = new_df[new_df['p-value'] < nwp_val]
new_df = new_df.set_index('index')
#writer_new_df = ExcelWriter('nmf_summary_table_others.xlsx')
#new_df.to_excel(writer_new_df, 'final_summary_table', index=True)
#writer_new_df.save()

print('Data summary table generation done!')

# select induced data from single microbe
new_mi_df = new_mi_df[new_mi_df['Dominant'] == 1]
#writer_new_mi_df_whole = ExcelWriter('whole_single_microbe_AUC.xlsx')
#new_mi_df.to_excel(writer_new_mi_df_whole, 'single_microbe_AUC_p_val', index=True)
#writer_new_mi_df_whole.save()

# select induced data from clustered data
new_df = new_df[new_df['Dominant'] == 1]

# box plot generate script
script_body = '''file_list <- Sys.glob(file.path("*.data"))
df <- data.frame()
library(stringr)
for(f in file_list){data <- read.table(f)
'''
level_list = ''''''
for gr_idx in range(len(groupl)):
    if gr_idx == 0:
        level_list += '''"kk_%s","mi_%s"''' %(groupl[gr_idx], groupl[gr_idx])
        command = '''
        if(str_detect(f,"%s")){
        if(str_detect(f,"kk")){
        a <- data.frame(ID = data[,1], AUC = sort(as.numeric(data[,-1]),decreasing = T), BodySite="%s", Method="NMF", index=as.character(1:(ncol(data)-1)), Loc="I")}
        else if(str_detect(f,"mi")){
        a <- data.frame(ID = data[,1], AUC = sort(as.numeric(data[,-1]),decreasing = T), BodySite="%s", Method="Microbe", index=as.character(1:(ncol(data)-1)), Loc="I")}}''' %(groupl[gr_idx],groupl[gr_idx],groupl[gr_idx])
        script_body += command
    else:
        level_list += ''',"kk_%s","mi_%s"''' %(groupl[gr_idx], groupl[gr_idx])
        command ='''
        else if(str_detect(f,"%s")){if( str_detect(f,"kk") ){
        a <- data.frame(ID = data[,1], AUC = sort(as.numeric(data[,-1]),decreasing = T), BodySite="%s", Method="NMF", index=as.character(1:(ncol(data)-1)), Loc="C")}
        else if( str_detect(f,"mi")){
        a <- data.frame(ID = data[,1], AUC = sort(as.numeric(data[,-1]),decreasing = T), BodySite="%s", Method="Microbe", index=as.character(1:(ncol(data)-1)), Loc="C")}}''' %(groupl[gr_idx],groupl[gr_idx],groupl[gr_idx])
        script_body += command
script_tail1 = '''
df <- rbind(df,a)}
library(ggplot2)
df$ID <- factor(df$ID, levels=c('''
script_tail3 ='''))
tiff("NMF_microbe_AUC_comparison.tiff", width = 12, height = 8, units = 'in', res = 350)
ggplot(df, aes(x=Method,  y=AUC)) + 
geom_boxplot(aes(fill=Method), alpha = 0.5) + 
geom_line(aes(group = interaction(index,BodySite)), alpha = 0.5, colour = "darkgrey") +
facet_wrap( ~ BodySite,ncol=5) + theme_classic() + 
theme(strip.text.x = element_text(size=30,face="bold"), axis.text=element_text(size=20,face="bold"),text = element_text(size=30), panel.background = element_rect(colour = "black", size=1.5),legend.position="none")
dev.off()
'''

rscript_f_name = 'AUC_comparison_rscript.R'
rscript_full = script_body + script_tail1 + level_list + script_tail3
rscript_write = open(rscript_f_name ,'w')
rscript_write.write(rscript_full)
rscript_write.close()
rscript_cmd = 'Rscript %s' %(rscript_f_name)
os.system(rscript_cmd)
print('Drawing boxplot for AUC comparison done!')

# extract auc mi
convert_to_str = []
for one_c in list(new_mi_df['AUC']):
    str_one = '%f' %one_c
    convert_to_str.append(str_one)
mi_auc_l = '\t'.join(convert_to_str)
# single_mi_ filter
new_mi_df_pv_cut = new_mi_df[new_mi_df['p-value'] < nwp_val]
#writer_new_df = ExcelWriter('single_microbe_AUC_pval_table.xlsx')
#new_mi_df_pv_cut.to_excel(writer_new_df, 'single_microbe', index=True)
#writer_new_df.save()

    
#extract auc nmf
convert_to_str = []
for one_c in list(new_df['AUC']):
    str_one = '%f' %one_c
    convert_to_str.append(str_one)
nmf_auc_l = '\t'.join(convert_to_str)
new_df.to_csv('NMF_result_summary.txt', header=True, index=True, sep='\t')

#AUC data writing
auc_list_string = "microb\t%s\nNMF\t%s\n" %(mi_auc_l, nmf_auc_l)
handle_auc = open('AUC_value_list.txt','w')
handle_auc.write(auc_list_string)
handle_auc.close()

#make a list for highest auc each species under condition
status = list(set(new_df['status']))
microbe_list = []
microbe_con = []
microbe_auc = []
for condi in status:
    if condi != 'nan':
        true_data = new_df[new_df['status'] == condi]
        td_idx = true_data.index
        for tdid in td_idx:
            sll = true_data.loc[tdid]['Species_list']
            auc_v = true_data.loc[tdid]['AUC']
            spl = sll.split(',')
            for sp in spl:
                if sp not in microbe_list:
                    temp_l = []
                    temp_l2 = []
                    microbe_list.append(sp)
                    temp_l.append(condi)
                    temp_l2.append(auc_v)
                    microbe_con.append(temp_l)
                    microbe_auc.append(temp_l2)
                    del temp_l
                    del temp_l2
                else:
                    if condi not in microbe_con[microbe_list.index(sp)]:
                        microbe_con[microbe_list.index(sp)].append(condi)
                        microbe_auc[microbe_list.index(sp)].append(auc_v)
                    else:
                        if microbe_auc[microbe_con[microbe_list.index(sp)].index(condi)][0] < auc_v:
                            microbe_auc[microbe_con[microbe_list.index(sp)].index(condi)][0] = auc_v                
pv_tag = str(nwp_val).replace('0','').replace('.','')
print('writing AUC value is done!')

# make new dataframe from
kcluster_stat_l = new_df.index
k_clster_species_list_df = pd.DataFrame()
for s_index in kcluster_stat_l:
    kcluster_info_df = new_df[new_df.index == s_index]
    if str(kcluster_info_df['status'].to_numpy()[0]).startswith('H-'):
        k_clster_species_list_df = pd.concat([k_clster_species_list_df, kcluster_info_df])

# sorting dictionary list ordered by values
def dict_sort(dict_str):
    each_sp = dict_str.strip().split(',')
    spl = []
    values = []
    for sp in each_sp:
        data = sp.strip("'").split(':')
        spl.append(data[0])
        fval = float(data[1])
        values.append(fval)
    spv_dict = dict(zip(spl, values))
    values_set = list(set(values))
    values_set.sort(reverse=True)
    newl = []
    for vl in values_set:
        for spn, val in spv_dict.items():
            if val == vl:
                out = spn+':'+str(val)
                newl.append(out)
    outstr = ','.join(newl)
    return outstr

# MMFs_microbe_contiribution value ordering
handle_cont = open('MMF_microbe_list.txt','r')
rcont = handle_cont.read()
rc_ls = rcont.strip().replace('g__','').replace("NA : NA ','\n","").replace(' ','').replace(",'\n",",'").replace(",'@@","@@").replace('\n\n','\n').split('\n')
han_out = open('MMF_microbe_list_ordered.txt','w')
header = 'index\tmicrobe_list\n'
han_out.write(header)
kk_name_list = []
kk_sp_contribution_l = []

for line in rc_ls:
    spl = line.split('@@')
    if spl[0] != '':
        cl_i = spl[1]
        cl_num_i = cl_i.split('-')
        cluster_name = 'K-%s_K%s' %(cl_num_i[0], cl_num_i[1])
        ordered_info = dict_sort(spl[0])
        kk_name_list.append(cluster_name)
        kk_sp_contribution_l.append(ordered_info)
        out_str = '%s\t%s\n' %(cluster_name, ordered_info)
        han_out.write(out_str)
han_out.close()

sp_cont_dict = {"cluster":kk_name_list,"mi_cont":kk_sp_contribution_l}
sp_contribution_df = pd.DataFrame(sp_cont_dict)
sp_contribution_df = sp_contribution_df.set_index('cluster')
print('MMF microbe file writing is done!')

## co-occurrence networks cluster relation btw kk clustering
data_list = []
clustering_out = open('MMF_high_AUC_and_list.data','w')
cl_header = 'index\thighest_microbe_AUC\tmicrobe_name\n'
clustering_out.write(cl_header)

## co-occurrence networks cluster relation btw kk clustering
## gathering highest single microbe AUC
## list for contribution and auc correlation
kk_nam_l = []
auc_val_l = []
contribution_l = []
def index_shift(input_query, order_num):
    s_or_f = 'no'
    try:
        new_mi_df_pv_cut[new_mi_df_pv_cut.index==input_query].AUC.to_numpy()[0]
        s_or_f = 'get'
    except IndexError:
        order_num += 1
    return s_or_f, order_num
print('mapping of AUC of each contribution microbes')

for k8_index in k_clster_species_list_df.index:
    one_row = k_clster_species_list_df[k_clster_species_list_df.index == k8_index]
    kk_nam = k8_index.split(':')[0]
    stat = k8_index.split(':')[1]
    # sp_contribution_df is dataframe which has species and contribution values
    sp_cont_l = sp_contribution_df[sp_contribution_df.index == kk_nam]['mi_cont'].to_numpy()[0].split(',')
    # if want to gather intersection with cluster species list
    new_sp_l = []
    for sp_n_cont in sp_cont_l:
        sp_n = sp_n_cont.split(':')[0]
        sp_cont = float(sp_n_cont.split(':')[1])
        query_input = '%s:%s' %(sp_n, stat)
        sp_auc_v = new_mi_df_pv_cut[new_mi_df_pv_cut.index==query_input].AUC.to_numpy()
        if sp_cont > 0:
            if len(sp_auc_v) == 1:
                out_info = "%s:%f" %(sp_n, float(sp_auc_v[0]))
                new_sp_l.append(out_info)
                print_str = '%s\t%s\t%s\t%s' %(k8_index, sp_n, sp_auc_v[0], sp_cont)
                kk_nam_l.append(k8_index+':'+sp_n)
                auc_val_l.append(sp_auc_v[0])
                contribution_l.append(sp_cont)
            else:
                out_info = "%s:0" %(sp_n)
                new_sp_l.append(out_info)
    # sort by AUC value
    if len(new_sp_l) != 0:
        inter_sec = dict_sort(','.join(new_sp_l))
    else:
        inter_sec = 'none'
    
    # highest microbe AUC selection only for statistically available
    get_it = 'no'
    inter_sec_n = ''
    highest_auc = 0.5
    species_order = 0
    while get_it != 'get':
        if len(inter_sec.split(',')) == 1:
            sp_name = inter_sec[0].split(':')[0]
            query_input_2 = '%s:%s' %(sp_name, stat)
            get_it, species_order = index_shift(query_input_2, species_order)
            #this code for avoid eternal loop when the species length is 1 
            if species_order > len(inter_sec.split(',')):
                get_it = 'get'
        elif species_order <= len(inter_sec.split(','))-1:
            sp_name = inter_sec.split(',')[species_order].split(':')[0]
            query_input_2 = '%s:%s' %(sp_name, stat)
            get_it, species_order = index_shift(query_input_2, species_order)
            if get_it == 'get':
                query_input_2 = '%s:%s' %(sp_name, stat)
                highest_auc = new_mi_df_pv_cut[new_mi_df_pv_cut.index==query_input_2].AUC.to_numpy()[0]
        else:
            get_it = 'get'
    out_strings = '%s\t%s\t%s\n' %(k8_index, highest_auc, inter_sec)
    clustering_out.write(out_strings)
clustering_out.close()
print('data gen for correlation contribution and AUC')

# datagen for correlation between contribution and auc
data_l2 = [kk_nam_l, auc_val_l, contribution_l]
newData = np.transpose(data_l2)
con_df = pd.DataFrame(newData, columns=['MMF_microbe','AUC','Contribution'])
con_df_out = con_df.sort_values('AUC', ascending=False)
con_df_out = con_df_out.set_index('MMF_microbe')

# code for making cont_auc_txt file
df1 = pd.read_table('NMF_result_summary.txt', sep='\t')
df2 = pd.read_table('MMF_microbe_list_ordered.txt', sep='\t')
df3 = pd.read_table('MMF_high_AUC_and_list.data', sep='\t')

#removing NA data in status column from summary table 
df1 = df1.dropna(subset=['status'])

#gathering only Highly induced data
df1 = df1[df1['status'].str.startswith('H-')]

#remove duplicate
df3 = df3.drop_duplicates(subset=['index'], keep='first')

#extract only index and AUC
df1 = df1.loc[:,['index','AUC']]
results = df1.merge(df3, on= 'index')
results['cluster'] = results['index'].str.split(':').str.get(0)
results['class'] = results['index'].str.split(':').str.get(1)
results = results.sort_values('class')
cont_auc = pd.merge(results, df2, how='left', left_on=['cluster'], right_on=['index'])
cont_auc = cont_auc.drop_duplicates(subset=['index_x'], keep='first')

#extract only index and AUC
final_cont_auc = cont_auc.loc[:,['cluster','class','AUC','highest_microbe_AUC','microbe_name','microbe_list']]
final_cont_auc = final_cont_auc.sort_values(['class','cluster'])

#reset index remove previous index
final_cont_auc = final_cont_auc.reset_index(drop=True)
final_cont_auc.columns = ['Cluster', 'Group category','AUC','Highest microbe AUC','Microbe list and AUC','Microbe list and contribution']
writer_cont_auc = ExcelWriter('MMF_microbe_AUC_data_summary.xlsx')
final_cont_auc.to_excel(writer_cont_auc, 'MMF_AUC_and_microbe', index=True)
writer_cont_auc.save()

# ECDF
def ecdf(data_list, cut_off):
    raw_data = np.array(data_list)
    # create a sorted series of unique data
    cdfx = np.sort(np.unique(data_list))
    # x-data for the ECDF: evenly spaced sequence of the uniques
    x_values = np.linspace(start=min(cdfx),stop=max(cdfx),num=len(cdfx))
    # size of the x_values
    size_data = raw_data.size
    cut_value = 0
    for i in x_values:
        # all the values in raw data less than the ith value in x_values
        temp = raw_data[raw_data <= i]
        # fraction of that value with respect to the size of the x_values
        value = temp.size / size_data
        if value >= cut_off:
            cut_value = i
            break
    return cut_value

#gathering features which have difference between kk AUC and highest species
## define symmf
def class_label(x):
    if x>0: 
        return 'SYMMF'
    else: 
        return 'MMF'
def color_label(x):
    if x>0: 
        return 'Red'
    else: 
        return 'Blue'
def condition_color_label(x):
    return cmap_l[groupl.index(x)]

final_cont_auc['Diff'] = final_cont_auc['AUC'] - final_cont_auc['Highest microbe AUC']        
final_cont_auc['Class'] = final_cont_auc['Diff'].apply(class_label)
final_cont_auc['Color'] = final_cont_auc['Diff'].apply(color_label)
final_cont_auc['Group_color'] = final_cont_auc['Group category'].apply(condition_color_label)

ecdf_cut_off_AUC_l = []
condition_cutoff_value_l = []
final_cont_auc_cutoff = pd.DataFrame()

for part_n in groupl:
    temp_p_df = final_cont_auc[final_cont_auc['Group category']==part_n]
    cutv = ecdf(list(temp_p_df['AUC']), ecdf_cut_v)
    ecdf_cut_off_AUC_l.append(cutv)
    out_str = '%s : %s' %(part_n, cutv)
    condition_cutoff_value_l.append(out_str)
    print(out_str)
    temp_p_cut_df = temp_p_df[temp_p_df['AUC'] >= cutv]
    final_cont_auc_cutoff = pd.concat([final_cont_auc_cutoff, temp_p_cut_df], axis = 0)
    
#High AUC SYMMF##
symmf_cutoff = final_cont_auc_cutoff.copy()
symmf_cutoff = symmf_cutoff[symmf_cutoff['Diff'] > dif_auc_val].sort_values(by='Diff', axis=0, ascending=False)
symmf_cutoff.to_excel("SYMMF_summary_table.xlsx", sheet_name="SYMMFs_data")
minmum_auc = min(ecdf_cut_off_AUC_l)
print('Minimum auc : ',minmum_auc)
    
### network_generation
#cluster_sp_list = open('SYMMF_microbe_list.txt','w')
cluster_sp = ''''''
heatmap_sp_list = []
high_auc_kcluster_r_coefmat = pd.DataFrame()

df_for_heatmap = pd.DataFrame()
symmf_groupl = list(set(list(symmf_cutoff['Group category'])))
for one_c in symmf_groupl:
    temp_l_for_group = []
    temp_l_for_group.append(one_c)
    temp_l_for_val = []
    temp_l_for_sp = []
    stat_df = symmf_cutoff[symmf_cutoff['Group category'] == one_c].sort_values('AUC', ascending=False)
    if len(stat_df.index) == 0:
        print_str = '%s has no network data' %(one_c)
        
    #extract top AUC kkcluster list
    top1_diff_series = stat_df.iloc[0]
    top2_diff_l = list(stat_df['Cluster'])
	
	#gathering top_auc species_l
    spe_diff_c_l = top1_diff_series['Microbe list and AUC'].split(',')
    sp_l_diff = []
    for sp_diff_n in spe_diff_c_l:
        sp_l_diff.append(str(sp_diff_n.split(':')[0]))
    for kk_name in top2_diff_l:
        one_df = rcoef_df.loc[kk_name,]
        if len(high_auc_kcluster_r_coefmat.index) == 0:
            high_auc_kcluster_r_coefmat = one_df
        else:
            high_auc_kcluster_r_coefmat = pd.concat([high_auc_kcluster_r_coefmat, one_df], axis=1)
            
    #network file auc cut off
    sp_pair_l_ol = []
    sp_pair_l_olv = []
    sp_pair_l_olv_stack = []
    total_species_l = []
    total_conv_auc_a = []
    total_category_l = []
    for index, tdl in stat_df.iterrows():
        class_name = tdl[1]
        auc_v = float(tdl[2])
        sp_con_l = tdl[4]
        sp_name_l = []
        sp_adj_vl = []
        temp_set = {}
        if one_c == class_name:
            for one_sp in sp_con_l.strip().split(','):
                info = one_sp.split(':')
                #filtering microbe species contribution value over 1% 
                if float(info[1]) > microbe_contribution_cuttoff:
                    sp_name_l.append(info[0])
                    adj_v = auc_v * float(info[1]) * 100
                    sp_adj_vl.append(adj_v)
                    if info[0] not in total_species_l:
                        total_species_l.append(info[0])
                        total_conv_auc_a.append(adj_v)
                        temp_l = []
                        temp_l.append(class_name)
                        total_category_l.append(temp_l)
                    else:
                        total_conv_auc_a[total_species_l.index(info[0])] += adj_v
                        if class_name not in total_category_l[total_species_l.index(info[0])]:
                            total_category_l[total_species_l.index(info[0])].append(class_name)
                            
            for x in range(len(sp_name_l)):
                for y in range(len(sp_name_l)):
                    if x > y:
                        sp_pair = []
                        sp_pair.append(sp_name_l[x])
                        sp_pair.append(sp_name_l[y])
                        sp_pair.sort()
                        spp = ' i '.join(sp_pair)                      
                        if spp not in sp_pair_l_ol:
                            sp_pair_l_ol.append(spp)
                            sp_pair_l_olv.append(auc_v)
                            sp_pair_l_olv_stack.append(auc_v)
                        else:
                            sp_pair_l_olv_stack[sp_pair_l_ol.index(spp)] += auc_v
                            if sp_pair_l_olv[sp_pair_l_ol.index(spp)] < auc_v:
                                sp_pair_l_olv[sp_pair_l_ol.index(spp)] = auc_v
                                
    sp_polv_cut = ecdf(sp_pair_l_olv, network_cut_off)
    sp_pair_l_ol_fil = []
    sp_pair_l_olv_fil = []
    sp_pair_l_olv_stack_fil = []
    for idx in range(len(sp_pair_l_ol)):
        pn = sp_pair_l_ol[idx]
        val = sp_pair_l_olv[idx]
        val2 = sp_pair_l_olv_stack[idx]
        if val >= sp_polv_cut:
            sp_pair_l_ol_fil.append(pn)
            sp_pair_l_olv_fil.append(val)
            sp_pair_l_olv_stack_fil.append(val2)
    for spp in sp_pair_l_ol_fil:
        spl = spp.split(' i ')
        for sn in spl:
            if sn not in heatmap_sp_list:
                heatmap_sp_list.append(sn)
    f_name_1 = '%s_pairs_list.sif' %(one_c)
    #f_name_5 = '%s_pairs_edge_value_AUC_thk.attrs' %(one_c)
    #f_name_3 = '%s_microbe_node_value.attrs' %(one_c)
    #f_name_4 = '%s_microbe_cate.attrs' %(one_c)
    #f_name_6 = '%s_highest_auc_microbe.attrs' %(one_c)
    handle_pair_out = open(f_name_1,'w')
    #handle_pair_v_stack_out = open(f_name_5,'w')
    #handle_pair_v_stack_out.write('interaction_auc_stack_value\n')
    for num in range(len(sp_pair_l_ol_fil)):
        pair = sp_pair_l_ol_fil[num]
        pvalue = sp_pair_l_olv_fil[num]
        pvalue2 = sp_pair_l_olv_stack_fil[num]
        outstrings2 = '%s\n' %(pair)
        pair = sp_pair_l_ol_fil[num].replace(" i "," (i) ")
        outstrings1 = '%s = %0.5f\n' %(pair, pvalue)
        outstrings3 = '%s = %0.5f\n' %(pair, pvalue2)
        handle_pair_out.write(outstrings2)
        #handle_pair_v_stack_out.write(outstrings3)
    handle_pair_out.close()
    #handle_pair_v_stack_out.close()
    
    #writing species related inforation
    #handle_pair_v_highest_auc = open(f_name_6,'w')
    #handle_one = open(f_name_3,'w')
    #handle_pair_cat_out = open(f_name_4,'w')
    #handle_pair_v_highest_auc.write('highest_auc\n')
    #handle_pair_cat_out.write('microbe_cat\n')
    #handle_one.write('microbe_n_value\n')
    for num2 in range(len(total_species_l)):
        species_n = total_species_l[num2]
        if species_n in sp_l_diff:
            string_high = '%s = highest_auc\n' %(species_n)
            #handle_pair_v_highest_auc.write(string_high)
        species_node_val = total_conv_auc_a[num2]
        string_out = '%s = %0.5f\n' %(species_n, species_node_val)
        ##line for data frame generation
        temp_l_for_val.append(species_node_val)
        temp_l_for_sp.append(species_n)
        cat_l = total_category_l[num2]
        cat_l.sort()
        cat_str = ','.join(cat_l)
        string_out2 = '%s = "%s"\n' %(species_n, cat_str)
        #handle_pair_cat_out.write(string_out2)
        #handle_one.write(string_out)
    #handle_pair_v_highest_auc.close()
    #handle_one.close()
    #handle_pair_cat_out.close()
        
    # heatmap
    templ_val = []
    templ_val.append(temp_l_for_val)
    df_data = dict(zip(temp_l_for_group,templ_val))
    temp_df = pd.DataFrame(data=df_data, columns=temp_l_for_group, index=temp_l_for_sp)
    df_for_heatmap = df_for_heatmap.join(temp_df, how='outer')

    # searching  clustered microbelist
    mic_p_list = []
    for idx in range(len(sp_pair_l_ol_fil)):
        pair_l = sp_pair_l_ol_fil[idx].split(' i ')
        mic_p_list.append(set(pair_l))
    
    #gathering clusters
    for x in range(len(mic_p_list)):
        for y in range(len(mic_p_list)):
            if x > y:
                ins_c = mic_p_list[x].intersection(mic_p_list[y])
                if len(ins_c) != 0:
                    mic_p_list[x] = mic_p_list[x].union(mic_p_list[y])
                    mic_p_list[y] = {}

    #remove null sets
    cluster_l = []
    for cls in mic_p_list:
        if len(cls) != 0:
            cluster_l.append(cls)
    for cls in cluster_l:
        mi_single_AUC_list = []
        cls_p_list = []
        cls_p_v_list = []
        for one_pair in sp_pair_l_ol_fil:
            pair_i = one_pair.split(' i ')
            if pair_i[0] in cls:
                if pair_i[1] in cls:
                    cls_p_list.append(one_pair)
                    cls_p_v_list.append(sp_pair_l_olv_fil[sp_pair_l_ol_fil.index(one_pair)])
        for microbe in list(cls):
            auc_val = mcm_ks_df[one_c]['AUC'][microbe]
            mi_single_AUC_list.append(auc_val)
        cluster_sp += 'label : %s\nmicrobe_list : %s\nmax_edge_AUC : %0.5f\nAverage_edge_AUC : %0.5f\nmax_single_microbe_AUC : %0.3f\nAverage_microbe_AUC : %0.3f\n\n' %(one_c, cls, max(cls_p_v_list), sum(cls_p_v_list)/len(cls_p_v_list), max(mi_single_AUC_list), sum(mi_single_AUC_list)/len(mi_single_AUC_list))
#cluster_sp_list.write(cluster_sp)
#cluster_sp_list.close()

#heamap drawing
plt.clf()
plt.figure(figsize=(8,23))
hmp_no = sns.set()
hmp_no = sns.heatmap(df_for_heatmap, xticklabels=1, yticklabels='auto')

# options recording
optin_out = '''iteration count of AUC calculation = %d
SYMMF-net edge p-value cut-off = %f
minmum AUC = %f
AUC p-value = %f
''' %(val_it, network_cut_off, minmum_auc, nwp_val)

min_auc_str = '%0.5f' %(minmum_auc)
min_auc_str = min_auc_str.replace('.','')
c_pv_str = '%0.5f' %(nwp_val)
cut_pval_str = c_pv_str.replace('.','')

option_out = open('options_record.txt','w')
option_out.write(optin_out)
option_out.close()
con_df_out = con_df.sort_values('AUC', ascending=False)
con_df_out = con_df_out.set_index('MMF_microbe')

data_for_plot_all_mmf = final_cont_auc[final_cont_auc['Highest microbe AUC'] != 0]
data_for_plot_all_mmf.columns = ["Cluster","Group category","AUC of MMF","Highest microbe AUC","Microbe list and AUC","Microbe list and contribution","Diff","Class","Color","Group color"]

max_x = float(data_for_plot_all_mmf['AUC of MMF'].max()) + 0.01
min_x = float(data_for_plot_all_mmf['AUC of MMF'].min()) - 0.01
max_y = float(data_for_plot_all_mmf['Highest microbe AUC'].max()) + 0.01
min_y = float(data_for_plot_all_mmf['Highest microbe AUC'].min()) - 0.01
sns.set(style='ticks')

g = sns.lmplot(x="AUC of MMF", y="Highest microbe AUC", 
               data=data_for_plot_all_mmf,
               palette=dict(zip(groupl,cmap_l)),
               fit_reg=False,
               truncate=True,
               hue='Group category',
               legend=False,
               scatter_kws = {"s":20,"alpha":0.5})
g.set(xlim=(min_x, max_x), ylim=(min_y, max_y))
g.fig.set_figwidth(7)
g.fig.set_figheight(7)

# Move the legend to an empty part of the plot
plt.plot([0,1],'-k', alpha = 0.3)
plt.legend(loc='lower right')
g.savefig('all_SYMMF_MMF_scatter_plot.jpg')

data_for_plot = final_cont_auc_cutoff[final_cont_auc_cutoff['Highest microbe AUC']!=0]
data_for_plot.columns = ["Cluster","Group category","AUC of MMF","Highest microbe AUC","Microbe list and AUC","Microbe list and contribution","Diff","Class","Color","Group color"]

max_x = float(data_for_plot['AUC of MMF'].max()) + 0.01
min_x = float(data_for_plot['AUC of MMF'].min()) - 0.01
max_y = float(data_for_plot['Highest microbe AUC'].max()) + 0.01
min_y = float(data_for_plot['Highest microbe AUC'].min()) - 0.01
sns.set(style='ticks')

g = sns.lmplot(x="AUC of MMF", y="Highest microbe AUC", 
               data=data_for_plot,
               palette=dict(zip(groupl,cmap_l)),
               fit_reg=False,
               truncate=True,
               hue='Group category',
               legend=False,
               scatter_kws = {"s":20,"alpha":0.5})
g.set(xlim=(min_x, max_x), ylim=(min_y, max_y))
g.fig.set_figwidth(7)
g.fig.set_figheight(7)

# Move the legend to an empty part of the plot
plt.plot([0,1],'-k', alpha = 0.3)
plt.legend(loc='lower right')
smsplot_fig = 'SYMMF_MMF_scatter_plot_ECDF%s.jpg' %(ecdf_cut_str)
g.savefig(smsplot_fig)

g2 = sns.jointplot(x='AUC of MMF', 
                    y='Highest microbe AUC', 
                    kind='hex',
                    data=data_for_plot,
                   marginal_kws=dict(bins=30, rug=True),
                  # annot_kws=dict(stat='r')
		  )

g2.ax_joint.set_xlim(min_x, max_x)
g2.ax_joint.set_ylim(min_y, max_y)
g2.ax_joint.plot([0,1.3], [0, 1.3], '-k', alpha = 0.3)
mash_fig = "MMF_AUC_save_hex_ECDF%s.png" %(ecdf_cut_str)
g2.savefig(mash_fig)

g2 = sns.jointplot(x='AUC of MMF', 
                    y='Highest microbe AUC', 
                    data=data_for_plot,
                   kind='kde',
                   space=0,
                    color="k",)
g2.ax_joint.set_xlim(min_x, max_x)
g2.ax_joint.set_ylim(min_y, max_y)
g2.ax_joint.plot([0,1.3], [0, 1.3], '-k', alpha = 0.3)
mask_fig = "MMF_AUC_save_kde_ECDF%s.png" %(ecdf_cut_str)
g2.savefig(mask_fig)

# generate table for plot
data_for_plot = data_for_plot.sort_values(by='AUC of MMF', ascending=False).reset_index()
table_1 = data_for_plot[data_for_plot['AUC of MMF']>minmum_auc]
table_2 = table_1[table_1['Diff']> 0]

# generate symmf bmap
class_name = list(set(table_2['Group category'].values))
k_list = []
feat_list = []
for cls_c in class_name:
    class_k_l = []
    class_f_l = []
    class_t = table_2[table_2['Group category']==cls_c]
    for k_name in class_t['Cluster']:
        name_info = k_name.split('_')
        k_name = name_info[0].replace('K-','')
        f_name = name_info[1]
        if k_name not in k_list:
            class_k_l.append(k_name)
            temp_l = []
            temp_l.append(f_name)
            class_f_l.append(temp_l)
        else:
            kl_idx = class_k_l.index(k_name)
            class_f_l[kl_idx].append(f_name)
    k_list.append(class_k_l)
    feat_list.append(class_f_l)

## gen symmf basismap
bmap_col_name = []
total_bmap_table = pd.DataFrame()
count = 0
for class_idx in range(len(class_name)):
    class_k_l = k_list[class_idx]
    class_f_l = feat_list[class_idx]
    for idx in range(len(class_k_l)):
        k_nam = class_k_l[idx]
        f_nam = class_f_l[idx]
        read_f_name = 'raw_bmap_mat_%s.txt' %(k_nam)
        read_table = pd.read_table(read_f_name, sep=' ')
        for k_idx in f_nam:
            count += 1
            full_nam_of_f = 'K%s_%s' %(k_nam, k_idx)
            bmap_col_name.append(full_nam_of_f)
            if count == 1:
                total_bmap_table = read_table[k_idx]
            else:
                total_bmap_table = pd.concat([total_bmap_table, read_table[k_idx]], axis=1)

total_bmap_table.columns = bmap_col_name
total_bmap_table_name = 'SYMMF_feature_bmap_selected.txt'
total_bmap_table.to_csv(total_bmap_table_name, sep='\t')

## generation of node value with condtribution value in SYMMF
total_bmap_raw_sum = total_bmap_table.sum(axis=1)
sum_t = total_bmap_raw_sum * 10
selected_bmap_node_contribution_name = 'whole_network_nodes_size_p%s_auc%s_ECDF%s.attrs' %(cut_pval_str, min_auc_str, ecdf_cut_str)
write_bmap_node = open(selected_bmap_node_contribution_name, 'w')
write_bmap_node.write('node_size\n')
for idx_num in range(len(sum_t.index)):
    node_name = sum_t.index[idx_num].replace('g__','')
    node_value = "%0.5f" %(float(sum_t[idx_num]))
    w_string = '%s = %s\n' %(node_name, node_value)
    write_bmap_node.write(w_string)
write_bmap_node.close()

## gen symmf coefmat
coef_col_name_l = []
total_coef_table = pd.DataFrame()
count = 0
handle_for_col = open('raw_coefmat_2.txt','r')
col_line = handle_for_col.readline()
class_info = col_line.strip().split(' ')
for idx in range(len(k_list)):
    k_nam_l = k_list[idx]
    f_nam_l = feat_list[idx]
    for idx in range(len(k_nam_l)):
        k_nam = k_nam_l[idx]
        f_nam = f_nam_l[idx]
        read_f_name = 'raw_coefmat_%s.txt' %(k_nam)
        read_table = pd.read_table(read_f_name, sep=' ')
        t_table = read_table.T
        for k_idx in f_nam:
            count += 1
            full_nam_of_f = 'K-%s_%s' %(k_nam, k_idx)
            coef_col_name_l.append(full_nam_of_f)
            if count == 1:
                total_coef_table = t_table[k_idx]
            else:
                total_coef_table = pd.concat([total_coef_table, t_table[k_idx]], axis=1)
total_coef_table.index = class_info
total_coef_table.columns = coef_col_name_l
total_coef_table_name = 'SYMMF_feature_coef_selected.txt'
total_coef_table.to_csv(total_coef_table_name, sep='\t')

orderby_auc_l = list(final_cont_auc_cutoff.sort_values(['AUC'],ascending=False)['Cluster'])

# total coef group average table
total_coef_table_idx = total_coef_table.reset_index()
total_coef_group_table = total_coef_table_idx.groupby(['index'], as_index = False).mean()
total_coef_group_table = total_coef_group_table.set_index("index")
row_g_colors = []
for class_name in total_coef_group_table.index:
    for condition_1 in groupl:
        if condition_1 == class_name:
            idx_clr = cmap_l[groupl.index(condition_1)]
            row_g_colors.append(idx_clr)
colname_idx = list(total_coef_table.columns)

new_col_l = []
for x in orderby_auc_l:
    if x in colname_idx:
        new_col_l.append(x)

total_coef_group_table_reord = total_coef_group_table.loc[:,new_col_l]
hmap_mean_coef_reord = sns.clustermap(total_coef_group_table_reord, 
                           row_colors= row_g_colors, 
                           row_cluster = False,
                           col_cluster = False,
                           cmap="YlGnBu",
                           figsize=(20,2))
ax = hmap_mean_coef_reord.ax_heatmap
scmre_fig = "SYMMFs_coefmap_mean_reord_save_ECDF%s.png" %(ecdf_cut_str)
hmap_mean_coef_reord.savefig(scmre_fig)

# drawing heatmap of concatenated coefmap
row_colors = []
for class_name in total_coef_table.index:
    for condition_1 in groupl:
        if condition_1 == class_name:
            idx_clr = cmap_l[groupl.index(condition_1)]
            row_colors.append(idx_clr)
hmap_coef = sns.clustermap(total_coef_table.loc[:,new_col_l], 
                           row_colors= row_colors, 
                           row_cluster = False,
                           col_cluster = False,
                           cmap="YlGnBu",
                           figsize=(20,40))
ax = hmap_coef.ax_heatmap
h_coef_fig = "SYMMFs_coefmap_save_ECDF%s.png" %(ecdf_cut_str)
hmap_coef.savefig(h_coef_fig)

# total basis table W_symmf
norm_total_bmap_table = preprocessing.normalize(total_bmap_table)
norm_total_bmap_table_pd = pd.DataFrame(norm_total_bmap_table)
norm_total_bmap_table_pd.index = total_bmap_table.index
col_name_list_before_modified = list(total_bmap_table.columns)
modified_string = '\t'.join(col_name_list_before_modified)
modified_string_final = modified_string.replace('K','K-').replace('_K-','_K')
norm_total_bmap_table_pd.columns = modified_string_final.split('\t')
idx_col_name_l = norm_total_bmap_table_pd.columns
norm_total_bmap_table_pd_reord = norm_total_bmap_table_pd.loc[:,new_col_l]

bmap_basis = sns.clustermap(norm_total_bmap_table_pd_reord,  
                           row_cluster = True,
                           col_cluster = False,
                           cmap="YlGnBu",
                           figsize=(15,25))

w_basis_fig = "SYMMFs_basismap_reord_save_ECDF%s.png" %(ecdf_cut_str)
bmap_basis.savefig(w_basis_fig)

#calculation of correlation value
file_list = os.listdir("./")
cor_gen_R_script = '''options<-commandArgs(trailingOnly=T)
condition_name <- options[1]
r_file_name <- paste("SYMMF_feature_bmap_selected.txt", sep="")
sp_info_f_name <- paste("SYMMF_feature_selected_microbe_list.txt", sep="")
count_mat_name <- paste("SYMMF_feature_selected_microbe_count_mat.txt", sep="")
cor_file_name <- paste("cor_SYMMF_feature_bmap_selected.txt", sep="")

micro_abu <- read.table(r_file_name, stringsAsFactors=FALSE)
species_info <- rownames(micro_abu)
disease_class <- colnames(micro_abu)
micro_abu<-as.matrix(sapply(micro_abu, as.numeric))
rownames(micro_abu) <- species_info
colnames(micro_abu) <- disease_class

#row_rem_idx <- c()
#for(i in c(1:nrow(micro_abu)) ){	# for each microbiome
#	if( length(which(micro_abu[i,]>0)) < length(micro_abu[i,])*0.15 ){	# at least 15% species have contribution values	
#	row_rem_idx <- c(row_rem_idx,i)
#	}
#}
#if(length(row_rem_idx)>1){ 
#	micro_abu <- micro_abu[-row_rem_idx,]
#}

##############################
#remove all row data sum is 0
##############################
all_0_row_id <- which(apply(micro_abu, 1, sum)==0)
if (length(all_0_row_id) != 0){
	micro_abu <- micro_abu[-all_0_row_id,]
}
write.table(species_info,file=sp_info_f_name,quote=F,row.names=F,col.names=F)
write.table(micro_abu, file=count_mat_name, quote=F, col.names=T,row.names=T)
disease_class <- colnames(micro_abu)
species_name <- rownames(micro_abu)

nrow_micro = nrow(micro_abu)
for(i in 1:nrow_micro){
	mi_name <- rownames(micro_abu)[i]
	for(j in i:nrow_micro){
	if(j == i) next
	mj_name <- rownames(micro_abu)[j]
	try({cor.test(micro_abu[i,],micro_abu[j,])
	pval <- round( (cor.test(micro_abu[i,],micro_abu[j,]))$p.value , 7)
    cor_val <- round( (cor.test(micro_abu[i,],micro_abu[j,]))$estimate , 4)
	if(pval < 0.001){
	cat(paste(i,mi_name,j,mj_name,pval,cor_val,"\n"), file=cor_file_name, append=TRUE)}}
	, silent = TRUE)
	}} 
'''
handle_cor_rscript = open('SYMMF_feature_cor_gen.R','w')
handle_cor_rscript.write(cor_gen_R_script)
handle_cor_rscript.close()

file_list = os.listdir("./")
aval_condname_l = []
for one_f_name in file_list:
    if one_f_name.endswith('_microbe_node_value.attrs'):
        av_condition_name = one_f_name.replace('_microbe_node_value.attrs','')
        aval_condname_l.append(av_condition_name)

##calculation of correlation between each species for correlation networks
cor_gen_script_command = 'Rscript SYMMF_feature_cor_gen.R'
if os.path.isfile("cor_SYMMF_feature_bmap_selected.txt"):
    os.remove('cor_SYMMF_feature_bmap_selected.txt')
    os.system(cor_gen_script_command)
else:
    os.system(cor_gen_script_command)

#SYMMF_cutoff
fca_table = symmf_cutoff
f_in_symmf = list(set(fca_table['Group category']))

whole_sp_l_in_symmf = []
whole_sp_l_in_symmf_count = []
whole_sp_l_in_symmf_contribution = []
whole_sp_l_in_symmf_each_part = []
whole_sp_l_in_symmf_each_part_count = []
whole_sp_l_in_symmf_each_part_contribution = []

for pn_num in range(len(f_in_symmf)):
    temp_l = []
    temp_c_l = []
    temp_con_l = []
    #sp_l in each parts 
    whole_sp_l_in_symmf_each_part.append(temp_l)
    whole_sp_l_in_symmf_each_part_count.append(temp_c_l)
    whole_sp_l_in_symmf_each_part_contribution.append(temp_con_l)
    for x in fca_table[fca_table['Group category'] == f_in_symmf[pn_num]]['Microbe list and contribution']:
        sp1 = x.strip().split(',')
        for y in sp1:
            sname = y.split(':')[0]
            s_cont_v = y.split(':')[1]
            if sname in whole_sp_l_in_symmf:
                whole_sp_l_in_symmf_count[whole_sp_l_in_symmf.index(sname)] += 1
                whole_sp_l_in_symmf_contribution[whole_sp_l_in_symmf.index(sname)] += float(s_cont_v)
            else:
                whole_sp_l_in_symmf.append(sname)
                whole_sp_l_in_symmf_count.append(1)
                whole_sp_l_in_symmf_contribution.append(float(s_cont_v))
                
            if sname in whole_sp_l_in_symmf_each_part[pn_num]:
                whole_sp_l_in_symmf_each_part_count[pn_num][whole_sp_l_in_symmf_each_part[pn_num].index(sname)] += 1
                whole_sp_l_in_symmf_each_part_contribution[pn_num][whole_sp_l_in_symmf_each_part[pn_num].index(sname)] += float(s_cont_v)
            else:
                whole_sp_l_in_symmf_each_part[pn_num].append(sname)
                whole_sp_l_in_symmf_each_part_count[pn_num].append(1)
                whole_sp_l_in_symmf_each_part_contribution[pn_num].append(float(s_cont_v))

pval_cutoff = whole_network_pval_cutoff_l[0]
pval_cutoff_str = str(pval_cutoff).replace('.','').replace("-","")
ratio_file_name = 'whole_network_node_piechart_count_ratio_cor%s_p%s_auc%s_ECDF%s.table' %(pval_cutoff_str, cut_pval_str, min_auc_str, ecdf_cut_str)
ratio_write_f = open(ratio_file_name , 'w')
cont_ratio_file_name = 'whole_network_node_piechart_contribution_ratio_cor%s_p%s_auc%s_ECDF%s.table' %(pval_cutoff_str, cut_pval_str, min_auc_str, ecdf_cut_str)
cont_ratio_write_f = open(cont_ratio_file_name , 'w')
header = '\t'.join(f_in_symmf)
header = 'name\t%s\n' %(header)
ratio_write_f.write(header)
cont_ratio_write_f.write(header)
whole_sp_l_in_symmf_count_info = []
whole_sp_l_in_symmf_contrib_info = []
for s_num in range(len(whole_sp_l_in_symmf)):
    ss_nam = whole_sp_l_in_symmf[s_num]
    ss_t_count = whole_sp_l_in_symmf_count[s_num]
    ss_t_contribution = whole_sp_l_in_symmf_contribution[s_num]
    if ss_t_count != 0:
        info_string = ''
        info_string2 = ''
        temp_count_l = []
        temp_contrib_l = []
        for p_idx_num in range(len(f_in_symmf)):
            bp_name = f_in_symmf[p_idx_num]
            bp_sp_l = whole_sp_l_in_symmf_each_part[p_idx_num]
            bp_sp_cl = whole_sp_l_in_symmf_each_part_count[p_idx_num]
            bp_sp_contrib = whole_sp_l_in_symmf_each_part_contribution[p_idx_num]
            if ss_nam in bp_sp_l:
                bp_count = bp_sp_cl[bp_sp_l.index(ss_nam)]
                temp_count_l.append(bp_count)
                out_string = '%s:%d ' %(bp_name, bp_count)
                info_string += out_string
                bp_contrib = bp_sp_contrib[bp_sp_l.index(ss_nam)]
                temp_contrib_l.append(bp_contrib)
                out_string2 = '%s:%f ' %(bp_name, bp_contrib)
                info_string2 += out_string2
            else:
                temp_count_l.append(0)
                temp_contrib_l.append(0)
          
        total_count = float(whole_sp_l_in_symmf_count[whole_sp_l_in_symmf.index(ss_nam)])
        total_contrib = float(whole_sp_l_in_symmf_contribution[whole_sp_l_in_symmf.index(ss_nam)])
        temp_ratio_l = []
        temp_cont_ratio_l = []
        #for the count ratio
        for count_n in temp_count_l:
            if count_n != 0:
                ratio = (count_n/total_count) * 100
                temp_ratio_l.append(ratio)
            else:
                temp_ratio_l.append(0)
        #for the contribution
        for contrib_n in temp_contrib_l:
            if contrib_n != 0:
                cont_ratio = (contrib_n/total_contrib) * 100
                temp_cont_ratio_l.append(cont_ratio)
            else:
                temp_cont_ratio_l.append(0)
                
        whole_sp_l_in_symmf_count_info.append(info_string)
        whole_sp_l_in_symmf_contrib_info.append(info_string2)
        temp_ratio_l = list(map(str, temp_ratio_l))
        temp_cont_ratio_l = list(map(str, temp_cont_ratio_l))
        ratio_str = '\t'.join(temp_ratio_l)
        out_str = '%s\t%s\n' %(ss_nam, ratio_str)
        ratio_write_f.write(out_str)
        cont_ratio_str = '\t'.join(temp_cont_ratio_l)
        cont_out_str = '%s\t%s\n' %(ss_nam, cont_ratio_str)
        cont_ratio_write_f.write(cont_out_str)
    else:
        whole_sp_l_in_symmf_count_info.append('None')
        whole_sp_l_in_symmf_contrib_info.append('None')
ratio_write_f.close()
cont_ratio_write_f.close()


#while correaltion network generation
# pval_cutoff_l = [1e-07,1e-08,1e-09,1e-10,1e-11,1e-12,1e-13]
int_l = []
cov_l = []
neg_int_l = []
neg_cov_l = []
whole_l = []
whole_con_l = []
int_chr_l = []
file_name = 'cor_SYMMF_feature_bmap_selected.txt'
handle_f = open(file_name, 'r')
f_read = handle_f.read()
st_r = f_read.strip().replace('g__','').split('\n')
for one_i in st_r:
    i_data = one_i.split(' ')
    sp1 = i_data[1]
    sp2 = i_data[3]
    pair_l = []
    pair_l.append(sp1)
    pair_l.append(sp2)
    intr = list(set(pair_l))
    itr_pair = ' i '.join(intr)
    pval_12 = float(i_data[4])
    cor_val_12 = float(i_data[5])
    if pval_12 <= pval_cutoff:
        if cor_val_12 > 0:
            if itr_pair not in int_l:
                int_l.append(itr_pair)
                cov_l.append(cor_val_12)
                whole_l.append(itr_pair)
                whole_con_l.append(cor_val_12)
                int_chr_l.append('pos')
            else:
                one_it_idx = int_l.index(itr_pair)
                cov_l[one_it_idx] += cor_val_12
                whole_it_idx = whole_l.index(itr_pair)
                whole_con_l[whole_it_idx] += cor_val_12
        else:
            if itr_pair not in neg_int_l:
                neg_int_l.append(itr_pair)
                neg_cov_l.append(abs(cor_val_12))
            else:
                one_it_idx = neg_int_l.index(itr_pair)
                neg_cov_l[one_it_idx] += abs(cor_val_12)

whole_sif_f_name = 'whole_network_cor%s_p%s_auc%s_ECDF%s.sif' %(pval_cutoff_str, cut_pval_str, min_auc_str,ecdf_cut_str)
whole_attr_f_name = 'whole_network_edges_cor%s_p%s_auc%s_ECDF%s.attrs' %(pval_cutoff_str, cut_pval_str, min_auc_str,ecdf_cut_str)        

print('writing %s and %s is done!' %(whole_sif_f_name, whole_attr_f_name))

#writing whole network
whole_cor_write = open(whole_sif_f_name,'w')
whole_cor_write_val = open(whole_attr_f_name,'w')
whole_cor_write_val.write('interactions\n')
for idx_cro in range(len(whole_l)):
    i_name = whole_l[idx_cro]
    isp_name_l = i_name.split(' i ')
    cnum = 0
    for isp in isp_name_l:
        if isp in whole_sp_l_in_symmf:
            cnum+=1
    if cnum == 2:
        iname2 = i_name.replace(' i ',' (i) ')
        i_val = float(whole_con_l[idx_cro])
        i_char = str(int_chr_l[idx_cro])
        whole_cor_write_str = '%s\n' %(i_name)
        whole_cor_write.write(whole_cor_write_str)
        whole_cor_write_val_str = '%s = %0.5f\n'%(iname2, i_val)
        whole_cor_write_val.write(whole_cor_write_val_str)
        whole_cor_write_chr_str = '%s = %s\n'%(iname2, i_char)
whole_cor_write.close()
whole_cor_write_val.close()
