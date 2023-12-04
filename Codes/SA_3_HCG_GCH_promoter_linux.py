import random
import pandas as pd
import numpy as np
import csv
import matplotlib.pyplot as plt

print('parse mod was the 1st fun used')

def parse_intersect_result(inter_out_promo):
    inter_dict = dict()
    seen_ids = set()

    with open(inter_out_promo) as f:
        reader = csv.reader(f, delimiter = '\t')
        for line_s in reader:    
            chrom = line_s[0]            
            promo_start = int(line_s[1])
            promo_end = int(line_s[2])
            refid = line_s[3]
            gene_name = line_s[4]
            TSS = int(line_s[5])
            TES = int(line_s[6])
            strand = line_s[7]

            meth_start = int(line_s[9])
            meth_end = int(line_s[10])
            meth_rate = float(line_s[12])
            nt = line_s[13]

            #diff transcript variants, e.g. gene NM_001105522 with '89644578-89653576' and '89657231-89666229' and same methylation pattern, if key is refid, loose one entry
            dict_id = refid + "-" + str(promo_start)
            
            if strand == "+":
                promo_abs_pos = meth_start - promo_start
            if strand == "-":
                promo_abs_pos = promo_end - meth_end
        
            rel_pos = promo_abs_pos - 2000
        
            if dict_id not in seen_ids:
                seen_ids.add(dict_id)
                inter_dict[dict_id] = dict()
            
            inter_dict[dict_id][meth_start] = dict()
            
            inter_dict[dict_id][meth_start]["chrom"] = chrom
            inter_dict[dict_id][meth_start]["promoter_start"] = promo_start
            inter_dict[dict_id][meth_start]["promoter_end"] = promo_end
            inter_dict[dict_id][meth_start]["gene_name"] = gene_name
            inter_dict[dict_id][meth_start]["TSS"] = TSS
            inter_dict[dict_id][meth_start]["TES"] = TES
            inter_dict[dict_id][meth_start]["strand"] = strand
            
            inter_dict[dict_id][meth_start]["meth_start_genome"] = meth_start
            inter_dict[dict_id][meth_start]["meth_end_genome"] = meth_end
            inter_dict[dict_id][meth_start]["meth_pos_promo_abs"] = promo_abs_pos
            inter_dict[dict_id][meth_start]["meth_pos_promo_rel"] = rel_pos
                        
            inter_dict[dict_id][meth_start]["meth_rate"] = meth_rate
            inter_dict[dict_id][meth_start]["nt"] = nt

            chrom, promo_start, promo_end, refid, gene_name, TSS, TES, strand, meth_start, meth_end, meth_rate, nt = ['ERROR']*12
                
    return inter_dict

def parse_intersect_result_mod(inter_out_promo):
    inter_dict = dict()
    seen_ids = set()

    with open(inter_out_promo) as f:
        reader = csv.reader(f, delimiter = '\t')
        read = False
        for line_s in reader:
            if line_s[0] == '':
                meth_start = int(line_s[2])
                meth_end = int(line_s[3])
                meth_rate = float(line_s[5])
                nt = line_s[6]
                read = True
            else:
                chrom = line_s[0]            
                promo_start = int(line_s[1])
                promo_end = int(line_s[2])
                refid = line_s[3]
                gene_name = line_s[4]
                TSS = int(line_s[5])
                TES = int(line_s[6])
                strand = line_s[7]

            if read:
                read = False

                #diff transcript variants, e.g. gene NM_001105522 with '89644578-89653576' and '89657231-89666229' and same methylation pattern, if key is refid, loose one entry
                dict_id = refid + "-" + str(promo_start)
                
                if strand == "+":
                    promo_abs_pos = meth_start - promo_start
                if strand == "-":
                    promo_abs_pos = promo_end - meth_end
            
                rel_pos = promo_abs_pos - 2000
            
                if dict_id not in seen_ids:
                    seen_ids.add(dict_id)
                    inter_dict[dict_id] = dict()
                
                inter_dict[dict_id][meth_start] = dict()
                
                inter_dict[dict_id][meth_start]["chrom"] = chrom
                inter_dict[dict_id][meth_start]["promoter_start"] = promo_start
                inter_dict[dict_id][meth_start]["promoter_end"] = promo_end
                inter_dict[dict_id][meth_start]["gene_name"] = gene_name
                inter_dict[dict_id][meth_start]["TSS"] = TSS
                inter_dict[dict_id][meth_start]["TES"] = TES
                inter_dict[dict_id][meth_start]["strand"] = strand
                
                inter_dict[dict_id][meth_start]["meth_start_genome"] = meth_start
                inter_dict[dict_id][meth_start]["meth_end_genome"] = meth_end
                inter_dict[dict_id][meth_start]["meth_pos_promo_abs"] = promo_abs_pos
                inter_dict[dict_id][meth_start]["meth_pos_promo_rel"] = rel_pos
                            
                inter_dict[dict_id][meth_start]["meth_rate"] = meth_rate
                inter_dict[dict_id][meth_start]["nt"] = nt

                chrom, promo_start, promo_end, refid, gene_name, TSS, TES, strand, meth_start, meth_end, meth_rate, nt = ['ERROR']*12
                
    return inter_dict

def make_pandas_df(inter_out_promo, outfile):
    column_names = ['trans_id', 'refid', 'gene_name', 'chrom', 'promoter_start', 'promoter_end', 
                    'TSS', 'TES', 'strand', 'meth_start_genome', 'meth_end_genome', 'meth_pos_promo_abs', 'meth_pos_promo_rel', 'nt', 'meth_rate']
    
    info_dict = dict()
    for col in column_names:
        info_dict[col] = []
    
    inter_dict = parse_intersect_result(inter_out_promo)
    
    # take only genes with both NOME and WGBS info
    
    for dict_id in inter_dict.keys():
        refid = dict_id.split("-")[0]
        for promo_start in inter_dict[dict_id].keys():
            info_dict["trans_id"].append(dict_id) # refid+promo_start to have various transcript ids and do not miss some
            info_dict["refid"].append(refid)
            info_dict["gene_name"].append(inter_dict[dict_id][promo_start]["gene_name"])
            info_dict["chrom"].append(inter_dict[dict_id][promo_start]["chrom"])
            info_dict["promoter_start"].append(inter_dict[dict_id][promo_start]["promoter_start"])
            info_dict["promoter_end"].append(inter_dict[dict_id][promo_start]["promoter_end"])
            info_dict["TSS"].append(inter_dict[dict_id][promo_start]["TSS"])
            info_dict["TES"].append(inter_dict[dict_id][promo_start]["TES"])
            info_dict["strand"].append(inter_dict[dict_id][promo_start]["strand"])
            
            info_dict["meth_start_genome"].append(inter_dict[dict_id][promo_start]["meth_start_genome"])
            info_dict["meth_end_genome"].append(inter_dict[dict_id][promo_start]["meth_end_genome"])
            info_dict["meth_pos_promo_abs"].append(inter_dict[dict_id][promo_start]["meth_pos_promo_abs"])
            info_dict["meth_pos_promo_rel"].append(inter_dict[dict_id][promo_start]["meth_pos_promo_rel"])
            
            info_dict["nt"].append(inter_dict[dict_id][promo_start]["nt"])
            info_dict["meth_rate"].append(inter_dict[dict_id][promo_start]["meth_rate"])
            
    print('created dictionary---')
    
    # Built dataframe
    df = pd.DataFrame(0, index = np.arange(len(info_dict["refid"])), columns = column_names)
    for feat in column_names:
        df[feat] = info_dict[feat]  
        
    df = df.sort_values(by = ['chrom', 'promoter_start'], ascending = [True, True])
      
    df.to_csv(outfile)

def randomize_genome_methylation(df_promo_nuc_WGBS_pandas, inter_out_promo_WGBS, outfile_random):    
    all_meth_rates = list(df_promo_nuc_WGBS_pandas['meth_rate'])    # [0:10]
    # Return a k length list of unique elements chosen from the population sequence. Used for random sampling without replacement.
    all_meth_rates_random = random.sample(all_meth_rates, len(all_meth_rates))   
    
    # read csv file and add methylation rate from new random list
    """
        There is some error here
    """
    line_pos = 0
    with open(outfile_random, 'w') as fout:
        with open(inter_out_promo_WGBS, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for line_s in reader:
                line_s[12] = str(all_meth_rates_random[line_pos%len(all_meth_rates)])
                line_pos += 1
                fout.write('\t'.join(line_s) + '\n')

def plot_NOMe_WGBS_in_promoter(df_inter_promo_NOMe, dataset_name, plot_out):
    all_trans_ids = set(list(df_inter_promo_NOMe['trans_id']))
    all_meth_pos_promo_rel = list(df_inter_promo_NOMe['meth_pos_promo_rel'])
    all_meth_rate = list(df_inter_promo_NOMe['meth_rate'])
    
    # with open(plot_out.replace('.pdf','.txt'), 'w') as txt_info:
    #     txt_info.write("Average methylation per DNA position \n Nbr. genes/promoters/transcripts: " + 
    #                    str(len(all_trans_ids)) + " Nbr. methylation rates/promoter positions: " + 
    #                    str(len(all_meth_rate))+"\n ")
    
    av_dict = dict()
    seen_pos = set()
    for p in range(len(all_meth_pos_promo_rel)):
        rel_pos_x = all_meth_pos_promo_rel[p]
        meth_rate = all_meth_rate[p]
        
        #AVERAGE
        if rel_pos_x not in seen_pos:
            seen_pos.add(rel_pos_x)
            av_dict[rel_pos_x] = []
            
        if dataset_name == "NOMe":
            av_dict[rel_pos_x].append(100-meth_rate)

        if dataset_name == "WGBS":
            av_dict[rel_pos_x].append(meth_rate)
    
    '''
    AVERAGE
    '''
    inter_start = -2000
    inter_end = 1000
    step = 200
    inter = range(inter_start,inter_end+1,step)
    
    size_plot = 14
    axis_font = {'size':str(size_plot)}
    
    plt.figure(figsize=(20, 7), facecolor='w', edgecolor='k')
    
    x = []
    y = []
    
    for rel_pos_x in sorted(av_dict.keys()):
        if rel_pos_x >= inter_start and rel_pos_x <= inter_end:
            x.append(rel_pos_x)
            y.append(np.mean(av_dict[rel_pos_x]))

    plt.plot(x, y, "-", color = "grey")
        
    plt.xlabel("DNA position [bp]", **axis_font)
    
    if dataset_name == "NOMe":
        ylab = "100-GpC methylation level"

    if dataset_name == "WGBS":
        ylab = "CpG methylation level"

        # plt.ylim(10,60)
        
    plt.ylabel(ylab, **axis_font)
        
    # plt.xticks(inter)
    plt.axvline(x=0,linestyle='dashed',linewidth=1,color='grey')
    
    plt.savefig(plot_out, bbox_inches = 'tight')
    #plt.show()
    plt.close()

def plot_WGBS_CpG_distances_in_promoter(df_inter_promo_WGBS, plot_out):
    max_dist = 20

    dist_list = []
    chroms = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']
    for chrom in chroms:
        #### ALL methylated pos
        all_meth_pos = list(df_inter_promo_WGBS.loc[(df_inter_promo_WGBS["nt"] == "C") & 
                                                    (df_inter_promo_WGBS["meth_rate"] != 0) & 
                                                    (df_inter_promo_WGBS["chrom"] == chrom)]["meth_start_genome"])
        all_meth_pos_sorted = sorted(list(set(all_meth_pos))) # since mapped to promoter, CpG can be twice if promoter overlap, make set

        for pos in range(len(all_meth_pos_sorted) - 1):
            start = all_meth_pos_sorted[pos]
            next_start = all_meth_pos_sorted[pos+1]
            dist = next_start-start - 1 #is number of bases inbetween
            
            #if dist == 1:
            #    print start,next_start
            #    break
            
            if dist <= max_dist:
                dist_list.append(dist)

        
    print(len(dist_list))
    
    with open(plot_out.replace('.pdf','.txt'), 'w') as txt_info:
        txt_info.write("Number of cytosines with: meth. rate > 0, max. dist <= " + 
                       str(max_dist) + ", promoter regions: " + str(len(dist_list)) + " CpG positions\n\n")
        txt_info.write("Nbr of cytosines for every distance r <= " + str(max_dist) + "\n")
    
        y = []
        for r in range(max_dist+1):
            txt_info.write(str(r)+": "+str(dist_list.count(r))+"\n")
            c = float(dist_list.count(r))/float(len(dist_list))*100
            y.append(c)
    
    ####### PLOT
    size_plot = 18
    axis_font = {'size':str(size_plot)}
    
    plt.figure(figsize=(10, 8), facecolor='w', edgecolor='k')
    ax = plt.subplot(1,1,1)
    
    zed = [tick.label.set_fontsize(size_plot) for tick in ax.xaxis.get_major_ticks()]
    zed = [tick.label.set_fontsize(size_plot) for tick in ax.yaxis.get_major_ticks()]
    
    x = np.arange(len(y))
    c  = "grey" if "RANDOM" in plot_out else "#333232"
    plt.bar(x, y, align='center', color=c)
    
    plt.xlabel("Distance to next 5meC [bp]",**axis_font)
    plt.ylabel("Percentage [%]",**axis_font)
    
    ax.set_xlim(1,max_dist+1)
    plt.xticks(x,range(max_dist+1))
    
    #plt.xlim(0,max_dist+1)
    ax.set_ylim(0,10)
    
    plt.savefig(plot_out, bbox_inches='tight')
    plt.show()
    plt.close()

def plot_NOMe_WGBS_RANDOM_in_promoter(df_inter_promo_NOMe, dataset_name, ds):
    all_trans_ids = set(list(df_inter_promo_NOMe["trans_id"]))
    all_meth_pos_promo_rel = list(df_inter_promo_NOMe["meth_pos_promo_rel"])
    all_meth_rate = list(df_inter_promo_NOMe["meth_rate"])
    
    av_dict = dict()
    seen_pos = set()
    for p in range(len(all_meth_pos_promo_rel)):
        rel_pos_x = all_meth_pos_promo_rel[p]
        meth_rate = all_meth_rate[p]
        
        #AVERAGE
        if rel_pos_x not in seen_pos:
            seen_pos.add(rel_pos_x)
            av_dict[rel_pos_x] = []
            
        if dataset_name == "NOMe":
            av_dict[rel_pos_x].append(100-meth_rate)
        if dataset_name == "WGBS":
            av_dict[rel_pos_x].append(meth_rate)
    
    '''
    AVERAGE
    '''
    inter_start = -2000
    inter_end = 1000
    step = 200
    inter = range(inter_start,inter_end+1,step)
    
    size_plot = 14
    axis_font = {'size':str(size_plot)}
    
    
    
    zed = [tick.label.set_fontsize(size_plot) for tick in ax.xaxis.get_major_ticks()]
    zed = [tick.label.set_fontsize(size_plot) for tick in ax.yaxis.get_major_ticks()]
    
    
    x = []
    y = []
    for rel_pos_x in sorted(av_dict.keys()):
        if rel_pos_x >= inter_start and rel_pos_x<=inter_end:
            x.append(rel_pos_x)
            y.append(np.mean(av_dict[rel_pos_x]))

    c  = "grey" if ds=="Random" else "#333232"
    plt.plot(x,y,"-",color=c,linewidth=3,label=ds)
        
    plt.xlabel("DNA position [bp]",**axis_font)
    
    if dataset_name == "NOMe":
        ylab = "100-GpC methylation level"
    if dataset_name == "WGBS":
        ylab = "CpG methylation level"
        # plt.ylim(10,60)
        
    plt.ylabel(ylab,**axis_font)
        
    plt.xticks(inter)
    
    plt.axvline(x=0,linestyle='dashed',linewidth=1,color="black")
    plt.legend(frameon=1,prop={'size':14})


main_path = '/home/kevin/DNA-Methylation-patterns/'
data_path = main_path + 'Datasetb37_bismark/'
inpath_processing = main_path + 'Datasetb37_bismark/downstream/'
outpath_processing = main_path + 'Datasetb37_bismark/downstream/'

inter_out_promo_NOMe = inpath_processing + 'NOMe_intersect.csv'
inter_out_promo_WGBS = inpath_processing + 'WGBS_intersect.csv'
# inter_out_promo_NOMe = inpath_processing + 'NOMe_intersect_10k.csv'
# inter_out_promo_WGBS = inpath_processing + 'WGBS_intersect_10k.csv'

#pandas intersect outfiles
outfile_pandas_NOMe = outpath_processing + 'pandas_promoter_NOMe.csv'
outfile_pandas_WGBS = outpath_processing + 'pandas_promoter_WGBS.csv'
# outfile_pandas_NOMe = outpath_processing + 'pandas_promoter_NOMe_10m.csv'
# outfile_pandas_WGBS = outpath_processing + 'pandas_promoter_WGBS_10k.csv'

#random files for sliding window
inter_out_promo_WGBS_RANDOM = outpath_processing + 'intersectBed_UCSC_complete_promoter_WGBS_RANDOM.csv'
outfile_pandas_WGBS_RANDOM = outpath_processing + 'pandas_promoter_WGBS_RANDOM.csv'

'''
Parse files and make pandas df
'''
if False:
    # print('NOMe')
    # make_pandas_df(inter_out_promo_NOMe,outfile_pandas_NOMe)
    
    print('WGBS')
    make_pandas_df(inter_out_promo_WGBS,outfile_pandas_WGBS)

'''
Randomize WGBS file for sliding window approach
Randomize it here s.t. NORs and NDRs based on same randomized WGBS file
'''
if False:
    print('randomize promoter methylation')
    df_promo_nuc_WGBS_pandas = pd.read_csv(outfile_pandas_WGBS, index_col = 0)
    
    randomize_genome_methylation(df_promo_nuc_WGBS_pandas, inter_out_promo_WGBS, inter_out_promo_WGBS_RANDOM)    
    make_pandas_df(inter_out_promo_WGBS_RANDOM, outfile_pandas_WGBS_RANDOM)

'''
Analyze GCH and HCG in promoter
'''

#NOMe
if False:
    #plot average 100-GpC methylation around promoter
    df_inter_promo_NOMe = pd.read_csv(outfile_pandas_NOMe, index_col = 0)
    plot_out = outpath_processing + 'NOMe_promoter_average.png'
    plot_NOMe_WGBS_in_promoter(df_inter_promo_NOMe, "NOMe", plot_out)

#WGBS  
if True:
    df_inter_promo_WGBS = pd.read_csv(outfile_pandas_WGBS, index_col = 0)
    
    #calc distances from C to C with methylation rate != 0  
    plot_out = outpath_processing + "WGBS_promoter_distances_CpGs.png"
    plot_WGBS_CpG_distances_in_promoter(df_inter_promo_WGBS,plot_out)
    
    #calc 100-CgG methylation
    plot_out = outpath_processing + "WGBS_promoter_average_CpG_meth.png"
    plot_NOMe_WGBS_in_promoter(df_inter_promo_WGBS, "WGBS", plot_out)

#RANDOM
if True: 
    df_inter_promo_WGBS_RANDOM = pd.read_csv(outfile_pandas_WGBS_RANDOM, index_col=0)
    
    plot_out = outpath_processing + "WGBS_RANDOM_distances_CpGs.png"
    plot_WGBS_CpG_distances_in_promoter(df_inter_promo_WGBS_RANDOM,plot_out)
    
    plot_out = outpath_processing + "WGBS_RANDOM_promoter_average_CpG_meth.png"
    plot_NOMe_WGBS_in_promoter(df_inter_promo_WGBS_RANDOM,"WGBS",plot_out)

#RANDOM + WGBS
if True:
    df_inter_promo_WGBS = pd.read_csv(outfile_pandas_WGBS, index_col=0)
    df_inter_promo_WGBS_RANDOM = pd.read_csv(outfile_pandas_WGBS_RANDOM, index_col=0)
    plot_out = outpath_processing + "WGBS_AND_RANDOM_promoter_average_CpG_meth.png"
    
    plt.figure(figsize=(20, 7), facecolor='w', edgecolor='k')
    ax = plt.subplot(1,1,1)
    
    plot_NOMe_WGBS_RANDOM_in_promoter(df_inter_promo_WGBS_RANDOM,"WGBS","Random")   #random first to be in background
    plot_NOMe_WGBS_RANDOM_in_promoter(df_inter_promo_WGBS,"WGBS","Experimental")
    

    plt.savefig(plot_out,bbox_inches='tight')
    plt.show()
    plt.close() 
