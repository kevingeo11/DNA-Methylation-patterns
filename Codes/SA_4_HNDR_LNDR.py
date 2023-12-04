import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def convert_to_karlHMM_format(infile, outfile):
    itr = 0 

    with open(infile, 'r') as fin:
        with open(outfile, 'w') as fout:
            fout.write('track name=KarlInp type=bedDetail' + '\n')
            for line in fin:
                itr += 1

                line_s = line.strip().split('\t')
                chrom = line_s[0].replace('chr','')
                start = line_s[1]
                end = line_s[2]
                methyl_rate = line_s[4]
                coverage = line_s[6]
                strand = line_s[3]

                fout.write('\t'.join([str(s) for s in [chrom, start, end, methyl_rate, coverage, strand]]) + '\n')
            
def count(infile, check):
    c = 0
    with open(infile, 'r') as fin:
        for line in fin:
            if line.startswith(check):
                c += 1

    print('count =', c, 'for check: ', check)

def parse_NDR_karl(karl_file, outfile_NDR):
    seen_chroms = []
    with open(karl_file, 'r') as fin:
        with open(outfile_NDR, 'w') as fout:
            for line in fin:
                line_s = line.strip().split('\t')
                chrom = 'chr' + line_s[0]
                start = line_s[1]
                end = line_s[2]

                if line_s[0] == 'chr':
                    continue

                fout.write('\t'.join([chrom, start, end]) + '\n')
                if chrom not in seen_chroms:
                    seen_chroms.append(chrom)

    print(seen_chroms)

def make_pandas_promo_nuc_file(inter_out_promo, outfile):
    column_names = ["refid", "gene_name", "chrom", "promoter_start", "promoter_end", "TSS", "TES",
                     "strand", "nuc_region_start_genome", "nuc_region_end_genome", "nuc_start_promo_abs", 
                     "nuc_end_promo_abs", "nuc_start_promo_rel", "nuc_end_promo_rel", "region_length"]
    
    info_dict = dict()
    for col in column_names:
        info_dict[col] = []
    
    with open(inter_out_promo) as f:
        reader = csv.reader(f, delimiter ='\t')
        read = False
        for line_s in reader:
            chrom = line_s[0]            
            promo_start = int(line_s[1])
            promo_end = int(line_s[2])
            refid = line_s[3]
            gene_name = line_s[4]
            TSS = int(line_s[5])
            TES = int(line_s[6])
            strand = line_s[7]

            nuc_region_start = int(line_s[9])
            nuc_region_end = int(line_s[10])
            

            if nuc_region_start >= promo_start and nuc_region_end <= promo_end:    
                if strand == "+":
                    promo_abs_start = nuc_region_start - promo_start
                    promo_abs_end = nuc_region_end - promo_start
                    
                if strand == "-":
                    promo_abs_end = promo_end - nuc_region_start
                    promo_abs_start = promo_end - nuc_region_end

                rel_start = promo_abs_start - 2000
                rel_end = promo_abs_end - 2000

                info_dict["chrom"].append(chrom)
                info_dict["promoter_start"].append(promo_start)
                info_dict["promoter_end"].append(promo_end)
                info_dict["refid"].append(refid)
                info_dict["gene_name"].append(gene_name)
                info_dict["TSS"].append(TSS)
                info_dict["TES"].append(TES)
                info_dict["strand"].append(strand)
                
                info_dict["nuc_region_start_genome"].append(nuc_region_start)
                info_dict["nuc_region_end_genome"].append(nuc_region_end)
                
                info_dict["nuc_start_promo_abs"].append(promo_abs_start)
                info_dict["nuc_end_promo_abs"].append(promo_abs_end)
                
                info_dict["nuc_start_promo_rel"].append(rel_start)
                info_dict["nuc_end_promo_rel"].append(rel_end)
                
                info_dict["region_length"].append(rel_end-rel_start)

            
    # Built dataframe
    df = pd.DataFrame(0, index = np.arange(len(info_dict["refid"])), columns = column_names)
    for feat in column_names:
        df[feat] = info_dict[feat]  
        
    df = df.sort_values(by = ['chrom', 'promoter_start'], ascending=[True, True])
      
    df.to_csv(outfile) 

def plot_nuc_pos_distribution(df_inter, region, outname):
    with open(outname + "_nuc_promo_info.txt", 'w') as txt_info:
    
        rel_starts = list(df_inter["nuc_start_promo_rel"])
        rel_ends = list(df_inter["nuc_end_promo_rel"])
        region_length = list(df_inter["region_length"])
        
        refids = set(list(df_inter["refid"]))
        
        txt_info.write(", ".join(["number of refids (promoters)",str(len(refids))])+"\n")
        txt_info.write(", ".join(["number of nucleosome positions in promoters",str(len(rel_starts))])+"\n\n")
        
        txt_info.write("region length:\n")
        txt_info.write(", ".join(["mean",str(round(np.mean(region_length)))])+"\n")
        txt_info.write(", ".join(["std",str(round(np.std(region_length)))])+"\n")
        txt_info.write(", ".join(["median",str(round(np.median(region_length)))])+"\n")
        txt_info.write(", ".join(["max",str(np.max(region_length))])+"\n")
        txt_info.write(", ".join(["min",str(np.min(region_length))])+"\n")
    
    ######### PLOT 1: distribution
    
    size_plot = 18
    axis_font = {'fontname':'Arial', 'size':str(size_plot)}

    inter_start = -2000
    inter_end = 1000
    step = 200
    inter = range(inter_start,inter_end+1,step)
    
    plt.figure(num=None, figsize=(10, 8), facecolor='w', edgecolor='k')
    ax = plt.subplot(1,1,1)
    
    zed = [tick.label.set_fontsize(size_plot) for tick in ax.yaxis.get_major_ticks()]
    zed = [tick.label.set_fontsize(size_plot) for tick in ax.xaxis.get_major_ticks()]
    
    sns.distplot(rel_starts,color="#D35400",kde_kws={"linewidth":3},label = "Start positions")   #bins=50    region.replace("s","") +" starts"
    sns.distplot(rel_ends,color="#2874A6",kde_kws={"linewidth":3},label = "End positions")  #region.replace("s","")+" ends"
    
    plt.xlabel("DNA position [bp]",**axis_font)
    plt.ylabel("Density",**axis_font)
    plt.legend(frameon=1,prop={'size':14})
    
    plt.xlim(-2000,1000)
    #plt.xticks(inter)
    plt.xticks(range(-2000,1001,500))
    plt.axvline(x=0,linestyle='dashed',linewidth=2,color='grey')
    
    
    plt.savefig(outname + "nuc_pos_distr_promo.png",bbox_inches='tight')
    plt.show()
    plt.close()    
    

    ######### PLOT 2: nuc length
    plt.figure(num=None, figsize=(10, 8), facecolor='w', edgecolor='k')
    ax = plt.subplot(1,1,1)
    
    zed = [tick.label.set_fontsize(size_plot) for tick in ax.yaxis.get_major_ticks()]
    zed = [tick.label.set_fontsize(size_plot) for tick in ax.xaxis.get_major_ticks()]
    
    sns.distplot(region_length,kde_kws={"linewidth":3},color="grey")   #bins=50
    
    plt.xlabel("Length [bp]",**axis_font)
    plt.ylabel("Density",**axis_font)
    
    plt.xlim(0,2000)
    plt.xticks(range(0,2000,200))
    
    plt.savefig(outname + "nuc_length_promo.png",bbox_inches='tight')
    plt.close()   

main_path = '/home/kevin/DNA-Methylation-patterns/'
inpath_processing = main_path + 'Datasetb37_bismark/downstream/'
outpath_processing = main_path + 'Datasetb37_bismark/downstream/'

'''
HCG, GCH
'''
NOMe_infile = inpath_processing + 'NOMe_filtered_bed_sorted.csv'
NOMe_outfile = outpath_processing + 'NOMe_in_karl.csv'

## NOMe sorted file erased
# convert_to_karlHMM_format(NOMe_infile, NOMe_outfile)

# count(NOMe_outfile, 'chrX\t')

if True:
    regions = ["NORs", "NDRs"]
    
    for region in regions:
        print(region)
        
        ####### MAKE NOR and NDR files
        
        outfile_NDR = outpath_processing + "NOMe_NDR.bed"      # bed file from karl file
        outfile_NOR = outpath_processing + "NOMe_NOR.bed"
        
        if False and region == 'NDRs':
            karl_file = outpath_processing + "NOMe_out_karl.bed"
            parse_NDR_karl(karl_file, outfile_NDR)
            break
        
        #take complement of NDR to get dense regions
        # if region == "NORs":
        #     if False:
        #         file_genome = outpath_processing + "1_hg19_chrom_sizes.txt"
        #         make_complement(outfile_NDR, file_genome, outfile_NOR)

        ####### MAP PROMOTER to NOR/NDR and make PANDAS table
        outfile_UCSC_promo = main_path + 'refSeq h19\\' + '0-refGene_complete_hg19.csv'
        output_file_nuc_regions = outfile_NDR if region=="NDRs" else outfile_NOR
        output_file_inter = outpath_processing + 'NOMe_{}_intersect.bed'.format(region[:-1])
        outfile_inter_pandas = outpath_processing + 'pandas_NOMe_{}_intersectBed_promoter_NucRegions.csv'.format(region[:-1])

        if False:
            print("map:promoter - nuc position")
            make_intersect(outfile_UCSC_promo, output_file_nuc_regions, output_file_inter)

        if False:
            print("parse map:promoter - nuc position")
            make_pandas_promo_nuc_file(output_file_inter, outfile_inter_pandas)

        if True:
            print("plot nucleosome position in the promoter region")
            df_inter = pd.read_csv(outfile_inter_pandas, index_col = 0)
            outname = outpath_processing + "plot_" + region[:-1] + "_"
            plot_nuc_pos_distribution(df_inter, region, outname)


