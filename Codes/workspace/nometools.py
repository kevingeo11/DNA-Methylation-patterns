import numpy as np
import pandas as pd
import subprocess
import tqdm
import time
import pickle
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency

def filter_bed_files(infile, outfile, min_cov=3):
    chroms_seen = set()
    non_std_chroms = set()
    chroms = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']

    n = 0
    with open(infile, 'r') as fin:
        n = sum(1 for _ in fin)

    progress_check = int(n/10)
    itr = 0
    with open(outfile, 'w') as fout:
        with open(infile, 'r') as fin:
            print('start traversing bed file ', infile, 'n = ', n)
            curr_time = time.time()
            for line in fin:
                if itr != 0 and itr%progress_check == 0:
                    print('progress : {}% and time elapsed {} min'.format(round(itr*100/n,2), round((time.time()-curr_time)/60,2)))

                itr += 1
                if not line.startswith('track'):
                    line_s = line.strip().split('\t')

                    assert len(line_s) == 11
                    
                    chrom = 'chr' + line_s[0]
                    if chrom in chroms:
                        chroms_seen.add(chrom)

                        start = int(line_s[1])
                        end = int(line_s[2])
                        methyl_rate = float(line_s[3])
                        coverage = int(line_s[4])
                        strand = line_s[5]

                        if coverage >= min_cov:
                            fout.write('\t'.join([str(x) for x in [chrom, start, end, strand, methyl_rate, coverage]]) + '\n')
                    else:
                        if chrom not in non_std_chroms:
                            non_std_chroms.add(chrom)

    assert len(chroms_seen) == 24

    print('filered bed file with min coverage {} and saved to {}'.format(min_cov, outfile))
    print('Non standard chrs seen : ', non_std_chroms)

def sort_bed(infile, outfile):
    '''
        sort -k1,1V -k2,2n /loc/GCH.filtered.bed > /loc/sorted.bed
    '''
    res = subprocess.run(['sort', '-k1,1V', '-k2,2n', f'{infile}'], capture_output=True, text=True)

    if res.returncode == 0:
        with open(outfile, 'w') as fout:
            fout.write(res.stdout)
    else:
        print('error in sorting')

def get_num_reads(file):
    n = 0
    with open(file, 'r') as fin:
        n = sum(1 for _ in fin)

    return n

def get_promoters_refGene(infile, outfile):
    seen_coords = set()
    itr = 0
    with open(outfile, 'w') as fout:
        with open(infile, 'r') as fin:
            for line in fin:
                if not line.startswith('#'):
                    line_s = line.strip().split('\t')

                    refid = line_s[1]
                    chrom = line_s[2]
                    strand = line_s[3]
                    txStart = int(line_s[4])
                    txEnd = int(line_s[5])
                    cdsStart = int(line_s[6])
                    cdsEnd = int(line_s[7])
                    geneName = ''.join(c for c in line_s[12] if c.isalnum())

                    coords = chrom + str(txStart) + str(txEnd) + strand
                    if coords not in seen_coords:
                        seen_coords.add(coords)

                        promo_start = 2000
                        promo_end = 1000

                        if strand == "+":
                            promoter_start = txStart - promo_start
                            promoter_end = txStart + promo_end
                        else:
                            promoter_start = txEnd - promo_end
                            promoter_end = txEnd + promo_start
                        
                        if promoter_start < 0:
                                promoter_start = 0

                        if cdsStart != cdsEnd and geneName[0:3] != 'MIR' and geneName[0:3] != 'SNO' and '_' not in chrom:
                            fout.write('\t'.join([str(s) for s in [chrom, promoter_start, promoter_end, refid, geneName, txStart, txEnd, strand]]) + '\n')
                            itr += 1

    print('no of promoters defined :', itr)
    print('promoters saved to {}'.format(outfile))

def intersect_bed(file1, file2, outfile):
    '''
        bedtools intersect -a fa -b fb -wa -wb -sorted
    '''
    res = subprocess.run(['bedtools', 'intersect', '-a', f'{file1}', '-b', f'{file2}', '-wa', '-wb', '-sorted'],
                          capture_output=True, text=True)

    if res.returncode == 0:
        with open(outfile, 'w') as fout:
            fout.write(res.stdout)
    else:
        print('error in intersect')

def complement_bed(infile, sizes, outfile):
    '''
        bedtools complement -i infile.bed -g sizes.genome > outfile.bed
    '''
    res = subprocess.run(['bedtools', 'complement', '-i', f'{infile}', '-g', f'{sizes}'],
                          capture_output=True, text=True)

    if res.returncode == 0:
        with open(outfile, 'w') as fout:
            fout.write(res.stdout)
    else:
        print('error in complement')

def filter_by_chr(infile, chrs = [], fstrand = None, SILENT=False):
    if len(chrs) == 0:
        print('input filter list')
        return
    
    n = 0
    with open(infile, 'r') as fin:
        n = sum(1 for _ in fin)
    
    progress_check = int(n/10)
    res = []
    loc = 0
    itr = 0
    with open(infile, 'r') as fin:
        if not SILENT:
            print('start traversing bed file ', infile, 'n = ', n)
        curr_time = time.time()
        for line in fin:
            if itr != 0 and itr%progress_check == 0 and not SILENT:
                print('progress : {}% and time elapsed {} min'.format(round(itr*100/n,2), round((time.time()-curr_time)/60,2)))
            
            itr += 1
            line_s = line.strip().split('\t')
            chr = line_s[0]
            if chr in chrs:
                start = int(line_s[1])
                end = int(line_s[2])
                strand = line_s[3]
                methyl_rate = float(line_s[4])
                coverage = int(line_s[5])
                
                if fstrand is None:
                    res.append([chr, start, end, strand, methyl_rate, coverage])
                else:
                    if fstrand in ['+', '-']:
                        if strand == fstrand:
                            res.append([chr, start, end, strand, methyl_rate, coverage])
                    else:
                        print('Incorrect strand given to be filtered')
                        return None
                
                loc += 1

    return res

def get_promoter_methylation(infile):
    inter_dict = dict()
    seen_ids = set()

    with open(infile) as f:
        for line in f:
            line_s = line.strip().split('\t')

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
            coverage = int(line_s[13])

            dict_id = refid + "-" + str(promo_start)
            
            ## We are checking + and -ve strands based on promoter
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
            inter_dict[dict_id][meth_start]["coverage"] = coverage

    column_names = ['trans_id', 'refid', 'gene_name', 'chrom', 'promoter_start', 'promoter_end', 
                    'TSS', 'TES', 'strand', 'meth_start_genome', 'meth_end_genome', 'meth_pos_promo_abs', 'meth_pos_promo_rel', 'meth_rate', 'coverage']
    
    info_dict = dict()
    for col in column_names:
        info_dict[col] = []

    ## Can create a filter for no of GCH or HCG in a promoter
    for dict_id in inter_dict.keys():
        refid = dict_id.split("-")[0]
        for promo_start in inter_dict[dict_id].keys():
            info_dict["trans_id"].append(dict_id)
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
            
            info_dict["meth_rate"].append(inter_dict[dict_id][promo_start]["meth_rate"])
            info_dict["coverage"].append(inter_dict[dict_id][promo_start]["coverage"])

    # Built dataframe
    df = pd.DataFrame(0, index = np.arange(len(info_dict["refid"])), columns = column_names)
    for feat in column_names:
        df[feat] = info_dict[feat]  
        
    df = df.sort_values(by = ['chrom', 'promoter_start'], ascending = [True, True])

    return inter_dict, df

def plot_avg_methylation_levels(df, context, fig = True, c = 'grey', label = None):
    all_meth_pos_promo_rel = list(df['meth_pos_promo_rel'])
    all_meth_rate = list(df['meth_rate'])

    av_dict = dict()
    seen_pos = set()
    for p in range(len(all_meth_pos_promo_rel)):
        rel_pos_x = all_meth_pos_promo_rel[p]
        meth_rate = all_meth_rate[p]
        
        #AVERAGE
        if rel_pos_x not in seen_pos:
            seen_pos.add(rel_pos_x)
            av_dict[rel_pos_x] = []
        
        if context == 'GCH':
            av_dict[rel_pos_x].append(100-meth_rate)
        elif context == 'HCG':
            av_dict[rel_pos_x].append(meth_rate)
        else:
            print('Check context')
            return

    inter_start = -2000
    inter_end = 1000

    if fig:
        plt.figure(figsize=(15, 7), facecolor='w', edgecolor='k')

    x = []
    y = []

    for rel_pos_x in sorted(av_dict.keys()):
        if rel_pos_x >= inter_start and rel_pos_x <= inter_end:
            x.append(rel_pos_x)
            y.append(np.mean(av_dict[rel_pos_x]))

    if label is None:
        plt.plot(x, y, "-", color=c)
    else:
        plt.plot(x, y, "-", color=c, label=label, alpha=0.8)
        plt.legend()
        
    plt.xlabel("DNA position [bp]")

    if context == 'GCH':
        ylab = "100-GpC methylation level"
    if context == 'HCG':
        ylab = "CpG methylation level"
    plt.ylabel(ylab)

def parse_peaks(infile, outfile, minq = 0.1):
    seen_chroms = []
    with open(infile, 'r') as fin:
        with open(outfile, 'w') as fout:
            for line in fin:
                line_s = line.strip().split('\t')

                if line_s[0] == 'chr':
                    print('unidentified chr')
                    continue

                chrom = 'chr' + line_s[0]
                start = line_s[1]
                end = line_s[2]
                q_value = float(line_s[8])
                
                if q_value >= minq:
                    fout.write('\t'.join([chrom, start, end]) + '\n')
                if chrom not in seen_chroms:
                    seen_chroms.append(chrom)

    print('chrs seen : ', seen_chroms)

def get_nuc_positions(infile):
    column_names = ["trans_id", "refid", "gene_name", "chrom", "promoter_start", "promoter_end", "TSS", "TES",
                     "strand", "nuc_region_start_genome", "nuc_region_end_genome", "nuc_start_promo_abs", 
                     "nuc_end_promo_abs", "nuc_start_promo_rel", "nuc_end_promo_rel", "region_length"]
    
    info_dict = dict()
    for col in column_names:
        info_dict[col] = []
    
    with open(infile) as fin:
        for line in fin:
            line_s = line.strip().split('\t')

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

                info_dict["trans_id"].append(refid+'-'+str(promo_start))
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
                
                info_dict["region_length"].append(rel_end - rel_start)
            
            else:
                # print('please check intersection file')
                pass

            
    # Built dataframe
    df = pd.DataFrame(0, index = np.arange(len(info_dict["refid"])), columns = column_names)
    for feat in column_names:
        df[feat] = info_dict[feat]  
        
    df = df.sort_values(by = ['chrom', 'promoter_start'], ascending=[True, True])
      
    return df

def getGCHcount_in_window(window, res, pos):
    '''
        element in res = ['chr20', 60110, 60111, '+', 0.0, 3]
    '''
    start, end = window
    freq = []
    meth_c_count = 0
    net_c_count = 0
    for itr in range(pos, len(res)):
        if start > res[itr][1]:
            ## increase the start of the search such that we start near the window
            pos += 1
        elif start <= res[itr][1] and end >= res[itr][1]:
            freq.append((res[itr][1], res[itr][4], res[itr][5], round(res[itr][4]*res[itr][5]/100)))
            meth_c_count += round(res[itr][4]*res[itr][5]/100)
            net_c_count += res[itr][5]
        elif end < res[itr][1]:
            ## the rest of sites are away from the window
            break

    unmeth_c_count = net_c_count - meth_c_count
    return meth_c_count, unmeth_c_count, freq, pos

def findNDR(filtered_res, chr, win_len = 100, multiprocess=False, temp_loc=None):
    if multiprocess:
        if temp_loc is None:
            raise Exception('Input temp loc missing')

    start = filtered_res[0][1] - 90
    end = start + win_len
    step = 20
    pos_window = 0
    pos_bg_left = 0
    pos_bg_right = 0

    n = (filtered_res[-1][1] - filtered_res[0][1])/step
    progress_check = int(n/20)

    print(f'Finding NDR windows for {chr} with window len {win_len}')

    regions = []
    itr_count = 0
    curr_time = time.time()
    while (start < filtered_res[-1][1]):
        if itr_count != 0 and itr_count%progress_check == 0:
            print('progress for {} : {}% and time elapsed {} min'.format(chr, round(itr_count*100/n,2), round((time.time()-curr_time)/60,2)))
        
        window = (start, end)
        meth_c_count, unmeth_c_count, freq, pos_window = getGCHcount_in_window(window, filtered_res, pos_window)

        if not meth_c_count > 0:
            start += step
            end += step
            itr_count += 1
            continue
        
        window_bf_left = (start-4000, start)
        meth_c_count_bg_left, unmeth_c_count_bg_left, freq_bg_left, pos_bg_left = getGCHcount_in_window(window_bf_left, filtered_res, pos_bg_left)
        
        window_bf_right = (end, end+4000)
        meth_c_count_bg_right, unmeth_c_count_bg_right, freq_bg_right, pos_bg_right = getGCHcount_in_window(window_bf_right, filtered_res, pos_bg_right)

        '''
                    Window(NDR)         background
            Meth    meth_c_count    meth_c_count_bg
            Unmeth  unmeth_c_count  unmeth_c_count_bg   
        '''

        meth_c_count_bg = meth_c_count_bg_left + meth_c_count_bg_right
        unmeth_c_count_bg = unmeth_c_count_bg_left + unmeth_c_count_bg_right

        ## same condition as above -- if condition for bg
        if not meth_c_count_bg > 0:
            start += step
            end += step
            itr_count += 1
            continue
        
        table = np.array([[meth_c_count, meth_c_count_bg],
                        [unmeth_c_count, unmeth_c_count_bg]])
        
        chi_test = chi2_contingency(table)

        if -np.log10(chi_test.pvalue) > 5:
            regions.append(window)
        
        start += step
        end += step
        itr_count += 1

    print(f'Merging found windows for {chr}')

    merged_regions = []

    n = len(regions)
    progress_check = int(n/10)

    i = 0
    curr_time = time.time()
    while i < n:
        # if i != 0 and i%progress_check == 0:
        #     print('progress : {}% and time elapsed {} min'.format(round(i*100/n,2), round((time.time()-curr_time)/60,2)))

        start, end = regions[i]

        if i+1 < n:
            for j in range(i+1, n):
                if regions[j][0] <= end:
                    end = regions[j][1]
                else:
                    i = j
                    break
        
            if j == n-1 and i != j:
                i = n

            merged_regions.append((start, end))
        else:
            i += 1
            assert i == n
            merged_regions.append((start, end))

    if multiprocess:
        file_loc = temp_loc + 'temp.NDR.' + chr + '.bed'
        with open(file_loc, 'w') as fout:
            for reg in merged_regions:
                temp = [chr, str(reg[0]), str(reg[1])]
                fout.write('\t'.join(temp) + '\n')
    else:
        return merged_regions

def plot_CpG_CpG_dist(df, max_dist=20, c='grey'):
    dist_list = []
    chroms = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']
    for chrom in chroms:
        #### ALL methylated pos
        all_meth_pos = list(df.loc[(df["meth_rate"] != 0) & 
                                                    (df["chrom"] == chrom)]["meth_start_genome"])
        all_meth_pos_sorted = sorted(list(set(all_meth_pos))) # since mapped to promoter, CpG can be twice if promoter overlap, make set

        for pos in range(len(all_meth_pos_sorted) - 1):
            start = all_meth_pos_sorted[pos]
            next_start = all_meth_pos_sorted[pos+1]
            dist = next_start - start - 1 # is number of bases inbetween
            
            if dist <= max_dist:
                dist_list.append(dist)

    y = []
    for r in range(max_dist+1):
        freq = float(dist_list.count(r))/float(len(dist_list))*100
        y.append(freq)


    plt.figure(figsize=(10, 8), facecolor='w', edgecolor='k')

    x = np.arange(len(y))
    plt.bar(x, y, align='center', color=c)

    plt.xlabel("Distance to next 5meC [bp]")
    plt.ylabel("Percentage [%]")
    plt.xlim(1, max_dist+1)
    plt.xticks(x, range(max_dist+1))
    plt.ylim(0, 10)

def get_nuc_pos_methylation(infile):
    column_names = ["trans_id", "refid", "gene_name", "chrom", "promoter_start", "promoter_end", "TSS", "TES", "strand",
                    "nuc_region_start_genome", "nuc_region_end_genome", "nuc_region_length",
                    "meth_start_genome", "meth_end_genome", "meth_rate"]
    
    info_dict = dict()
    for col in column_names:
        info_dict[col] = []
    
    with open(infile, 'r') as fin:
        for line in fin:
            line_s = line.strip().split('\t')

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
            
            meth_start = int(line_s[20])
            meth_end = int(line_s[21])
            meth_rate = float(line_s[23])
            
            if nuc_region_start >= promo_start and nuc_region_end <= promo_end:
                info_dict["trans_id"].append(refid + "-" + str(promo_start))
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
                info_dict["nuc_region_length"].append(nuc_region_end-nuc_region_start)
                
                info_dict["meth_start_genome"].append(meth_start)
                info_dict["meth_end_genome"].append(meth_end)
                info_dict["meth_rate"].append(meth_rate)

            
    #Built dataframe
    df = pd.DataFrame(0, index = np.arange(len(info_dict["refid"])), columns = column_names)
    for feat in column_names:
        df[feat] = info_dict[feat]  
        
    df = df.sort_values(by=['chrom', 'promoter_start'], ascending=[True, True])

    return df

def plot_nbr_clashs(steric_path, clash_dict_loc, info_dict_loc, save=False):
    with open(clash_dict_loc, 'rb') as fin:
        clash_dict = pickle.load(fin)
    with open(info_dict_loc, 'rb') as fin:
        info_nbr_dict = pickle.load(fin)

    nbr_dnmt_residues = info_nbr_dict["model_dnmt_nbr_residues"]
    nbr_dnmt_atoms = info_nbr_dict["model_dnmt_nbr_atoms"]
    dnmt_atoms_consider_clash = info_nbr_dict["dnmt_atoms_consider_clash"]
    
    x = sorted(clash_dict.keys())
    #x = [12]
    y = []

    #lab_txt = []
    for pos in x:
        #nbr_clash_residues = len(clash_dict[pos]['steric_clash_list'].keys())
        #perc_clash_res = (float(nbr_clash_residues)/float(nbr_dnmt_residues))*100
        #rmsd = round(clash_dict[pos]['rmsd'],2)
        all_clashes = 0
        for clash_res in clash_dict[pos]['steric_clash_list'].keys():
            atom_clashes = clash_dict[pos]['steric_clash_list'][clash_res]
            nbr_clash_atoms = 0
            seen_atoms = set()
            for ac in atom_clashes: #['X_TYR923_N', 'A_ARG42_CD', 2.0103261]
                dnmt_atom = ac[0]
                nuc_atom = ac[1]
                if dnmt_atom not in seen_atoms:
                    seen_atoms.add(dnmt_atom)
                    nbr_clash_atoms += 1
            all_clashes += nbr_clash_atoms

        perc_clash_atoms = (float(all_clashes)/float(dnmt_atoms_consider_clash))*100
        y.append(perc_clash_atoms)
    
    #norm clashed st acc is probability between 0 and 1
    max_clash = np.max(y)
    min_clash = np.min(y)
    print ("max_clash", max_clash)
    print ("min_clash", min_clash)
    diff = max_clash - min_clash
    y_access = []
    for y_val in y:
        norm_val = 100-(float(y_val-min_clash)/float(diff))*100 #100- because is accessibility
        y_access.append(norm_val)
    
    if save:
        with open(steric_path + 'x_y_perc_file.txt', 'w') as fout:
            x_y_acc_dict = dict() 
            x_y_dict = dict()
            for i in range(0,len(x)):
                fout.write(str(x[i]) + "----" + str(y[i]) + "\n")
                x_y_dict[x[i]] = y[i]
                x_y_acc_dict[x[i]] = y_access[i]

        with open(steric_path + 'x_y_dict', 'wb') as fout:
            pickle.dump(x_y_dict, fout, -1)
        with open(steric_path + 'x_y_acc_dict', 'wb') as fout:
            pickle.dump(x_y_acc_dict, fout, -1)
    
    '''
    % ATOMS CLASHED
    '''
    fig = plt.figure(figsize=(15, 7))
    ax = plt.subplot(111)
    
    #plt.title('Superposition of DNMT1 (3PTA) and the nucleosome (1KX5)')
    ax.plot(x,y,linestyle="-",marker="o",linewidth=2)
    for pos in range(0,len(x)):
        #ax.text(x[pos]+0.3,y[pos],rmsd_list[pos])
        ax.text(x[pos]+0.3,y[pos],x[pos],fontsize=12)
        
    plt.xlim(0,140)
    plt.xticks(range(0,len(x)+5,5))
    plt.xlabel("Nucleosome position [bp]")
    plt.ylabel("DNMT1 atoms clashing [%]")
    
    ax.axvline(x=74,c="grey",linewidth=2,linestyle="--")#4169e1
    
    
    '''
    ACCESSIBILITY: 100-clasehd atoms
    '''
    fig = plt.figure(figsize=(15, 7))
    ax = plt.subplot(111)
    plt.title('Superposition of DNMT1 (3PTA) and the nucleosome (1KX5)')
    ax.plot(x,y_access,linestyle="-",marker="o",linewidth=2)
    for pos in range(0,len(x)):
        ax.text(x[pos]+0.3,y_access[pos],x[pos])
        
    plt.xlim(0,140)
    plt.xticks(range(0,len(x)+5,5))
    plt.xlabel("Nucleosome position")
    plt.ylabel("Accessibility [%]")
    
    ax.axvline(x=74,c="black",linewidth=2,linestyle="-")#4169e1
    ax.axhline(y=95,c="black",linewidth=2,linestyle="--")#4169e1
    ax.text(-1,94.5,'95')



