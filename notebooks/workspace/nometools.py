import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import time
import pickle
import multiprocessing 
from scipy.stats import chi2_contingency

from . import bedtools as bedtools

def get_num_reads(file):
    '''
        return number of read in the file
    '''
    n = 0
    with open(file, 'r') as fin:
        n = sum(1 for _ in fin)

    return n

def filter_bed_files(infile, outfile, min_cov=3):
    '''
        Reads a BED file
        columns : chr start end methylation coverage strand 
        filters reads less than min_cov
    '''
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

                    # bed files form vcf2bed.pl contain 11 columns
                    assert len(line_s) == 11
                    
                    chrom = line_s[0]
                    if chrom in chroms:
                        chroms_seen.add(chrom)

                        start = int(line_s[1])
                        end = int(line_s[2])
                        methyl_rate = float(line_s[3])
                        coverage = int(line_s[4])
                        strand = line_s[5]
                        nt = "G" if strand == "-" else "C"

                        if coverage >= min_cov:
                            fout.write('\t'.join([str(x) for x in [chrom, start, end, strand, methyl_rate, coverage, nt]]) + '\n')
                    else:
                        if chrom not in non_std_chroms:
                            non_std_chroms.add(chrom)

    assert len(chroms_seen) == 24

    print('filered bed file with min coverage {} and saved to {}'.format(min_cov, outfile))
    print('Non standard chrs seen : ', non_std_chroms)

def filter_by_chr(infile, chrs = [], fstrand = None, SILENT=False):
    '''
        Filters the bed file by chromosome
        Input bed file : chr, start, end, strand, methylation, coverage
        We can also filter by strand with fstrand ('+' or '-')

        output : List[list]
                 [['chr1', 64863, 64864, '+', 20.0, 5],
                 .
                 .
                 ]
    '''
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

def findNDR_for_chr(filtered_res, chr, win_len = 100, multiprocess=False, temp_loc=None):
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
            assert (start, end) == window
            window_ = (start, end, -np.log10(chi_test.pvalue))
            regions.append(window_)
        
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

        start, end, pval = regions[i]

        if i+1 < n:
            for j in range(i+1, n):
                if regions[j][0] <= end:
                    end = regions[j][1]
                    pval = min(pval, regions[j][2])
                else:
                    i = j
                    break
        
            if j == n-1 and i != j:
                i = n

            merged_regions.append((start, end, pval))
        else:
            i += 1
            assert i == n
            merged_regions.append((start, end, pval))

    if multiprocess:
        file_loc = temp_loc / f'temp.NDR.{chr}.bed'
        with open(file_loc, 'w') as fout:
            for reg in merged_regions:
                temp = [chr, str(reg[0]), str(reg[1]), str(reg[2])]
                fout.write('\t'.join(temp) + '\n')
    else:
        return merged_regions
    
def findNDRs(infile, win_len=200):
    '''
        Input bed file : chr, start, end, strand, methylation, coverage
    '''
    chrs = ['chr' + str(c) for c in range(1, 23)] + ['chrX', 'chrY']
    myprocess= []
    win_len = 200
    for chr in chrs:
        print(f'Run : {chr}')
        res = filter_by_chr(infile, chrs=[chr], SILENT=True)
        print(f'no of reads for {chr} : {len(res)}')
        print(f'{chr} res check -', len(res), res[0], res[-1])

        myprocess.append(multiprocessing.Process(target=findNDR_for_chr, args=(res, chr, win_len, True, infile.parent, )))
        myprocess[-1].start()

    for p in myprocess:
        p.join()

    ndr_regions = []
    for chr in chrs:
        fpath = infile.parent / f'temp.NDR.{chr}.bed'
        with open(fpath, 'r') as fin:
            for line in fin:
                ndr_regions.append(line)

    print(f'no of NDR regions {len(ndr_regions)}')

    outfile = infile.parent / 'NDR.bed'
    with open(outfile, 'w') as fout:
        for reg in ndr_regions:
            fout.write(reg)

    for chr in chrs:
        fpath = infile.parent / f'temp.NDR.{chr}.bed'
        os.remove(fpath)

    print('sorting NDR file')
    bedtools.sort_bed(infile=outfile, outfile=outfile)

def get_methylation(infile, region='region'):
    inter_dict = dict()
    seen_ids = set()

    with open(infile) as f:
        for line in f:
            line_s = line.strip().split('\t')

            chrom = line_s[0]            
            region_start = int(line_s[1])
            region_end = int(line_s[2])
            refid = line_s[3]
            TSS = int(line_s[4])
            TES = int(line_s[5])
            strand = line_s[6]
            meth_start = int(line_s[8])
            meth_end = int(line_s[9])
            meth_rate = float(line_s[11])
            coverage = int(line_s[12])
            nt = line_s[13]

            dict_id = refid + "-" + str(region_start)
            
            ## We are checking + and -ve strands based on promoter
            if strand == "+":
                intron_abs_pos = meth_start - region_start
            if strand == "-":
                intron_abs_pos = region_end - meth_end
        
            rel_pos = intron_abs_pos - 2000
        
            if dict_id not in seen_ids:
                seen_ids.add(dict_id)
                inter_dict[dict_id] = dict()
            
            inter_dict[dict_id][meth_start] = dict()
            
            inter_dict[dict_id][meth_start]["chrom"] = chrom
            inter_dict[dict_id][meth_start][f"{region}_start"] = region_start
            inter_dict[dict_id][meth_start][f"{region}_end"] = region_end
            inter_dict[dict_id][meth_start]["TSS"] = TSS
            inter_dict[dict_id][meth_start]["TES"] = TES
            inter_dict[dict_id][meth_start]["strand"] = strand
            
            inter_dict[dict_id][meth_start]["meth_start_genome"] = meth_start
            inter_dict[dict_id][meth_start]["meth_end_genome"] = meth_end
            inter_dict[dict_id][meth_start]["meth_pos_abs"] = intron_abs_pos
            inter_dict[dict_id][meth_start]["meth_pos_rel"] = rel_pos
                        
            inter_dict[dict_id][meth_start]["meth_rate"] = meth_rate
            inter_dict[dict_id][meth_start]["coverage"] = coverage
            inter_dict[dict_id][meth_start]["nt"] = nt

    column_names = ['trans_id', 'refid', 'chrom', f'{region}_start', f'{region}_end', 
                    'TSS', 'TES', 'strand', 'meth_start_genome', 'meth_end_genome', 'meth_pos_abs',
                    'meth_pos_rel', 'meth_rate', 'coverage', 'nt']

    info_dict = dict()
    for col in column_names:
        info_dict[col] = []

    ## Can create a filter for no of GCH or HCG in a promoter
    for dict_id in inter_dict.keys():
        refid = dict_id.split("-")[0]
        for promo_start in inter_dict[dict_id].keys():
            info_dict["trans_id"].append(dict_id)
            info_dict["refid"].append(refid)
            info_dict["chrom"].append(inter_dict[dict_id][promo_start]["chrom"])
            info_dict[f"{region}_start"].append(inter_dict[dict_id][promo_start][f"{region}_start"])
            info_dict[f"{region}_end"].append(inter_dict[dict_id][promo_start][f"{region}_end"])
            info_dict["TSS"].append(inter_dict[dict_id][promo_start]["TSS"])
            info_dict["TES"].append(inter_dict[dict_id][promo_start]["TES"])
            info_dict["strand"].append(inter_dict[dict_id][promo_start]["strand"])
            
            info_dict["meth_start_genome"].append(inter_dict[dict_id][promo_start]["meth_start_genome"])
            info_dict["meth_end_genome"].append(inter_dict[dict_id][promo_start]["meth_end_genome"])
            info_dict["meth_pos_abs"].append(inter_dict[dict_id][promo_start]["meth_pos_abs"])
            info_dict["meth_pos_rel"].append(inter_dict[dict_id][promo_start]["meth_pos_rel"])
            
            info_dict["meth_rate"].append(inter_dict[dict_id][promo_start]["meth_rate"])
            info_dict["coverage"].append(inter_dict[dict_id][promo_start]["coverage"])
            info_dict["nt"].append(inter_dict[dict_id][promo_start]["nt"])

    # Built dataframe
    df = pd.DataFrame(0, index = np.arange(len(info_dict["refid"])), columns = column_names)
    for feat in column_names:
        df[feat] = info_dict[feat]  
        
    df = df.sort_values(by = ['chrom', f'{region}_start'], ascending = [True, True])

    return df

def get_nuc_positions(infile, region='region'):
    column_names = ["trans_id", "refid", "chrom", f"{region}_start", f"{region}_end",
                    "TSS", "TES", "strand", "nuc_region_start_genome", "nuc_region_end_genome",
                    "nuc_start_promo_abs", "nuc_end_promo_abs", "nuc_start_promo_rel",
                    "nuc_end_promo_rel", "region_length"]
    
    info_dict = dict()
    for col in column_names:
        info_dict[col] = []
    
    with open(infile) as fin:
        for line in fin:
            line_s = line.strip().split('\t')

            chrom = line_s[0]            
            region_start = int(line_s[1])
            region_end = int(line_s[2])
            refid = line_s[3]
            TSS = int(line_s[4])
            TES = int(line_s[5])
            strand = line_s[6]

            nuc_region_start = int(line_s[8])
            nuc_region_end = int(line_s[9])
            

            if nuc_region_start >= region_start and nuc_region_end <= region_end:    
                if strand == "+":
                    nuc_abs_start = nuc_region_start - region_start
                    nuc_abs_end = nuc_region_end - region_start
                    
                if strand == "-":
                    nuc_abs_end = region_end - nuc_region_start
                    nuc_abs_start = region_end - nuc_region_end

                rel_start = nuc_abs_start - 2000
                rel_end = nuc_abs_end - 2000

                info_dict["trans_id"].append(refid+'-'+str(region_start))
                info_dict["chrom"].append(chrom)
                info_dict[f"{region}_start"].append(region_start)
                info_dict[f"{region}_end"].append(region_end)
                info_dict["refid"].append(refid)
                info_dict["TSS"].append(TSS)
                info_dict["TES"].append(TES)
                info_dict["strand"].append(strand)
                
                info_dict["nuc_region_start_genome"].append(nuc_region_start)
                info_dict["nuc_region_end_genome"].append(nuc_region_end)
                
                info_dict["nuc_start_promo_abs"].append(nuc_abs_start)
                info_dict["nuc_end_promo_abs"].append(nuc_abs_end)
                
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
        
    df = df.sort_values(by = ['chrom', f'{region}_start'], ascending=[True, True])
      
    return df

def get_nuc_pos_methylation(infile, region='region'):
    column_names = ["trans_id", "refid", "chrom", f"{region}_start", f"{region}_end",
                    "TSS", "TES", "strand", "nuc_region_start_genome", "nuc_region_end_genome",
                    "nuc_region_length", "meth_start_genome", "meth_end_genome", "meth_rate", "nt"]
    
    info_dict = dict()
    for col in column_names:
        info_dict[col] = []
    
    with open(infile, 'r') as fin:
        for line in fin:
            line_s = line.strip().split('\t')

            chrom = line_s[0]            
            region_start = int(line_s[1])
            region_end = int(line_s[2])
            refid = line_s[3]
            
            TSS = int(line_s[4])
            TES = int(line_s[5])
            strand = line_s[6]
            
            nuc_region_start = int(line_s[8])
            nuc_region_end = int(line_s[9])
            
            meth_start = int(line_s[18])
            meth_end = int(line_s[19])
            meth_rate = float(line_s[21])
            nt = line_s[23]
            
            if nuc_region_start >= region_start and nuc_region_end <= region_end:
                info_dict["trans_id"].append(refid + "-" + str(region_start))
                info_dict["chrom"].append(chrom)
                info_dict[f"{region}_start"].append(region_start)
                info_dict[f"{region}_end"].append(region_end)
                info_dict["refid"].append(refid)
                info_dict["TSS"].append(TSS)
                info_dict["TES"].append(TES)
                info_dict["strand"].append(strand)
                
                info_dict["nuc_region_start_genome"].append(nuc_region_start)
                info_dict["nuc_region_end_genome"].append(nuc_region_end)                
                info_dict["nuc_region_length"].append(nuc_region_end-nuc_region_start)
                
                info_dict["meth_start_genome"].append(meth_start)
                info_dict["meth_end_genome"].append(meth_end)
                info_dict["meth_rate"].append(meth_rate)
                info_dict["nt"].append(nt)

            
    #Built dataframe
    df = pd.DataFrame(0, index = np.arange(len(info_dict["refid"])), columns = column_names)
    for feat in column_names:
        df[feat] = info_dict[feat]  
        
    df = df.sort_values(by=['chrom', f'{region}_start'], ascending=[True, True])

    return df

def make_sliding_windows_file(df_promo_nuc_WGBS, x_y_clash_dict_norm, mask=True, region='region'):
    if mask:
        pos_true = pd.DataFrame([df_promo_nuc_WGBS['meth_start_genome'] > df_promo_nuc_WGBS['nuc_region_start_genome'], 
                                 df_promo_nuc_WGBS['meth_start_genome'] < df_promo_nuc_WGBS['nuc_region_end_genome'], 
                                 df_promo_nuc_WGBS['strand'] == '+']).T.all(axis=1)
        neg_true = pd.DataFrame([df_promo_nuc_WGBS['meth_end_genome'] > df_promo_nuc_WGBS['nuc_region_start_genome'], 
                                 df_promo_nuc_WGBS['meth_end_genome'] < df_promo_nuc_WGBS['nuc_region_end_genome'], 
                                 df_promo_nuc_WGBS['strand'] == '-']).T.all(axis=1)
        df_promo_nuc_WGBS = df_promo_nuc_WGBS[np.logical_or(pos_true, neg_true)]


    column_names = ["trans_id", "refid", "NOR_nbr", "window_nbr", "nbr_meth_CpGs", 
                    "nuc_region_length", "nuc_rel_center", "meth_rates_window"] 
    info_dict = dict()
    for col in column_names:
        info_dict[col] = []
        
    all_nuc_pos = x_y_clash_dict_norm.keys()
    nbr_bases_nuc = len(all_nuc_pos)
    
    all_trans_ids = list(set(list(df_promo_nuc_WGBS["trans_id"])))
    
    c = 0
    for trans_id in tqdm(all_trans_ids):
        refid = trans_id.split("-")[0]
        c += 1
        
        df_WGBS_tmp = df_promo_nuc_WGBS.loc[df_promo_nuc_WGBS["trans_id"] == trans_id]
        
        NOR_number = 1
        nuc_region_starts = list(set(list(df_WGBS_tmp["nuc_region_start_genome"])))
        for NOR_start in nuc_region_starts:
            df_NOR_tmp = df_WGBS_tmp.loc[df_WGBS_tmp["nuc_region_start_genome"] == NOR_start]
            
            start = NOR_start
            win_nbr = 1
                    
            list_meth_positions = list(df_NOR_tmp["meth_start_genome"])
            list_meth_rates = list(df_NOR_tmp["meth_rate"])
            meth_dict_promoter = dict(zip(list_meth_positions, list_meth_rates))
    
            NOR_end = list(df_NOR_tmp["nuc_region_end_genome"])[0]

            assert df_NOR_tmp[f"{region}_start"].unique().shape[0] == 1
            region_start = df_NOR_tmp[f"{region}_start"].unique()[0]
            assert df_NOR_tmp[f"{region}_end"].unique().shape[0] == 1
            region_end = df_NOR_tmp[f"{region}_end"].unique()[0]
            assert df_NOR_tmp[f"strand"].unique().shape[0] == 1
            strand = df_NOR_tmp[f"strand"].unique()[0]

            if NOR_start >= region_start and NOR_end <= region_end:    
                if strand == "+":
                    nuc_abs_start = NOR_start - region_start
                    nuc_abs_end = NOR_end - region_start
                    
                if strand == "-":
                    nuc_abs_end = region_end - NOR_start
                    nuc_abs_start = region_end - NOR_end

                rel_start = nuc_abs_start - 2000
                rel_end = nuc_abs_end - 2000
            else:
                raise Exception("problem")
            
            while start + nbr_bases_nuc - 1 <= NOR_end:
                #get relative and absolute position+methylation
                window_start = start
                window_end = window_start+nbr_bases_nuc-1
                
                meth_in_window_tmp = dict()
                for meth_start_genome in meth_dict_promoter.keys():
                    if (window_start <= meth_start_genome) & (meth_start_genome <= window_end):
                        rel_pos = meth_start_genome-window_start+1
                        
                        meth_rate = meth_dict_promoter[meth_start_genome]
                        meth_in_window_tmp[rel_pos] = meth_rate
                
                if len(meth_in_window_tmp) != 0:
                    #exp
                    info_dict["trans_id"].append(trans_id)
                    info_dict["refid"].append(refid)
                    info_dict["NOR_nbr"].append(NOR_number)
                    info_dict["window_nbr"].append(win_nbr)
                    info_dict["nbr_meth_CpGs"].append(len(meth_in_window_tmp.keys()))
                    info_dict["nuc_region_length"].append(np.abs(NOR_end-NOR_start)+1)
                    info_dict["nuc_rel_center"].append((rel_start+rel_end)/2)
                    info_dict["meth_rates_window"].append(meth_in_window_tmp)
                
                start += 1
                win_nbr += 1
        
            NOR_number += 1
                
                
    #Built dataframe
    df = pd.DataFrame(0, index = np.arange(len(info_dict["trans_id"])),columns = column_names)
    for feat in column_names:
        df[feat] = info_dict[feat]

    return df