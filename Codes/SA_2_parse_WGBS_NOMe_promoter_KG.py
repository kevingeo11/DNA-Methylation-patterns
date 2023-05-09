import gzip
import csv

def filter_exp_files(infile, outfile):
    chroms_seen = set()
    non_std_chroms = set()
    itr = 0
    chroms = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']

    with gzip.open(infile, 'rt') as fin:
        with open(outfile, 'w', newline='') as fout:
            writer = csv.writer(fout)
            header = ['chrom', 'start', 'end', 'strand', 'Methyl rate', 'Nucleotide', 'Coverage', 'Is Methylated']
            writer.writerow(header)
            for line in fin:
                itr += 1
                if itr%1000000 == 0:
                    print('processed', itr, 'lines')

                if line.startswith('track'): print(line)

                if not line.startswith('track'):
                    line_s = line.strip().split('\t')
                    assert len(line_s) == 8

                    chrom = 'chr' + line_s[0]
                    if chrom in chroms:
                        chroms_seen.add(chrom)

                        start = int(line_s[1])
                        end = int(line_s[2])
                        methyl_rate = float(line_s[6])
                        is_methylated = 1 if methyl_rate > 0 else 0
                        coverage = int(line_s[7])
                        strand = line_s[5]
                        nt = 'G' if strand == '-' else 'C'

                        if coverage >= 3:
                            writer.writerow([chrom, start, end, strand, methyl_rate, nt, coverage, is_methylated])
                    else:
                        if chrom not in non_std_chroms:
                            print('Some non standard chromosome seen', chrom)
                            non_std_chroms.add(chrom)

    assert len(chroms_seen) == 24

def parse_UCSC_file(infile, outfile):
    seen_coords = set()
    itr = 0

    with gzip.open(infile, 'rt') as fin:
        csvfile = csv.reader(fin)
        with open(outfile, 'w', newline='') as fout:
            writer = csv.writer(fout)
            header = ['chrom', 'prom start', 'prom end', 'ID', 'name', 'tx start', 'tx end', 'strand']
            writer.writerow(header)
            for line in csvfile:
                itr += 1

                if not line[0].startswith('#'):

                    refid = line[1]
                    chrom = line[2]
                    strand = line[3]
                    txStart = int(line[4])
                    txEnd = int(line[5])
                    cdsStart = int(line[6])
                    cdsEnd = int(line[7])
                    geneName = ''.join(c for c in line[12] if c.isalnum())

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

                        if cdsStart != cdsEnd and geneName[0:3] != 'MIR' and geneName[0:3] != 'SNO' and '_' not in chrom:
                            writer.writerow([chrom, promoter_start, promoter_end, refid, geneName, txStart, txEnd, strand])


main_path = 'D:\\Work\\Helms-Lab\\DNA-Methylation-patterns\\'
data_path = main_path + 'Data\\'
outpath_processing = main_path + 'Data\\Filtered Data\\'

'''
HCG, GCH
'''
WGBS_outfile = outpath_processing + "WGBS_filtered.csv"
NOMe_outfile = outpath_processing + "NOMe_filtered.csv"

if False:
    print('WGBS')
    WGBS_file = data_path + 'Human Data\\' + 'GSM5695527_IMR90_bT_WGBS_rep1_final.b37.calmd.cytosine.filtered.sort.HCG.strand.6plus2.bed.gz'
    filter_exp_files(WGBS_file, WGBS_outfile)
    
if False:
    print('NOMe')
    NOMe_file = data_path + 'Human Data\\' + 'GSM5695527_IMR90_bT_WGBS_rep1_final.b37.calmd.cytosine.filtered.sort.GCH.strand.6plus2.bed.gz'
    filter_exp_files(NOMe_file, NOMe_outfile)
    
'''
Promoter
'''    
if False:
    UCSC_file = data_path + 'refSeq h19\\' + '0-refGene_complete_hg19.csv.gz'
    outfile_UCSC_promo = outpath_processing + '0-refGene_complete_hg19_promoter.csv'
    parse_UCSC_file(UCSC_file,outfile_UCSC_promo)

    
    
    