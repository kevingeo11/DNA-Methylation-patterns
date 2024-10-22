import numpy as np

def normalize_clash_dict(x_y_clash_dict):
    min_clash = min(x_y_clash_dict.values())
    max_clash = max(x_y_clash_dict.values())

    x_y_clash_dict_norm = {k: (float(v-min_clash)/float(max_clash-min_clash))*100 for k, v in x_y_clash_dict.items()}
    return x_y_clash_dict_norm

def get_promoters_refGene(infile, outfile):
    '''
        get promoter region based on transcript
    '''
    seen_coords = set()
    itr = 0
    with open(outfile, 'w') as fout:
        with open(infile, 'r') as fin:
            for line in fin:
                if not line.startswith('#'):
                    line_s = line.strip().split('\t')

                    refid = line_s[0]
                    chrom = line_s[1]
                    strand = line_s[2]
                    txStart = int(line_s[3])
                    txEnd = int(line_s[4])
                    cdsStart = int(line_s[5])
                    cdsEnd = int(line_s[6])

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

                        if cdsStart != cdsEnd and '_' not in chrom:
                            fout.write('\t'.join([str(s) for s in [chrom, promoter_start, promoter_end, refid, txStart, txEnd, strand]]) + '\n')
                            itr += 1

    print('no of promoters defined :', itr)
    print('promoters saved to {}'.format(outfile))

def get_introns_refGene(infile, outfile, anchor=1, pos='start'):
    '''
        get intron regions based on transcript
        anchor : the position of intron (int)
        pos : 'start' or 'end'
    '''
    seen_coords = set()
    itr = 0
    with open(outfile, 'w') as fout:
        with open(infile, 'r') as fin:
            for line in fin:
                if not line.startswith('#'):
                    line_s = line.strip().split('\t')

                    refid = line_s[0]
                    chrom = line_s[1]
                    strand = line_s[2]
                    txStart = int(line_s[3])
                    txEnd = int(line_s[4])
                    cdsStart = int(line_s[5])
                    cdsEnd = int(line_s[6])

                    coords = chrom + str(txStart) + str(txEnd) + strand
                    if coords not in seen_coords:
                        seen_coords.add(coords)

                        if cdsStart != cdsEnd and '_' not in chrom:
                            assert line_s[8].endswith(',')
                            assert line_s[9].endswith(',')

                            exon_starts = np.array(list(map(int, line_s[8][:-1].strip().split(','))))
                            exon_ends = np.array(list(map(int, line_s[9][:-1].strip().split(','))))

                            assert exon_starts.shape == exon_ends.shape

                            intron_starts = exon_ends[:-1].copy()
                            intron_ends = exon_starts[1:].copy()

                            assert np.all((intron_ends - intron_starts) >= 0)

                            for i, (start, end) in enumerate(zip(intron_starts, intron_ends)):
                                if i+1 == anchor:
                                    if end > start:
                                        start_offset = 2000
                                        end_offset = 1000

                                        if pos == 'start':
                                            if strand == "+":
                                                intron_start = start - start_offset
                                                intron_end = start + end_offset
                                            else:
                                                intron_start = end - end_offset
                                                intron_end = end + start_offset
                                        
                                        elif pos == 'end':
                                            if strand == "+":
                                                intron_start = end - start_offset
                                                intron_end = end + end_offset
                                            else:
                                                intron_start = start - end_offset
                                                intron_end = start + start_offset
                                            
                                        if intron_start < 0:
                                                intron_start = 0

                                        fout.write('\t'.join([str(s) for s in [chrom, intron_start, intron_end, refid, txStart, txEnd, strand]]) + '\n')
                                        itr += 1

                                    break

    print('no of introns defined :', itr)
    print('introns saved to {}'.format(outfile))