import subprocess

def sort_bed(infile, outfile=None):
    '''
        sort -k1,1V -k2,2n /loc/GCH.filtered.bed > /loc/sorted.bed
    '''
    res = subprocess.run(['sort', '-k1,1V', '-k2,2n', f'{infile}'], capture_output=True, text=True)

    if outfile is None:
        outfile = infile.with_stem(infile.stem + '.sorted')
        
    if res.returncode == 0:
        with open(outfile, 'w') as fout:
            fout.write(res.stdout)
    else:
        print('error in sorting')

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

def complement_bed(infile, genome_sizes, outfile):
    '''
        bedtools complement -i fin -g genome_size
    '''
    res = subprocess.run(['bedtools', 'complement', '-i', f'{infile}', '-g', f'{genome_sizes}'],
                          capture_output=True, text=True)

    if res.returncode == 0:
        with open(outfile, 'w') as fout:
            fout.write(res.stdout)
    else:
        print('error in complement')