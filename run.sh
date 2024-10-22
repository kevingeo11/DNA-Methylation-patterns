#!/bin/bash

usage() {
    echo "Usage: $0 --download"
    exit 1
}

DOWNLOAD=false
if [ "$1" == "--download" ]; then
    DOWNLOAD=true
fi

if [ "$DOWNLOAD" = true ]; then
    # Check if the credentials file exists
    if [ ! -f Data/credentials.json ]; then
        echo "Error: credentials.json file not found in the Data directory."
        exit 1
    fi
fi

# Create the output directory if it doesn't exist
if [ ! -d Data ]; then
    mkdir -p Data
fi

if [ "$DOWNLOAD" = true ]; then
    # Run the pyega3 fetch command
    # If you run into recurring md5sum error
    # Change according to https://github.com/EGA-archive/ega-download-client/issues/192#issuecomment-1692123299
    pyega3 -c 30 -cf Data/credentials.json fetch EGAF00001273727 --output-dir Data/ --max-retries -1
    pyega3 -c 30 -cf Data/credentials.json fetch EGAF00001273728 --output-dir Data/ --max-retries -1
fi

# Trim galore version 0.6.10
# https://bio-protocol.org/exchange/minidetail?type=30&id=9632690&utm_source=miniprotocol
if [ ! -d Data/trim ]; then
    mkdir -p Data/trim
fi
# trim_galore --cores 4 --clip_R1 9 --clip_R2 9 --paired -o Data/trim/ Data/EGAF00001273727/01_HepG2_LiHG_Ct1_NOMe_S_1_ACAGTG_SN7001355_0058_AC2L6FACXX_L007_R1.fastq.gz Data/EGAF00001273728/01_HepG2_LiHG_Ct1_NOMe_S_1_ACAGTG_SN7001355_0058_AC2L6FACXX_L007_R2.fastq.gz

# Bowtie 2 version 2.5.1
# Bismark Version: v0.24.2
# packages/Bismark-0.24.2/bismark_genome_preparation --parallel 16 --verbose Data/GRCh38genome
if [ ! -d Data/bismark ]; then
    mkdir -p Data/bismark
fi
# packages/Bismark-0.24.2/bismark --fastq --parallel 6 --output_dir Data/bismark --genome Data/GRCh38genome -1 Data/trim/01_HepG2_LiHG_Ct1_NOMe_S_1_ACAGTG_SN7001355_0058_AC2L6FACXX_L007_R1_val_1.fq.gz -2 Data/trim/01_HepG2_LiHG_Ct1_NOMe_S_1_ACAGTG_SN7001355_0058_AC2L6FACXX_L007_R2_val_2.fq.gz

# Adding Read Groups using picard tools
if [ ! -d Data/bissnp ]; then
    mkdir -p Data/bissnp
fi
# packages/jdk-17.0.8/bin/java -jar packages/picard.jar AddOrReplaceReadGroups \
# -I Data/bismark/01_HepG2_LiHG_Ct1_NOMe_S_1_ACAGTG_SN7001355_0058_AC2L6FACXX_L007_R1_val_1_bismark_bt2_pe.bam \
# -O Data/bissnp/01_HepG2_LiHG_Ct1_NOMe_S_1_ACAGTG_SN7001355_0058_AC2L6FACXX_L007_R1_val_1_bismark_bt2_pe.withRG.bam \
# --RGPL illumina --RGLB lib --RGPU run --RGSM name --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --SORT_ORDER coordinate

# Mark Duplicates using picard tools
# packages/jdk-17.0.8/bin/java -jar packages/picard.jar MarkDuplicates \
# I=Data/bissnp/01_HepG2_LiHG_Ct1_NOMe_S_1_ACAGTG_SN7001355_0058_AC2L6FACXX_L007_R1_val_1_bismark_bt2_pe.withRG.bam \
# O=Data/bissnp/01_HepG2_LiHG_Ct1_NOMe_S_1_ACAGTG_SN7001355_0058_AC2L6FACXX_L007_R1_val_1_bismark_bt2_pe.withRG.mdups.bam \
# M=Data/bissnp/01_HepG2_LiHG_Ct1_NOMe_S_1_ACAGTG_SN7001355_0058_AC2L6FACXX_L007_R1_val_1_bismark_bt2_pe.withRG.metric.txt \
# CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT


# packages/jdk-17.0.8/bin/java -jar packages/BisSNP-0.82.2.jar -R Data/GRCh38genome/GRCh38.p14.genome.fa \
# -I Data/bissnp/01_HepG2_LiHG_Ct1_NOMe_S_1_ACAGTG_SN7001355_0058_AC2L6FACXX_L007_R1_val_1_bismark_bt2_pe.withRG.mdups.bam \
# -T BisulfiteCountCovariates -knownSites Data/dbsnp/GCF_000001405.40.4.1.annotate.vcf \
# -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -recalFile Data/bissnp/recalFile.csv -nt 16

# packages/jdk-17.0.8/bin/java -jar packages/BisSNP-0.82.2.jar -R Data/GRCh38genome/GRCh38.p14.genome.fa \
# -I Data/bissnp/01_HepG2_LiHG_Ct1_NOMe_S_1_ACAGTG_SN7001355_0058_AC2L6FACXX_L007_R1_val_1_bismark_bt2_pe.withRG.mdups.bam \
# -o Data/bissnp/01_HepG2_LiHG_Ct1_NOMe_S_1_ACAGTG_SN7001355_0058_AC2L6FACXX_L007_R1_val_1_bismark_bt2_pe.withRG.mdups.recal.bam \
# -T BisulfiteTableRecalibration -recalFile Data/bissnp/recalFile.csv -maxQ 40

# packages/jdk-17.0.8/bin/java -jar packages/BisSNP-0.82.2.jar -R Data/GRCh38genome/GRCh38.p14.genome.fa -T BisulfiteGenotyper \
# -I Data/bissnp/01_HepG2_LiHG_Ct1_NOMe_S_1_ACAGTG_SN7001355_0058_AC2L6FACXX_L007_R1_val_1_bismark_bt2_pe.withRG.mdups.recal.bam \
# -D Data/dbsnp/GCF_000001405.40.4.1.annotate.vcf -vfn1 Data/bissnp/cpg.raw.vcf -vfn2 Data/bissnp/snp.raw.vcf \
# -stand_call_conf 20 -stand_emit_conf 0 -mmq 30 -mbq 0 -sm GM -out_modes EMIT_VARIANT_AND_CYTOSINES -nt 8

# perl packages/sortByRefAndCor.pl --k 1 --c 2 Data/bissnp/cpg.raw.vcf Data/GRCh38genome/GRCh38.p14.genome.fa.fai > Data/bissnp/cpg.raw.sort.vcf
# perl packages/sortByRefAndCor.pl --k 1 --c 2 Data/bissnp/snp.raw.vcf Data/GRCh38genome/GRCh38.p14.genome.fa.fai > Data/bissnp/snp.raw.sort.vcf

# perl packages/vcf2bed.pl Data/bissnp/cpg.raw.sort.vcf GCH
# perl packages/vcf2bed.pl Data/bissnp/cpg.raw.sort.vcf HCG