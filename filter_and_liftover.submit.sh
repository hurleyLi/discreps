#!/bin/bash

#PBS -N filter_and_liftover
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -q analysis
#PBS -d /discreps
#PBS -e /discreps/log/filter.err
#PBS -o /discreps/log/filter.out
#PBS -t 1-300

# This script is used to submit jobs for variant filtering and variant lift-over
# The project-level .bcf files are divided into 300 chunks to enable parallele computing
# The code has been refactored. Please submit an issue or contact the author if there is problem or question

baseDir=/discreps
max_allele_INDEL=3
max_allele_SNV=4
total_sample=1572
version=$1 # specify version of the reference [hg19|hg38]
assignMissing=/utils/assignMissing.py

if [[ "$version" == "hg19" ]]; then
    reference=/reference/hs37d5/hs37d5.fa # NOT included in github submission
    reference_other=/reference/GRCh38_1000Genomes/reference.fa # NOT included in github submission
    chain_file=/reference/hg37ToHg38.over.chain # downloaded from UCSC genome browser
    exon_capture=/reference/HG19_vcrome2.1_with_PKv2.100bp.bed.gz
    chunkfile=/reference/reference_hg19.auto.chunks
elif [[ "$version" == "hg38" ]]; then
    reference=/reference/GRCh38_1000Genomes/reference.fa # NOT included in github submission
    reference_other=/reference/hs37d5/hs37d5.fa # NOT included in github submission
    chain_file=/reference/hg38ToHg37.over.chain # downloaded from UCSC genome browser
    exon_capture=/reference/HG38_vcrome2.1_with_PKv2.100bp.bed.gz
    chunkfile=/reference_hg38.auto.chunks
else
    echo "Specify either hg19 or hg38 as version"; exit 1
fi

# read chromosome position of each chunk 
chunk=$PBS_ARRAYID
line=`sed -n "${PBS_ARRAYID}p" $chunkfile`
line=($line)
chr=${line[0]}
start=${line[1]}
end=${line[2]}

# choose filtering criteria
dps=(15) # minimum depth
ratios=(0.25) # minimum minor allelic ratio
genos=(0.85) # minimum genotyping rate across all samples

# We also tried other different combinations of QC criteria
# dps=(5 10 15 20 25)
# ratios=(0.15 0.2 0.25 0.3 0.35)
# genos=(0.75 0.8 0.85 0.9 0.95)

for type in SNV INDEL; do
    prefix=$version.$type
    input=/inputs/$prefix.bcf # obtain from dbGaP
    workingDir=/by_chunk/$chunk
    chr_vcf=$workingDir/$prefix.vcf.gz # vcf file separated by each chr chunk
    mkdir -p $workingDir/ID $workingDir/assignMissing
    cd $workingDir
    echo -e "chunk: $chunk, $chr:$start-$end" > chunk
    
    # different filtering criteria for SNV and INDEL for multi-allelic variants
    if [[ $type == "SNV" ]]; then
        max_allele=$max_allele_SNV
    else
        max_allele=$max_allele_INDEL
    fi

    # 1. separate the file into chunk
    bcftools view -r $chr:$start-$end -M $max_allele -f '.' $input \
    | awk -v start=$start -v end=$end '$1~/^#/ || ($2>=start && $2<=end)' \
    | bgzip -c > $chr_vcf
    
    varCount=`bcftools stats $chr_vcf | grep "number of records:" | awk -v FS="\t" '{print $4}'`
    echo "varCount: $varCount" >> chunk
    
    # determine whether there is variant within this chunk
    if [[ $varCount -eq 0 ]]; then touch warning.noVariants; fi
    if [[ ! -f warning.noVariants ]]; then
        tabix -p vcf $chr_vcf
        bcftools view -R $exon_capture $chr_vcf \
        | bgzip -c > $chr_vcf.2
        mv $chr_vcf.2 $chr_vcf
        rm $chr_vcf.tbi
        
        varCount_exon=`bcftools stats $chr_vcf | grep "number of records:" | awk -v FS="\t" '{print $4}'`
        echo "varCount_exon: $varCount_exon" >> chunk
        
        # determine whether there is variant in exon capture regions within this chunk
        if [[ $varCount_exon -eq 0 ]]; then touch warning.noExonVariants; fi
        if [[ ! -f warning.noExonVariants ]]; then
            tabix -p vcf $chr_vcf
        fi
    fi

    if [[ ! -f warning.noVariants ]] && [[ ! -f warning.noExonVariants ]]; then
        # 2. liftover
        java -Xmx7500m -Xms7500m \
        -jar picard.jar LiftoverVcf \
        I=$chr_vcf \
        O=$prefix.lifted.vcf \
        CHAIN=$chain_file \
        REJECT=$prefix.lifted.reject.vcf \
        RECOVER_SWAPPED_REF_ALT=true \
        TAGS_TO_DROP=AC \
        R=$reference_other
        
        bcftools plugin fill-tags $prefix.lifted.vcf -- -t AC \
        | vcf-sort -c | bgzip -c > $prefix.lifted.vcf.gz
        tabix -p vcf $prefix.lifted.vcf.gz
        rm $prefix.lifted.vcf $prefix.lifted.vcf.idx
        
        rejectCount=`bcftools stats $prefix.lifted.reject.vcf | grep "number of records:" | awk -v FS="\t" '{print $4}'`
        if [[ $rejectCount -eq 0 ]]; then rm $prefix.lifted.reject.vcf; fi
        
        # 3. assign missing variants based on chosen criteria
        # different QC criteria were tested, but only use was finally chosen and used in this study
        for depth_thres in ${dps[@]}; do for allelic_ratio_thres in ${ratios[@]}; do for genotypingRateThres in ${genos[@]}; do
            total_AN=`echo "$total_sample * 2 * $genotypingRateThres" | bc -l`
            python $assignMissing -i $chr_vcf -d $depth_thres -r $allelic_ratio_thres \
            | bcftools plugin fill-tags -- -t AN,AF,AC \
            | bcftools norm -c wx -f $reference -m- -Ou \
            | bcftools norm -d none -Ou \
            | bcftools view -i "AN>$total_AN & AC>0" -Ou \
            | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' \
            | bgzip -c > $prefix.filtered.vcf.gz
            tabix -p vcf $prefix.filtered.vcf.gz
        done; done; done
    fi
done

