#PBS -l walltime=999:00:00
#PBS -l nodes=2:ppn=20,mem=32gb
#PBS -q cv2
#PBS -j oe
    /home/jia_wenqi/anaconda3/bin/fastp -i DP8480008987TL_L01_7268_1.part_001.fq.gz -I DP8480008987TL_L01_7268_2.part_001.fq.gz -o DP8480008987TL_L01_7268_1.fp.fq.gz -O DP8480008987TL_L01_7268_2.fp.fq.gz -g --thread=10 --length_required=10 --n_ba
#mapping reads
##去除read长度低于10，大于700bp，精确比对，双端必须匹配，不保留unmapping reads
align(){
sambamba=/home/yuan_jianwen/anaconda3/bin/sambamba
spe=mm10
prefix=`pwd`
out1=$prefix/bam
out2=$prefix/macs2
out3=$prefix/bw
genome=/mnt/1/genome/mm10/bt2_index/mm10

root=`basename $1`
echo $root
    bf=`echo $root | sed -r 's/\_1.fastq.fp.fq.gz//g' | sed -r 's/.gz//g'`
echo $bf
    p1=`echo $root`
echo $p1
    p2=`echo $root | sed 's/\_1.fastq.fp.fq.gz/\_2.fastq.fp.fq.gz/g'`
echo $p2
echo PE...mapping $bf".sam" $bf"_q30.bam" $bf"rmdup_q30.bam"



###ATAC###
bowtie2 \
-p 30 \
--very-sensitive \
--end-to-end \
--no-unal \
-x $genome \
-1 $prefix/$p1 \
-2 $prefix/$p2 \
--phred33 \
2> $out1/$bf".mapping_summary" \
-S $out1/$bf".sam"


###CUTTag###
bowtie2 \
-p 10 \
--local \
--very-sensitive-local \
--no-mixed \
--no-discordant \
--no-unal \
-x $genome \
-1 $prefix/$p1 \
-2 $prefix/$p2 \
-I 10 \
-X 700 \
--phred33 \
2> $out1/$bf".mapping_summary" \
-S $out1/$bf".sam"

#生成bam文件统计信息
#samtools flagstat $out1/$bf".sam" >> $out1/$bf".mapping_summary"
#去除低质量mapping reads
samtools view -u -bS -F 0x04 -q 30 $out1/$bf".sam" > $out1/$bf"_q30.bam"
#去重
$sambamba markdup -t 10 -r $out1/$bf"_q30.bam" $out1/$bf"rmdup_q30.bam" 2>>$out1/$bf".mapping_summary"
#bam文件重排，用于后续分析
samtools sort $out1/$bf"rmdup_q30.bam" > $out1/$bf"rmdupchrm_q30_sorted.bam"
samtools index $out1/$bf"rmdupchrm_q30_sorted.bam"
#删除中间文件
rm $out1/$bf".sam" $out1/$bf"_q30.bam" $out1/$bf"rmdup_q30.bam"
}



##去除read长度低于10，大于700bp，精确比对，双端必须匹配，不保留unmapping reads
align() {
sambamba=/home/yuan_jianwen/anaconda3/bin/sambamba
spe=mm10
prefix=`pwd`
out1=$prefix/bam
out2=$prefix/macs2
out3=$prefix/bw
genome=/mnt/1/genome/mm10/bt2_index/mm10


root=`basename $1`
echo $root
    bf=`echo $root | sed -r 's/\_1.fp.fq//g' | sed -r 's/.gz//g'`
echo $bf
    p1=`echo $root`
echo $p1
    p2=`echo $root | sed 's/\_1/\_2/g'`
echo $p2
echo PE...mapping

bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -p 20 --phred33 -I 10 -X 700 -x $genome -1  $prefix/$p1 -2  $prefix/$p2 2> $out1/$bf".mapping_summary"  | samtools view -bS -q 30 - > $out1/$bf"_q30.bam"
$sambamba markdup -t 10 -r $out1/$bf"_q30.bam" $out1/$bf"rmdup_q30.bam" 2>>$out1/$bf".mapping_summary"
samtools sort $out1/$bf"rmdup_q30.bam" > $out1/$bf"rmdup_sorted_q30.bam"
samtools index $out1/$bf"rmdup_sorted_q30.bam"
#删除中间文件
rm $out1/$bf"_q30.bam" $out1/$bf"rmdup_q30.bam"
}

export -f align
ls ./*_1.fp.fq.gz | cat | parallel --tmpdir `pwd` --jobs 2 align


#统计fragment分部情况
docker run --tty --interactive --mount src=/mnt,target=/mnt,type=bind --workdir="${PWD}" broadinstitute/gatk bash
for f in $out1/*rmdup_q30.bam
do
    root=`basename $f`
    echo $root"...Count_Fragment"
    bf=`echo $root | sed -r 's/rmdup_q30.bam//g'`
gatk CollectInsertSizeMetrics INPUT=${out1}/${bf}rmdup_q30.bam OUTPUT=${out1}/${bf}rmdup_q30.bam.fragment.txt Histogram_FILE=${out1}/${bf}rmdup_q30.bam.pic METRIC_ACCUMULATION_LEVEL=ALL_READS
done


#BAM生成BW
docker run --tty --interactive --mount src=/mnt,target=/mnt,type=bind --workdir="${PWD}" ubuntu:18.04_atac bash

for f in $out1/*rmdupchrm_q30_sorted.bam
for f in DP8480008980TL_L01_7268_rmdup_sorted_q30.bam
do
root=`basename $f`
/home/yuan_jianwen/anaconda3/bin/bamCoverage -b $f \
--normalizeUsing RPKM \
--effectiveGenomeSize 2652783500 \
--minMappingQuality 30 \
--numberOfProcessors 5 \
--ignoreDuplicates \
-o /mnt/12/wu_xinyu/20241002_4CL_8CLC_primed_H3K79me3/bam/$root"cpm.bw" &
done 

--effectiveGenomeSize 2652783500 \

#call peak
for f in $out1/*rmdupchrm_q30_sorted.bam
for f in *sorted.bam
do
root=`basename $f`
echo $root"...Calling peak"
bf=`echo $root | sed -r 's/.bam/sort_all/g'`
echo $bf

macs2 callpeak \
-n $bf \
-t $f \
-f BAMPE \
-B -g mm \
--outdir $out2 \
--SPMR \
--nomodel --keep-dup 1 \
--extsize 200 \
-q 0.00001 &
done


#去除peak中blacklist region
for f in *peaks.narrowPeak
do
root=`basename $f`
bf=`echo $root | sed -r 's/_peaks.narrowPeak//g'`
bedtools intersect -a $f -b /mnt/6/WXY/mm10.blacklist.bed -v > $bf"_deblack.peak"
done

for f in *cleanchr.bed
do
root=`basename $f`
bf=`echo $root | sed -r 's/.bed//g'`
echo $bf
bedtools intersect -a $f -b /mnt/6/WXY/mm10.blacklist.bed -v > $bf"_deblack.peak"
done


#生成logLR.bw
for f in $out2/*_treat_pileup.bdg
do
root=`basename $f`
echo $root"...Doing bdgcmp"
bf=`echo $root | sed -r 's/_treat_pileup.bdg//g'`
FE=`echo $bf'_FE.bdg'`
logLR=`echo $bf'_logLR.bdg'`
control_lambda=`echo $bf'_control_lambda.bdg' `
echo $bf $control_lambda
macs2 bdgcmp \
-t $f \
-c $out2/$control_lambda \
--outdir $out3 \
-o $logLR \
-m logLR \
-p 0.01 &
done

wait
chrsize=/mnt/6/WXY/mm10_chromInfo.txt

for f in $out3/*_logLR.bdg
do
    root=`basename $f`
    echo $root"...Doing bdg2bw"
    bf=`echo $root | sed -r 's/_FE.bdg//g' | sed -r 's/_logLR.bdg//g'`
    bash /mnt/11/WXY_classic_data/download/bin/bdg2bw.sh $f $chrsize 
done
wait
rm $out3/*_logLR.bdg


