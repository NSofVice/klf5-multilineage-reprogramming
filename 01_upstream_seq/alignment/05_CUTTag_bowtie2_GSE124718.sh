#PBS -l walltime=999:00:00
#PBS -l nodes=2:ppn=20,mem=32gb
#PBS -q cv2
#PBS -j oe

#mapping reads
##ȥ��read���ȵ���10������700bp����ȷ�ȶԣ�˫�˱���ƥ�䣬������unmapping reads
align(){
sambamba=/home/yuan_jianwen/anaconda3/bin/sambamba
prefix=`pwd`
out1=$prefix/bam
out2=$prefix/macs2
out3=$prefix/bw
genome=/mnt/1/genome/mm10/bt2_index/mm10

mkdir $out1 $out2 $out3


root=`basename $1`
echo $root
    bf=`echo $root | sed -r 's/\_1.fp.fq//g' | sed -r 's/.gz//g'`
echo $bf
    p1=`echo $root`
echo $p1
    p2=`echo $root | sed 's/\_1.fp.fq/\_2.fp.fq/g'`
echo $p2
echo PE...mapping $bf".sam" $bf"_q30.bam" $bf"rmdup_q30.bam"


bowtie2 \
-p 10 \
-N 1 \
-L 25 \
--end-to-end \
--very-sensitive \
--no-unal \
-x $genome \
-U $prefix/$root \
-I 10 \
-X 1000 \
--phred33 \
2> $out1/$bf".mapping_summary" | samtools view -u -bS -F 0x04 -q 20 - > $out1/$bf"_q20.bam"


#����bam�ļ�ͳ����Ϣ
#samtools flagstat $out1/$bf".sam" >> $out1/$bf".mapping_summary"
#ȥ��������mapping reads
#ȥ��
$sambamba markdup -t 10 -r $out1/$bf"_q20.bam" $out1/$bf"rmdup_q20.bam" 2>>$out1/$bf".mapping_summary"
#bam�ļ����ţ����ں�������
samtools sort $out1/$bf"rmdup_q20.bam" > $out1/$bf"rmdupchrm_q20_sorted.bam"
samtools index $out1/$bf"rmdupchrm_q20_sorted.bam"
#ɾ���м��ļ�
rm $out1/$bf"_q20.bam" $out1/$bf"rmdup_q20.bam"
}

bowtie2 \
-p 10 \
-N 1 \
-L 25 \
--end-to-end \
--very-sensitive \
--no-mixed \
--no-discordant \
-x $genome \
-1 $prefix/$p1 \
-2 $prefix/$p2 \
-I 10 \
-X 1000 \
--phred33 \
2> $out1/$bf".mapping_summary" | samtools view -u -bS -F 0x04 -q 20 - > $out1/$bf"_q20.bam"


##ȥ��read���ȵ���10������700bp����ȷ�ȶԣ�˫�˱���ƥ�䣬������unmapping reads
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
#ɾ���м��ļ�
rm $out1/$bf"_q30.bam" $out1/$bf"rmdup_q30.bam"
}

export -f align
ls ./*.fp.fq.gz | cat | parallel --tmpdir `pwd` --jobs 5 align


#ͳ��fragment�ֲ����
docker run --tty --interactive --mount src=/mnt,target=/mnt,type=bind --workdir="${PWD}" broadinstitute/gatk bash
for f in $out1/*rmdup_q30.bam
do
    root=`basename $f`
    echo $root"...Count_Fragment"
    bf=`echo $root | sed -r 's/rmdup_q30.bam//g'`
gatk CollectInsertSizeMetrics INPUT=${out1}/${bf}rmdup_q30.bam OUTPUT=${out1}/${bf}rmdup_q30.bam.fragment.txt Histogram_FILE=${out1}/${bf}rmdup_q30.bam.pic METRIC_ACCUMULATION_LEVEL=ALL_READS
done


#BAM����BW
docker run --tty --interactive --mount src=/mnt,target=/mnt,type=bind --workdir="${PWD}" ubuntu:18.04_atac bash

for f in DP8480008980TL_L01_7268_rmdup_sorted_q30.bam
for f in *.bam
do
root=`basename $f`
/home/yuan_jianwen/anaconda3/bin/bamCoverage -b $f \
--normalizeUsing RPKM \
--effectiveGenomeSize 2652783500 \
--numberOfProcessors 30 \
--ignoreDuplicates \
-o $out3/$root"RPKM.bw" 
done 

--effectiveGenomeSize 2652783500 \

#call peak
out1=/mnt/2/wuliang/05.Project/00.TGFB_project/00.Dataset/HistoneModification/bam/
out2=/mnt/6/WXY/GSE124718_human_early_embryo_H3K4me3_K27me3_K27ac/8CLC_histone_modification/peak/
for f in *H3K27me3*.bam
do
root=`basename $f`
echo $root"...Calling peak"
bf=`echo $root | sed -r 's/.fastq.fp.fq//g'`
echo $bf

macs2 callpeak \
-n $bf \
-t $f \
-f BAM \
-g hs \
--outdir $out2 \
--SPMR \
--broad --broad-cutoff 0.05 \
--nomodel --nolambda &
done

for f in $out1/*H3K27ac*.bam
do
root=`basename $f`
echo $root"...Calling peak"
bf=`echo $root | sed -r 's/.fastq.fp.fq//g'`
echo $bf

macs2 callpeak \
-n $bf \
-t $f \
-f BAM \
-g hs \
--outdir $out2 \
--SPMR \
-q 0.05 \
--nomodel --nolambda &
done


#ȥ��peak��blacklist region
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


#����logLR.bw
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


#MAnorm2

in=/mnt/6/WXY/GSE124718_human_early_embryo_H3K4me3_K27me3_K27ac/8CLC_histone_modification/merge_bam/
peak=/mnt/6/WXY/GSE124718_human_early_embryo_H3K4me3_K27me3_K27ac/8CLC_histone_modification/peak/
manorm \
--p1 $peak/8CLC-sorting_H3K27me3_merged.bam_peaks.broadPeak \
--p2 $peak/4CL_H3K27me3_merged.bam_peaks.broadPeak \
--pf broadpeak \
--r1 $in/8CLC-sorting_H3K27me3_merged.bam \
--r2 $in/4CL_H3K27me3_merged.bam \
--n1 8CLC-sorting_H3K27me3 \
--n2 4CL_H3K27me3 \
--rf bam \
--pe \
-w 1000 \
-p 0.05 \
-m 0.5849625 \
-o $peak/MAnorm &



#featureCounts


out1=/mnt/2/wuliang/05.Project/00.TGFB_project/00.Dataset/HistoneModification/bam/
out1=/mnt/6/WXY/GSE124718_human_early_embryo_H3K4me3_K27me3_K27ac/8CLC_histone_modification/merge_bam
out2=/mnt/6/WXY/GSE124718_human_early_embryo_H3K4me3_K27me3_K27ac/8CLC_histone_modification/peak/
counts=/mnt/6/WXY/GSE124718_human_early_embryo_H3K4me3_K27me3_K27ac/8CLC_histone_modification/counts/


for bed in $out2/*merged*.broadPeak
do
root=`basename $bed`
bf=`echo $root | sed -r 's/_peaks.broadPeak//g'`
annotatePeaks.pl $bed mm10 -annStats $counts/$bf".count.txt" > $counts/$bf".annotate.txt"
awk '{OFS="\t"; print $1,$2,$3,$4,"."}' $counts/$bf".annotate.txt" > $counts/$bf".annotate.saf"
sed -i '1d' $counts/$bf".annotate.saf"

featureCounts -O -F SAF -T 30 -p -B --minOverlap 26  \
-a $counts/$bf".annotate.saf" -o $counts/$bf".count" $out1/$bf

done


#PE
for bam in $out1
do
g=`basename $bam | sed 's/_VSI.bam/.count/g'`
featureCounts -O -F SAF -T 4 -p -B --minOverlap 26  \
-a $anno -o $counts/$g $bam 
done
#SE
for bam in $bamin1 $bamin2
do
g=`basename $bam | sed 's/_VSI.bam/.count/g'`
featureCounts -O -F SAF -T 4 --minOverlap 26  \
-a $anno -o $counts/$g $bam 
done

list=`ls $counts/*.count`
for f in $list
do
g=`basename $f | sed 's/.count/.c/g'`
cut -f 7 $f > $counts/$g
cut -f 1 $f > $counts/peak_id
done

paste $counts/peak_id $counts/*.c > $counts/cut-tag.countmatrix

