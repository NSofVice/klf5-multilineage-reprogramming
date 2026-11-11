
out=/mnt/6/WXY/20210308_compare_TSC_TF/bam

in=/mnt/6/WXY/TSC-GSE110950_mm9_TSC_TF_ChIPseq_ATACseq/bam



for f in $in/*_sorted.bam
do
root=`basename $f`
bf=`echo $root | sed -r 's/_sorted.bam//g'`
samtools view -u -bS -q 20 -F 0x04 $f | samtools sort > $in/$bf"_VSI.bam"
samtools index -b $in/$bf"_VSI.bam"
done