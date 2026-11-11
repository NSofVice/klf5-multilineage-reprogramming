
out=/mnt/6/WXY/20210308_compare_TSC_TF/bw_center_reads

in=/mnt/6/WXY/20210308_compare_TSC_TF/bam



for f in $in/*_sorted.bam
do
root=`basename $f`
bamCoverage -b $f \
--normalizeUsing RPKM \
--extendReads 200 \
--binSize 10 \
--centerReads \
--effectiveGenomeSize 2652783500 \
-o $out/$root".bw" &
done 