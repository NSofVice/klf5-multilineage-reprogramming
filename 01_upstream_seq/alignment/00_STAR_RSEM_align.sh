## STAR + RSEM, parallel mode (miguel-cluster)
## Corresponds to manuscript step_01: Bulk RNA-seq alignment & quantification
##   STAR v2.7.3a + RSEM v1.3.1, reference Ensembl/UCSC mm10 (mm10_gencode_M25_rsem_star)
## Merge each sequencing run of a sample beforehand (cat R1 files, cat R2 files).
##
## Usage (inside the docker image):
##   ls ./*.fastq.gz | cat | parallel --jobs 4 align
##
## Notes on edits made to the archived original (all marked "# FIX"):
##   - removed an undefined `name=`basename $f`` line ($f was never set)
##   - `g=` now derives from the BAM name instead of $1 (originally basename of the
##     input FASTQ, which produced a wrong RSEM output prefix)
##   - RSEM output prefix given an explicit rsem_out/ directory
##   - quoted variables so paths containing spaces cannot break the call

docker run --tty --interactive --mount src=/mnt,target=/mnt,type=bind --workdir="${PWD}" ubuntu:18.04_fastp_rsem_star_deep bash

align() {

	r1=$1
	out=$(echo "$r1" | sed 's/.fastq.gz//g')
	bam="${out}Aligned.toTranscriptome.out.bam"
	echo "$r1" "$out" "$bam"

	genome=/mnt/1/genome/mm10/mm10_gencode_M25_rsem_star

	mkdir -p ./rsem_out

	STAR  \
	--genomeDir $genome  --outSAMunmapped Within  \
	--outFilterType BySJout  --outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20  \
	--outFilterMismatchNmax 999  --outFilterMismatchNoverLmax 0.04  --alignIntronMin 20  \
	--alignIntronMax 1000000  --alignMatesGapMax 1000000  --alignSJoverhangMin 8  --alignSJDBoverhangMin 1  \
	--sjdbScore 1  --runThreadN 5 --genomeLoad LoadAndKeep  --outSAMtype BAM Unsorted  \
	--quantMode TranscriptomeSAM  --outSAMheaderHD \@HD VN:1.4 SO:unsorted  \
	--outFileNamePrefix $out  \
	--readFilesCommand zcat  \
	--readFilesIn $r1

	# FIX: derive the RSEM sample name from the BAM, not from the FASTQ argument
	g=$(basename "$bam" | sed 's/Aligned.toTranscriptome.out.bam//g')

	rsem-calculate-expression \
	--bam "$bam" \
	-p 10 \
	/mnt/1/genome/mm10/mm10_gencode_M25_rsem_star/index ./rsem_out/"$g"

	wait
	rm "$bam" "$out"Aligned.out.bam
}
export -f align

ls ./*.fastq.gz | cat | parallel --jobs 4 align