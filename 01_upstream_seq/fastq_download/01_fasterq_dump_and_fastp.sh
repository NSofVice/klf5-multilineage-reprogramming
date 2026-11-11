#PBS -l walltime=999:00:00
#PBS -l nodes=1:ppn=15,mem=10gb
#PBS -q cv2
#PBS -j oe






cat download.txt | while read id
do
prefetch $id -X 500G -O `pwd` & 
done





in=$1
out=`dirname $in`/fastp_output
mkdir $out


Fastq_dump(){
#fasterq-dump --split-3 $1 -e 2 -O `pwd`
/mnt/11/WXY_classic_data/download/sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump --split-3 $1 -O `pwd`  --gzip
}  #qstat -Qf可以查看有哪些能提交的队列

export -f Fastq_dump

ls *.sra | cat | parallel --tmpdir `pwd` --jobs 5 Fastq_dump


#fasterq-dump --split-3 $f -e 10 -O $in 

echo "Fasterq-dump is done"
echo "Doning fastp"
##fastp is a new data quality control software, also see https://www.jianshu.com/p/6f492058da5b and https://github.com/OpenGene/fastp

Fastp(){
in=`pwd`
out=`pwd`
name=`basename $1`
#fastq-dump --split-3 --gzip $f -O $out
#fasterq-dump --split-3 $f -e 10 -O $in 
root=`echo $name`
if [[ $root == *_1.fq.gz ]]
then
echo PE...fastp
echo $root
bf=`echo $root | sed -r 's/\_1.fq//g' | sed -r 's/.gz//g'`
echo $bf
p1=`echo $root`
echo $p1
p2=`echo $root | sed 's/\_1/\_2/g'`
echo $p2
fastp -i $in/$p1 -o $out/$bf"_1.fp.fq.gz" \
-I $in/$p2 -O $out/$bf"_2.fp.fq.gz" \
--thread=5 --compression=10 --correction
# not a PE read
else 
if [[ $root != *_2.fq.gz ]]
then
echo SE...fastp
bf=`echo $root | sed 's/.fq.gz//g'`
echo $bf
fastp -i $in/$root -o $out/$bf".fp.fq.gz" \
--thread=5 --compression=10 
fi
fi
}
export -f Fastp

ls  SRR*794*.fastq.gz SRR*793*.fastq.gz SRR*792*.fastq.gz SRR*788*.fastq.gz SRR*791*.fastq.gz | cat | parallel --tmpdir `pwd` --jobs 5 Fastp &
ls  *fq.gz | cat | parallel --tmpdir `pwd` --jobs 4 Fastp 