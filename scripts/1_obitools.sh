#############################################
###Sequence data processing using obitools###
###   Yue Shi, University of Washington   ###  
#############################################

#home directory .
#pair-end raw sequencing for all samples are stored in the directory ./data/all
mkdir ./data/merged
mkdir ./data/aligned
mkdir ./data/noSingleton
mkdir ./data/addSampleTag

##############################
# Step 1: merge paired reads #
##############################

# merge with illuminapairedend command and filter out reads with alignment score less than 200; 
cd ./data/all
for R1 in *_R1_001.fastq
do
	R2="$(basename $R1 _R1_001.fastq)_R2_001.fastq"
	results_file="$(basename $R1 _R1_001.fastq)_merged.fastq"
	illuminapairedend --score-min=200 -r $R2 $R1 > ../merged/$results_file
done

# count raw sequence reads
for f in ./data/all/*.fastq
do
echo $(basename $f) >> file.names
obicount $f >> file.count
paste file.names file.count > sum.raw
done
rm file.count

cd ../..
# remove nonaligned reads
for f in ./data/merged/*.fastq
do 
	results_file=$(basename $f .fastq)_aligned.fastq
	obigrep -p 'mode!="joined"' $f > ./data/aligned/$results_file
done

################################
# Step 2: dereplicate globally #
################################

# add sample tag
cd ./data/aligned
for f in *.fastq
do
	tag=`echo "$f" | cut -d"_" -f1`
	results_file=$(basename $f .fastq)_sampletagadded.fastq
	obiannotate -S sample:$tag $f > ../addSampleTag/$results_file
done
cd ..
# merge all *.sampletagadded.fastq files and dereplicate
#compare all reads, group identical reads, output the sequence for each group and its count, remove duplicate reads. 
cat ./addSampleTag/*.fastq > ./dereplicated/canid.diet.fastq
cd ./dereplicated/
obiuniq -m sample canid.diet.fastq > canid.diet.uniq.fasta 
#-m sample: keep information of the samples of origin for each unique sequence.

# get total count per sequence
obiannotate -k count -k merged_sample canid.diet.uniq.fasta > $$ ; mv $$ canid.diet.uniq.shortheader.fasta 
#-k count -k merged-sample: to only keep per sample count and total count;
obistat -c count canid.diet.uniq.shortheader.fasta | sort -nk1 > canid_diet.count.txt 
#3 columns: occurence, number of unique sequences, total reads.

###################
# Step 3: denoise #
###################

#remove reads that are shorter than 80 bp and less than 400 copies with obigrep command;
obigrep -l 80 -p 'count>=400' canid.diet.uniq.shortheader.fasta > canid.diet.uniq.shortheader.c400.l80.fasta
#remove PCR and sequencing errors with obiclean command; 
obiclean -s merged_sample -r 0.05 -H canid.diet.uniq.shortheader.c400.l80.fasta > canid.diet.uniq.shortheader.c400.l80.clean.fasta
#-H: keep the head sequence; 
#-r 0.05: no variants with a count greater than 5% of their own count; 

##########################
# Step 4: sort and count #
##########################
#sequences can be sorted by decreasing order of count;
obisort -k count -r canid.diet.uniq.shortheader.c400.l80.clean.fasta > canid.diet.uniq.shortheader.c400.l80.clean.sort.fasta
#generate a count table and deduce sequence count of each sample
obitab -o canid.diet.uniq.shortheader.c400.l80.clean.sort.fasta > canid.diet.obitools.tab


