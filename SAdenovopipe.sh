#!/bin/bash

	### FOR SRST2, dependent on older versions of bowtie ###
# set PATH to include srst2 if it exists #
if [ -d "$HOME/gen-soft/srst2" ] ; then
	PATH="$HOME/gen-soft/srst2/scripts:$PATH"
fi
# set PATH to include bowtie2 for srst2 if it exists #
if [ -d "$HOME/gen-soft/bowtie2-2.2.9" ] ; then
	PATH="$HOME/gen-soft/bowtie2-2.2.9/bowtie2-2.2.9:$PATH"
fi
# set PATH to include samtools (old version) for srst2 if it exist #
if [ -d "$HOME/gen-soft/samtools" ] ; then
	PATH="$HOME/gen-soft/samtools:$PATH"
fi

### CREATE all necessary output folders ###
folder="SA2022"
bacteria="Staphylococcus_aureus"
mkdir -p ~/Documents/"$folder"/{corr,seeker,denovo,readqc,SPAdeslog,MLST/out,serotype/out,AMR/{amrtxt,amrjson},metaoutput} &&

### START OF ANALYSIS ###
for X in LIST_OF_SAMPLES ; do
		cd ~/Documents/"$folder"
		cp ./fastq/"$X"_R1.fastq.gz . &&
		cp ./fastq/"$X"_R2.fastq.gz . &&
		
		#trim R1 reads to 25k (100k lines) for strainseeker#
		gzip -d "$X"_R1.fastq.gz &&
		head -n 100000 "$X"_R1.fastq >"$X"_R1-25k.fastq &&
		gzip "$X"_R1.fastq &&
		
		#strainseeker#
		perl /home/user/gen-soft/strainseeker/seeker.pl -i "$X"_R1-25k.fastq -d /home/user/gen-soft/strainseeker/ss_db_w32_4324 -o "$X"_ss.txt &&
		echo "	*** done with Strainseeker for $X"
		### START SPADES de novo assembly; saves corrected reads for use in subsequent steps ###
		spades.py -t 28 -1 "$X"_R1.fastq.gz -2 "$X"_R2.fastq.gz -o "$X"_spades-de-novo  &&
			#cleanup
		cd "$X"_spades-de-novo/corrected/ &&
		mv "$X"_R1.fastq.00.0_0.cor.fastq.gz ../../corr/"$X"_R1-corr.fastq.gz &&
		mv "$X"_R2.fastq.00.0_0.cor.fastq.gz ../../corr/"$X"_R2-corr.fastq.gz &&
		cd .. &&
		mv contigs.fasta ../"$X"_contigs.fasta &&
		mv spades.log ../SPAdeslog/"$X"_spades.log &&
		cd .. &&
		echo "	*** done with SPAdes for $X"
		### RENAME contigs ###
		awk '/^>/{print ">""'"$X"'""_contigs-" ++i; next}{print}' < "$X"_contigs.fasta > "$X"_contigs2.fasta &&
		
		### Remove contigs <4K ###
		perl /home/user/gen-soft/removesmalls.pl 1000 "$X"_contigs2.fasta >"$X"_contigs-1K.fasta &&
		
		#Raspberry for qc on corrected reads#
		raspberry ./corr/"$X"* >"$X"_readqc.txt &&
		rm ./corr/*.rlen
		echo "	*** done with read_quality for $X"
			#cleanup
		mv "$X"_ss.txt ./seeker/ &&
		mv "$X"_contigs-1K.fasta ./denovo/ &&
		gzip ./denovo/"$X"_contigs-1K.fasta
		mv "$X"_readqc.txt ./readqc/ &&
		rm -r "$X"_spades-de-novo &&
		rm *.fastq &&
		rm *.fastq.gz &&
		rm "$X"_contigs2.fasta &&
		rm "$X"_contigs.fasta &&
		
		echo "	*** Starting analysis on gene seq"
		cd ~/Documents/GBS_reference_files/ &&
		
		### MLST ###
		# Requires correct MLST input files for bacterial species >Staphylococcus_aureus.fasta and >>>>staph galactiae.txt
		cp ../"$folder"/fastq/"$X"_R1.fastq.gz ./"$X"_1.fastq.gz &&
		cp ../"$folder"/fastq/"$X"_R2.fastq.gz ./"$X"_2.fastq.gz &&
			# srst2 command
		srst2.py --threads 20 --output "$X" --input_pe "$X"_*.fastq.gz --mlst_db Staphylococcus_aureus.fasta --mlst_definitions !!!!!Staphylococcus_aureus!!!!!!.txt --mlst_delimiter '_' &&
			#cleanup
		mv "$X"__mlst__Staphylococcus_aureus__results.txt ../"$folder"/MLST/out/"$X"_mlst.txt &&
		rm *.pileup &&
		rm *.bam &&
		echo "	*** Done with MLST for $X"
		
		### SEROTYPE ###
		srst2.py --input_pe "$X"_1.fastq.gz "$X"_2.fastq.gz --gene_db SA!!!!!!_serotype.fa --gene_max_mismatch 10 --min_coverage 98 --min_depth 10 --max_divergence 1 --output "$X" --threads 20 &&
		#mismatch changed from 1 to 10, depth from 15 to 10, divergence from 0 to 1#
			#cleanup
		mv "$X"__fullgenes__SA_serotype__results.txt ../"$folder"/serotype/out/"$X"_serotype.txt &&
		rm "$X"__*.bam &&
		rm "$X"__*.pileup &&
		rm "$X"__genes*.txt &&
		rm "$X"_1.fastq.gz &&
		rm "$X"_2.fastq.gz &&
		echo "	*** Done with serotyping for $X"
		
		### AMR calling ###
		cp ../"$folder"/denovo/"$X"_contigs-1K.fasta.gz . &&
		gunzip "$X"_contigs-1K.fasta.gz &&
		rgi main --input_sequence "$X"_contigs-1K.fasta --output_file "$X" --input_type contig --num_threads 20 --local --clean &&
			#cleanup
		rm "$X"_contigs-1K.fasta &&
		mv "$X".txt ../"$folder"/AMR/amrtxt/ &&
		mv "$X".json ../"$folder"/AMR/amrjson/ &&
		echo "	*** done with AMR calling for $X"
done

cd ~/Documents/"$folder"/MLST/ &&
head -1 ./out/header_mlst.txt > ../metaoutput/all_mlst.txt && tail -n +2 -q ./out/*.txt >> ../metaoutput/all_mlst.txt &&
cd ../serotype
cat ./out/*.txt > ../metaoutput/all_serotype.txt
cd ../AMR/ &&
rgi heatmap -i ./amrjson -clus both -f &&
mv ./*.csv ../metaoutput/ &&
cd ../seeker/ &&
cat *.txt > ../metaoutput/all_seeker.txt &&
echo "Done combining MLST, serotype and AMR data"
