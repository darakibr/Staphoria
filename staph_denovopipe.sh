#!/bin/bash

	### FOR SRST2 ###
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

#### CHANGE TO APPROPRIATE NAME ####
folder="Staphoria"
bacteria="Staphylococcus_aureus"

### CREATE all necessary output folders ###
mkdir corr &&
mkdir seeker &&
mkdir readqc &&
mkdir SPAdeslog &&
mkdir -p MLST/out &&
mkdir -p serotype/out &&
mkdir -p AMR/amrtxt &&
mkdir -p AMR/amrjson &&
mkdir metaoutput &&

### START OF ANALYSIS ###
for X in 180N12 247N6F 235N3F 94A12 196N12 167N12 45N3 243A3 258NI 84N12 288N3 205N6F 180N9F 49NI 114N12 62N12 243A6 63N12F 248N3 140NI 158N12 70NI 146N12 47N12 179N12 120N9 250A9 78NI 78N12F 208N12 66NI 73N12 248A6F 21N12 18N9 97N12 176N6 235A3F 144N9 68NI 194N12F 148N9 207N9F 14N9 60NI 236N3 194N9 73N6 148N12 73N12F 114N9 176N12 233N6F 84N12F 52N12F 99NI 261N3 191N6F 68A12F 57N12 143N12 132N12 102N9 99N9F 17N9 78N12 68A12 58N9 129N6F 94A12F 208N6F 98N12 246N3 71N12F 136N9 129N9 57N12F 153N12 210N12 178N12 242N6F 210N9F 124N6 218N6F 208N9F 207A12 235N6F 132N9 73A9 98N9F 250A3 97A12 65NIF 194A12 150A12 62N12F 77A9 52N12F 125N6 24AIF 140N9 160N12 139NF 64AIF 71N12 197N12 83N9 246N6F 163N12 4A9 62NI 219N6F 150N6F 47N9F 107N6F 247N3 176N9 152N12 233N3 83NI 180N6 269AI 256N3 20N9F 100N6 80N9F 96N12F 73A6 43N9F 81NI 181N12 84N9 78N9F 60N9 244N3 233N9 22N12 58NI 248N6F 270N3 160A12 100N12 15N9 148N6 ; do
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
