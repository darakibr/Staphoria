### This program will prepare demultiplexed fastq files for analysis
### Program used to demultiplex is bcl-convert from Illumina data (MiSeq or NextSeq)
### This will rename the fastq files and prepare a sample text file

### load libraries
import os
import re

### Paths saved to variables
cwd = os.getcwd()
fastqfolder = cwd+"/fastq/"

def prepilluminafiles(fastqfolder):
	### Prepare lists to hold files
	samples = []
	other_files = []
	counter = 0
	for file in os.listdir(fastqfolder):
		try:
			if file.endswith("R1_001.fastq.gz"):
				file_rn = str(file)
				file_rn = re.sub('_S.+?R','_R', file_rn)
				file_rn = re.sub('_001','',file_rn)
				os.rename(fastqfolder+file, fastqfolder+file_rn)
				samples.append(re.sub('_R1.fastq.gz','', file_rn))
				counter += 1
			elif file.endswith("R2_001.fastq.gz"):
				file_rn=str(file)
				file_rn = re.sub('_S.+?R','_R', file_rn)
				file_rn = re.sub('_001','',file_rn)
				os.rename(fastqfolder+file,fastqfolder+file_rn)
				counter += 1
			elif file.endswith("_R1.fastq.gz"):
				file_rn=str(file)
				samples.append(re.sub('_R1.fastq.gz','',file_rn))
				counter += 1
			elif file.endswith("_R2.fastq.gz"):
				counter += 1
			else:
				other_files.append(str(file))
				counter += 1
		except Exception as e:
			raise e
			print("No files found...")
	print("Number of files renamed = ",counter)
	print("Number of unique samples found = ",len(samples))
	samplesinput = ' '.join(map(str, samples))
	with open("samples.txt", "w") as file:
		file.write(samplesinput)
	print("================ FINISHED RENAMING FASTQ TO PROPER INPUT NAMES ================")

### Use function on the fastqfolder ###
prepilluminafiles(fastqfolder)


