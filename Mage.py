#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Nexus import Nexus
import os,sys,csv
import subprocess
import argparse
import logging

"""
A program to take ab1 files from a range of target amplicons and trim + align forward and reverse reads. Then generate a consensus for each pair, align all sequences from each target, trim alignments, concatenate targets and generate
a phylogenetic tree using RAxML. This allows a quick and easy Multilocus sequence alignment and analysis for the user.
"""
# Formatting definitions for help menu ASCII art
END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
CYAN = '\033[36m'
YELLOW="\033[33m"

def green(text):
    return GREEN + text + END_FORMATTING

def bold_green(text):
    return GREEN + BOLD + text + END_FORMATTING

def magenta(text):
    return MAGENTA + text + END_FORMATTING

def bold(text):
    return BOLD + text + END_FORMATTING

def bold_yellow(text):
    return YELLOW + BOLD + text + END_FORMATTING

def bold_cyan(text):
    return CYAN + BOLD + text + END_FORMATTING

def remove_formatting(text):
    return re.sub('\033.*?m', '', text)

def get_ascii_art():
    ascii_art = (bold_green(" __  __    _    ____      \n") +
                 bold_green("|  \/  |  / \  / ___| ___ \n")+
                 bold_green("| |\/| | / _ \| |  _ / _ \      ")+
                 magenta(" /\ \n") +
                 bold_green("| |  | |/ ___ \ |_| |  __/     ")+
                 magenta("------\n") +
                 bold_green("|_|  |_/_/   \_\____|\___|    ")+
                 bold_cyan("(∩｀-｀)⊃━☆")+bold_yellow("ﾟ.*･｡ﾟ \n"))
    return(ascii_art)

MAGe=bold_yellow(magenta("M")+"ultilocus " + magenta("A")+"lignment " + magenta("Ge") +"nerator\n")
full_description=MAGe+"\n"+get_ascii_art()

# Calling argparser and defining the arguments
parser = argparse.ArgumentParser(description=full_description,formatter_class=argparse.RawTextHelpFormatter,add_help=False)

required = parser.add_argument_group('Required Arguments')
required.add_argument('-p', '--path', type=str, required=True, help="folder of abi files")
required.add_argument('-c', '--csv', type=str, required=True, help="spreadsheet of forward read, reverse read, sample ID and target region")

optional = parser.add_argument_group('Optional Arguments')
optional.add_argument("-f", "--force",default=False, action='store_true',required=False, help="Force overwrite of previous run")
optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

args=parser.parse_args()
path=os.path.normpath(args.path)
csv_input=args.csv
force=args.force

logging.basicConfig(level=logging.INFO, format='%(name)s - %(levelname)s - %(message)s')


if force:
	print("Forcing overwrite of previous analysis")

# STEP 1: format AB1 files as FASTQ
def reformat(dir):
	print("Starting Step 1: Reformat of AB1 to FASTQ")
	files=os.listdir(dir)
	fastq_output=os.path.join(path,"1_fastq_files")
	trimmed_output=os.path.join(path,"2_trimmed_reads")
	if not os.path.exists(fastq_output):
		os.mkdir(fastq_output)
	elif force:
		pass
	else:
		sys.exit("Output files already exist, please remove files to continue")
	for file in files:
		if os.path.isfile(os.path.join(path,file)): 
			(shortname, extension) = os.path.splitext(file)
			if extension == ".ab1":
				with open(os.path.join(path, file), "rb") as input_handle:
					with open(os.path.join(path,"1_fastq_files",shortname)+".fastq", "w") as output_handle:
						sequences=SeqIO.parse(input_handle, "abi")
						SeqIO.write(sequences, output_handle, "fastq")
	return(os.listdir(fastq_output))

# STEP 2: Trims low quality bases from the end of fastq files
def trim(fastq):
	print("Starting Step 2: Trimming")
	fastq_output=os.path.join(path,"1_fastq_files")
	trimmed_output=os.path.join(path,"2_trimmed_reads")
	if not os.path.exists(trimmed_output):
		os.mkdir(trimmed_output)
	elif force:
		pass
	else:
		sys.exit("Output files already exist, please remove files to continue")
	for file in fastq:
		(shortname, extension) = os.path.splitext(file)
		subprocess.run(["cutadapt","-q","40,40",os.path.join(fastq_output,file),"-o",os.path.join(trimmed_output,shortname+".trimmed"+extension)], stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

# STEP 3: Assemble forward and reverse reads
def assemble():
	print("Starting Step 3: Assembling forward and reverse sequences")
	fastq_input=os.path.join(path,"2_trimmed_reads")
	aln_output=os.path.join(path,"3_aligned_reads")
	missing_files=False
	if not os.path.exists(aln_output):
		os.mkdir(aln_output)
	elif force:
		for file in os.listdir(aln_output):
			subprocess.run(["rm",os.path.join(aln_output,file)])
	else:
		sys.exit("Output files already exist, please remove files to continue")
	with open(csv_input) as csv_file:
		csv_reader=csv.reader(csv_file,delimiter=',')
		for row in csv_reader:
			# If forward sample is missing from csv
			if not row[0]:
				logging.error("Forward sequence missing in csv for target '" + row[3] + "' of sample '"+row[2] + "'. Please ensure all samples have sequences for all targets. Exiting" )
				sys.exit()
			else:
				# If reverse sample is missing from csv
				if not row[1]:
					logging.error("Reverse sequence missing in csv for target '" + row[3] + "' of sample '"+row[2] + "'. Please ensure all samples have sequences for all targets. Exiting" )
					sys.exit()
				# Else if both forward or reverse are not missing
				else:
					if not row[0] == "Reference":
						pair=[]
						### BREAKS HERE IF FILE MISSING
						# Parses in the forward sample and appends it to 'pair'
						if os.path.exists(os.path.join(fastq_input,row[0])+".trimmed.fastq") and os.path.exists(os.path.join(fastq_input,row[1])+".trimmed.fastq"):
							with open(os.path.join(fastq_input,row[0])+".trimmed.fastq", "rU") as input_handle_F:
								for seq_record in SeqIO.parse(input_handle_F, "fastq"):
									pair.append(seq_record)
							# Parse in the reverse sample and appends it to 'pair'
							with open(os.path.join(fastq_input,row[1])+".trimmed.fastq", "rU") as input_handle_R:
								for seq_record in SeqIO.parse(input_handle_R, "fastq"):
									pair.append(seq_record)
							# Writes out the pair as a multifasta to use in mafft
							with open(os.path.join(fastq_input,row[2])+"mafft_merged"+row[3]+".fasta", "w") as output_handle:
								SeqIO.write(pair, output_handle, "fasta")
							# runs mafft on the newly created fasta file and redirects the stdout to a file
							with open(os.path.join(aln_output,str(row[2]).replace(" ","_"))+"mafft_merged"+row[3]+".aln","w") as aln:
								subprocess.run(['mafft','--adjustdirection',os.path.join(fastq_input,row[2])+"mafft_merged"+row[3]+".fasta"],stdout=aln,stderr=subprocess.DEVNULL)
						else:
							print("File/s for sample '"+row[2]+"' missing. Please check all files that are specified in the csv are present")
							missing_files=True
	if missing_files:
		print("Missing files. Exiting")
		sys.exit()
	else:
		pass

# STEP 4: Generate a consensus from forward and reverse reads
def consensus():
	print("Starting Step 4: Generating consensus of forward and reverse")
	alignment_input=os.path.join(path,"3_aligned_reads")
	consensus_output=os.path.join(path,"4_sanger_consensus")
	if not os.path.exists(consensus_output):
		os.mkdir(consensus_output)
	elif force:
		pass
	else:
		sys.exit("Output files already exist, please remove files to continue")
	for aln in os.listdir(alignment_input):
		(shortname,ext)=os.path.splitext(aln)
		sample_id=shortname.split("mafft_merged")
		subprocess.run(['cons', '-sequence',os.path.join(alignment_input,aln),"-name",sample_id[0],'-outseq',os.path.join(consensus_output,sample_id[0]+"_"+sample_id[1])+".consensus.fasta"],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

def reference_format():
	alignment_folder=os.path.join(path,"3_aligned_reads")
	consensus_output=os.path.join(path,"4_sanger_consensus")
	for file in os.listdir(path):
		(shortname, extension) = os.path.splitext(file)
		if extension == ".fasta" or extension == ".fna" or extension == ".fa":
			subprocess.run(["cp",os.path.join(path,file),os.path.join(alignment_folder,shortname)+".fasta"])
	with open(csv_input) as csv_file:
		csv_reader=csv.reader(csv_file,delimiter=',')
		for row in csv_reader:
			if str.lower(row[0]) == "reference":
				with open(os.path.join(alignment_folder, row[1])+".fasta", "rU") as input_handle:
					with open(os.path.join(consensus_output,row[1])+"_"+row[3]+".header_renamed.fasta", "w") as output_handle:
						for seq_record in SeqIO.parse(input_handle,"fasta"):
							seq_record.id=str(row[2]).replace(" ", "_")
							seq_record.description=""
							seq_record.name=""
							SeqIO.write(seq_record,output_handle,"fasta")

# STEP 5: Generate alignments for each target region
def alignment():
	print("Starting Step 5: Generating target alignments")
	consensus_input=os.path.join(path,"4_sanger_consensus")
	alignments_output=os.path.join(path,"5_alignments")
	if not os.path.exists(alignments_output):
		os.mkdir(alignments_output)
	elif force:
		pass
	else:
		sys.exit("Output files already exist, please remove files to continue")
	# Get list of MLSA targets from csv
	targets=[]
	with open(csv_input) as csv_file:
		csv_reader=csv.reader(csv_file,delimiter=',')
		for row in csv_reader:
			if row[3] not in targets:
				targets.append(row[3])

	# Take those targets and create a multifasta for mafft
	for target in targets:
		multifasta=[]
		for file in os.listdir(consensus_input):
			if target in file:
				with open(os.path.join(consensus_input,file),"rU") as input_handle:
					for seq_record in SeqIO.parse(input_handle,"fasta"):
						multifasta.append(seq_record)
		with open(os.path.join(alignments_output,target)+'.mafft_input.fasta',"w") as output_handle:
			SeqIO.write(multifasta,output_handle,"fasta")

		for file in os.listdir(alignments_output):
			(shortname, ext)= os.path.splitext(file)
			if ext == ".fasta":
				if target in file:
					with open(os.path.join(alignments_output,target)+".aln","w") as consensus:
						subprocess.run(['mafft','--adjustdirection',os.path.join(alignments_output,file)], stdout=consensus,stderr=subprocess.DEVNULL)

def trimal():
	print("                 Trimming alignments")
	alignment_input=os.path.join(path,"5_alignments")
	for file in os.listdir(alignment_input):
		(shortname,ext)=os.path.splitext(file)
		if ext == ".aln":
			with open(os.path.join(alignment_input,file), "r") as alignment_handle:
				with open(os.path.join(alignment_input,shortname)+"renamed"+ext,"w") as alignment_out_handle:
					data=alignment_handle.read()
					changed=data.replace("_R_","")
					alignment_out_handle.write(changed)
					subprocess.run(['trimal','-in',os.path.join(alignment_input,shortname)+"renamed"+ext,'-out',os.path.join(alignment_input,shortname+'.nex'),'-nexus','-nogaps'],stdout=subprocess.DEVNULL)


#STEP 6: Creating a supermatrix from the aligned regions
def supermatrix(): 
	print("Starting Step 6: Preparation of supermatrix")
	alignment_input=os.path.join(path,"5_alignments")
	supermatrix_output=os.path.join(path,"6_supermatrix")
	if not os.path.exists(supermatrix_output):
		os.mkdir(supermatrix_output)
	elif force:
		pass
	else:
		sys.exit("Output files already exist, please remove files to continue")
	matrix=[]
	for file in os.listdir(alignment_input):
		(shortname,ext)=os.path.splitext(file)
		if ext == ".nex":
			matrix.append(os.path.join(alignment_input,file))
	nexi = [(fname, Nexus.Nexus(fname)) for fname in matrix]
	combined = Nexus.combine(nexi)
	combined.write_nexus_data(filename=open(os.path.join(supermatrix_output,"supermatrix.nex"), "w"))

#STEP 7: Generating a bootstrap replicate tree from the supermatrix
def RAxML():
	print("Starting Step 7: Generating raxml tree")
	supermatrix_input=os.path.join(path,"6_supermatrix")
	tree_output=os.path.join(path,"7_FinalTree")
	if not os.path.exists(tree_output):
		os.mkdir(tree_output)
	elif force:
		pass
	else:
		sys.exit("Output files already exist, please remove files to continue")
	with open(os.path.join(supermatrix_input, "supermatrix.nex"), "rU") as input_handle:
		with open(os.path.join(supermatrix_input,"supermatrix.aln"), "w") as output_handle:
			sequences=SeqIO.parse(input_handle, "nexus")
			SeqIO.write(sequences, output_handle, "fasta")
	subprocess.run(['raxmlHPC-PTHREADS','-#','20','-T','8','-m','GTRGAMMA','-p','1337','-s',os.path.join(supermatrix_input,"supermatrix.aln"),'-w',os.path.abspath(tree_output),"-n","initial_tree"],stdout=subprocess.DEVNULL)
	subprocess.run(['raxmlHPC-PTHREADS','-T','8','-#','100','-p','1337','-b','1337',"-n","bootstrap","-s",os.path.join(supermatrix_input,"supermatrix.aln"),"-m","GTRGAMMA","-w",os.path.abspath(tree_output)],stdout=subprocess.DEVNULL)
	subprocess.run(['raxmlHPC-PTHREADS','-T','8',"-m","GTRGAMMA",'-p','1337',"-f","b","-t",os.path.join(tree_output,"RAxML_bestTree.initial_tree"),"-z",os.path.join(tree_output,"RAxML_bootstrap.bootstrap"),"-n","final","-w",os.path.abspath(tree_output)],stdout=subprocess.DEVNULL)
	out_file=open(os.path.join(tree_output,"MAGE_tree_output.final"),"w")
	with open(os.path.join(tree_output,"RAxML_bipartitionsBranchLabels.final"), "r") as tree_handle:
		with open(os.path.join(tree_output,"MAGE_tree_output.final"),"w") as mage_handle:
			data=tree_handle.read()
			mage_handle.write(data)

def main():
	trim(reformat(path))
	assemble()
	consensus()
	reference_format()
	alignment()
	trimal()
	supermatrix()
	RAxML()
	print(green('Done!'))
main()	
