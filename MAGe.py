from Bio import SeqIO
from Bio.Nexus import Nexus
import os
import sys
import csv
import subprocess
import argparse
from gooey import Gooey

class Mage():
	def __init__(self,path,csv,force):
		self.path = path
		self.csv = csv
		self.force = force
		self.abifiles = os.listdir(path)
		self.fastq_dir1=os.path.join(path,"1_fastq_files")
		self.trimmed_dir2=os.path.join(path,"2_trimmed_reads")
		self.aligned_dir3=os.path.join(path,"3_aligned_reads")
		self.consensus_dir4=os.path.join(path,"4_sanger_consensus")
		self.alignments_dir5=os.path.join(path,"5_alignments")
		self.supermatrix_dir6=os.path.join(path,"6_supermatrix")
		self.tree_dir7=os.path.join(path,"7_FinalTree")

	def reformat(self, path):
		print("Starting Step 1: Reformat of AB1 to FASTQ")
		if not os.path.exists(self.fastq_dir1):
			os.mkdir(self.fastq_dir1)
		elif self.force:
			pass
		else:
			sys.exit("Output files already exist, please remove files to continue")
		for file in self.abifiles:
			if os.path.isfile(os.path.join(self.path,file)): 
				(shortname, extension) = os.path.splitext(file)
				if extension == ".ab1":
					with open(os.path.join(self.path, file), "rb") as input_handle:
						with open(os.path.join(self.fastq_dir1,shortname)+".fastq", "w") as output_handle:
							sequences=SeqIO.parse(input_handle, "abi")
							SeqIO.write(sequences, output_handle, "fastq")
		return(os.listdir(self.fastq_dir1))
	
	def trim(self, fastq):
		print("Starting Step 2: Trimming")
		if not os.path.exists(self.trimmed_dir2):
			os.mkdir(self.trimmed_dir2)
		elif self.force:
			pass
		else:
			sys.exit("Output files already exist, please remove files to continue")
		for file in fastq:
			(shortname, extension) = os.path.splitext(file)
			subprocess.run(["cutadapt","-q","40,40",os.path.join(self.fastq_dir1,file),"-o",os.path.join(self.trimmed_dir2,shortname+".trimmed"+extension)], stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

	def assemble(self):
		print("Starting Step 3: Assembling forward and reverse sequences")
		missing_files=False
		if not os.path.exists(self.aligned_dir3):
			os.mkdir(self.aligned_dir3)
		elif force:
			for file in os.listdir(self.aligned_dir3):
				subprocess.run(["rm",os.path.join(self.aligned_dir3,file)])
		else:
			sys.exit("Output files already exist, please remove files to continue")
		with open(self.csv) as csv_file:
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
							# Parses in the forward sample and appends it to 'pair'
							if os.path.exists(os.path.join(self.trimmed_dir2,row[0])+".trimmed.fastq") and os.path.exists(os.path.join(self.trimmed_dir2,row[1])+".trimmed.fastq"):
								with open(os.path.join(self.trimmed_dir2,row[0])+".trimmed.fastq", "rU") as input_handle_F:
									for seq_record in SeqIO.parse(input_handle_F, "fastq"):
										pair.append(seq_record)
								# Parse in the reverse sample and appends it to 'pair'
								with open(os.path.join(self.trimmed_dir2,row[1])+".trimmed.fastq", "rU") as input_handle_R:
									for seq_record in SeqIO.parse(input_handle_R, "fastq"):
										pair.append(seq_record)
								# Writes out the pair as a multifasta to use in mafft
								with open(os.path.join(self.trimmed_dir2,row[2])+"mafft_merged"+row[3]+".fasta", "w") as output_handle:
									SeqIO.write(pair, output_handle, "fasta")
								# runs mafft on the newly created fasta file and redirects the stdout to a file
								with open(os.path.join(self.aligned_dir3,str(row[2]).replace(" ","_"))+"mafft_merged"+row[3]+".aln","w") as aln:
									subprocess.run(['mafft','--adjustdirection',os.path.join(self.trimmed_dir2,row[2])+"mafft_merged"+row[3]+".fasta"],stdout=aln,stderr=subprocess.DEVNULL)
							else:
								print("File/s for sample '"+row[2]+"' missing. Please check all files that are specified in the csv are present")
								missing_files=True
		if missing_files:
			print("Missing files. Exiting")
			sys.exit()
		else:
			pass

	def consensus(self):
		print("Starting Step 4: Generating consensus of forward and reverse")
		if not os.path.exists(self.consensus_dir4):
			os.mkdir(self.consensus_dir4)
		elif force:
			pass
		else:
			sys.exit("Output files already exist, please remove files to continue")
		for aln in os.listdir(self.aligned_dir3):
			(shortname,ext)=os.path.splitext(aln)
			sample_id=shortname.split("mafft_merged")
			subprocess.run(['cons', '-sequence',os.path.join(self.aligned_dir3,aln),"-name",sample_id[0],'-outseq',os.path.join(self.consensus_dir4,sample_id[0]+"_"+sample_id[1])+".consensus.fasta"],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

	def reference_format(self):
		for file in os.listdir(self.path):
			(shortname, extension) = os.path.splitext(file)
			if extension == ".fasta" or extension == ".fna" or extension == ".fa":
				subprocess.run(["cp",os.path.join(self.path,file),os.path.join(self.aligned_dir3,shortname)+".fasta"])
		with open(self.csv) as csv_file:
			csv_reader=csv.reader(csv_file,delimiter=',')
			for row in csv_reader:
				if str.lower(row[0]) == "reference":
					try:
						with open(os.path.join(self.aligned_dir3, row[1])+".fasta", "rU") as input_handle:
							with open(os.path.join(self.consensus_dir4,row[1])+"_"+row[3]+".header_renamed.fasta", "w") as output_handle:
								for seq_record in SeqIO.parse(input_handle,"fasta"):
									seq_record.id=str(row[2]).replace(" ", "_")
									seq_record.description=""
									seq_record.name=""
									SeqIO.write(seq_record,output_handle,"fasta")
					except Exception as error:
						print("File '" + row[2] + "' in csv, does not exist in directory" + self.path)

	def alignment(self):
		print("Starting Step 5: Generating target alignments")
		if not os.path.exists(self.alignments_dir5):
			os.mkdir(self.alignments_dir5)
		elif force:
			pass
		else:
			sys.exit("Output files already exist, please remove files to continue")
		# Get list of MLSA targets from csv
		targets=[]
		with open(self.csv) as csv_file:
			csv_reader=csv.reader(csv_file,delimiter=',')
			for row in csv_reader:
				if row[3] not in targets:
					targets.append(row[3])

		# Take those targets and create a multifasta for mafft
		for target in targets:
			multifasta=[]
			for file in os.listdir(self.consensus_dir4):
				if target in file:
					with open(os.path.join(self.consensus_dir4,file),"rU") as input_handle:
						for seq_record in SeqIO.parse(input_handle,"fasta"):
							multifasta.append(seq_record)
			with open(os.path.join(self.alignments_dir5,target)+'.mafft_input.fasta',"w") as output_handle:
				SeqIO.write(multifasta,output_handle,"fasta")

			for file in os.listdir(self.alignments_dir5):
				(shortname, ext)= os.path.splitext(file)
				if ext == ".fasta":
					if target in file:
						with open(os.path.join(self.alignments_dir5,target)+".aln","w") as consensus:
							subprocess.run(['mafft','--adjustdirection',os.path.join(self.alignments_dir5,file)], stdout=consensus,stderr=subprocess.DEVNULL)

	def trimal(self):
		print("                 Trimming alignments")
		for file in os.listdir(self.alignments_dir5):
			(shortname,ext)=os.path.splitext(file)
			if ext == ".aln":
				with open(os.path.join(self.alignments_dir5,file), "r") as alignment_handle:
					with open(os.path.join(self.alignments_dir5,shortname)+"renamed"+ext,"w") as alignment_out_handle:
						data=alignment_handle.read()
						changed=data.replace("_R_","")
						alignment_out_handle.write(changed)
						subprocess.run(['trimal','-in',os.path.join(self.alignments_dir5,shortname)+"renamed"+ext,'-out',os.path.join(self.alignments_dir5,shortname+'.nex'),'-nexus','-nogaps'],stdout=subprocess.DEVNULL)

	def supermatrix(self): 
		print("Starting Step 6: Preparation of supermatrix")
		if not os.path.exists(self.supermatrix_dir6):
			os.mkdir(self.supermatrix_dir6)
		elif force:
			pass
		else:
			sys.exit("Output files already exist, please remove files to continue")
		matrix=[]
		for file in os.listdir(self.alignments_dir5):
			(shortname,ext)=os.path.splitext(file)
			if ext == ".nex":
				matrix.append(os.path.join(self.alignments_dir5,file))
		nexi = [(fname, Nexus.Nexus(fname)) for fname in matrix]
		combined = Nexus.combine(nexi)
		combined.write_nexus_data(filename=open(os.path.join(self.supermatrix_dir6,"supermatrix.nex"), "w"))

	def RAxML(self):
		print("Starting Step 7: Generating raxml tree")
		if not os.path.exists(self.tree_dir7):
			os.mkdir(self.tree_dir7)
		elif force:
			pass
		else:
			sys.exit("Output files already exist, please remove files to continue")
		with open(os.path.join(self.supermatrix_dir6, "supermatrix.nex"), "rU") as input_handle:
			with open(os.path.join(self.supermatrix_dir6,"supermatrix.aln"), "w") as output_handle:
				sequences=SeqIO.parse(input_handle, "nexus")
				SeqIO.write(sequences, output_handle, "fasta")
		subprocess.run(['raxmlHPC-PTHREADS','-#','20','-T','8','-m','GTRGAMMA','-p','1337','-s',os.path.join(self.supermatrix_dir6,"supermatrix.aln"),'-w',os.path.abspath(self.tree_dir7),"-n","initial_tree"],stdout=subprocess.DEVNULL)
		subprocess.run(['raxmlHPC-PTHREADS','-T','8','-#','100','-p','1337','-b','1337',"-n","bootstrap","-s",os.path.join(self.supermatrix_dir6,"supermatrix.aln"),"-m","GTRGAMMA","-w",os.path.abspath(self.tree_dir7)],stdout=subprocess.DEVNULL)
		subprocess.run(['raxmlHPC-PTHREADS','-T','8',"-m","GTRGAMMA",'-p','1337',"-f","b","-t",os.path.join(self.tree_dir7,"RAxML_bestTree.initial_tree"),"-z",os.path.join(self.tree_dir7,"RAxML_bootstrap.bootstrap"),"-n","final","-w",os.path.abspath(self.tree_dir7)],stdout=subprocess.DEVNULL)
		out_file=open(os.path.join(self.tree_dir7,"MAGE_tree_output.final"),"w")
		with open(os.path.join(self.tree_dir7,"RAxML_bipartitionsBranchLabels.final"), "r") as tree_handle:
			with open(os.path.join(self.tree_dir7,"MAGE_tree_output.final"),"w") as mage_handle:
				data=tree_handle.read()
				mage_handle.write(data)

	def run(self):
		self.trim(self.reformat(self.path))
		self.assemble()
		self.consensus()
		self.reference_format()
		self.alignment()
		self.trimal()
		self.supermatrix()
		self.RAxML()

@Gooey
def main():
	parser = argparse.ArgumentParser(description="MLSA Generator", formatter_class=argparse.RawTextHelpFormatter, add_help=False)
	optional = parser.add_argument_group('Optional Arguments')
	required = parser.add_argument_group('Required Arguments')

	required.add_argument('-p', '--path', type=str, required=True, help="folder of abi files")
	required.add_argument('-c', '--csv', type=str, required=True, help="spreadsheet of forward read, reverse read, sample ID and target region")
	optional.add_argument("-f", "--force",default=False, action='store_true',required=False, help="Force overwrite of previous run")
	optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

	args=parser.parse_args()
	path=os.path.normpath(args.path)
	csv_input=args.csv
	force=args.force

	job=Mage(path, csv_input, force)
	job.run()

if __name__ == '__main__':
	main()
