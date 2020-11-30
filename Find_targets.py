#!/usr/bin/env python3

import os
import sys
import subprocess
from subprocess import check_output
import argparse
from gooey import Gooey


class MLSA_Blast():
	def __init__(self,ref_genomes,ref_targets,output):
		self.genomes = os.listdir(ref_genomes)
		self.targets = os.listdir(ref_targets)
		self.targets_folder = ref_targets
		self.genomes_folder = ref_genomes
		self.outdir = output

	def blast(self,file):
		(shortname,ext)=os.path.splitext(file)
		for genome in self.genomes:
			with open (os.path.join(self.outdir,shortname)+"_" + genome, "w") as output_handle:
				header=">" + file + "_" + genome + "\n"
				output_handle.write(header)
				seq=subprocess.Popen(["blastn","-query",os.path.join(self.targets_folder,file),"-subject",os.path.join(self.genomes_folder,genome),"-outfmt","6 sseq"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
				stdout,stderr=seq.communicate()
				output_handle.write(stdout.decode('utf-8'))

	def run(self):
		for file in self.targets:
			self.blast(file)

@Gooey
def main():
	parser = argparse.ArgumentParser(description="MAGE gene screen ", add_help=False)

	required = parser.add_argument_group('Required Arguments')
	required.add_argument('-i', '--input_refs', type=str, required=True, help="folder of reference genome files")
	required.add_argument('-t', '--targets', type=str, required=True, help="folder of targets in fasta format")
	required.add_argument('-o', '--output_dir', type=str, required=True, help="output folder")

	optional = parser.add_argument_group('Optional Arguments')
	optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

	args=parser.parse_args()
	input=os.path.normpath(args.input_refs)
	targets_folder=args.targets
	outdir=args.output_dir
	
	job = MLSA_Blast(input, targets_folder, outdir)
	job.run()


if __name__ == '__main__':
	main()
