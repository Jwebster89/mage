from Bio import SeqIO
import os
import sys
import csv
import subprocess
import argparse
import shutil
import logging


class Mage():
    def __init__(self, path, csv, force):
        self.path = path
        self.csv = csv
        self.force = force
        self.abifiles = os.listdir(path)

        self.fastq_dir1 = os.path.join(path, "1_fastq_files")
        self.trimmed_dir2 = os.path.join(path, "2_trimmed_reads")
        self.aligned_dir3 = os.path.join(path, "3_aligned_reads")
        self.consensus_dir4 = os.path.join(path, "4_sanger_consensus")
        self.alignments_dir5 = os.path.join(path, "5_alignments")
        self.supermatrix_dir6 = os.path.join(path, "6_supermatrix")
        self.tree_dir7 = os.path.join(path, "7_FinalTree")

    def force_removal(self, force):
        if self.force is True:
            for d in [
                self.fastq_dir1,
                self.trimmed_dir2,
                self.aligned_dir3,
                self.consensus_dir4,
                self.alignments_dir5,
                self.supermatrix_dir6,
                self.tree_dir7,
            ]:
                if os.path.exists(d):
                    shutil.rmtree(d)

    def reformat(self, path):
        print("Starting Step 1: Reformat of AB1 to FASTQ")
        if not os.path.exists(self.fastq_dir1):
            os.mkdir(self.fastq_dir1)
        elif self.force:
            pass
        else:
            sys.exit("Output files already exist, please remove files to continue")

        for file in self.abifiles:
            if os.path.isfile(os.path.join(self.path, file)):
                (shortname, extension) = os.path.splitext(file)
                if extension == ".ab1":
                    with open(os.path.join(self.path, file), "rb") as input_handle:
                        with open(os.path.join(self.fastq_dir1, shortname) + ".fastq", "w") as output_handle:
                            sequences = SeqIO.parse(input_handle, "abi")
                            SeqIO.write(sequences, output_handle, "fastq")
        return os.listdir(self.fastq_dir1)

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
            subprocess.run(
                ["cutadapt", "-q", "40,40", os.path.join(self.fastq_dir1, file),
                 "-o", os.path.join(self.trimmed_dir2, shortname + ".trimmed" + extension)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )

    def assemble(self):
        print("Starting Step 3: Assembling forward and reverse sequences")
        missing_files = False

        if not os.path.exists(self.aligned_dir3):
            os.mkdir(self.aligned_dir3)
        elif self.force:
            pass
        else:
            sys.exit("Output files already exist, please remove files to continue")

        with open(self.csv) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                if not row[0]:
                    logging.error("Forward sequence missing in csv for target '" + row[3] + "' of sample '" + row[2] + "'. Exiting")
                    sys.exit()
                if not row[1]:
                    logging.error("Reverse sequence missing in csv for target '" + row[3] + "' of sample '" + row[2] + "'. Exiting")
                    sys.exit()

                if row[0] != "Reference":
                    pair = []
                    f_trim = os.path.join(self.trimmed_dir2, row[0]) + ".trimmed.fastq"
                    r_trim = os.path.join(self.trimmed_dir2, row[1]) + ".trimmed.fastq"

                    if os.path.exists(f_trim) and os.path.exists(r_trim):
                        with open(f_trim, "r") as input_handle_F:
                            for seq_record in SeqIO.parse(input_handle_F, "fastq"):
                                pair.append(seq_record)

                        with open(r_trim, "r") as input_handle_R:
                            for seq_record in SeqIO.parse(input_handle_R, "fastq"):
                                pair.append(seq_record)

                        merged_fa = os.path.join(self.trimmed_dir2, row[2]) + "mafft_merged" + row[3] + ".fasta"
                        with open(merged_fa, "w") as output_handle:
                            SeqIO.write(pair, output_handle, "fasta")

                        aln_out = os.path.join(self.aligned_dir3, str(row[2]).replace(" ", "_")) + "mafft_merged" + row[3] + ".aln"
                        with open(aln_out, "w") as aln:
                            subprocess.run(
                                ['mafft', '--adjustdirection', merged_fa],
                                stdout=aln,
                                stderr=subprocess.DEVNULL
                            )
                    else:
                        print(f"File/s for sample '{row[2]}' missing. Please check all files specified in the csv are present")
                        missing_files = True

        if missing_files:
            print("Missing files. Exiting")
            sys.exit()

    def consensus(self):
        print("Starting Step 4: Generating consensus of forward and reverse")
        if not os.path.exists(self.consensus_dir4):
            os.mkdir(self.consensus_dir4)
        elif self.force:
            pass
        else:
            sys.exit("Output files already exist, please remove files to continue")

        for aln in os.listdir(self.aligned_dir3):
            (shortname, ext) = os.path.splitext(aln)
            sample_id = shortname.split("mafft_merged")
            subprocess.run(
                ['cons', '-sequence', os.path.join(self.aligned_dir3, aln),
                 "-name", sample_id[0],
                 '-outseq', os.path.join(self.consensus_dir4, sample_id[0] + "_" + sample_id[1]) + ".consensus.fasta"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )

    def reference_format(self):
        # Copy reference files exactly as-is into 3_aligned_reads
        for file in os.listdir(self.path):
            (shortname, extension) = os.path.splitext(file)
            if extension.lower() in [".fasta", ".fna", ".fa"]:
                src = os.path.join(self.path, file)
                dst = os.path.join(self.aligned_dir3, file)
                shutil.copy2(src, dst)

        # Rename headers for reference sequences into 4_sanger_consensus
        with open(self.csv) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                if str(row[0]).lower() == "reference":
                    try:
                        ref_path = os.path.join(self.aligned_dir3, row[1])
                        with open(ref_path, "r") as input_handle:
                            out_path = os.path.join(self.consensus_dir4, row[1] + "_" + row[3] + ".header_renamed.fasta")
                            with open(out_path, "w") as output_handle:
                                for seq_record in SeqIO.parse(input_handle, "fasta"):
                                    seq_record.id = str(row[2]).replace(" ", "_")
                                    seq_record.description = ""
                                    seq_record.name = ""
                                    SeqIO.write(seq_record, output_handle, "fasta")
                    except Exception:
                        print(f"Reference file '{row[1]}' in csv does not exist in directory {self.path}")

    def alignment(self):
        print("Starting Step 5: Generating target alignments")
        if not os.path.exists(self.alignments_dir5):
            os.mkdir(self.alignments_dir5)
        elif self.force:
            pass
        else:
            sys.exit("Output files already exist, please remove files to continue")

        # Get unique targets
        targets = []
        with open(self.csv) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                if row[3] not in targets:
                    targets.append(row[3])

        # For each target, build multifasta and run mafft
        for target in targets:
            multifasta = []
            for file in os.listdir(self.consensus_dir4):
                if target in file:
                    with open(os.path.join(self.consensus_dir4, file), "r") as input_handle:
                        for seq_record in SeqIO.parse(input_handle, "fasta"):
                            multifasta.append(seq_record)

            in_fa = os.path.join(self.alignments_dir5, target) + '.mafft_input.fasta'
            with open(in_fa, "w") as output_handle:
                SeqIO.write(multifasta, output_handle, "fasta")

            aln_out = os.path.join(self.alignments_dir5, target) + ".aln"
            with open(aln_out, "w") as consensus:
                subprocess.run(
                    ['mafft', '--adjustdirection', in_fa],
                    stdout=consensus,
                    stderr=subprocess.DEVNULL
                )

    def trimal(self):
        """
        Output trimmed alignments as FASTA (not NEXUS),
        avoiding Biopython Nexus parsing issues later.
        """
        print("                 Trimming alignments")
        for file in os.listdir(self.alignments_dir5):
            (shortname, ext) = os.path.splitext(file)
            if ext == ".aln":
                in_path = os.path.join(self.alignments_dir5, file)

                renamed_path = os.path.join(self.alignments_dir5, shortname + "_renamed.aln")
                trimmed_fa = os.path.join(self.alignments_dir5, shortname + ".trimmed.fasta")

                # remove "_R_" then close file before trimal
                with open(in_path, "r") as alignment_handle:
                    data = alignment_handle.read()

                changed = data.replace("_R_", "")

                with open(renamed_path, "w") as out_handle:
                    out_handle.write(changed)

                subprocess.run(
                    ['trimal', '-in', renamed_path, '-out', trimmed_fa, '-fasta', '-nogaps'],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )

    def supermatrix(self):
        """
        Concatenate per-target trimmed FASTA alignments into a single supermatrix FASTA.
        Missing taxa in a target are padded with gaps of that target's alignment length.
        """
        print("Starting Step 6: Preparation of supermatrix")
        if not os.path.exists(self.supermatrix_dir6):
            os.mkdir(self.supermatrix_dir6)
        elif self.force:
            pass
        else:
            sys.exit("Output files already exist, please remove files to continue")

        aln_files = []
        for f in os.listdir(self.alignments_dir5):
            if f.endswith(".trimmed.fasta"):
                aln_files.append(os.path.join(self.alignments_dir5, f))

        if not aln_files:
            sys.exit("No trimmed fasta alignments were produced (Step 5/trimAl failed). Cannot build supermatrix.")

        # Read all alignments, track taxa and alignment lengths
        all_taxa = set()
        per_alignment = []  # list of (aln_len, {taxon: seq})
        locus_stats = []    # list of dicts for summary output

        for fa in sorted(aln_files):
            seqs = {}
            aln_len = None
            taxa_present = 0

            for rec in SeqIO.parse(fa, "fasta"):
                taxa_present += 1
                tax = rec.id
                s = str(rec.seq).upper()
                if aln_len is None:
                    aln_len = len(s)
                elif len(s) != aln_len:
                    sys.exit(f"Alignment length mismatch inside file {fa}: {tax} has {len(s)} but expected {aln_len}")
                seqs[tax] = s
                all_taxa.add(tax)

            if aln_len is None:
                sys.exit(f"Trimmed alignment is empty: {fa}")

            per_alignment.append((aln_len, seqs))

            locus_name = os.path.basename(fa).replace(".trimmed.fasta", "")
            locus_stats.append({
                "locus": locus_name,
                "aln_len": aln_len,
                "taxa_present": taxa_present
            })

        # Concatenate
        superseq = {tax: "" for tax in sorted(all_taxa)}
        for aln_len, seqs in per_alignment:
            gap = "-" * aln_len
            for tax in superseq:
                superseq[tax] += seqs.get(tax, gap)

        # Write supermatrix FASTA
        out_fa = os.path.join(self.supermatrix_dir6, "supermatrix.fasta")
        with open(out_fa, "w") as handle:
            for tax, s in superseq.items():
                handle.write(f">{tax}\n{s}\n")

        # Write summary TSV
        summary_path = os.path.join(self.supermatrix_dir6, "supermatrix_summary.tsv")
        total_taxa = len(superseq)
        total_len = len(next(iter(superseq.values()))) if total_taxa > 0 else 0

        with open(summary_path, "w") as out:
            out.write("metric\tvalue\n")
            out.write(f"targets\t{len(locus_stats)}\n")
            out.write(f"taxa_total\t{total_taxa}\n")
            out.write(f"supermatrix_length_bp\t{total_len}\n")
            out.write("\n")
            out.write("locus\tlength_bp\ttaxa_present\ttaxa_missing\n")
            for st in sorted(locus_stats, key=lambda x: x["locus"]):
                missing = total_taxa - st["taxa_present"]
                out.write(f"{st['locus']}\t{st['aln_len']}\t{st['taxa_present']}\t{missing}\n")

        print(f"Step 6 summary written to: {summary_path}")

    def RAxML(self):
        """
        Run RAxML directly on supermatrix FASTA.
        """
        print("Starting Step 7: Generating raxml tree")
        if not os.path.exists(self.tree_dir7):
            os.mkdir(self.tree_dir7)
        elif self.force:
            pass
        else:
            sys.exit("Output files already exist, please remove files to continue")

        aln = os.path.join(self.supermatrix_dir6, "supermatrix.fasta")
        if not os.path.exists(aln):
            sys.exit("supermatrix.fasta not found. Step 6 must succeed before RAxML can run.")

        subprocess.run(
            ['raxmlHPC-PTHREADS', '-#', '20', '-T', '8', '-m', 'GTRGAMMA', '-p', '1337',
             '-s', aln,
             '-w', os.path.abspath(self.tree_dir7), "-n", "initial_tree"],
            stdout=subprocess.DEVNULL
        )

        subprocess.run(
            ['raxmlHPC-PTHREADS', '-T', '8', '-#', '100', '-p', '1337', '-b', '1337',
             "-n", "bootstrap",
             "-s", aln,
             "-m", "GTRGAMMA",
             "-w", os.path.abspath(self.tree_dir7)],
            stdout=subprocess.DEVNULL
        )

        subprocess.run(
            ['raxmlHPC-PTHREADS', '-T', '8', "-m", "GTRGAMMA", '-p', '1337', "-f", "b",
             "-t", os.path.join(self.tree_dir7, "RAxML_bestTree.initial_tree"),
             "-z", os.path.join(self.tree_dir7, "RAxML_bootstrap.bootstrap"),
             "-n", "final",
             "-w", os.path.abspath(self.tree_dir7)],
            stdout=subprocess.DEVNULL
        )

        with open(os.path.join(self.tree_dir7, "RAxML_bipartitionsBranchLabels.final"), "r") as tree_handle:
            with open(os.path.join(self.tree_dir7, "MAGE_tree_output.final"), "w") as mage_handle:
                mage_handle.write(tree_handle.read())

    def run(self):
        print(self.csv)
        self.force_removal(self.force)
        self.trim(self.reformat(self.path))
        self.assemble()
        self.consensus()
        self.reference_format()
        self.alignment()
        self.trimal()
        self.supermatrix()
        self.RAxML()


def main():
    parser = argparse.ArgumentParser(description="MLSA Generator", formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    optional = parser.add_argument_group('Optional Arguments')
    required = parser.add_argument_group('Required Arguments')

    required.add_argument('-p', '--path', type=str, required=True, help="folder of abi files")
    required.add_argument('-c', '--csv', type=str, required=True, help="spreadsheet of forward read, reverse read, sample ID and target region")
    optional.add_argument("-f", "--force", default=False, action='store_true', required=False, help="Force overwrite of previous run")
    optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

    args = parser.parse_args()
    path = os.path.normpath(args.path)
    csv_input = args.csv
    force = args.force

    job = Mage(path, csv_input, force)
    job.run()


if __name__ == '__main__':
    main()
