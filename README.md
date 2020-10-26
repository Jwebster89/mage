# MAGe
[M]ultilocus [A]lignment [Ge]nerator

## Introduction
MAGe takes a set of forward and reverse Sanger sequencing reads with the .ab1 extension, a set of reference fasta files and a csv file to generate a multilocus sequence alignment and produce a RAxML tree from the concatenated supermatrix.

## Quick Usage
All Sanger reads and reference fasta files are place in a single directory and specified with the `-p`/`--path` option. MAGe will produce all output files in subdirectories of this path. The csv file specified with `-c`/`--csv` names all files and target regions that MAGe uses to build the final alignment.

### CSV format
The format of the csv is 4 columns with forward sequence, reverse sequence, sample name and target region. Sequence names should **match the filename exactly without the extension**. E.g. in the example below, "1" relates to file 1.ab1 and "2" relates to file 2.ab1.

For reference files, column two should again include the filename minus the extension, e.g. 'Ref 1 BT1' for filename 'Ref 1 BT1.fasta'. Note that filenames may have whitespace characters. Column one should just contain the word "Reference"

Columns 3 and 4 contain sample names and target regions respectively. Sample names must remain consistent or they will be treated as separate samples, e.g. "Sample 1" and "Sample.1" will be treated separately. Sample names may also contain whitespace. Target names must also remain consistent or will be treated as separate amplicons. Target names are not currently supported with whitespace or -. Best practice currently is to include only alphanumeric characters and underscores, as below.
|           |            |              |              |
|-----------|------------|--------------|--------------|
| 1         | 2          | Sample 1     | Beta_Tubulin |
| 3         | 4          | Sample 2     | Beta_Tubulin |
| 5         | 6          | Sample 1     | TEF          |
| 7         | 8          | Sample 2     | TEF          |
| Reference | Ref 1 BT   | Reference 1  | Beta_Tubulin |
| Reference | Ref 1 TEF  | Reference 1  | TEF          |
| Reference | Ref 2 BT   | Reference 2  | Beta_Tubulin |
| Reference | Ref 2 TEF  | Reference 2  | TEF          |

## Output
Seven subfolders of intermediate files are produced in the `--path` directory.
- 1_fastq_files: Contains ab1 files converted to fastq
- 2_trimmed_reads: Quality filtered fastq reads
- 3_alignmed_reads: Mafft alignment of forward and reverse reads
- 4_sanger_consensus: Consensus files from the previous Mafft alignment. Reference files are also moved to this folder and their headers renamed.
- 5_alignments: A multi fasta of each target region for input into Mafft (.fasta), a Mafft alignment (.aln) and a trimmed alignment (.nex)
- 6_supermatrix: A concatenation of the nexus files in the previous step as well as a conversion of the nexus file back to an alignment (.aln).
- 7_final_tree: The output files from the RAxML steps and a copy of the tree renamed to "MAGE_tree_output.final" as a final output file.


## Options and Usage
```
usage: MAGe.py -p PATH -c CSV [-f] [-h]

Multilocus Alignment Generator

 __  __    _    ____
|  \/  |  / \  / ___| ___
| |\/| | / _ \| |  _ / _ \       /\
| |  | |/ ___ \ |_| |  __/     ------
|_|  |_/_/   \_\____|\___|    (∩｀-｀)⊃━☆ﾟ.*･｡ﾟ


Required Arguments:
  -p PATH, --path PATH  folder of abi files
  -c CSV, --csv CSV     spreadsheet of forward read, reverse read, sample ID and target region

Optional Arguments:
  -f, --force           Force overwrite of previous run
  -h, --help            show this help message and exit
```
