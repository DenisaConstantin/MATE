# MATE (Mapper and Annotator of Transposable Elements)

**MATE** is a resource-efficient, bash-based bioinformatics pipeline designed for the high-resolution mapping and annotation of transposable element (TE) insertions. 

The pipeline runs TE queries against *de novo* assemblies via `BLASTN`. It unifies adjacent signals originating from the same insertion and extracts the insertion along with variable-sized genomic flanking sequences (up to 3000 bp upstream and downstream). 

## ⚙️ How It Works
By mapping these junction segments back to a reference genome, MATE categorizes the insertions into three main groups:
1. **Genotype-Specific Insertions:** Putative *de novo* or genotype-specific TEs present only in the target assembly (absent in the reference).
2. **Conserved and Annotated Insertions:** TEs present in both the target and the reference genome, matching existing TE annotations.
3. **Conserved but Unannotated Insertions:** TEs present in both genomes but lacking prior annotation in the reference database.

## 🛠️ Dependencies
Ensure the following tools are installed and available in your system's `$PATH`:
* `bash` & `awk`
* `blast+` (specifically `blastn`)
* `seqkit`

* ## 💻 Usage
The main execution script is `mate.sh`. It requires paths to your genotype assembly, reference genome, TE sequences, and a coordinate file.

./mate.sh -genotype <assembly.fasta> \
          -reference <reference.fasta> \
          -te <te_sequences.fasta> \
          -te_coord <te_coord.txt> \
          -out <output_folder_name>

### Arguments:
* `-genotype`: The *de novo* genome assembly FASTA file.
* `-reference`: The reference genome FASTA file used for comparison.
* `-te`: A FASTA file containing the sequences of the Transposable Elements of interest.
* `-te_coord`: A tab-separated reference file containing TE coordinates/annotations.
* `-out`: Name of the directory where results will be saved.
