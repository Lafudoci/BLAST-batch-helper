# BLAST-batch-helper
A script for better handle BLAST job on large fasta query.

Highlight:
* Automatically resume last unfinished BLAST job (reconized by output filename).
* Periodical report the last blast progress.
* Predict BLAST finish time (Work in progress).

## Requirements
* Working [BLAST+ Command Line Applications](https://www.ncbi.nlm.nih.gov/books/NBK279671/), environment variable PATH needs to be set.
* [BLAST database](ftp://ftp.ncbi.nlm.nih.gov/blast/db/), environment variable BLASTDB needs to be set.
* Python >= 3.5

## Usage
```
usage: blast_batch_helper.py  aln_prog -db DB -num_threads NUM_THREADS -query QUERY -out OUT [-others OTHERS] [-h]

required arguments:
  aln_prog            Aln program: eg. blasx, blastn, blastp.
  -db                 BLAST database name, eg. nr, swissprot.
  -num_threads        Number of thread used to BLAST.
  -query              Path to fasta file.
  -out                Path to BLAST result output.

optional arguments:
  -others OTHERS      Pass other BLAST args.
  -h, --help          Show this help message and exit
```
## Example
```
python blast_batch_helper.py \
blastx \
-db nr \
-num_threads 10 \
-query fasta_all.fasta \
-out fasta_all_blastx_nr.txt \
-others "-task blastx-fast -max_target_seqs 1 -evalue 1e-3"
```
The example script will pass following command to BLAST:
`blastx -db nr -num_threads 10 -query fasta_all.fasta -out fasta_all_blastx_nr.txt -outfmt 6 -task blastx-fast -max_target_seqs 1 -evalue 1e-3 
`
## Output
Output formate uses BLAST argument `-outfmt 6`
```
TRINITY_DN89960_c0_g1_i1	F8VQB6.1	86.667	75	10	0	1	225	490	564	1.21e-027	106
TRINITY_DN89978_c0_g1_i1	P52179.2	59.140	93	38	0	2	280	1214	1306	3.36e-033	123
TRINITY_DN53731_c0_g1_i1	Q4KTA7.1	81.690	71	13	0	2	214	639	709	2.84e-036	131
TRINITY_DN76613_c0_g1_i1	Q99NI3.1	39.683	63	38	0	191	3	792	854	2.87e-009	54.3
TRINITY_DN58502_c0_g1_i1	Q5E9R3.1	86.923	520	68	0	421	1980	12	531	0.0	951
TRINITY_DN21354_c0_g1_i1	Q9YGK2.1	71.946	221	57	2	74	721	1	221	1.28e-075	232
```
