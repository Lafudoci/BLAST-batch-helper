# BLAST-batch-helper
A script for better handle BLAST job on large fasta query.

Highlight:
* Automatically resume last unfinished BLAST job.
* Periodical report the last blast progress.
* Predict BLAST finish time.
* [NEW] Compatible with GNU parallel.

## Requirements
* Working [BLAST+ Command Line Applications](https://www.ncbi.nlm.nih.gov/books/NBK279671/), environment variable PATH needs to be set.
* [BLAST database](ftp://ftp.ncbi.nlm.nih.gov/blast/db/), environment variable BLASTDB needs to be set. (Not necessary when using remote BLAST)
* Python >= 3.5
* [GNU parallel](https://www.gnu.org/software/parallel/) (optional) 

## Usage
```
usage: blast_batch_helper.py  aln_prog -db DB -query QUERY -out OUT [-others OTHERS] [-h]

required arguments:
  aln_prog            Aln program: eg. blasx, blastn, blastp.
  -db                 BLAST database name, eg. nr, swissprot.
  -query              Path to fasta file.
  -out                Path to BLAST result output.

optional arguments:
  -others             Pass other BLAST args.
  -gnu_parallel       Use GNU parallel.
  -gnu_parallel_b     Set GNU parallel block size (Default is 100k).
  -gnu_parallel_j     Set GNU parallel job number.
  -no_rollback        Ignore previous queries w/o hit in a resuming job.  
  -h, --help          Show this help message and exit
```
Please note that `-outfmt 6` is hard coded in script, so do NOT include `-outfmt` in `-others`.

## Example
### Local BLAST
```
python blast_batch_helper.py \
blastx \
-db swissprot \
-query fasta_all.fasta \
-out fasta_all_blastx_swissprot.fmt6 \
-others "-task blastx-fast -evalue 1e-6 -num_threads 10"
```
The example script will pass following command to local BLAST+:
`blastx -db swissprot -query fasta_all.fasta -out fasta_all_blastx_swissprot.txt -outfmt 6 -task blastx-fast -evalue 1e-6 -num_threads 10`

### Remote BLAST
```
python blast_batch_helper.py \
blastx \
-db swissprot \
-query fasta_all.fasta \
-out fasta_all_blastx_swissprot.fmt6 \
-others "-task blastx-fast -evalue 1e-6 -remote"
```
If you like to run BLAST on NCBI sever instead of your local computer, just replace the `-num_threads` with `-remote`.

### Resume an unfinished BLAST
When you want to continue an unfinished BLAST job, just run the script with the same arguments again. The script will look for the same output filename from `-out` argument. If the output already exists, script parses the last BLAST hit out. Then new subfasta file will be made for continuing BLAST job. All the results will be extracted and save back into the original output file.

### Use GNU parallel (Optional)
```
python3 blast_batch_helper.py \
blastx \
-db swissprot \
-query fasta_all.fasta \
-out fasta_all_blastx_nr.txt \
-others "-task blastx-fast -evalue 1e-6"
-gnu_parallel \
-gnu_parallel_b 10k \
-gnu_parallel_j 10
```
If you like to use [GNU parallel](https://www.gnu.org/software/parallel/) to speed up the job (ONLY for local BLAST), just add `-gnu_parallel`. You could also set `gnu_parallel_b` and `gnu_parallel_j` to fit your condition.

## Output
Shows information in blasting
```
[2018/10/08 17:37:48]
Check BLAST status...
Last hit: TRINITY_DN2431_c1_g1_i2
Finished percentage: 0.88 % (83/9468)
Current BLASTing speed: 3596 seqs/hour.
Finish time is predicted: 2018/10/08 20:14:06
BLASTing...update in 20 secs.
```
Output formate uses BLAST argument `-outfmt 6`
```
TRINITY_DN89960_c0_g1_i1	F8VQB6.1	86.667	75	10	0	1	225	490	564	1.21e-027	106
TRINITY_DN89978_c0_g1_i1	P52179.2	59.140	93	38	0	2	280	1214	1306	3.36e-033	123
TRINITY_DN53731_c0_g1_i1	Q4KTA7.1	81.690	71	13	0	2	214	639	709	2.84e-036	131
TRINITY_DN76613_c0_g1_i1	Q99NI3.1	39.683	63	38	0	191	3	792	854	2.87e-009	54.3
TRINITY_DN58502_c0_g1_i1	Q5E9R3.1	86.923	520	68	0	421	1980	12	531	0.0	951
TRINITY_DN21354_c0_g1_i1	Q9YGK2.1	71.946	221	57	2	74	721	1	221	1.28e-075	232
```
