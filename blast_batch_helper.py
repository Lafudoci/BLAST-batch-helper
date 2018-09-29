import os
import subprocess
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('aln_prog', help="Aln program: blasx, blastn, blastp.")
parser.add_argument('-db', help="BLAST database name.", required=True)
parser.add_argument('-query', help="Path to fasta file.", required=True)
parser.add_argument('-out', help="Path to blast result output.", required=True)
parser.add_argument('-gnu_parallel', help="Use GNU parallel.", action='store_true', required=False)
parser.add_argument('-gnu_parallel_b', help="Set GNU parallel block size.", default= '100k', required=False)
parser.add_argument('-gnu_parallel_j', help="Set GNU parallel job number.",required=False)
parser.add_argument('-others', help="Pass other blast args.")
args = parser.parse_args()


fasta_file = args.query
blast_output = args.out

fasta_ids = []

def parse_fasta():
	with open (fasta_file, 'r') as f:
		for line in f.readlines():
			if line.startswith('>'):
				fasta_ids.append(line.split(' ')[0][1:])
	print('Total fasta counts: '+ str(len(fasta_ids)))
	return fasta_ids

def parse_blast_id():
	blast_ids = []
	if os.path.exists(blast_output):
		with open (blast_output, 'r') as f:
			for line in f.readlines():
				blast_ids.append(str(line).split('\t')[0])
	else:
		print('No BLAST output yet. Skipping parse.')

	return blast_ids

def blast_last_result():
	"""
	Check last blast result from parse_blast_id() and return last hit id
	Return -1 if there is no result yet
	"""
	blast_ids = parse_blast_id()
	if len(blast_ids) > 0 :
		print('Total hits: '+ str(len(blast_ids)))
		print('Last hit: '+blast_ids[-1])
		finished_fasta = fasta_ids.index(blast_ids[-1])+1
	else:
		finished_fasta = 0
		return -1
	print('Finished fasta: '+ str(finished_fasta))
	print('Finished percentage: %.02f %% (%d/%d)' % (finished_fasta/len(fasta_ids)*100, finished_fasta, len(fasta_ids)))
	return blast_ids[-1]

def prepare_subfasta():
	"""
	Prepare sub fasta file from unfinished work for continuing blasting job
	Return the starting id of fasta file (str starts from id_)
	Return origin fasta file path if there is no result yet
	"""
	last_hit = blast_last_result()
	if last_hit != -1:
		start_id = fasta_ids[fasta_ids.index(last_hit)+1]
		subfasta_file = fasta_file+'.from_'+start_id+'.fasta'
		if os.path.exists(subfasta_file):
			os.remove(subfasta_file)
			print('Remove old duplicate subfasta file.')
		
		print('BLAST+ will start from: '+ start_id + '\nGenerating subfasta file... It may take a while.')
		subfile = False
		with open (fasta_file, 'r') as f:
			for line in f.readlines():
				if line.startswith('>'+start_id):
					subfile = True
				if subfile == True:
					file = open(subfasta_file, 'a')
					file.write(line)
					file.close()
		return 'id_' + start_id
	else:
		print('BLAST+ will start from origin fasta file')
		return fasta_file

def write_ok_mark(file_name):
	file = open(file_name+'.ok', 'w')
	file.close()

def write_timing_mark(file_name):
	file = open(file_name+'.timing', 'w')
	file.write(str(time.time()))
	file.close()

def clean_tmp_file(file_name):
	if os.path.exists(file_name+'.tmp'):
		os.remove(file_name+'.tmp')
	if os.path.exists(file_name+'.timing'):
		os.remove(file_name+'.timing')

def blast_work(fasta_query):
	"""
	If input starts with id_ means it's a continuing work and a subfasta file will be made.
	"""
	clean_tmp_file(blast_output)
	if fasta_query.startswith('id_'):
		query_command = '-query %s.from_%s.fasta' % (fasta_file, fasta_query[3:])
	else:
		query_command = '-query %s' % (fasta_query)
	out_command = '-out %s.tmp' % (blast_output)
	if args.gnu_parallel == True:
		print('\nGNU parallel is enabled. The statement from authors:\nWhen using programs that use GNU Parallel to process data for publication please cite:\nO. Tange (2011): GNU Parallel - The Command-Line Power Tool,\n;login: The USENIX Magazine, February 2011:42-47.\n')
		gnu_parallel_args = "parallel --no-notice --recstart '>' --block %s --pipe" % (args.gnu_parallel_b)
		if args.gnu_parallel_j:
			gnu_parallel_args = '%s -j %s' % (gnu_parallel_args, args.gnu_parallel_j)
		blast_command = "cat %s | %s %s -db %s -outfmt 6 %s -query -> %s" % (query_command[7:], gnu_parallel_args, args.aln_prog, args.db, args.others, out_command[5:])
	else:
		blast_command = '%s -db %s %s %s -outfmt 6 %s' % (args.aln_prog, args.db, query_command, out_command, args.others)
	print('Executing: ' + blast_command)
	write_timing_mark(blast_output)
	blast_process = subprocess.Popen(blast_command, shell=True)
	while (blast_process.poll()==None):
		extract_blast_output()
		blast_last_result()
		predict_finish_time()
		print('BLASTing...\n')
		time.sleep(20)
	print('BLAST finished')
	extract_blast_output()	# extract the remaining hits after BLAST finish
	write_ok_mark(blast_output)
	clean_tmp_file(blast_output)

def predict_finish_time():
	tmp_hit = []
	timing_file = open(blast_output+'.timing', 'r')
	time_start = float(timing_file.read())
	timing_file.close()
	if os.path.exists(blast_output+'.tmp'):
		with open(blast_output+'.tmp', 'r') as f:
			for line in f.readlines():
				tmp_hit.append(str(line).split()[0])
		if len(tmp_hit) > 0:
			num_finished_fasta = fasta_ids.index(tmp_hit[-1]) - fasta_ids.index(tmp_hit[0]) + 1
			num_wait_fasta = len(fasta_ids) - fasta_ids.index(tmp_hit[-1])
			time_spent = time.time() - time_start
			blast_speed_per_sec = num_finished_fasta / time_spent
		else:
			num_finished_fasta = 0

		if num_finished_fasta == 0:
			print('Not enough hit for finish time prediction')
		else:
			time_remaining = num_wait_fasta / blast_speed_per_sec
			time_finish = time.time() + time_remaining
			print('Current BLASTing speed: %d seqs/hour.' % int(blast_speed_per_sec*3600))
			print('Finish time is predicted: '+ str(time.asctime(time.localtime(time_finish))))
	else:
		print('Not enough hit for finish time prediction')

def extract_blast_output():
	"""
	Extract blast result from tmp file to blast_output file.
	"""
	if os.path.exists(blast_output+'.tmp'):
		existing_id = parse_blast_id()
		with open (blast_output+'.tmp', 'r') as f:
			for line in f.readlines():
				result = []
				result = line.split()
				if (result[0] not in existing_id):
					file = open(blast_output, 'a')
					file.write(line)
					file.close()
					# print('New hit: '+ result[0] + ' ' + result[1])
	else:
		print('No Blast output yet. Skipping extraction.')

def main():
	if os.path.exists(blast_output):
		print(blast_output + ' exists.')
		if os.path.exists(blast_output+'.ok'):
			print('Blast+ work was already finished.')
		else:
			print('Last Blast+ work was not finished.\nContinuing Blast+ work')
			parse_fasta()
			blast_work(prepare_subfasta())
	else:
		print(blast_output + ' does not exist.\nStarting BLAST+ work')
		parse_fasta()
		blast_work(fasta_file)


if __name__ == '__main__':
	main()