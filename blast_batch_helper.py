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

version = '0.7.1'

fasta_file = args.query
blast_output = args.out
check_interval = 20
nowtime = time.strftime('[%Y/%m/%d %H:%M:%S]')

fasta_ids = []

def parse_fasta_id():
	"""
	Parse fasta and keep query ids in global fasta_ids
	"""
	global fasta_ids
	with open (fasta_file, 'r') as f:
		for line in f.readlines():
			if line.startswith('>'):
				fasta_ids.append(line.split()[0][1:])
	print('Total fasta counts: '+ str(len(fasta_ids)))

def parse_blast_id():
	"""
	Parse blast and output hits ids list
	"""
	blast_ids = []
	if os.path.exists(blast_output):
		with open (blast_output, 'r') as f:
			for line in f.readlines():
				hit_line = str(line).strip().split('\t')
				if len(hit_line) == 12:
					blast_ids.append(str(line).split('\t')[0])
				else:
					print('\nERROR: Invalid hit in ouputfile:\n'+ str(line))
					print('Please check if the ouputfile is in default fmt6 format. Or it contains incomplete hits.')
					raise SystemExit(0)
	else:
		print('No BLAST output yet. Skipping parse.')
	return blast_ids

def finish_and_unfinished_id():
	"""
	Output finish and unfinished query ids list
	"""
	finish_ids = []
	re_check_ids = []
	no_hit_ids = []
	unfinished_ids = []
	blast_ids = parse_blast_id()
	blast_ids_set = set(blast_ids)	# use set to speed up
	for query in fasta_ids:
		if query in blast_ids_set:
			finish_ids.append(query)
		else:
			no_hit_ids.append(query)
	
	if len(finish_ids)>0:
		last_hit_index = fasta_ids.index(finish_ids[-1])
	else:
		last_hit_index = 0
	
	for no_hit in no_hit_ids:
		if fasta_ids.index(no_hit)<last_hit_index:
			re_check_ids.append(no_hit)
		else:
			unfinished_ids.append(no_hit)

	# print('finish_ids:%d\nre_check_ids:%d\nunfinished_ids:%d'% (len(finish_ids), len(re_check_ids), len(unfinished_ids)))
	return finish_ids, re_check_ids, unfinished_ids

def parse_tmp_id():
	"""
	Parse blast tmp file and output only complete hit's id in a list
	"""
	tmp_ids = []
	if os.path.exists(blast_output+'.tmp'):
		with open (blast_output+'.tmp', 'r') as f:
			for line in f.readlines():
				hit_line = str(line).strip().split('\t')
				if len(hit_line) == 12:
					tmp_ids.append(hit_line[0])
	else:
		print('No BLAST tmp yet. Skipping parse.')
	return tmp_ids

def last_blast_result():
	"""
	Check last blast hit from parse_blast_id().
	"""
	finish_ids = parse_blast_id()
	if len(finish_ids) > 0 :
		finished_fasta = fasta_ids.index(finish_ids[-1])+1
		print('Last hit: %s'% (finish_ids[-1]))
		print('Finished percentage: %.02f %% (%d/%d)' % (finished_fasta/len(fasta_ids)*100, finished_fasta, len(fasta_ids)))
	else:
		finished_fasta = 0

def prepare_subfasta():
	"""
	Prepare sub fasta file from unfinished work & query w/o any hit for continuing blasting job
	Return the starting id of fasta file (str starts from id_)
	Return origin fasta file path if there is no result yet
	"""
	finish_ids, re_check_ids, unfinished_ids = finish_and_unfinished_id()
	print('\nQuery w/  hit counts: '+ str(len(finish_ids)))
	print('Query w/o hit counts: '+ str(len(re_check_ids)+len(unfinished_ids)))
	
	if len(finish_ids) > 0:
		last_hit = finish_ids[-1]
		if len(unfinished_ids) > 0:
			start_id = unfinished_ids[0]
		else:
			start_id = last_hit
		subfasta_file = fasta_file+'.from_'+start_id+'.fasta'
		if os.path.exists(subfasta_file):
			os.remove(subfasta_file)
			print('Remove old duplicate subfasta file.')
		
		print('\nLast BLAST+ stopped at '+ start_id + '\nBut ALL query w/o hit before it will run BLAST again.\nGenerating subfasta file... It may take a while.\n')
		subfile = False
		with open (fasta_file, 'r') as f:
			for line in f.readlines():
				if line[0] == '>':
					if line.split()[0][1:] in re_check_ids or line.split()[0][1:] in unfinished_ids:
						subfile = True
					else:
						subfile = False
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
		print('\n%s\nCheck BLAST status...'% nowtime)
		extract_blast_output('tmp')
		last_blast_result()
		predict_finish_time()
		print('BLASTing...update in %d secs.'% check_interval)
		time.sleep(check_interval)
	print('\nBLAST finished. Extracting final results...')
	extract_blast_output('final')	# extract the remaining hits after BLAST finish
	write_ok_mark(blast_output)
	clean_tmp_file(blast_output)
	print('\n%s\nBLAST job is done.'% nowtime)

def predict_finish_time():
	timing_file = open(blast_output+'.timing', 'r')
	time_start = float(timing_file.read())
	timing_file.close()
	if os.path.exists(blast_output+'.tmp'):
		tmp_hit = parse_tmp_id()
		tmp_hit_set = set(tmp_hit)
		if len(tmp_hit_set) > 0:
			tmp_index_set = set()
			for hit in tmp_hit_set:
				tmp_index_set.add(fasta_ids.index(hit))
			num_finished_fasta = max(tmp_index_set)-min(tmp_index_set)+1
			num_wait_fasta = len(fasta_ids) - max(tmp_index_set)
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
			print('Finish time is predicted: '+ str(time.strftime('%Y/%m/%d %H:%M:%S',(time.localtime(time_finish)))))
	else:
		print('Not enough hit for finish time prediction')

def extract_blast_output(arg):
	"""
	Extract blast result from tmp file to blast_output file.
	"""
	existing_id = parse_blast_id()
	existing_id_set = set(existing_id)	# Use set to speed up
	tmp_existing_id = parse_tmp_id()	# Build query id set from tmp file
	if len(tmp_existing_id) > 0:
		tmp_last_id = tmp_existing_id[-1]
		with open (blast_output+'.tmp', 'r') as f:	# Extract hit from tmp file
			for line in f.readlines():
				result = line.split('\t')
				if arg == 'tmp':
					if (result[0] not in existing_id_set):
						if result[0] == tmp_last_id:	# Stop at current query to avoid incomplete BLAST hit
							break
						else:
							file = open(blast_output, 'a')
							file.write(line)
							file.close()
				elif arg == 'final':
					if (result[0] not in existing_id_set):	# Extract all the new query
						file = open(blast_output, 'a')
						file.write(line)
						file.close()
				else:
					print('ERROR: unknown arg of extract_blast_output: '+arg)
	else:
		print('No BLAST tmp yet. Skipping extraction.')

def main():
	print('\nBLAST-batch-helper v%s\n'%(version))
	if os.path.exists(blast_output):
		print(blast_output + ' exists.')
		if os.path.exists(blast_output+'.ok'):
			print('BLAST+ work was already finished.')
		else:
			print('Last BLAST+ work was not finished.\nContinuing BLAST+ work')
			parse_fasta_id()
			blast_work(prepare_subfasta())
	else:
		print(blast_output + ' does not exist.\nStarting BLAST+ work')
		parse_fasta_id()
		blast_work(fasta_file)


if __name__ == '__main__':
	main()