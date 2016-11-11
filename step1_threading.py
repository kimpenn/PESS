#==================================================================================================
# Sarah Middleton
# Junhyong Kim Lab
# University of Pennsylvania
# Last update: Nov 2016
#==================================================================================================
# step1_threading.py
# Threads protein sequences (provided in FASTA format) against a set of 1,814 representative 
# templates using the CNFalign_lite module of RaptorX. Outputs a table of scores (.scoremat)
# 
#
# Notes:
# - IMPORTANT: This script must be run from within the directory that holds the RaptorX executables, 
#   otherwise those executables won't run properly.
# - Individual fasta files are created for each input sequence in a folder called /indiv_fasta/
# - Feature profiles are built, if they don't already exist, in /tgt_files/
# - These folders are created in same folder as the fasta file unless otherwise specified with --out
# - If the feature files (phase I) already exist, the script can start at threading (phase II). You
#   will be prompted to choose whether to do this when the script detects that /tgt_files/ exists.
# 
# Usage:
#    python step1_threading.py FASTA [options]
#
# Examples:
#    python ~/pess/step1_threading.py ~/pess/demo/demo.fa
#    python ~/pess/step1_threading.py ~/pess/demo/demo.fa --cpu=32 --out="~/pess/demo/demo_results"
#==================================================================================================
import subprocess, sys, os, time
from multiprocessing import Pool, cpu_count
from optparse import OptionParser


#------------------------------------------------------------------
# Score an individual seq against all templates. Calls RaptorX.
# Parameters: list with format [seqID, TPL_DIR, tgtOut, [tplList]]
# Returns: list with format [seqID, {resultDictionary}]
#------------------------------------------------------------------
def score_seq(params):
	seqID = params[0]
	tplDir = params[1]
	tgtDir = params[2]
	tplList = params[3]
	outDir = params[4]
	
	results = {}
	for tplID in tplList:
		command = "./CNFalign_lite -t %s -q %s -l %s -g %s" % (tplID, seqID, tplDir, tgtDir)
		job = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		jobOutput = []
		for line in job.stdout:
			jobOutput.append(line)
		result = job.wait()
		
		if result != 0:
			outStr = ">>Error returned by CNFalign (non-zero return status).\n"
			outStr += "  Command:" + command + "\n"
			outStr += str(jobOutput)
			print(outStr)
		
		#print jobOutput
		things = jobOutput[0].split()
		score = float(things[2])
		results[tplID] = score
	
	# save results to output file in case script gets terminated early
	outFile = "%s%s.scores" % (outDir, seqID)
	outStr = ""
	for tplID in results:
		outStr += "%s\t%s\n" % (tplID, results[tplID])
	outs = open(outFile, 'w')
	outs.write(outStr)
	outs.close()
	
	print("Finished %s" % seqID)
	
	return [seqID, results]
	

#------------------------------------------------------------------
# Build the tgt file for a sequence. Calls RaptorX and BLAST.
# Parameters: list with format [seqFile, tgtFile]
# Returns: result code
#------------------------------------------------------------------
def build_tgt(params):
	seqFile = params[0]
	tgtFile = params[1]
		
	command = "./buildFeature -i %s -o %s -c 1" % (seqFile, tgtFile)
	job = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	jobOutput = job.stdout.read()
	result = job.wait()
	
	if result != 0:
		outStr = ">>Error returned by buildFeature (non-zero return status).\n"
		outStr += "  Command:" + command + "\n"
		outStr += jobOutput
		print(outStr)
		
	print("Finished %s" % seqFile)
	
	return result	

	
#------------------------------------------------------------------
# Read a fasta format file into a dictionary (dict[header] = sequence)
# Parameters: path to fasta file
# Returns: dictionary & whether there was an error (boolean)
#------------------------------------------------------------------
def read_fasta(fileName):
	seqs = {}
	error = False
	
	try:
		ins = open(fileName, 'r')
	except IOError:
		error = True
		print("Error: in read_fasta(): could not open", fileName)
	else:
		id = ""
		for line in ins:
			line = line.rstrip('\r\n')
			if (len(line) > 0) and (">" == line[0]):
				id = line[1:]
				if id in seqs:
					print("Warning: in read_fasta(): repeat id (%s) in file. Overwriting." % id)
				seqs[id] = ""
			else:
				seqs[id] += line.upper()
	
		ins.close()
		
	return (seqs, error)	

#------------------------------------------------------------------
# Main - set things up, launch jobs, compile results.
#------------------------------------------------------------------
if __name__ == '__main__':
	
	# parse args
	usageMsg = "Usage: %prog FASTA [options]"
	parser = OptionParser(usage=usageMsg)
	parser.add_option("--cpus", action="store", type='int', default=1, dest="MAX_CPU", help="Maximum number of CPUs to use. Default is [%default].")
	parser.add_option("--out", action="store", default=None, dest="OUT_DIR", help="Path to desired output directory. Default is the directory of the input fasta.")
	
	# read/process args
	(opts, args) = parser.parse_args()
	if len(args) == 1:
		SEQ_FILE = args[0]
	else:
		print(">> Missing input file. Use -h for help. Exiting.")
		sys.exit()
	if opts.OUT_DIR == None:
		workingDir = os.path.dirname(SEQ_FILE)
	else:
		workingDir = opts.OUT_DIR
	
	
	# set up additional file names
	scriptDir = os.path.dirname(sys.argv[0]) #gives relative path to script, w/o last '/'
	filename = os.path.basename(SEQ_FILE)
	fnParts = filename.split(".")
	seqdbName = fnParts[0]
	indivOut = "%s/indiv_fasta/" % (workingDir)
	tgtOut = "%s/tgt_files/" % (workingDir)
	scoreOut = "%s/score_files/" % (workingDir)
	scorematOut = "%s/%s.scoremat" % (workingDir, seqdbName)
	skipTgt = False
	
	# templates
	TPL_LIST = "%s/templates/reference_tpl_list" % scriptDir
	TPL_DIR = "%s/templates/CAL_TPL/" % scriptDir
	
	# check that everything exists
	if not os.path.exists(SEQ_FILE):
		print(">> Error: could not find indicated fasta file.")
		print("   (tried: %s)" % SEQ_FILE)
		print(">> Exiting.")
		sys.exit()
	if not os.path.exists(TPL_LIST):
		print(">> Error: could not find TPL_LIST.")
		print("   (tried: %s)" % TPL_LIST)
		print(">> Exiting.")
		sys.exit()
	if not os.path.exists(TPL_DIR):
		print(">> Error: could not find TPL_DIR.")
		print("   (tried: %s)" % TPL_DIR)
		print(">> Exiting.")
		sys.exit()
	
	print("")
	print("Files/paths to be used:\n")
	print("    output directory | %s" % workingDir)
	print("       sequence file | %s" % SEQ_FILE)
	print("            tpl list | %s" % TPL_LIST)
	print("")
	
	
	# create output directory if necessary
	if not os.path.exists(workingDir):
		print("Output directory %s does not exist, creating." % workingDir)
		os.makedirs(workingDir)
		
	if os.path.exists(indivOut):
		print("")
		print("Output directory %s already exists." % indivOut)
		response = raw_input("Ok to overwrite existing output files? (y/n) ")
		if response != "y":
			print("Exiting.")
			sys.exit()
	else:
		print("Output directory %s does not exist, creating." % indivOut)
		os.makedirs(indivOut)
		
	if os.path.exists(tgtOut):
		print("")
		print("Output directory %s already exists." % tgtOut)
		response = raw_input("Skip feature profile (.tgt) file generation? (y/n) ")
		if response == "y":
			skipTgt = True
	else:
		print("Output directory %s does not exist, creating." % tgtOut)
		os.makedirs(tgtOut)
		
	if os.path.exists(scoreOut):
		print("")
		print("Output directory %s already exists." % scoreOut)
		response = raw_input("Ok to overwrite existing output files? (y/n) ")
		if response != "y":
			print("Exiting.")
			sys.exit()
	else:
		print("Output directory %s does not exist, creating." % scoreOut)
		os.makedirs(scoreOut)
	
	
	start = time.time()
	
	# read in fasta
	(seqs, error) = read_fasta(SEQ_FILE)
	if error:
		print(">> Error reading fasta sequence file. Exiting.")
		sys.exit()
	
	# print each seq to a separate file
	idList = []
	for id in seqs:
		if id in idList:
			print(">> Warning: id already in list:", id)
			print("   This sequence will be overwritten.")
		else:
			idList.append(id)
		
		outFile = "%s%s.fa" % (indivOut, id)
		outs = open(outFile, 'w')
		outStr = ">%s\n%s" % (id, seqs[id])
		outs.write(outStr)
		outs.close()
	
	# read list of templates
	tplList = []
	ins = open(TPL_LIST, 'r')
	for line in ins:
		line = line.rstrip('\r\n')
		tplList.append(line)
	ins.close
	tplList.sort()
	
	print("")
	print("Read in %s sequences." % len(idList))
	print("Read in %s templates." % len(tplList))
	print("")
	print("Starting searches. This can take several minutes per sequence.")
	print("")
	
	# create .tgt files
	print("")
	print("======================================================")
	print("")
	print("Phase I: Creating feature files...")
	tgtFinishedCount = 0
	tgtStart = time.time()
	if skipTgt == True:
		print(">> Skipping tgt file creation.")
	else:
		jobList = []
		for seqID in idList:
			seqFile = "%s%s.fa" % (indivOut, seqID)
			tgtFile = "%s%s.tgt" % (tgtOut, seqID)
			params = [seqFile, tgtFile]
			jobList.append(params)
		
		# create a pool of processes to run scoring in parallel
		pool1 = Pool(processes = opts.MAX_CPU)
		result1 = pool1.map_async(build_tgt, jobList)
		result1.wait()		
			
	
	tgtElapsed = time.time() - tgtStart
	tgtElapsedHrs = (float(tgtElapsed) / 60) / 60
	print("")
	print("Time elapsed: %.2f s (%.2f hr)" % (tgtElapsed, tgtElapsedHrs))
	print("")
	print("======================================================")
	
	# create param package for each seq [seqID, tplDir, tgtDir, [tplIDs]]
	execList = []
	for seqID in idList:	
		tgtFile = "%s%s.tgt" % (tgtOut, seqID)
		if os.path.exists(tgtFile):
			params = [seqID, TPL_DIR, tgtOut, tplList, scoreOut]
			execList.append(params)
		else:
			print("Skipping", seqID, "-- no .tgt file.")
	
	print("")
	print("Phase II: RaptorX threading process.")
	print("%s processes will be created." % opts.MAX_CPU)
	print("%s jobs will be assigned to the process pool." % len(execList))
	print("")
	
	# create a pool of processes to run scoring in parallel
	pool = Pool(processes = opts.MAX_CPU)
	scoreStart = time.time()
	result = pool.map_async(score_seq, execList)
	result.wait()
	elapsedTime = time.time() - scoreStart
	elapsedHrs = (float(elapsedTime) / 60) / 60
	
	print("")
	print("Time elapsed: %.2f s (%.2f hr)" % (elapsedTime, elapsedHrs))
	print("")
	print("======================================================")
	
	# compile scores
	results = result.get()
	scoreMat = {}
	for entry in results:
		seqID = entry[0]
		hash = entry[1]
		scoreMat[seqID] = {}
		for tplID in hash:
			scoreMat[seqID][tplID] = hash[tplID]
	
	# print scoremat file
	outs1 = open(scorematOut, 'w')
	header = "\t".join(tplList)
	outs1.write(header + "\n")
	for id in sorted(scoreMat):
		outStr = id
		for tplID in tplList:
			outStr += "\t%.2f" % scoreMat[id][tplID]
		outs1.write(outStr + "\n")
	outs1.close()
	
	totalTime = time.time() - start
	totalHrs = (float(totalTime) / 60) / 60
	
	print("")
	print("Finished. Total time: %.2f s (%.2f hr)" % (totalTime, totalHrs))
	print("")

