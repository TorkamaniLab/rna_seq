import subprocess as sp
import os
import sys
import time

## RNASEQ pipe for paired end alignment
## NOTE - single functionality to be added

class RNASEQ(object):

	def __init__(self,DataDir,OutDir,GroupsFile,ComparisonsFile,CellType=None):

		## Groups file
		## group	read1	read2
		## Comparison File
		## group1	group2

		self.ddir  = DataDir
		self.odir  = OutDir
		
		GFP = os.path.join(DataDir,GroupsFile)
		CFP = os.path.join(DataDir,ComparisonsFile)
		self.groups, self.comp = self.read_samplesheet(GFP,CFP)
		if CellType == None:
			self.cell = ''
		else:
			self.cell  = CellType
	
	def main_smarter(self):
		
		Groups = self.groups
		Comps  = self.comp
		CellT  = self.cell
		DDirec = self.ddir
		OutDir = self.odir

		CountFiles = {}
		CountJobs  = ''
		for group in Groups:
			
			CountFiles[group] = []
			i = 0	
			while i < len(Groups[group]['R1']):
		
				read1 = Groups[group]['R1'][i]
				read2 = Groups[group]['R2'][i]

				print "Submitting jobs for: "+read1+', '+read2

				## 1. Trim Adapters:
				trimjobID, readtF1, readtF2 = self.trim(DDirec,OutDir,read1,read2)
				
				## 2. Remove Smarter Randomer and Poly C/G
				cutjobID, readF1, readF2 = self.cutadapt(OutDir,OutDir,readtF1,readtF2,trimJobs=trimjobID)

				## 2. Run Aligner:
				alignJobID, alignF = self.aligner(OutDir,readF1,readF2,trimJob=cutjobID)
		
				## 3. Count reads	
				countJobID, countF = self.HTSeq_count(OutDir,alignF,alignJob=alignJobID)
				
				CountFiles[group].append(countF)
				CountJobs += ':'+countJobID

				i += 1
		
		#self.run_deseq(OutDir,Comps,CountFiles,countJobs = CountJobs)
	
	def main(self):
		
		Groups = self.groups
		Comps  = self.comp
		CellT  = self.cell
		DDirec = self.ddir
		OutDir = self.odir

		CountFiles = {}
		CountJobs  = ''
		for group in Groups:
			
			CountFiles[group] = []
			i = 0	
			while i < len(Groups[group]['R1']):
		
				read1 = Groups[group]['R1'][i]
				read2 = Groups[group]['R2'][i]

				print "Submitting jobs for: "+read1+', '+read2

				## 1. Trim Adapters:
				trimjobID, readF1, readF2 = self.trim(DDirec,OutDir,read1,read2)

				## 2. Run Aligner:
				alignJobID, alignF = self.aligner(OutDir,readF1,readF2,trimJob=trimjobID)
		
				## 3. Count reads	
				countJobID, countF = self.HTSeq_count(OutDir,alignF,alignJob=alignJobID)
				
				CountFiles[group].append(countF)
				CountJobs += ':'+countJobID

				i += 1
					
		#self.run_deseq(OutDir,Comps,CountFiles,countJobs = CountJobs)
	
	def submit_job(self,args):

		# Function to make sure a job has been submitted properly
		retry = 0
		args = ['qsub']+args
					
		while retry<5:
			p1     = sp.Popen(args,stdout=sp.PIPE)
			StdOut = p1.communicate()
			jobID  = StdOut[0].replace('\n','') 
			print jobID 

			try:
				job_number = int(jobID[0:jobID.index('.')]) #Make sure you have a proper job ID
				retry = 100 #Exit while loop
			except:
				retry += 1
				time.sleep(retry*60) #Sleep for some time and then retry
		
		return jobID
	
	def read_samplesheet(self,sampleSheetPath,compFilePath):
		
		##	Group	Read1	Read2
		groups = {}
		comps  = []

		with open(sampleSheetPath) as inp1:
			
			line = inp1.readline()
			while line:
				line = line.replace('\n','').split('\t')
				
				group = line[0]
				read1 = line[1]
				read2 = line[2]

				if group in groups:
					groups[group]['R1'].append(read1)
					groups[group]['R2'].append(read2)
				
				else:
					groups[group] = {'R1':[],'R2':[]}
					groups[group]['R1'].append(read1)
					groups[group]['R2'].append(read2)
				
				line = inp1.readline()
		
		with open(compFilePath) as inp2:
			
			line = inp2.readline()
                        while line:
				line = line.replace('\n','').split('\t')
				comps.append([line[0],line[1]])
				line = inp2.readline()
		
		return groups,comps

	def cutadapt(self,datadir,outdir,reads1,reads2,nbases='7',trimJobs=None):
		
		## NOTE - Use this instead of trimmomatic if you want to remove fixed n bases
		## reads1/2 - full filenames with the path
		input_1  = os.path.join(datadir,reads1)
		input_2  = os.path.join(datadir,reads2)
		output_1 = os.path.join(outdir,reads1.replace('.fastq.gz','.trim.fastq'))
		output_2 = os.path.join(outdir,reads2.replace('.fastq.gz','.trim.fastq'))
		
		with open('cut.job','w') as outp:
			outp.write('#!/bin/bash\n')
			outp.write('cutadapt --cut '+nbases+' -o '+output_1+'.temp '+input_1+' \n')
			outp.write('cutadapt --cut '+nbases+' -o '+output_2+'.temp '+input_2+' \n')
			outp.write('mv '+output_1+'.temp '+output_1+' \n')
			outp.write('mv '+output_2+'.temp '+output_2+' \n')

		if trimJobs == None:
			cutJobArgs = ['cut.job']
		else:
			cutJobArgs = ['-W','depend=afterok:'+trimJobs,'cut.job']

		cutJobID   = self.submit_job(cutJobArgs)

		return cutJobID, output_1, output_2
	
	def write_trim(self,arguments):

		## ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
		paired_output_1 = arguments[2] 
		paired_output_2 = arguments[4]
		
		arguments = (' ').join(arguments)

		with open('trim.job','w') as outp:
			outp.write('#!/bin/bash\n')
			outp.write("java -jar /gpfs/home/bhuvan/Programs/Trimmomatic-0.32/trimmomatic-0.32.jar PE "+arguments+" ILLUMINACLIP:/gpfs/home/bhuvan/Programs/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10:2:true LEADING:3 TRAILING:3 \n")
			outp.write('gzip -d '+paired_output_1+'\n')
			outp.write('gzip -d '+paired_output_2+'\n')
		
		return 'trim.job'

	def trim(self,datadir,outdir,reads1,reads2):
	
		## reads1/2 - full filenames with the path
		input_1 = os.path.join(datadir,reads1)
		input_2 = os.path.join(datadir,reads2)
		paired_output_1   = os.path.join(outdir,reads1.replace('.fastq.gz','.trim.fastq.gz'))
		unpaired_output_1 = os.path.join(outdir,reads1.replace('.fastq.gz','.unpair.fastq.gz'))
		paired_output_2   = os.path.join(outdir,reads2.replace('.fastq.gz','.trim.fastq.gz'))
		unpaired_output_2 = os.path.join(outdir,reads2.replace('.fastq.gz','.unpair.fastq.gz'))

		arguments = [input_1,input_2,paired_output_1,unpaired_output_1,paired_output_2,unpaired_output_2]

		trimJobFile = self.write_trim(arguments)
		trimJobArgs = [trimJobFile]
		trimJobID   = self.submit_job(trimJobArgs)
	
		paired_output_1 = os.path.join(outdir,reads1.replace('.fastq.gz','.trim.fastq'))
		paired_output_2 = os.path.join(outdir,reads2.replace('.fastq.gz','.trim.fastq'))
			
		return trimJobID, paired_output_1, paired_output_2


	def write_paired_aligner(self,readF1,readF2,OutDir,bowtie_ref):

		# Function to  write the file containing commands to run alignment using bowtie2
		# fdir = directory where the files are

		sep = os.sep
		if readF1.endswith('_R1.trim.fastq'):
			outfname = readF1.split(sep)[-1].replace('_R1.trim.fastq','.sam') 
			outfpath = os.path.join(OutDir,outfname)
		else:
			outfname = readF1.split(sep)[-1].replace('.trim.fastq','.sam') 
			outfpath = os.path.join(OutDir,outfname)

		with open('rna_bowtie2.job','w') as outp:
			outp.write('#!/bin/bash\n')
			outp.write('module load bowtie/2.2.3\n')
			outp.write('bowtie2 --very-sensitive -N 1 -p 8 -x '+bowtie_ref+' -q -1 '+readF1+' -2 '+readF2+" -S "+outfpath+'\n')
		
		return 'rna_bowtie2.job', outfname 

	def write_single_aligner(self,readF,OutDir,bowtie_ref):
		
		# Function to  write the file containing commands to run alignment using bowtie2
	
		sep = os.sep
		outfname = readF.split(sep)[-1].replace('_R1.trim.fastq','.sam')
		outfpath = os.path.join(OutDir,outfname)

		with open('rna_bowtie2.job','w') as outp:
			outp.write('#!/bin/bash\n')
			outp.write('module load bowtie/2.2.3\n')
			outp.write('cd '+DataDir+'\n')
			outp.write('bowtie2 --very-sensitive -N 1 -p 8 -x '+bowtie_ref+' -q '+readF+" -S "+outfpath+'\n')
				
		return 'rna_bowtie2.job', outfname

	def aligner(self,OutDir,readF1,readF2=None,trimJob=None,bowtie_ref="/gpfs/group/torkamani/bhuvan/Data/hg19/bowtie2_reference/hg19_ucsc"):

		# Function to run bowtie2 alignment for a list of files
		# NOTE - Requires full paths for the read files
		
		if readF2 == None:
			alignjob_file, align_file = self.write_single_aligner(readF1,OutDir,bowtie_ref)
		else:
			alignjob_file, align_file = self.write_paired_aligner(readF1,readF2,OutDir,bowtie_ref)

		if trimJob == None:
			args = ['-l', 'ncpus=8',alignjob_file]
		else:
			args = ['-W','depend=afterok:'+trimJob,'-l','ncpus=8',alignjob_file]
	
		jobID = self.submit_job(args)
	
		return jobID, align_file

	def write_count_job(self,datadir,fin,fout):
	
		gtfFile = '/gpfs/home/atorkama/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2011-08-30-21-45-18/Genes/genes.gtf'
		with open('htseq_count.job','w') as outp:
			outp.write('#!/bin/bash\n')
			outp.write('module load python\n')
			outp.write('cd '+datadir+'\n')
			outp.write('htseq-count '+fin+' '+gtfFile+' > '+fout+'.temp\n')

			## Modifying the count file to remove the last few lines of stats
			outp.write('head --lines -5 '+fout+'.temp > '+fout+'\n')
			#outp.write('rm '+fout+'.temp\n')
		
		return 'htseq_count.job'
	
	def HTSeq_count(self,datadir,fname,alignJob=None):
	
		outfile = fname.replace('.sam','.counts')
		count_job_file = self.write_count_job(datadir,fname,outfile)
		if alignJob == None:
			args = [count_job_file]
		else:
			args = ['-W','depend=afterok:'+alignJob,count_job_file]

		jobID = self.submit_job(args)
	
		return jobID, outfile

	def write_deseq_R(self,DataDir,testGrp,ctrlGrp,test_samples,ctrl_samples):
		
		## Function to write the R script file

		cellType   = self.cell 
		
		samples = []
		for i in test_samples:	
			samples.append(testGrp+'_'+str(test_samples.index(i)+1))
		for i in ctrl_samples:
			samples.append(ctrlGrp+'_'+str(ctrl_samples.index(i)+1))	

		samples = str(tuple(samples))
		Rscript = cellType+'_'+testGrp+'_'+ctrlGrp+'.R'
		scriptF = open(Rscript,'w')

		countfiles  = str(tuple(test_samples+ctrl_samples))
		test_number = str(len(test_samples))
		ctrl_number = str(len(ctrl_samples))	

		pdf_output  = os.path.join(DataDir,cellType+'_'+testGrp+'_'+ctrlGrp+'.pdf')
		diffExp_out = os.path.join(DataDir,cellType+'_'+testGrp+'_'+ctrlGrp+'_diff_exp.csv')
		
		Rcode  = ['library(DESeq)\n']
		Rcode += ['sampleData = data.frame(sample = c'+samples+', fnames = c'+countfiles+', condition = factor(c(rep("'+testGrp+'",'+test_number+'),rep("'+ctrlGrp+'",'+ctrl_number+')))) \n']
		Rcode += ['cds = newCountDataSetFromHTSeqCount(sampleData, directory = "'+DataDir+'")\n']
		Rcode += ['cds = estimateSizeFactors(cds)\n']
		Rcode += ['sizeFactors(cds)\n']
		Rcode += ['cds = estimateDispersions(cds)\n']
		Rcode += ['pdf("'+pdf_output+'")\n']
		Rcode += ['plotDispEsts(cds)\n']
		Rcode += ['res = nbinomTest( cds, "'+testGrp+'","'+ctrlGrp+'")\n']
		Rcode += ['plotMA(res)\nhist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")\ndev.off()\n']
		Rcode += ['write.csv( res, file="'+diffExp_out+'" )\n']
		
		scriptF.writelines(Rcode)
		scriptF.close()

		return Rscript

	def write_deseq_job(self,DataDir,testGrp,ctrlGrp,test_samples,ctrl_samples):
		
		Rscript = self.write_deseq_R(DataDir,testGrp,ctrlGrp,test_samples,ctrl_samples)
		with open('deseq.job','w') as outp:
			outp.write('#!/bin/bash\n')
			outp.write('module load R\n')
			outp.write('cd /gpfs/home/bhuvan/Chromatin_Remodelling/scripts\n')
			outp.write('echo Comparing '+testGrp+' '+ctrlGrp+'\n')
			outp.write('Rscript '+Rscript+'\n')
		
		return 'deseq.job'
							
	def run_deseq(self,DataDir,Comparisons,FileGroups,countJobs=None):

		jobList = ''
		for compar in Comparisons:
			
			test = compar[0]
			ctrl = compar[1]
			
			test_sam  = FileGroups[test]
			ctrl_sam  = FileGroups[ctrl]
			deseq_job = self.write_deseq_job(DataDir,test,ctrl,test_sam,ctrl_sam)
		
			if countJobs == None:
				args = [deseq_job]
			else:
	 			args = ['-W','depend=afterok'+countJobs,deseq_job]
	
			jobID = self.submit_job(args)
			jobList += ':'+jobID
			
		return jobList

		
if __name__ == '__main__':
	
	datadir = sys.argv[1]
	outdir  = sys.argv[2]
	groupsF = sys.argv[3]
	comparF = sys.argv[4] 		
	
	try:
		cell = sys.argv[5]	
	except IndexError:
		cell = None

	#datadir = "/gpfs/home/atorkama/cec/round2/"
	#outdir  = "/gpfs/home/bhuvan/Chromatin_Remodelling/scripts/test/"
	#cell_line = 'Nick'

	inst = RNASEQ(datadir,outdir,groupsF,comparF,cell)
	inst.main()




