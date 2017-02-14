import optparse
import os
import tempfile

from readbiasfunctions import *
from readbiasplots import *

from matplotlib.backends import backend_pdf

########
# MAIN #
########

if __name__ == "__main__":
	#################
	# Read Parameters
	usage="usage: %prog [options] vcf_file bam_file reference_sequence_file filtered_vcf_file"
	
	parser = optparse.OptionParser(usage)
	
	# Optional Arguments
	parser.add_option('--tempFolder', action='store', type='string', dest='temp_folder', help='Path to the folder where temporary files are stored [default: %default]', default=tempfile.gettempdir())
	parser.add_option('--mapq', action='store', type='int', dest='mapq', help='Minimal mapping quality of a read to be considered for error count calculation [default: %default]', default=30)
	parser.add_option('--baseq', action='store', type='int', dest='baseq', help='Minimal base quality to be considered for error count calculation [default: %default]', default=13)
	parser.add_option('--qualityScheme', action='store', type='string', dest='quality_scheme', help='Scheme of quality score used in sequencing (illumina or phred) [default: %default]', default="phred")
	parser.add_option('--pValThres', action='store', type='float', dest='p_val_thres', help='P-value threshold of binomial test for read bias [default: %default]', default=0.01)
	parser.add_option('--biasRatioMin', action='store', type='float', dest='bias_ratio_min', help='Minimal bias ratio for a variant to be considered as weakly biased [default: %default]', default=0.53)
	parser.add_option('--biasRatioMax', action='store', type='float', dest='bias_ratio_max', help='Minimal bias ratio for a variant to be considered as strongly biased [default: %default]', default=0.63)
	parser.add_option('--nReadsMin', action='store', type='int', dest='n_reads_min', help='Minimal number of reads observed for a certain variant to be considered for weak bias calculation [default: %default]', default=20)
	parser.add_option('--nMutMin', action='store', type='int', dest='n_mut_min', help='Minimal number of mutations observed for a certain variant to be considered for bias calculation [default: %default]', default=4)
	parser.add_option('--maxOpReadsPcrWeak', action='store', type='int', dest='max_op_reads_pcr_weak', help='Maximal number of reads observed on the opposite strand to flag a variant as being weakly pcr biased [default: %default]', default=0)
	parser.add_option('--maxOpReadsPcrStrong', action='store', type='int', dest='max_op_reads_pcr_strong', help='Maximal number of reads observed on the opposite strand to flag a variant as being strongly pcr biased [default: %default]', default=1)
	parser.add_option('--maxOpReadsSeqWeak', action='store', type='int', dest='max_op_reads_seq_weak', help='Maximal number of reads observed on the opposite strand to flag a variant as being weakly sequencing biased [default: %default]', default=0)
	parser.add_option('--maxOpReadsSeqStrong', action='store', type='int', dest='max_op_reads_seq_strong', help='Maximal number of reads observed on the opposite strand to flag a variant as being strongly sequencing biased [default: %default]', default=1)
	parser.add_option('--maxOpRatioPcr', action='store', type='float', dest='max_op_ratio_pcr', help='Maximal ratio of reads from opposite strand to flag a variant as pcr biased [default: %default]', default=0.1)
	parser.add_option('--maxOpRatioSeq', action='store', type='float', dest='max_op_ratio_seq', help='Maximal ratio of reads from opposite strand to flag a variant as pcr biased [default: %default]', default=0.1)
	parser.add_option('--filterCycles', action='store', type='int', dest='filter_cycles', help='Number of filtering cycles. If number of cycles is 0, then the vcf file is only annotated with ACGTNacgtn<PLUS | MINUS> entries in the INFO field, and bias plots are created before filtering [default: %default]', default=2)
	parser.add_option('-q', '--writeQC', action='store_true', dest='write_qc', help='Write quality control? If true, then a folder is created within the same folder as the filtered vcf file storing bias plots and qc files')
	
	(options,args) = parser.parse_args()	
	
	# Positional (required) arguments
	vcf_filename=os.path.abspath(args[0])
	bam_filename=os.path.abspath(args[1])
	reference_sequence_filename=os.path.abspath(args[2])
	filtered_vcf_filename=os.path.abspath(args[3])


	##################
	# Perform Analysis
	pdf_pages_pcr = None
	pdf_pages_seq = None
	pdf_pages_pcr_simplified = None
	pdf_pages_seq_simplified = None

	qc_folder = None
	plot_file_prefix = None
	pcr_matrix_file_prefix = None
	seq_matrix_file_prefix = None
	mut_matrix_file_prefix = None
	if(options.write_qc):
		# Set QC folder
		split_vcf_path = os.path.split(filtered_vcf_filename)
		# cut file suffix extension from filtered vcf filename
		filtered_vcf_file_prefix = ".".join(split_vcf_path[1].split(".")[:-1])
		qc_folder = (os.sep).join([split_vcf_path[0], filtered_vcf_file_prefix+"_qcSummary"])

		# Create QC folder, if it does not exist
		if(not(os.path.exists(qc_folder))):
			os.mkdir(qc_folder)


		plot_file_prefix = (os.sep).join([qc_folder, "plots", filtered_vcf_file_prefix])
		pcr_matrix_file_prefix = (os.sep).join([qc_folder, "pcr_matrices", filtered_vcf_file_prefix])
		seq_matrix_file_prefix = (os.sep).join([qc_folder, "seq_matrices", filtered_vcf_file_prefix])
		mut_matrix_file_prefix = (os.sep).join([qc_folder, "mut_matrices", filtered_vcf_file_prefix])

		# Create QC subfolders, if they do not exist yet
		if(not(os.path.exists(os.path.dirname(plot_file_prefix)))):
			os.mkdir(os.path.dirname(plot_file_prefix))
		if(not(os.path.exists(os.path.dirname(pcr_matrix_file_prefix)))):
			os.mkdir(os.path.dirname(pcr_matrix_file_prefix))
		if(not(os.path.exists(os.path.dirname(seq_matrix_file_prefix)))):
			os.mkdir(os.path.dirname(seq_matrix_file_prefix))
		if(not(os.path.exists(os.path.dirname(mut_matrix_file_prefix)))):
			os.mkdir(os.path.dirname(mut_matrix_file_prefix))

		# Create PdfPages object for pcr and sequencing bias plots
		pdf_pages_pcr = backend_pdf.PdfPages(plot_file_prefix+"_pcr_bias.pdf")
		pdf_pages_seq = backend_pdf.PdfPages(plot_file_prefix+"_seq_bias.pdf")
		pdf_pages_pcr_simplified = backend_pdf.PdfPages(plot_file_prefix+"_pcr_bias_simplified.pdf")
		pdf_pages_seq_simplified = backend_pdf.PdfPages(plot_file_prefix+"_seq_bias_simplified.pdf")

	# Annotate VCF file with ACGTNacgtn<PLUS | MINUS> fields and get first round of error matrices
	vcf_filename_temp = options.temp_folder+os.sep+os.path.basename(vcf_filename)+".tmp"
	error_matrix_pcr, error_matrix_sequencing, mutation_count_matrix = calculateErrorMatrix(vcf_filename, vcf_filename_temp, reference_sequence_filename, bam_filename, options.mapq, options.baseq, options.quality_scheme)

	bias_matrix_pcr = None
	bias_matrix_seq = None
	plot_str_pcr = "PCR bias before filtering"
	plot_str_seq = "Sequencing bias before filtering"
	matrix_cycle_string = "before_filtering.csv"
	for i in range(options.filter_cycles):
		# Create bias matrices
		bias_matrix_pcr = calculateBiasMatrix(options.p_val_thres, options.bias_ratio_min, options.bias_ratio_max, options.n_reads_min, options.n_mut_min, error_matrix_pcr, mutation_count_matrix)
		bias_matrix_seq = calculateBiasMatrix(options.p_val_thres, options.bias_ratio_min, options.bias_ratio_max, options.n_reads_min, options.n_mut_min, error_matrix_sequencing, mutation_count_matrix)

		if(options.write_qc):
			# Plot results
			plotErrorMatrix(error_matrix_pcr, mutation_count_matrix, plot_str_pcr, pdf_pages_pcr)
			plotErrorMatrix(error_matrix_sequencing, mutation_count_matrix, plot_str_seq, pdf_pages_seq)
			plotErrorMatrix(bias_matrix_pcr, mutation_count_matrix, plot_str_pcr, pdf_pages_pcr_simplified, is_bias=True)
			plotErrorMatrix(bias_matrix_seq, mutation_count_matrix, plot_str_seq, pdf_pages_seq_simplified, is_bias=True)

			# Write matrices
			error_matrix_pcr_filename = pcr_matrix_file_prefix+"_pcr_error_matrix_"+matrix_cycle_string
			error_matrix_seq_filename = seq_matrix_file_prefix+"_seq_error_matrix_"+matrix_cycle_string
			bias_matrix_pcr_filename = pcr_matrix_file_prefix+"_pcr_bias_matrix_"+matrix_cycle_string
			bias_matrix_seq_filename = seq_matrix_file_prefix+"_seq_bias_matrix_"+matrix_cycle_string
			mut_matrix_filename = mut_matrix_file_prefix+"_mut_matrix_"+matrix_cycle_string
			writeMatrix(error_matrix_pcr, error_matrix_pcr_filename)
			writeMatrix(error_matrix_sequencing, error_matrix_seq_filename)
			writeMatrix(bias_matrix_pcr, bias_matrix_pcr_filename, is_bias=True)
			writeMatrix(bias_matrix_seq, bias_matrix_seq_filename, is_bias=True)
			writeMatrix(mutation_count_matrix, mut_matrix_filename, is_bias=True)
	
		# Update plot titles
		plot_str_pcr = "PCR bias after "+str(i+1)+" rounds of filtering"
		plot_str_seq = "Sequencing bias after "+str(i+1)+" rounds of filtering"
		matrix_cycle_string = str(i+1)+"rounds_of_filtering.csv"
	
		# Flag biased variants
		error_matrix_pcr, error_matrix_sequencing, mutation_count_matrix = flagBiasedMutations(vcf_filename_temp, filtered_vcf_filename, reference_sequence_filename, bias_matrix_pcr, bias_matrix_seq, options.max_op_reads_pcr_weak, options.max_op_reads_pcr_strong, options.max_op_reads_seq_weak, options.max_op_reads_seq_strong, options.max_op_ratio_pcr, options.max_op_ratio_seq)

		# Move filtered_vcf_filename to vcf_filename_temp
		os.rename(filtered_vcf_filename, vcf_filename_temp)

	# Create bias matrices
	bias_matrix_pcr = calculateBiasMatrix(options.p_val_thres, options.bias_ratio_min, options.bias_ratio_max, options.n_reads_min, options.n_mut_min, error_matrix_pcr, mutation_count_matrix)
	bias_matrix_seq = calculateBiasMatrix(options.p_val_thres, options.bias_ratio_min, options.bias_ratio_max, options.n_reads_min, options.n_mut_min, error_matrix_sequencing, mutation_count_matrix)

	if(options.write_qc):
		# Plot results
		plotErrorMatrix(error_matrix_pcr, mutation_count_matrix, plot_str_pcr, pdf_pages_pcr)
		plotErrorMatrix(error_matrix_sequencing, mutation_count_matrix, plot_str_seq, pdf_pages_seq)
		plotErrorMatrix(bias_matrix_pcr, mutation_count_matrix, plot_str_pcr, pdf_pages_pcr_simplified, is_bias=True)
		plotErrorMatrix(bias_matrix_seq, mutation_count_matrix, plot_str_seq, pdf_pages_seq_simplified, is_bias=True)

		# Write matrices
		error_matrix_pcr_filename = pcr_matrix_file_prefix+"_pcr_error_matrix_"+matrix_cycle_string
		error_matrix_seq_filename = seq_matrix_file_prefix+"_seq_error_matrix_"+matrix_cycle_string
		bias_matrix_pcr_filename = pcr_matrix_file_prefix+"_pcr_bias_matrix_"+matrix_cycle_string
		bias_matrix_seq_filename = seq_matrix_file_prefix+"_seq_bias_matrix_"+matrix_cycle_string
		mut_matrix_filename = mut_matrix_file_prefix+"_mut_matrix_"+matrix_cycle_string
		writeMatrix(error_matrix_pcr, error_matrix_pcr_filename)
		writeMatrix(error_matrix_sequencing, error_matrix_seq_filename)
		writeMatrix(bias_matrix_pcr, bias_matrix_pcr_filename, is_bias=True)
		writeMatrix(bias_matrix_seq, bias_matrix_seq_filename, is_bias=True)
		writeMatrix(mutation_count_matrix, mut_matrix_filename, is_bias=True)
	
		# Close PDF files containing plots
		pdf_pages_pcr.close()
		pdf_pages_seq.close()
		pdf_pages_pcr_simplified.close()
		pdf_pages_seq_simplified.close()

	# Move temporary vcf file to its final destination
	os.rename(vcf_filename_temp, filtered_vcf_filename)
