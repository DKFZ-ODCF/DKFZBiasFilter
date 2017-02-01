"""
    This module provides access to functions used for the calculation
    of sequence read biases in DNA mutations detected from NGS data.
"""
import sys
import pysam
import numpy as np
import math
from scipy.stats import binom

#################
# Help Routines #
#################

def qualFromASCII(ch):
    """
    Return an integer converted base quality.

    Args
    ----
      ch: string
        1 letter string

    Value
    -----
      int
        Integer value corresponding to ASCII character
    """
    return(ord(ch) - qualScoreOffset)

def transformQualStr(s):
    """
    Return an integer converted base quality list derived from a string.

    Args
    ----
      s: string

    Value
    -----
      list
        List of integer values realting to integer converted ASCII characters.
    """
    return map(qualFromASCII,s)

def getIndexACGTNacgtn(is_reverse, is_read1, base):
    """
    Return index of a base in ACGTNacgtn list based on read strand information.

    Args
    ----
      is_reverse: bool
        Is read reverse?
      is_read1: bool
        Is read the first in sequencing?
      base: string
        1 letter string (A | C | C | T | N | a | c | g | t | n)

    Value
    -----
      int
        index of base in ACGTNacgtn list
    """
    if (is_reverse):
        if(is_read1):
            if(base == "a"):
                return ["minus", 5]
            elif(base == "c"):
                return ["minus", 6]
            elif(base == "g"):
                return ["minus", 7]
            elif(base == "t"):
                return ["minus", 8]
            elif(base == "n"):
                return ["minus", 9]
        else:
            if(base == "a"):
                return ["plus", 5]
            elif(base == "c"):
                return ["plus", 6]
            elif(base == "g"):
                return ["plus", 7]
            elif(base == "t"):
                return ["plus", 8]
            elif(base == "n"):
                return ["plus", 9]
    else:
        if(is_read1):
            if(base == "a"):
                return ["plus", 0]
            elif(base == "c"):
                return ["plus", 1]
            elif(base == "g"):
                return ["plus", 2]
            elif(base == "t"):
                return ["plus", 3]
            elif(base == "n"):
                return ["plus", 4]
        else:
            if(base == "a"):
                return ["minus", 0]
            elif(base == "c"):
                return ["minus", 1]
            elif(base == "g"):
                return ["minus", 2]
            elif(base == "t"):
                return ["minus", 3]
            elif(base == "n"):
                return ["minus", 4]

def complement(base):
    """
        Return complement of a base.

        Args
        ----
          base: string
            1 letter string

        Value
        -----
          string
            1 letter string
    """
    if(base == "A"):
        return "T"
    elif(base == "C"):
        return "G"
    elif(base == "G"):
        return "C"
    elif(base == "T"):
        return "A"
    elif(base == "a"):
        return "t"
    elif(base == "c"):
        return "g"
    elif(base == "g"):
        return "c"
    elif(base == "t"):
        return "a"
    elif(base == "n"):
        return "n"
    elif(base == "N"):
        return "N"

def splitMetaInfoString(meta_info_string):
    meta_info_string_split_raw = meta_info_string.split(",")
    meta_info_string_split = []
    open_field = False
    for element in meta_info_string_split_raw:
        if('"' in element and not open_field):
            open_field = True
            meta_info_string_split += [element]
            if(element[-1] == '"'):
                open_field = False
        elif(open_field):
            meta_info_string_split[-1] = ",".join([meta_info_string_split[-1], element])
            if(element[-1] == '"'):
                open_field = False
        else:
            meta_info_string_split += [element]
    return meta_info_string_split

def createMetaInfoDict(meta_info_string):
    """
    Return a dictionary based on a meta information line in a vcf file.
    Args
    ----
      meta_info_string: string
    Value
    -----
      dictionary
        keys: string
          keys of meta information string
        values: string
          values of meta information string
    """
    meta_info_dict = {}
    if(meta_info_string[0] == "<" and meta_info_string[-1] == ">"):
        for tupel in splitMetaInfoString(meta_info_string[1:-1]):
            split_tupel = [tupel.split("=")[0], "=".join(tupel.split("=")[1:])]
            meta_info_dict[split_tupel[0]] = split_tupel[1]
    return meta_info_dict

def calculateACGTNacgtnFields(bamFile, chromosome, position, mapq, baseq, qual_score_offset):
    """
    Return ACGTNacgtn<PLUS | MINUS> fields, given a bam file and a genomic position.

    Args
    ----
      bamFile: pysam.Samfile instance
      chromosome: string
      position: int
      mapq: float
        Minimal mapping quality of a read to be considered
      baseq: float
        Minimal base quality to be considered
      qual_score_offset: int
        Quality score offset used to convert as ASCII character into an integer

    Value
    -----
      string
    """
    global qualScoreOffset
    qualScoreOffset = qual_score_offset

    ACGTNacgtn1 = [0]*10
    ACGTNacgtn2 = [0]*10

    readNameHash = {}

    for pileupcolumn in bamFile.pileup(chromosome, (position-1), position):
        if pileupcolumn.pos == (position-1):
            for pileupread in pileupcolumn.pileups:
                # TODO: Check if read is dupmarked!
                if pileupread.alignment.mapq >= mapq:
                    # Define Positions of base and base quality within read
					pos_base = pileupread.qpos
                    pos_qual = pos_base
                    cigar_tuples = pileupread.alignment.cigar
                    if(cigar_tuples[0][0] == 4):
                        pos_qual_offset = cigar_tuples[0][1]
                        pos_qual -= pos_qual_offset

                    # Only consider bases which are not! soft-clipped
                    if(pos_qual >= len(pileupread.alignment.qqual) or pos_qual < 0):
                        # We are looking at soft-clipped bases at the end of the read and skip!
                        continue

                    if transformQualStr(pileupread.alignment.qqual[pos_qual])[0] >= baseq:
                        try:
                            readNameHash[pileupread.alignment.qname] += 1
                        except KeyError:
                            readNameHash[pileupread.alignment.qname] = 1
                            # If read name was not seen include count to ACGTNacgtn list
                            is_reverse = pileupread.alignment.is_reverse
                            is_read1 = pileupread.alignment.is_read1
                            base = pileupread.alignment.seq[pos_base].lower()
                            ACGTNacgtn_index = getIndexACGTNacgtn(is_reverse, is_read1, base)
                            if(ACGTNacgtn_index[0] == "plus"):
                                ACGTNacgtn1[ACGTNacgtn_index[1]] += 1
                            else:
                                ACGTNacgtn2[ACGTNacgtn_index[1]] += 1

    ACGTNacgtn1_string = "ACGTNacgtnPLUS="+",".join([str(i) for i in ACGTNacgtn1])
    ACGTNacgtn2_string = "ACGTNacgtnMINUS="+",".join([str(i) for i in ACGTNacgtn2])

    return ACGTNacgtn1_string+";"+ACGTNacgtn2_string

def writeMatrix(matrix, output_filename, is_bias=False):
    """
    Write bias, or error matrix to a file.

    Args
    ----
      matrix: dictionary
        Dictionary containing for each possible mutation, and triplet context a number or a list of numbers
      output_filename: string
        filepath to the output file
      is_bias: bool
        Is matrix[mut][base_before][base_after] an integer (is_bias==True) or a list of integer (is_bias==False)
    """
    output_file = open(output_filename, "w")
	possible_mutations = ["CA", "CG", "CT", "TA", "TC", "TG"]
    possible_bases = ["A", "C", "G", "T"]

    for mut in possible_mutations:
        output_file.write(">"+mut[0]+"->"+mut[1]+"\n")
        output_file.write("\t".join([""]+possible_bases)+"\n")
        for base_before in possible_bases:
            entries=[]
            for base_after in possible_bases:
                if(is_bias):
                    entries += [str(matrix[mut][base_before][base_after])]
                else:
                    entries += [str(matrix[mut][base_before][base_after][0])+";"+str(matrix[mut][base_before][base_after][1])]
            output_file.write("\t".join([base_before]+entries)+"\n")

    output_file.close()

#####################################
# Main Routines for bias annotation #
#####################################

def calculateErrorMatrix(vcfFilename, vcf_filename_temp, referenceFilename, bamFilename, mapq, baseq, qualityScore):
    """
        Return read count matrices for plus and minus stranded PCR template-, and sequencing reads, and a mutation count matrix. 
        Write ACGTNacgtn<PLUS | MINUS> entries to a newly created vcf file.

        Args
        ----
          vcfFilename: string
            Filepath to the input vcf file
          vcf_filename: string
            Filepath to the output vcf file
          referenceFilename: string
            Filepath to the reference sequence (fasta format). A fasta index must exist in the same directory
          bamFilename: string
            Filepath to a bam file. A bam index file must exist in the same directory
          mapq: float
            Minimal mapping quality of a read to be considered
          baseq: float
            Minimal base quality to be considered
          qualityScore: string
            Quality scoring scheme for base qualities used in the bam file (values: "illumina" | "phred")

        Value
        -----
          dictionary
            Read counts PCR strands
          dictionary
            Read counts sequencing strands
          dictionary
            Mutation counts
    """
    if qualityScore == 'illumina': qualScoreOffset = 64
    elif qualityScore == 'phred': qualScoreOffset = 33

    # Open files
	vcfFile = open(vcfFilename, "r")
    reference = pysam.Fastafile(referenceFilename)
    bamFile = pysam.Samfile(bamFilename)
    vcf_file_temp=open(vcf_filename_temp, "w")

	possible_mutations = ["CA", "CG", "CT", "TA", "TC", "TG"]
    possible_bases_clean = ["A", "C", "G", "T"]

    # Initialize Error and Mutation Count Matrix
    error_matrix_pcr = {}
    error_matrix_sequencing = {}
    mutation_count_matrix = {}
    for mutation in possible_mutations:
        error_matrix_pcr[mutation] = {}
        error_matrix_sequencing[mutation] = {}
        mutation_count_matrix[mutation] = {}
        possible_bases = ["A", "C", "G", "T"]
        for base_before in possible_bases:
            error_matrix_pcr[mutation][base_before] = {}
            error_matrix_sequencing[mutation][base_before] = {}
            mutation_count_matrix[mutation][base_before] = {}
            for base_after in possible_bases:
                error_matrix_pcr[mutation][base_before][base_after] = [1, 1] # Initialize with pseudo counts
                error_matrix_sequencing[mutation][base_before][base_after] = [1, 1] # Initialize with pseudo counts
                mutation_count_matrix[mutation][base_before][base_after] = 0

    # Initialize header list for later getting indices of column names
    has_header = False
    header_written=False
    header = None
    header_lines = []
    for line in vcfFile:
        if(line[:2] == "##"):
            # Save meta information lines for later writing
            # Check if ACGTNacgtnPLUS, and ACGTNacgtnMINUS are already contained in the metainformation (FORMAT), and remove them if so.
            split_meta_info = line.rstrip().split("=")
            if(split_meta_info[0] == "##INFO"):
                meta_info_dict = createMetaInfoDict(line.rstrip()[7:])
                if(not(meta_info_dict["ID"] == "ACGTNacgtnPLUS" or meta_info_dict["ID"] == "ACGTNacgtnMINUS")):
                    header_lines += [line]
            else:
                header_lines += [line]
        elif(line[:1] == "#"):
            has_header = True
            header_lines += [line]
            header=line.rstrip().split("\t")
        else:
            # Exit, if vcf file does not contain a proper header line
            if(not(has_header)):
                exit(vcfFilename+" does not include a vcf conform header (\"#CHROM\tPOS\t...)\"!")

            # Write header
            if(not(header_written)):
                ACGTNacgtnPLUS_string="##INFO=<ID=ACGTNacgtnPLUS,Number=10,Type=Integer,Description=\"The first five numbers correspond to the number of bases on forward reads found to be A, C, G, T, or N, while the last five numbers correspond to bases on reverse reads found to be a, c, g, t, or n on plus stranded PCR templates (only reads with a mapping quality greater or equal to "+str(mapq)+", and bases with a base quality greater or equal to "+str(baseq)+" were considered).\">\n"
                ACGTNacgtnMINUS_string="##INFO=<ID=ACGTNacgtnMINUS,Number=10,Type=Integer,Description=\"The first five numbers correspond to the number of bases on forward reads found to be A, C, G, T, or N, while the last five numbers correspond to bases on reverse reads found to be a, c, g, t, or n on minus stranded PCR templates (only reads with a mapping quality greater or equal to "+str(mapq)+", and bases with a base quality greater or equal to "+str(baseq)+" were considered).\">\n"
                header_line=header_lines[-1]
                meta_information_lines=header_lines[:-1]
                meta_information_lines += [ACGTNacgtnPLUS_string, ACGTNacgtnMINUS_string, header_line]
                for l in meta_information_lines:
                    vcf_file_temp.write(l)
                header_written=True
                

            split_line = line.rstrip().split("\t")
            # Skip entry if it was already flagged as being biased
            flagged=False
            for filter_flag in split_line[header.index("FILTER")].split(";"):
                if(filter_flag == "bPcr" or filter_flag == "bSeq"):
                    flagged=True
            if(flagged):
                vcf_file_temp.write(line)
                continue

            chrom = split_line[header.index("#CHROM")]
            pos = int(split_line[header.index("POS")])
            context = reference.fetch(chrom, pos-2, pos+1)

            ref = split_line[header.index("REF")].split(",")
            alt = split_line[header.index("ALT")].split(",")

            # Define current mutation as string "<REF><ALT>". If two alternative bases exist, mutation is defined
            # as "<ALT1><ALT2>"
            current_mutation = ""
            if(len(alt) == 1):
                current_mutation = ref[0]+alt[0]
            else:
                current_mutation = alt[0]+alt[1]
            
            base_before = None
            base_after = None
            try:
                base_before = context[0].upper()
                base_after = context[2].upper()
            except IndexError:
                print "No reference sequence for SNV at: "+chrom+" "+str(pos)
                sys.exit()

            acgtn_fields = calculateACGTNacgtnFields(bamFile, chrom, pos, mapq, baseq, qualScoreOffset)

            # Append atcgn_fields entry to INFO fields and write split_line to vcf_file_temp
            # Remove ACGTNacgtn entries from INFO field, if they already exist
            info_field = split_line[header.index("INFO")].split(";")
            cleaned_info_field = []
            for element in info_field:
                split_element = element.split("=")
                if(not(split_element[0] == "ACGTNacgtnPLUS" or split_element[0] == "ACGTNacgtnMINUS")):
                    cleaned_info_field += [element]
            
            split_line[header.index("INFO")] = ";".join(cleaned_info_field+[acgtn_fields])
            vcf_file_temp.write("\t".join(split_line)+"\n")

            info_list = [i.split("=") for i in acgtn_fields.split(";")]

            # Get strand specific counts
            ACGTNacgtnPLUS = []
            ACGTNacgtnMINUS = []

            for element in info_list:
                if(element[0] == "ACGTNacgtnPLUS"):
                    ACGTNacgtnPLUS = [int(i) for i in element[1].split(",")]
                elif(element[0] == "ACGTNacgtnMINUS"):
                    ACGTNacgtnMINUS = [int(i) for i in element[1].split(",")]

            # Count number of alternative bases
            possible_bases = ["A", "C", "G", "T", "N", "a", "c", "g", "t", "n"]
            read1_nr = ACGTNacgtnPLUS[possible_bases.index(current_mutation[1])]
            read1_r = ACGTNacgtnMINUS[possible_bases.index(current_mutation[1].lower())]
            read2_nr = ACGTNacgtnMINUS[possible_bases.index(current_mutation[1])]
            read2_r = ACGTNacgtnPLUS[possible_bases.index(current_mutation[1].lower())]

            # Check if current_mutation is in set of possible mutations. If not, reverse complement current_mutation
            # as well as base_before and base_after
            reverse_mutation=False
            try:
                mutation_index = possible_mutations.index(current_mutation)
            except ValueError:
                current_mutation = complement(current_mutation[0])+complement(current_mutation[1])
                base_before_reverse_complement = complement(base_after)
                base_after_reverse_complement = complement(base_before)

                base_before = base_before_reverse_complement
                base_after = base_after_reverse_complement
                reverse_mutation=True

            # Extend Error Matrix for PCR errors
            PCR_plus = 0
            PCR_minus = 0
            if(not(reverse_mutation)):
                PCR_plus = read1_nr + read2_r
                PCR_minus = read2_nr + read1_r

            else:
                PCR_plus = read1_r + read2_nr
                PCR_minus = read1_nr + read2_r

            if(current_mutation in possible_mutations and base_before in possible_bases_clean and base_after in possible_bases_clean):
                error_matrix_pcr[current_mutation][base_before][base_after][0] += PCR_plus
                error_matrix_pcr[current_mutation][base_before][base_after][1] += PCR_minus

            # Extend Error Matrix for sequencing errors
            SEQ_plus = 0
            SEQ_minus = 0
            if(not(reverse_mutation)):
                SEQ_plus = read1_nr + read2_nr
                SEQ_minus = read1_r + read2_r

            else:
                SEQ_plus = read1_r + read2_r
                SEQ_minus = read1_nr + read2_nr

            if(current_mutation in possible_mutations and base_before in possible_bases_clean and base_after in possible_bases_clean):
                error_matrix_sequencing[current_mutation][base_before][base_after][0] += SEQ_plus
                error_matrix_sequencing[current_mutation][base_before][base_after][1] += SEQ_minus
                mutation_count_matrix[current_mutation][base_before][base_after] += 1

    # Close files
    vcf_file_temp.close()
    vcfFile.close()

    return error_matrix_pcr, error_matrix_sequencing, mutation_count_matrix

def calculateBiasMatrix(p_val_threshold, bias_ratio_min, bias_ratio_max, n_reads_min, n_muts_min, error_matrix, mutation_count_matrix):
    """
        Return bias matrix for all possible mutations.

        Args
        ----
          p_val_threshold: float
            Significance threshold of binomial test for bias testing
          bias_ratio_min: float
            Minimal ratio of reads from strand with major read count to consider a mutation for weak bias
          bias_ratio_max: float
            Minimal ratio of reads from strand with major read count to consider a mutation for strong bias
          n_reads_min: int
            Minimal number of reads found for a mutation to consider it for bias calculation
          n_muts_min: int
            Minimal number of mutations found for a certain mutation to consider it for bias calculation
          error_matrix: dictionary
            Read count matrix
          mutation_count_matrix: dictionary
            Mutation count matrix

        Value
        -----
          dictionary
            Dictionary, containing the bias information for all possible mutations in all possible triplet contexts
    """
	possible_mutations = ["CA", "CG", "CT", "TA", "TC", "TG"]
    possible_bases = ["A", "C", "G", "T"]

    bias_matrix = {}
    for mut in possible_mutations:
        bias_matrix[mut] = {}
        for base_before in possible_bases:
            bias_matrix[mut][base_before] = {}
            for base_after in possible_bases:
                bias_matrix[mut][base_before][base_after] = 0

                # Variables needed for bias calculation
                n_reads_plus = error_matrix[mut][base_before][base_after][0]
                n_reads_minus = error_matrix[mut][base_before][base_after][1]
                minor_read_count = min([n_reads_plus, n_reads_minus])
                n_reads = n_reads_plus + n_reads_minus

                # Catch zero division errors
                frac_plus = 0.5
                frac_minus = 0.5
                if(n_reads > 0):
                    frac_plus = float(n_reads_plus)/float(n_reads)
                    frac_minus = float(n_reads_minus)/float(n_reads)

                n_muts = mutation_count_matrix[mut][base_before][base_after]
                bias = None

                # If there are more plus stranded reads than minus stranded reads
                if(n_reads_plus >= n_reads_minus):
                    if(binom.cdf(minor_read_count, n_reads, 0.5) <= p_val_threshold and frac_plus >= bias_ratio_min and n_reads >= n_reads_min and n_muts >= n_muts_min):
                        bias = 1

                        # Strong bias if the fraction of plus stranded reads exceeds a given upper threshold
                        if(frac_plus >= bias_ratio_max):
                            bias = 2

                # If there are more minus stranded reads than plus stranded reads 
                else:
                    if(binom.cdf(minor_read_count, n_reads, 0.5) <= p_val_threshold and frac_minus >= bias_ratio_min and n_reads >= n_reads_min and n_muts >= n_muts_min):
                        bias = -1

                        # Strong bias if the fraction of minus stranded reads exceeds a given upper threshold
                        if(frac_minus >= bias_ratio_max):
                            bias = -2

                # If there is a bias write it to bias_matrix
                if(bias is not None):
                    bias_matrix[mut][base_before][base_after] = bias

    return bias_matrix
                    
def flagBiasedMutations(vcf_filename, vcf_filename_flagged, reference_filename, bias_matrix_pcr, bias_matrix_seq, max_num_opposite_reads_pcr_weak, max_num_opposite_reads_pcr_strong, max_num_opposite_reads_seq_weak, max_num_opposite_reads_seq_strong, max_opposite_ratio_pcr, max_opposite_ratio_seq):
    """
        Flag vcf file for read biases and return read count matrices for plus and minus stranded PCR template-, and sequencing reads, after filtering mutations showing a read bias. Furthermore, return a mutation count matrix after filtering for read biases.

        Args
        ----
          vcf_filename: string
            Filepath to vcf file being flagged for read biases
          vcf_filename_flagged: string
            Filepath to resulting (flagged) vcf file
          reference_filename: string
            Filepath to reference sequence (fasta format). A fasta file index must exist in the same directory
          bias_matrix_pcr: dictionary
            Bias matrix for pcr biases
          bias_matrix_seq: dictionary
            Bias matrix for sequencing biases
          max_num_opposite_reads_pcr_weak: int
            Maximal number of reads from opposite strand allowed to flag a mutation as weakly pcr biased
          max_num_opposite_reads_pcr_strong: int
            Maximal number of reads from opposite strand allowed to flag a mutation as strongly pcr biased
          max_num_opposite_reads_seq_weak: int
            Maximal number of reads from opposite strand allowed to flag a mutation as weakly sequencing biased
          max_num_opposite_reads_seq_strong: int
            Maximal number of reads from opposite strand allowed to flag a mutation as strongly sequencing biased
          max_opposite_ratio_pcr: float
            Maximal ratio of reads from opposite strand allowed to flag a mutation as pcr biased
          max_opposite_ratio_seq: float
            Maximal ratio of reads from opposite strand allowed to flag a mutation as sequencing biased

        Value
        -----
          dictionary
            Read counts pcr strands
          dictionary
            Read counts sequencing strands
          dictionary
            Mutation counts
    """
    # General data
	possible_mutations = ["CA", "CG", "CT", "TA", "TC", "TG"]
    possible_bases = ["A", "C", "G", "T", "N", "a", "c", "g", "t", "n"]
    possible_bases_clean = ["A", "C", "G", "T"]

    # Initialize Error Matrices for PCR and SEQ Error, as well as mutation counts matrix
    error_matrix_pcr = {}
    error_matrix_seq = {}
    mutation_count_matrix = {}
    for mut in possible_mutations:
        error_matrix_pcr[mut] = {}
        error_matrix_seq[mut] = {}
        mutation_count_matrix[mut] = {}
        for base_before in possible_bases_clean:
            error_matrix_pcr[mut][base_before] = {}
            error_matrix_seq[mut][base_before] = {}
            mutation_count_matrix[mut][base_before] = {}
            for base_after in possible_bases_clean:
                error_matrix_pcr[mut][base_before][base_after] = [1, 1]
                error_matrix_seq[mut][base_before][base_after] = [1, 1]
                mutation_count_matrix[mut][base_before][base_after] = 0

    # open files
	vcf_file = open(vcf_filename, "r")
    vcf_file_flagged = open(vcf_filename_flagged, "w")
    reference_fasta = pysam.Fastafile(reference_filename)

    # Flag biased Variants
    has_header = False
    header_written=False
    header = None
    header_lines = []
    for line in vcf_file:
        if(line[:2] == "##"):
            # Save meta information lines for later writing
            # Check if bPcr, and bSeq are already contained in the metainformation (FILTER), and remove them if so.
            split_meta_info = line.rstrip().split("=")
            if(split_meta_info[0] == "##FILTER"):
                meta_info_dict = createMetaInfoDict(line.rstrip()[9:])
                if(not(meta_info_dict["ID"] == "bPcr" or meta_info_dict["ID"] == "bSeq")):
                    header_lines += [line]
            else:
                header_lines += [line]
        elif(line[:1] == "#"):
            has_header = True
            header_lines += [line]
            header=line.rstrip().split("\t")
        else:
            # Exit, if vcf file does not contain a proper header line
            if(not(has_header)):
                exit(vcfFilename+" does not include a vcf conform header (\"#CHROM\tPOS\t...)\"!")

            # Write header
            if(not(header_written)):
                pcr_bias_filter_string="##FILTER=<ID=bPcr,Description=\"Variant allele shows a bias towards one PCR template strand.\">\n"
                seq_bias_filter_string="##FILTER=<ID=bSeq,Description=\"Variant allele shows a bias towards one sequencing strand.\">\n"

                header_line=header_lines[-1]
                meta_information_lines=header_lines[:-1]
                meta_information_lines += [pcr_bias_filter_string, seq_bias_filter_string, header_line]
                for l in meta_information_lines:
                    vcf_file_flagged.write(l)
                header_written=True

            split_line = line.rstrip().split("\t")

            # Skip entry if it was already flagged as being biased
            flagged=False
            for filter_flag in split_line[header.index("FILTER")].split(";"):
                if(filter_flag == "bPcr" or filter_flag == "bSeq"):
                    flagged=True
            if(flagged):
                vcf_file_flagged.write(line)
                continue

            # Extract information from vcf entry
            chrom=split_line[header.index("#CHROM")]
            pos=int(split_line[header.index("POS")])
            ref=split_line[header.index("REF")].split(",")
            alt=split_line[header.index("ALT")].split(",")
            context = reference_fasta.fetch(chrom, pos-2, pos+1)

            # Define current mutation as string "<REF><ALT>". If two alternative bases exist, mutation is defined
            # as "<ALT1><ALT2>"
            current_mutation = ""
            if(len(alt) == 1):
                current_mutation = ref[0]+alt[0]
            else:
                current_mutation = alt[0]+alt[1]
            base_before = context[0].upper()
            base_after = context[2].upper()

            # Get read count information from ACGTNacgtn fields
            split_info_entry = split_line[header.index("INFO")].split(";")
            ACGTNactgnPLUS=[]
            ACGTNacgtnMINUS=[]
            for element in split_info_entry:
                split_element=element.split("=")
                if(split_element[0] == "ACGTNacgtnPLUS"):
                    ACGTNacgtnPLUS = [int(i) for i in split_element[1].split(",")]
                elif(split_element[0] == "ACGTNacgtnMINUS"):
                    ACGTNacgtnMINUS = [int(i) for i in split_element[1].split(",")]
            read1_f = ACGTNacgtnPLUS[possible_bases.index(current_mutation[1])]
            read1_r = ACGTNacgtnMINUS[possible_bases.index(current_mutation[1].lower())]
            read2_f = ACGTNacgtnMINUS[possible_bases.index(current_mutation[1])]
            read2_r = ACGTNacgtnPLUS[possible_bases.index(current_mutation[1].lower())]

            n_reads = read1_f + read1_r + read2_f + read2_r

            # Check if there is a bias for the current mutation
            # Check if current_mutation is in set of possible mutations. If not, reverse complement current_mutation
            # as well as base_before and base_after
            reverse_mutation=False
            try:
                mutation_index = possible_mutations.index(current_mutation)
            except ValueError:
                current_mutation = complement(current_mutation[0])+complement(current_mutation[1])
                base_before_reverse_complement = complement(base_after)
                base_after_reverse_complement = complement(base_before)

                base_before = base_before_reverse_complement
                base_after = base_after_reverse_complement
                reverse_mutation=True

            # Extend Error Matrix for PCR errors
            n_plus_reads_pcr = 0
            n_minus_reads_pcr = 0
            if(not(reverse_mutation)):
                n_plus_reads_pcr = read1_f + read2_r
                n_minus_reads_pcr = read2_f + read1_r

            else:
                n_plus_reads_pcr = read1_r + read2_f
                n_minus_reads_pcr = read1_f + read2_r

            # Extend Error Matrix for sequencing errors
            n_plus_reads_seq = 0
            n_minus_reads_seq = 0
            if(not(reverse_mutation)):
                n_plus_reads_seq = read1_f + read2_f
                n_minus_reads_seq = read1_r + read2_r

            else:
                n_plus_reads_seq = read1_r + read2_r
                n_minus_reads_seq = read1_f + read2_f       

            frac_plus_reads_pcr = 0.5
            frac_minus_reads_pcr = 0.5
            frac_plus_reads_seq = 0.5
            frac_minus_reads_seq = 0.5

            if(n_reads > 0):
                frac_plus_reads_pcr = float(n_plus_reads_pcr)/float(n_reads)
                frac_minus_reads_pcr = float(n_minus_reads_pcr)/float(n_reads)
                frac_plus_reads_seq = float(n_plus_reads_seq)/float(n_reads)
                frac_minus_reads_seq = float(n_minus_reads_seq)/float(n_reads)

            pcr_bias=False
            seq_bias=False

            # Check if base context contains bases ACGT, and current mmutation is in set of possible mutations
            is_valid = base_before in possible_bases_clean and base_after in possible_bases_clean and current_mutation in possible_mutations
            if(is_valid):
                # Check for PCR bias
                bias_pcr = bias_matrix_pcr[current_mutation][base_before][base_after]
                if(bias_pcr > 0):
                    if(bias_pcr == 1 and n_minus_reads_pcr <= max_num_opposite_reads_pcr_weak and frac_minus_reads_pcr <= max_opposite_ratio_pcr):
                        pcr_bias = True
                    elif(bias_pcr == 2 and n_minus_reads_pcr <= max_num_opposite_reads_pcr_strong and frac_minus_reads_pcr <= max_opposite_ratio_pcr):
                        pcr_bias = True
                if(bias_pcr < 0):
                    if(bias_pcr == -1 and n_plus_reads_pcr <= max_num_opposite_reads_pcr_weak and frac_plus_reads_pcr <= max_opposite_ratio_pcr):
                        pcr_bias = True
                    elif(bias_pcr == -2 and n_plus_reads_pcr <= max_num_opposite_reads_pcr_strong and frac_plus_reads_pcr <= max_opposite_ratio_pcr):
                        pcr_bias = True
                # Check for Seq bias
                bias_seq = bias_matrix_seq[current_mutation][base_before][base_after]
                if(bias_seq > 0):
                    if(bias_seq == 1 and n_minus_reads_seq <= max_num_opposite_reads_seq_weak and frac_minus_reads_seq <= max_opposite_ratio_seq):
                        seq_bias = True
                    elif(bias_seq == 2 and n_minus_reads_seq <= max_num_opposite_reads_seq_strong and frac_minus_reads_seq <= max_opposite_ratio_seq):
                        seq_bias = True
                if(bias_seq < 0):
                    if(bias_seq == -1 and n_plus_reads_seq <= max_num_opposite_reads_seq_weak and frac_plus_reads_seq <= max_opposite_ratio_seq):
                        seq_bias = True
                    elif(bias_seq == -2 and n_plus_reads_seq <= max_num_opposite_reads_seq_strong and frac_plus_reads_seq <= max_opposite_ratio_seq):
                        seq_bias = True

            # Write flagged VCF
            filter_flags=[]
            current_filter_flag = split_line[header.index("FILTER")]
            if(not(current_filter_flag == "PASS" or current_filter_flag == ".")):
                filter_flags = current_filter_flag.split(";")
            if(pcr_bias):
                filter_flags += ["bPcr"]
            if(seq_bias):
                filter_flags += ["bSeq"]
            if(not(pcr_bias or seq_bias)):
                filter_flags = [current_filter_flag]
            split_line[header.index("FILTER")] = ";".join(filter_flags)
            vcf_file_flagged.write("\t".join(split_line)+"\n")

            # Update matrices for pcr and seq errors, as well as mutation counts
            if(not(pcr_bias or seq_bias) and is_valid):
                error_matrix_pcr[current_mutation][base_before][base_after][0] += n_plus_reads_pcr
                error_matrix_pcr[current_mutation][base_before][base_after][1] += n_minus_reads_pcr
                error_matrix_seq[current_mutation][base_before][base_after][0] += n_plus_reads_seq
                error_matrix_seq[current_mutation][base_before][base_after][1] += n_minus_reads_seq
                mutation_count_matrix[current_mutation][base_before][base_after] += 1

    # Close files
    vcf_file.close()
    vcf_file_flagged.close()

    return error_matrix_pcr, error_matrix_seq, mutation_count_matrix
