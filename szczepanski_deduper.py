#!/usr/bin/env python

# Deduper script
# Nikki Szczepanski

# Import necessary functions
import argparse
import re
import os

# Argparse function
def get_args():
    """Define the arguments needed for the script. 
    If `--paired` is specified, or if `--umi` in not specified, print error message and quit."""
    parser = argparse.ArgumentParser(description='Remove PCR duplicates in a SAM file, keeping the first occurence of read.\
        Outputs include the deduped SAM file named [input-filename]_deduped.sam, a file with all duplicate reads named \
        [input-filename]_dupes.sam, and a file with all unknown UMIs named [input-filename]_uknown_umi.sam. \
        Currently only works with single-end reads and specified UMIs.')
    parser.add_argument("-f", "--file", type=str, required=True, help='Absolute path for SAM file to be de-duplicated')
    parser.add_argument("-p", "--paired", required=False, action='store_true', help='Designates file is paired-end (unset if single-end)')
    parser.add_argument("-u", "--umi", type=str, required=False, help='Absolute path to file containing the list of UMIs (unset if randomers instead of UMIs)')
    args = parser.parse_args()
    # If `--paired` is specified, it will be set as `True` (action='store_true'), causing the program to quit.
    if args.paired is True:
        print("\n-----ERROR-----")
        print("\nSorry, the '--paired' feature is not available yet. Only single-end reads can be de-duplicated")
        print("\nProgram quitting.\n")
        print("---------------")
        exit()
    # If `--umi` is not specified, by default it is set to `None`, causing the program to quit.
    if args.umi is None:
        print("\n-----ERROR-----")
        print("\nSorry, the randomer-finder feature is not available yet. UMIs must be specified.")
        print("\nProgram quitting.\n")
        print("---------------")
        exit()
    return parser.parse_args()

parse_args = get_args()
# Rename variables from argparse arguments
sam_file = parse_args.file
paired = parse_args.paired
umi_file = parse_args.umi

# Define the files to be created (with the same path as the input SAM file)
deduped_file = sam_file.split(".sam")[0] + "_deduped.sam"
dupes_file = sam_file.split(".sam")[0] + "_dupes.sam"
unknown_umi_file = sam_file.split(".sam")[0] + "_unknown_umi.sam"

# Samtools sort
#   Sort the file by chromosome and position, and separate the SAM file into two files by forward and reverse strand.
#   The two files will be called 'fwd.sam' and 'rev.sam' and will be deleted at the end of the script.
os.system("samtools sort %s -o sorted.sam" % sam_file)
os.system("samtools view -F 16 sorted.sam > fwd.sam")
os.system("samtools view -f 16 sorted.sam > rev.sam")
split_files=["fwd.sam","rev.sam"]

# UMI dictionary function
def create_umi_dict():
    """Use the list of UMIs to create a dictionary where the key is the UMI and the value is the number of UMI appearances."""
    with open(umi_file, "r") as fh:
        umi_dict = {}
        for line in fh:
            line = line.strip("\n")
            umi_dict[line]=0
    return umi_dict

# UMI function
def find_umi(line):
    """Store the UMI listed at the end of the QNAME as the variable 'umi'."""
    line_list = line.split("\t")
    umi = line_list[0].split(":")[-1]
    return umi

# Chromosome function
    """Store the chromosome listed in Column 3 as the variable 'chrom'."""
def find_chrom(line):
    chrom = line.split("\t")[2]
    return chrom

# Position functions
#-------------------------------
# Forward strand position
def find_pos_fwd(line):
    """Use the CIGAR string (Column 6), position (Column 4), and strand (see above function)
    to determine the starting position, which is the leftmost position for forward strand reads,
    stored as the variable 'start_pos'."""
    pos = int(line.split("\t")[3])
    cigar = line.split("\t")[5]
    first_letter = re.search('[0-9]+([A-Z])', cigar).group(1)
    # The starting position must be adjusted for softclipping if the first letter of the CIGAR string is S.
    if first_letter == "S":
        adjust = int(cigar.split("S")[0])
        start_pos = pos - adjust
    # Otherwise, the position in column 4 automatically gives the leftmost starting position of the read.
    else:
        start_pos = pos
    return start_pos

# Reverse strand position
def find_pos_rev(line):
    """Use the CIGAR string (Column 6), position (Column 4), and strand (see above function)
    to determine the starting position, which is the rightmost position for reverse strand reads,
    stored as the variable 'start_pos'."""
    # Reads that align to the reverse strand of the reference sequence have their 5' end at the rightmost position.
    pos = int(line.split("\t")[3])
    cigar = line.split("\t")[5]
    first_letter = re.search('[0-9]+([A-Z])', cigar).group(1)
    if first_letter == "S":
        pairs_to_add = cigar.split("S")[1]
    else:
        pairs_to_add = cigar
    # Only include the numbers, not the letters, of the cigar string
    nums_to_add = re.findall('[0-9]+',pairs_to_add)
    # Convert numbers (currently characters) to integers
    for i in range(0, len(nums_to_add)):
        nums_to_add[i]=int(nums_to_add[i])
    # Add up all the integers to find the rightmost starting position
    adjust = sum(nums_to_add)
    start_pos = pos + adjust
    return start_pos

# Main function
def main_func():
    """Parse through the initial SAM file and the two split files ('fwd.sam' and 'rev.sam')
    and use the previously defined functions to store 3 key identifiers as a tuple:
    the chromosome, the position, and the UMI. Each line will result in a tuple that is added to a set of unique read identifiers."""
    # Copy the header section (where each line starts with '@') from the original SAM file to the new SAM file.
    with open(deduped_file,'w') as deduped, open(dupes_file,'w') as dupes, open(unknown_umi_file,'w') as unknown:
        with open(sam_file,'r') as sam:
            for line in sam:
                if line.startswith("@"):
                    deduped.write(line)
        # Initialize the set of unique read identifiers and set counters to 0.
        unique_reads = set()
        chrom = ''
        unique_count = 0
        dupes_count = 0
        unknown_umi_count = 0
        # Create a dictionary of UMIs from the given UMI file.
        umi_dict = create_umi_dict()
        for file in split_files:
            # Only one of the two forward or reverse SAM files is read at a time, 
            # since there cannot be PCR duplicates across the two files.
            with open(file,'r') as fh:
                # Read through each line of both the forward and reverse SAM files and store their tuple of identifiers.
                for line in fh:
                    umi = find_umi(line)
                    # Write lines for reads with unknown UMI to a separate file.
                    if umi not in umi_dict:
                        unknown_umi_count += 1
                        unknown.write(line)
                        continue
                    # Increase counter for reads known UMIs.
                    else:
                        umi_dict[umi]+=1
                    # Since the files are ordered by chromosome, once chromosome changes, the set can be cleared.
                    if chrom != find_chrom(line):
                        unique_reads = set()
                    chrom = find_chrom(line)
                    # Use the correct starting position function for the forward and reverse reads.
                    if file == "fwd.sam":
                        start_pos = find_pos_fwd(line)
                    elif file == "rev.sam":
                        start_pos = find_pos_rev(line)
                    # Store the key identifiers as a tuple.
                    identifiers = (chrom,start_pos,umi)
                    # If the tuple has already been seen, write to duplicates file.
                    if identifiers in unique_reads:
                        dupes_count +=1
                        dupes.write(line)
                    # If the tuple has not been seen, write to de-duplicated file and add to set of tuples.
                    else:
                        unique_count +=1
                        unique_reads.add(identifiers)
                        deduped.write(line)
    return dupes_count,unique_count,unknown_umi_count,umi_dict

dupes_count,unique_count,unknown_umi_count,umi_dict = main_func()

# Remove the temporarily made 'sorted.sam' file and split files ('fwd.sam' and 'rev.sam') 
# which were created in the inital samtools steps.
os.remove("fwd.sam")
os.remove("rev.sam")
os.remove("sorted.sam")

# Print out summary of results to standard out.
print("\nFinished.")
print("\nSummary:",unique_count,"unique reads,",dupes_count,"PCR duplicates and",unknown_umi_count,"unknown UMIs.\n")
print("Table of Known UMI Counts\n")
for umi,count in umi_dict.items():
    print(umi,' : ',count)
