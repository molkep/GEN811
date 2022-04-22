

#64 codons, 3 stop
#probability that you will see a stop codon: 3/64 -> 1/21(5% ish)
#how long is a gene? -> start stop should be longer than 21 codons, longer it is the more likely its an orf
#run BLAST within the program, do not outsource or it will make the code run forever

#import necessary modules
import argparse
import Bio.Seq
import os
import subprocess
import csv

#function to parse out fasta headers and create one contiguous sequence
def fileparse():
    header = ""
    seq = ""
    try:
        #open specified file from user
        with open("{0}".format(args.input), "r") as fasta_file:
            #for each line, look for a greater than sign. if its there,
            #continue to next line. if it doesnt (aka sequence data), add to a string until one long string of data is created
            for line in fasta_file:
                line=line.rstrip()
                if line.startswith(">"):
                    continue
                else:
                    seq += line
    except IOError:
        print("Error reading file")
    #return string of sequence data
    return seq

#function for detecting orfs in sequence data
def orffinder(seq):
    #iterate over each frame, 1-6
    for frame in range(1, 6):
        #make the sequence a seq instance
        sequence = Bio.Seq.Seq("{0}".format(seq))
        #create a reverse complement of sequence
        revcomp = sequence.reverse_complement()
        #if the frame is 1-3, use the normal sequence. if 4-6, use revcomp
        strand = sequence if frame < 4 else revcomp
        strand = strand.upper()
        #reference lists of all start and stop codons
        start_codon = ["ATG"]
        stop_codons = ["TAA", "TAG", "TGA"]
        #start scanning the sequence for the first start/stop codon and have it look at every 3 nucleotides
        for startorf in range(frame % 3, len(strand), 3):
           #define what a codon looks like
            begin = strand[startorf:startorf+3]
            #if the codon is in either of the reference lists, have that be the start of the orf
            if begin in start_codon or begin in stop_codons:
                coord1 = startorf
                #start scanning for the stop codon, aka the end of the orf, at the beginning of the start/stop codon
                for endorf in range(coord1, len(strand), 3):
                    end = strand[endorf:endorf+3]
#if a stop codon is hit, have that be the end of the orf
                    if end in stop_codons:
                        coord2 = endorf
                        #index the strand by the first and last codons found and save to a yeild function
                        yield (strand[coord1:coord2+3])
                        break
        return

#function to filter orfs by length
def filterorfs():
    count = 1
    try:
        #open file to write orfs into in fasta format
        with open("orf.fasta", "w") as orf_file:
            #if the orf found is longer than 300 bps (100 amino acids), it is most likely a true orf
            for orf in orffinder(fileparse()):
                if len(orf) >= 300:
                    if len(orf) <= 1200:
                       #a header for each orf is written followed by a number that increases with each orf added
                       orf_file.write("> seq{0}\n".format(count))
                       orf_file.write("{0}\n".format(orf))
                       count += 1
    except IOError:
        print("Error with file creation!")
    return orf_file

#function to perform the BLAST functionality
def BLAST():
    #because ron does not have any preloaded databases, the user has to specify what database they want to use from an outside platform. this is renamed toa generic term for easy use
    database_get = subprocess.run("wget -O database.fna.gz {0}".format(args.db), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #most databases come compressed so they need to be decompressed before use
    decompress = subprocess.run("gzip -d database.fna.gz", shell=True, stdout=subprocess.PIPE)
    #the database has to be created on the local directory
    db_create = subprocess.run("makeblastdb -in database.fna -parse_seqids -dbtype nucl -out db", shell=True, stdout=subprocess.PIPE)
    #the actual blast command. the output file is in format 6, which can be tailored to give all the information necessary fora gff file format
    blast_run = subprocess.run("blastn -query orf.fasta -db db -out orf_results.txt  -outfmt '6 qseqid sseqid sstart send score sframe sstrand salltitles'", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #remove fasta file, as is no longer needed
    os.remove("orf.fasta")
    return

#function to create gff file
def gffmaker():
    #empty dictionary for writting orf and information to
    orf_dict = {}
    try:
        #open file with user specified name to write to
        with open("{0}.gff".format(args.output), "w", newline="") as gff_file:
            writer = csv.writer(gff_file, delimiter="\t", quotechar="|")
            #open blast result document to parse
            with open("orf_results.txt", "r") as blast_results:
                #for each blast result, have it split on the tab so that each entry can be easily added to a dictionary
                for entry in blast_results:
                    entry = entry.rstrip()
                    entry = entry.split("\t")
                    # since the result file returns multiple blast results for each entry, only the top result is needed, so skip all others if it is already in the dictionary
                    if entry[0] in orf_dict.keys():
                        continue
                    #dict key is the sequence number from the fasta file, and the value is the information returned from the blast
                    orf_dict[entry[0]]=entry[1:]
                #iterating through the values, the data from the blast query can be parsed and a gff file entry can be constructed (tab delimited) from it
                for value in orf_dict.values():
                    writer.writerow(["{0}".format(value[0].split("|")[1]), "NCBI", "ORF", "{0}".format(value[1]), "{0}".format(value[2]), "{0}".format(value[3]), "{0}".format(value[5]), "{0}".format(value[4]), "{0}".format(value[6])])
    except IOError:
        print("Error handling files!")
    return

#positional argments for input file name, database URL, and output file name
parser = argparse.ArgumentParser()
parser.add_argument("input", help="input file name")
parser.add_argument("db", help="link to database to BLAST orfs against")
parser.add_argument("output", help="output file name")
args = parser.parse_args()

#running the functions with progress markers
print("Running sequence analysis...")
filterorfs()
print("ORFs found, BLAST query running...")
BLAST()
print("BLAST completed, compiling results...")
gffmaker()
print("File generation complete!")

#removing unnecessary files from the working directroy after completion of the  file creation
os.remove("orf_results.txt")
os.remove("database.fna")
remove = subprocess.run("rm db.*", shell=True, stdout=subprocess.PIPE)

 
