

#64 codons, 3 stop
#probability that you will see a stop codon: 3/64 -> 1/21(5% ish)
#how long is a gene? -> start stop should be longer than 21 codons, longer it is the more likely its an orf
#run BLAST within the program, do not outsource or it will make the code run forever

#import necessary modules
import argparse
import os
import subprocess
import csv

#function to parse out fasta headers and create one contiguous sequence
def fileparse():
    sequences = []
    descr = None
    try:
        #open specified file from user
        with open("{0}".format(args.input), "r") as fasta_file:
            line = fasta_file.readline()[:-1]  # always trim newline
            while line:
                if line[0] == '>':
                    if descr:  # any sequence found yet?
                        sequences.append((descr, seq))
                    descr = str(line[1:].split('>'))
                    seq = ''  # start a new sequence
                else:
                    seq += line
                line = fasta_file.readline()[:-1]
            sequences.append((descr, seq))
    except IOError:
        print("Error reading file")
    #return string of sequence data
    return sequences

#function for detecting orfs in sequence data
def orffinder(sequences):
    listOfOrf = list()
    for index, value in enumerate(sequences):  # looping over the fragments extracted
        frames = [] # storing the six frame translation that it zould be extacted from the fragments 
        dna = value[1]  # extract the fragment 
        description = value[0] #extact the desciption even were not use it, just for learning purpose
        reverseCdna = [] # storing the reverse compliments
    # create the positive frames
    # split the frames into codons for better performance
        frames.append([dna[i:i + 3] for i in range(0, len(dna), 3)])
        frames.append([dna[i:i + 3] for i in range(1, len(dna), 3)])
        frames.append([dna[i:i + 3] for i in range(2, len(dna), 3)])
    # reverse compliment of the fragment
        reverse = {"A": "T", "C": "G", "T": "A", "G": "C"}
        for i in range(len(dna)):
            reverseCdna.append(reverse[dna[-i - 1]]) if dna[-i - 1] in reverse.keys() else reverseCdna.append(dna[-i - 1])  # if any contamination found we keep it for further more check
        reverseCdna = ''.join(reverseCdna) # joining 
    # create the negative frames
        frames.append([reverseCdna[i:i + 3] for i in range(0, len(reverseCdna), 3)])
        frames.append([reverseCdna[i:i + 3] for i in range(1, len(reverseCdna), 3)])
        frames.append([reverseCdna[i:i + 3] for i in range(2, len(reverseCdna), 3)])
   
    for i in range(0,len(frames),1): #looping all the frames
        start=0
        while start <len(frames[i]): #looping each frame for start and stop codons 
            if frames[i][start]=="ATG":
                for stop in range(start+1,len(frames[i]),1):
                         if frames[i][stop]=="TAA" or  frames[i][stop]=="TAG" or  frames[i][stop]=="TGA" :
                                listOfOrf.append(frames[i][start:stop]) # retrieve the orf 
                                start=stop+1 # avoiding multiple start codons
                                break
            start+=1
    return listOfOrf

#function to filter orfs by length
def filterorfs():
    count = 1
    try:
        #open file to write orfs into in fasta format
        with open("orf.fasta", "w") as orf_file:
            #if the orf found is longer than 300 bps (100 amino acids), it is most likely a true orf
            for orf in orffinder(fileparse()):
                orf = "".join(orf)
                if len(orf) >= 300:
                    if len(orf) <= 1050:
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
    blast_run = subprocess.run("blastn -query orf.fasta -db db -out orf_results.txt  -outfmt '6 qseqid sseqid sstart send score qframe qstrand salltitles'", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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

 
