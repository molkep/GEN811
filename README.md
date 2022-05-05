#The point of this project is to expand and correct a previously started python code that locates ORFs from a FASTA file of DNA sequences of the users specifications and turning it itno a GFF file by comparing it to a reference genome and doing a local BLAST (database also specified by the usuer)
##
##
#The project was sucessful in finding ORFs in the test FASTA files and creating a GFF like file, but the orf finder function had some issues with going back into an orf and finding another orf that was not real
#the other issues is that orfs from different organisms were found and recorded after the blast
#these are some issues I would like to fix

![image](https://user-images.githubusercontent.com/83464534/166949877-c85dc343-ccd3-415f-92f6-6d526a015cff.png)
Above is a flowchart of how the script runs and takes user input and impliments it
