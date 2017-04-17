#!/bin/bash

: <<'END'
 *********************************************************************
 *   Consensus Finder 
 *
 *   Copyright (C) 2016 Bryan J. Jones (bryanjjones@gmail.com)
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ********************************************************************
END

usage()
{
cat << EOF
Usage: $0 [Option] [-q query file]		
	-h					This help text
	-s <number>				Set maximum sequences for BLAST search
	-e <number in scientific notation>	Set maximum e value for BLAST search
	-t <number>				Conservation shreshold for suggesting mutations
	-c					Do not use complete sequences, not only matched portions
	-i <number>				Iterations of Clustal Omega alingments
	-r <number>				CD-HIT redundancy threshold
	-v					verbose
	-l					turn on logging
	-k					keep all temp files for troubleshooting
	-q <FASTA file name>			REQUIRED protein sequence input file in FASTA format

Example usage:
$0 -s 1000 -e 1e-2 -t 0.9 -c -i 1 -r 0.9 -q FILE.FASTA

Consensus protein sequences are useful for numerous applications. Often, mutating a protein to be 
more like the consensus of homologs will often increase the stability of a protein, allowing it to  
function at higher temperatures, and have better soluble expression when expressed recombinatly in 
various hosts. "Consensus Finder" will help identify the consensus sequence and find potentially 
stabilizing mutations.

$0 will take given your given protein sequence, find similar sequences 
from the NCBI database, align them, remove redundant/highly similar sequences, trim alignments to 
the size of the original query, and analyze consensus. Output is trimmed alignment, consensus 
sequence, frequency and count tables for amino acids at each position, as well as a list of 
suggested mutations to consensus that may be stabilizing.
No flags required except -q with provided query file. Query file must be a protein sequence (not DNA)
text file in FASTA format.

BLAST parameters can be changed with "-s" to adjust maximum number of sequences and "-e" to adjust 
the maxumum e value. One or the other will limit the number of sequences returned. For both lower 
values will tend to return fewer more similar results. Defaults for "-s" 2000, with reasonable 
range would be 10-10,000, BLAST searchese will occasionally time out if too many sequences are 
requested. Default for "-e" 1e-3 (scientific notation), reasonable range would be 1e-30 to 1e-1.

Changing the threshold for "conserved residues" with "-t" will change how many mutation suggestions 
are returned. DEFAULT=0.6, reasonable range 0.05-0.99, with lower numbers returing more suggestions.

By default complete sequences will individually be downloaded from NCBI. This takes a lot of time, 
espicially with a lot of BLAST hits.  Adding "-c" flag uses only the partial sequences matching the 
target in the BLAST result, this can make the program run faster.

Increasing the number of Clustal Omega alingment iterations with "-i" will, in theory, will give better 
quality alignments, but it takes much longer and can give rise to other issues.  DEFAULT=1, 
reasonable options: integers from 1-5
The maximum threshold for eliminating redundant sequences with CD-HIT can be changed with "-r". The 
default is 0.9, which will remove any sequences with over 90% identity to prevent over sampling of 
over represented groups of proteins. Reasonable range would be 0.7-1.0

$0 should run on most 64bit linux systems with python and numpy installed.
Uses binaries (in "binaries" folder) for: blastp (2.4.0+), CD-HIT (4.6.4), & Clustal Omega (1.2.0).

Copyright 2016 Bryan J. Jones (bryanjjones@gmail.com)
$0 can be freely copied and distributed under the GNU General Public 
License version 2 (GPLv2) or later.

Citations:  
Weizhong Li, Lukasz Jaroszewski & Adam Godzik. "Clustering of highly homologous sequences to reduce 
  thesize of large protein database",  Bioinformatics, (2001) 17:282-283
Weizhong Li, Lukasz Jaroszewski & Adam Godzik. "Tolerating some redundancy significantly speeds up 
  clustering of large protein databases", Bioinformatics, (2002) 18:77-82
Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, 
  Thompson JD, Higgins DG. "Fast, scalable generation of high-quality protein multiple sequence 
  alignments using Clustal Omega." Mol Syst Biol. 2011 Oct 11;7:539. doi: 10.1038/msb.2011.75. 
  PMID: 21988835.
Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) 
  "BLAST+: architecture and applications." BMC Bioinformatics 10:421.

EOF
}
#Variables to change
MAXSEQ="2000" #Maximum number of sequences to return from BLAST search. Lower => fewer more similar results DEFAULT=2000, reasonable range would be 10-10,000 (lower numbers may run faster)
EVALUE="1e-3" #Maximum expect value for returned hits on BLAST search. Lower => fewer more similar results DEFAULT=1e-3 (scientific notation), reasonable range would be 1E-30 to 1e-1 (lower numbers may run faster)
THRESHOLD="0.6" #Threshold of conservation to suggest mutations, i.e. residue must be conserved higher than this to suggest mutation to consensus. DEFAULT=0.6, reasonable range 0.05-0.99
USE_COMPLETE_SEQS="1" #1" will use complete sequences of BLAST hits, while "0" will only use portions of the sequence that match. DEFAULT=1, reasonable options (BOOLEAN) 0 or 1. Setting to "0" will make the script run significantly faster.
ALIGNMENTITER="1" #Number of iterations of Clustal Omega alignment. In theory doing more iterations of the alignment will give better quality alignments, but it takes much longer and can give rise to other issues.  DEFAULT=1, reasonable options: integers from 1-5
ID_REDUNDANCY="0.9" #The maximum threshold for eliminating redundant sequences with CD-HIT, DEFAULT=0.9. Reasonable range would be 0.5-1.0
LOGGING="0" #default is logging off, if turned on log file will be saved in ./procissing/ folder
KEEPTEMP="0" #defaulet is to remove temporary files, only turn on for troubleshooting.
unset $SOURCE #make sure there is nothing stored as "SOURCE"

while getopts ":hs:e:t:ci:r:vklq:" opt; do
	case $opt in
		h)	
			usage
			exit 1
			;;
		s)
			MAXSEQ="$OPTARG"
			echo "Maximum BLAST sequences: $OPTARG" >&2
			;;
		e)
			EVALUE="$OPTARG"
			echo "Maximum BLAST e value: $OPTARG" >&2
			;;
		t)
			THRESHOLD="$OPTARG"
			echo "Minimum threshold for consensus: $OPTARG" >&2
			;;
		c)
			USE_COMPLETE_SEQS="0"
			echo "Use complete sequences: off" >&2
			;;
		i)
			ALIGNMENTITER="$OPTARG"
			echo "Clustal Omega alingment iterations: $OPTARG" >&2
			;;
		r)
			ID_REDUNDANCY="$OPTARG"
			echo "Maximum identity for redundant sequence cutoff: $OPTARG" >&2
			;;			
		v)	
			set -x #this turns on displaying every command for trouble shooting
			echo "verbose" >&2
			;;
		l)	
			LOGGING="1"
			echo "logging on" >&2
			;;
		k)	
			KEEPTEMP="1"
			echo "Temporary files will not be deleted" >&2
			;;
		q)
			SOURCE="$OPTARG"
			echo "query file: $OPTARG" >&2
			;;
		?)
			echo "Invalid option: -$OPTARG. Run Consensus Finder with -h flag for help." >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument. Run Consensus Finder with -h flag for help." >&2
			exit 1
			;;

	esac
done


#checking the settings

if [ "$LOGGING" = "1" ]; then
	exec 3>&1 4>&2 #Saves file descriptors 
	trap 'exec 2>&4 1>&3' 0 1 2 3 #Restore file descriptors for particular signals
	exec 1>./processing/${SOURCE}_log.out 2>&1 #Redirect stdout to file log.out then redirect stderr to stdout. 
fi
BLAST="./binaries/blastp" #change this if another version is installed locally
CDHIT="./binaries/cd-hit" #change this if another version is installed locally
CLUSTAL="./binaries/clustalo-1.2.0-Ubuntu-x86_64" #change this if another version is installed locally
#set -e #will cause the script to exit if any of the steps produce an error




#CHECKS TO MAKE SURE ALL THE FILES AND FOLDERS ARE WHERE THEY NEED TO BE

if [ ! -d "uploads" ]; then #make sure uploads directory exists
  mkdir uploads #make directory if it doesn't exsit
fi
if [ ! -d "processing" ]; then #make sure processing directory exists
  mkdir processing #make directory if it doesn't exsit
fi
if [ ! -d "completed" ]; then #make sure completed directory exists
  mkdir completed #make directory if it doesn't exsit
fi
if [ -z "$SOURCE" ] #make sure querry was given
then
   echo "NO QUERY FILE SPECIFIED"
   usage
   exit 1
fi
if test -f $SOURCE; then #check if the source file exists
	if sed 's: :_:g' ./$SOURCE | sed -r 's:>.+::g' | sed ':a;N;$!ba;s/\n//g' | grep -E [FLIMVSPYHQKDEWR]; then #see if the sequence contains amino acid letters in addition to ACGTN
		echo "$SOURCE looks like a protein sequence" 
	else 
		echo "$SOURCE does not look like a protein sequence, is it DNA? EXITING."
		exit 1
	fi
	mv $SOURCE ./processing/$SOURCE #move the source file to the processing directory
elif test -f ./uploads/$SOURCE; then
	if sed 's: :_:g' ./uploads/$SOURCE | sed -r 's:>.+::g' | sed ':a;N;$!ba;s/\n//g' | grep -E [FLIMVSPYHQKDEWR]; then #see if the sequence contains amino acid letters in addition to ACGTN
		echo "$SOURCE looks like a protein sequence" 
	else 
		echo "$SOURCE does not look like a protein sequence, is it DNA? EXITING."
		exit 1
	fi
	mv ./uploads/$SOURCE ./processing/$SOURCE
else 
	echo "COULD NOT FIND $SOURCE IN EITHER $PWD OR IN $PWD/uploads, ARE YOU SURE YOU SPECIFIED THE CORRECT FILE?"
	exit 1
fi 
if ! [ $(grep -c \> ./processing/$SOURCE) -eq 1 ]; then
	echo "$SOURCE LOOKS LIKE IT CONTAINS MORE THAN ONE SEQUENCE. PLEASE ENTER A FILE WITH ONLY ONE SEQUENCE."
	if ! [ "$KEEPTEMP" = "1" ]; then #if keeptemp is on, keep temp files
		rm -f  ./processing/$SOURCE #cleaning up the temporary files
	fi
	exit 1
fi
rm -f ./processing/1BLAST.out ./processing/2BLAST-FASTA.out ./processing/3CDHIT.out ./processing/3CDHIT.out.clstr ./processing/4QUERY_CDHIT.out ./processing/5CLUSTAL.out ./processing/6SEQSONLY.out ./processing/7VERSION.out #make sure temporary files don't exist yet	


#Now everything looks good, STARTING THE ACTUAL RUN
echo "starting run:" 1>&2
echo "looking for BLAST hits in the nr database (This may take a while)"
REDUCED_MAX_SEQS=0 #
timeout --signal=KILL 15m $BLAST -db nr -query ./processing/$SOURCE -out ./processing/1BLAST.out -evalue $EVALUE -max_target_seqs $MAXSEQ -outfmt "6 sacc sseq pident" -remote #BLAST ncbi nr database for up to MAXSEQ matches with expect values of EVALUE or less
if ! [ 0 -eq $? ] ; then 
	##NEED TO MAKE THIS TRY AGAIN W/ FEWER HITS REQUESTED
	echo "BLAST TIMED OUT, TRYING AGAIN WITH LOWER MAXSEQ. ONLY ASKING FOR 200 HITS NOW."
	REDUCED_MAX_SEQS=1
	timeout --signal=KILL 15m $BLAST -db nr -query ./processing/$SOURCE -out ./processing/1BLAST.out -evalue $EVALUE -max_target_seqs 200 -outfmt "6 sacc sseq pident" -remote #BLAST ncbi nr database for up to 200 matches with expect values of EVALUE or less	
	if ! [ 0 -eq $? ] ; then 
		echo "BLAST FAILED AGAIN, NOT SURE WHY. TRY AGAIN LATER. EXITING"
		if ! [ "$KEEPTEMP" = "1" ]; then #if keeptemp is on, keep temp files
			rm -f ./processing/1BLAST.out ./processing/2BLAST-FASTA.out ./processing/3CDHIT.out ./processing/3CDHIT.out.clstr ./processing/4QUERY_CDHIT.out ./processing/5CLUSTAL.out ./processing/6SEQSONLY.out ./processing/7VERSION.out ./processing/$SOURCE #cleaning up the temporary files
			rm -f ./processing/SEQ* #remove temp files
		fi
		exit 1
	fi
fi
COUNT=`wc -l ./processing/1BLAST.out | cut -f1 -d' '` #count lines of file
if ! ls ./processing/1BLAST.out ; then 
	echo "BLAST SEARCH FAILED, NO OUTPUT FOUND. EXITING"
	if ! [ "$KEEPTEMP" = "1" ]; then #if keeptemp is on, keep temp files
		rm -f ./processing/1BLAST.out ./processing/2BLAST-FASTA.out ./processing/3CDHIT.out ./processing/3CDHIT.out.clstr ./processing/4QUERY_CDHIT.out ./processing/5CLUSTAL.out ./processing/6SEQSONLY.out ./processing/7VERSION.out ./processing/$SOURCE #cleaning up the temporary files
		rm -f ./processing/SEQ* #remove temp files
	fi
	exit 1
elif [ "$COUNT" = "0" ]; then #if the BLAST file is empty
	echo "ERROR: BLAST SEARCH RETURNED ZERO RESULTS. EXITING"
	if ! [ "$KEEPTEMP" = "1" ]; then #if keeptemp is on, keep temp files
		rm -f ./processing/1BLAST.out ./processing/2BLAST-FASTA.out ./processing/3CDHIT.out ./processing/3CDHIT.out.clstr ./processing/4QUERY_CDHIT.out ./processing/5CLUSTAL.out ./processing/6SEQSONLY.out ./processing/7VERSION.out ./processing/$SOURCE #cleaning up the temporary files
		rm -f ./processing/SEQ* #remove temp files
	fi
	exit 1	
fi
COUNT=`wc -l ./processing/1BLAST.out | cut -f1 -d' '` #count lines of file
PCT=`tail -1 ././processing/1BLAST.out | grep -oE [0-9]+\.[0-9][0-9]\$`
echo "Found $COUNT hits with at $PCT% or greater identy to $SOURCE"
VERSIONS=( $(sed -r 's:'\(\\w\{1,\}\).*':'\\1':' ./processing/1BLAST.out) )
#OLD ONE USING GI NUMBERS: VERSIONS=( $(sed -r 's:'gi.\([0-9]\{1,\}\).*':'\\1':' ./processing/1BLAST.out) ) #store array of gi values from blast search results
if [ "$USE_COMPLETE_SEQS" = "1" ]; then
	echo "Downloading complete sequences from NCBI."
	DOWNLOADS=0
	for((i=0;i<${#VERSIONS[@]};i++)); do #fetch sequence
		DOWNLOADS=$((DOWNLOADS+1)) #COUNTER for counting simultanequsly downloading sequences.
		if [[ "$DOWNLOADS" -gt 10 ]]; then #up to 10 simultaneous curl commands, otherwise, wait until they are all done
			echo "waiting for current dowloads to finish"
			wait
			DOWNLOADS="0" #reset the counter
		fi 
		VERSION=">${VERSIONS[$i]}" #make name with ">" for fasta formatting
		eval VERSION$i=\$VERSION
		(echo "fetching sequence $VERSION"
		SEQ=$(curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&amp;id=${VERSIONS[i]}&amp;retmode=text&amp;rettype=fasta" | sed -r 's:'\>.+'::g' | sed ':a;N;$!ba;s/\n//g') #retrieve FASTA files from NCBI corrisponding to each GI value, sed commands extract only sequence info from files.
		if echo $SEQ | grep -q "Error"; then echo "Failed to return sequence, trying again" #try retrieving a few more times if it fails (retrieves data containing the word "Error"
			SEQ=$(curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&amp;id=${VERSIONS[i]}&amp;retmode=text&amp;rettype=fasta" | sed -r 's:'\>.+'::g' | sed ':a;N;$!ba;s/\n//g')
			if echo $SEQ | grep -q "Error"; then echo "Failed to return sequence, trying again"
				SEQ=$(curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&amp;id=${VERSIONS[i]}&amp;retmode=text&amp;rettype=fasta" | sed -r 's:'\>.+'::g' | sed ':a;N;$!ba;s/\n//g')
				if echo $SEQ | grep -q "Error"; then echo "Failed to return sequence, trying again"
					SEQ=$(curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&amp;id=${VERSIONS[i]}&amp;retmode=text&amp;rettype=fasta" | sed -r 's:'\>.+'::g' | sed ':a;N;$!ba;s/\n//g')
					if echo $SEQ | grep -q "Error"; then echo "Failed to return sequence with four tries, trying once more"
						SEQ=$(curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&amp;id=${VERSIONS[i]}&amp;retmode=text&amp;rettype=fasta" | sed -r 's:'\>.+'::g' | sed ':a;N;$!ba;s/\n//g')
						if echo $SEQ | grep -q "Error"; then echo "FIVE ATTEMPTS FAILED TO RETURN SEQUENCE FROM NCBI, EXITING"
							if ! [ "$KEEPTEMP" = "1" ]; then #if keeptemp is on, keep temp files
								rm -f ./processing/1BLAST.out ./processing/2BLAST-FASTA.out ./processing/3CDHIT.out ./processing/3CDHIT.out.clstr ./processing/4QUERY_CDHIT.out ./processing/5CLUSTAL.out ./processing/6SEQSONLY.out ./processing/7VERSION.out ./processing/$SOURCE #cleaning up the temporary files
								rm -f ./processing/SEQ* #remove temp files
							fi
							exit 1 #if we can't get the sequence, exit the substring, hopefully exit the parent string too
						fi
					fi
				fi
			fi
		fi
		echo $VERSION > ./processing/SEQ$i #add the gi number to a temporary file with the number
		echo $SEQ >> ./processing/SEQ$i ) & #add the sequence to the temp file. ") &" will allow this to run in parallel
	done
	wait #make sure any lagging substrings are done
	: > ./processing/2BLAST-FASTA.out #make sure the file exists and is empty
	for((i=0;i<${#VERSIONS[@]};i++)); do #combine files
		cat ./processing/SEQ$i>>./processing/2BLAST-FASTA.out
	done
else # i.e. if  "$USE_COMPLETE_SEQS" =/= "1" 
	sed -r 's:^:>:' ./processing/1BLAST.out > ./processing/2BLAST-FASTA.out #add ">" to name
	sed -r 's:\t:\n:' -i ./processing/2BLAST-FASTA.out #names and sequences to seperate lines (fasta format)
	sed -r 's:-::g' -i ./processing/2BLAST-FASTA.out #remove "-" from sequences to prepare for CD-HIT grouping
	sed -r 's:\t.*$::' -i ./processing/2BLAST-FASTA.out #remove % id from end of sequences
fi
$CDHIT -i ./processing/2BLAST-FASTA.out -o ./processing/3CDHIT.out -c $ID_REDUNDANCY -M 0 -T 0 #Remove redundant sequences
if ! [ 0 -eq $? ] ; then 
	echo "CD-HIT FAILED, EXITING"
	if ! [ "$KEEPTEMP" = "1" ]; then #if keeptemp is on, keep temp files
		rm -f ./processing/1BLAST.out ./processing/2BLAST-FASTA.out ./processing/3CDHIT.out ./processing/3CDHIT.out.clstr ./processing/4QUERY_CDHIT.out ./processing/5CLUSTAL.out ./processing/6SEQSONLY.out ./processing/7VERSION.out ./processing/$SOURCE #cleaning up the temporary files
		rm -f ./processing/SEQ* #remove temp files
	fi
	exit 1
fi
echo "Aligning those hits."
sed 's: :_:g' ./processing/$SOURCE > ./processing/4QUERY_CDHIT.out #add query sequence
cat ./processing/3CDHIT.out >> ./processing/4QUERY_CDHIT.out #add other sequences to query
$CLUSTAL --iter=$ALIGNMENTITER -i ./processing/4QUERY_CDHIT.out -o ./processing/5CLUSTAL.out --outfmt=fa -v -v --force #Aligning hits
if ! [ 0 -eq $? ] ; then 
	echo "CLUSTAL ALIGNMENT FAILED, EXITING"
	if ! [ "$KEEPTEMP" = "1" ]; then #if keeptemp is on, keep temp files
		rm -f ./processing/1BLAST.out ./processing/2BLAST-FASTA.out ./processing/3CDHIT.out ./processing/3CDHIT.out.clstr ./processing/4QUERY_CDHIT.out ./processing/5CLUSTAL.out ./processing/6SEQSONLY.out ./processing/7VERSION.out ./processing/$SOURCE #cleaning up the temporary files
		rm -f ./processing/SEQ* #remove temp files
	fi
	exit 1
fi
echo "Done Aligning"
SEQS=( $(sed -r 's:>.+:break:g' ./processing/5CLUSTAL.out | sed ':a;N;$!ba;s/\n//g' | sed -r 's:break:\n:g')) #store array of sequences
for SEQ in "${SEQS[@]}"; do echo $SEQ; done > ./processing/6SEQSONLY.out #save sequences for troubleshooting
grep '>' ./processing/5CLUSTAL.out > ./processing/7VERSION.out
echo "Trimming to $SOURCE sequence."
export VERSIONS
export THRESHOLD #export variable so python code can use it
export SOURCE #export file name so python code can use it
export REDUCED_MAX_SEQS
#python to calculate frequencies of residues
python <<END_OF_PYTHON
import sys #need sys to 
import numpy as np # need numpy for arrays and the like
import os # need os to be able to pull in variables like "SOURCE" and "THRESHOLD"
SOURCE = os.environ['SOURCE'] #pull source file name for naming output files
THRESHOLD = float(os.environ['THRESHOLD']) #pull thershold as float for suggesting mutation to consensus from bash script
REDUCED_MAX_SEQS = os.environ['REDUCED_MAX_SEQS'] #pull REDUCED_MAX_SEQS to know if the max sequences had to be reduced
VERSIONS = np.genfromtxt('./processing/7VERSION.out', dtype=str) #make a list to get GI numbers (initially contains sequences too)
AA = np.genfromtxt('./processing/6SEQSONLY.out', delimiter=1, dtype=str)#make array from file with each AA as an element
LENGTH = len(AA[0,:]) #how many positions do we need to look for gaps?
for index in range(len(AA[0,:])) :
	AA[0,index]=AA[0,index].upper()
for index in range(LENGTH):#at each position
	if AA[0,(LENGTH-1-index)] == "-":#LENGTH-1-index will start at the end, -1 to account for 0 based indexing, and find gaps
		AA=np.delete(AA, (LENGTH-1-index), 1)#delete the gaps
		print ''.join(AA[0,:])#print the target sequence to show progress
SEQS=[]
for index in range(len(AA[:,0])): #alternate ">[GI number]" and sequences
	SEQS=np.hstack((SEQS,VERSIONS[(index)]))
	SEQS=np.hstack((SEQS,(''.join(AA[(index),:]))))
np.savetxt(('./completed/'+SOURCE+'_trimmed_alignment.fst'),SEQS,delimiter="",fmt="%s") #save file with AA sequence of consensus sequence

COUNTS = np.zeros([22, len(AA[0,:])],dtype=object) #makes array the length of the alingment with 22 rows (20AAs + "-" + "other")
for index in range(len(AA[0,:])): # for each position along the alingment, count occourances of each AA
	COUNTS[0,index]=AA[:,index].tolist().count("G")
	COUNTS[1,index]=AA[:,index].tolist().count("P")
	COUNTS[2,index]=AA[:,index].tolist().count("A")
	COUNTS[3,index]=AA[:,index].tolist().count("V")
	COUNTS[4,index]=AA[:,index].tolist().count("L")
	COUNTS[5,index]=AA[:,index].tolist().count("I")
	COUNTS[6,index]=AA[:,index].tolist().count("M")
	COUNTS[7,index]=AA[:,index].tolist().count("C")
	COUNTS[8,index]=AA[:,index].tolist().count("F")
	COUNTS[9,index]=AA[:,index].tolist().count("Y")
	COUNTS[10,index]=AA[:,index].tolist().count("W")
	COUNTS[11,index]=AA[:,index].tolist().count("H")
	COUNTS[12,index]=AA[:,index].tolist().count("K")
	COUNTS[13,index]=AA[:,index].tolist().count("R")
	COUNTS[14,index]=AA[:,index].tolist().count("Q")
	COUNTS[15,index]=AA[:,index].tolist().count("N")
	COUNTS[16,index]=AA[:,index].tolist().count("E")
	COUNTS[17,index]=AA[:,index].tolist().count("D")
	COUNTS[18,index]=AA[:,index].tolist().count("S")
	COUNTS[19,index]=AA[:,index].tolist().count("T")
	COUNTS[20,index]=AA[:,index].tolist().count("-") #empty spaces
	COUNTS[21,index]=(len(AA[:,index]) - sum(COUNTS[:,index].tolist())) #other, not counted above
#make an array of AA codes as a key for the above counts
IDS = np.zeros([22,1],dtype=object)
IDS[0,0]="G"
IDS[1,0]="P"
IDS[2,0]="A"
IDS[3,0]="V"
IDS[4,0]="L"
IDS[5,0]="I"
IDS[6,0]="M"
IDS[7,0]="C"
IDS[8,0]="F"
IDS[9,0]="Y"
IDS[10,0]="W"
IDS[11,0]="H"
IDS[12,0]="K"
IDS[13,0]="R"
IDS[14,0]="Q"
IDS[15,0]="N"
IDS[16,0]="E"
IDS[17,0]="D"
IDS[18,0]="S"
IDS[19,0]="T"
IDS[20,0]="-"
IDS[21,0]="other" 

FREQS = np.zeros_like(COUNTS) #make an array for calculating frequencies of each AA
FREQS = np.float64(FREQS) #it needs to be numbers
for index in range(len(FREQS[0,:])): # calculate the frequencey of each AA as [occurrences]/[occurrences of all AAs], "-" and "other" not counted in total
	FREQS[0,index]=np.float64(COUNTS[0,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[1,index]=np.float64(COUNTS[1,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[2,index]=np.float64(COUNTS[2,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[3,index]=np.float64(COUNTS[3,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[4,index]=np.float64(COUNTS[4,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[5,index]=np.float64(COUNTS[5,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[6,index]=np.float64(COUNTS[6,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[7,index]=np.float64(COUNTS[7,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[8,index]=np.float64(COUNTS[8,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[9,index]=np.float64(COUNTS[9,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[10,index]=np.float64(COUNTS[10,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[11,index]=np.float64(COUNTS[11,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[12,index]=np.float64(COUNTS[12,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[13,index]=np.float64(COUNTS[13,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[14,index]=np.float64(COUNTS[14,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[15,index]=np.float64(COUNTS[15,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[16,index]=np.float64(COUNTS[16,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[17,index]=np.float64(COUNTS[17,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[18,index]=np.float64(COUNTS[18,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[19,index]=np.float64(COUNTS[19,index])/sum(np.float64(COUNTS[:20,index]))
	FREQS[20,index]=np.float64(COUNTS[20,index])/sum(np.float64(COUNTS[:,index])) #frequency of gaps "-" as fraction of all seqs
	FREQS[21,index]=np.float64(COUNTS[21,index])/sum(np.float64(COUNTS[:,index])) #frequency of gaps "other" as fraction of all seqs

IDCOUNTS = np.hstack((IDS,COUNTS)) #make list with AA counts and names of AAs
np.savetxt(('./completed/'+SOURCE+'_counts.csv'),IDCOUNTS,delimiter=",",fmt="%s") #save file with AA names and counts

IDFREQS = np.hstack((IDS,FREQS)) #make list with names and AA frequencies
np.savetxt(('./completed/'+SOURCE+'_frequencies.csv'),IDFREQS,delimiter=",",fmt="%s") #save file with AA names and frequencies

CONSENSUS_SEQ = np.zeros([1, len(FREQS[0,:])],dtype=object) #make an array to store consensus sequence
TYPES = [('wt', 'S1'), ('res', int), ('sug', 'S1'),('freq', float)] #make a list of data types for suggested mutation tuples
MUTATIONS_ARRAY=np.empty([0,])
MUTATIONS_ARRAY=np.array(MUTATIONS_ARRAY, dtype=TYPES)
for index in range(len(FREQS[0,:])): #for each AA position
	CONSENSUS_SEQ[0,index] = IDS[np.argmax(COUNTS[:20,index]),0] #find the largest value, and get the corrisponding AA from IDS, and add it to CONSENSUS_SEQ
	if THRESHOLD < max(FREQS[:20,index]): #if the consensus of a residue is greater than the threshold
		if IDS[np.argmax(FREQS[:20,index]),0] != AA[0,index]: #and the consens residue is different than the first sequence
			print "Residue number " + str(int(index) + 1)
			print str(int(100*max(FREQS[:20,index]))) + "% is greater than or equal to " + str(int(100*THRESHOLD)) + "%"
			SUGGESTION=np.array([(AA[0,index], int(index + 1), CONSENSUS_SEQ[0,index], (-1*max(FREQS[:20,index])))], dtype=TYPES) #entry with negative frequency to allow easy sorting.
			MUTATIONS_ARRAY = np.append(MUTATIONS_ARRAY,SUGGESTION, axis=0)#add new suggestion on to any existing "MUTATIONS_ARRAY"
MUTATIONS_ARRAY=np.sort(MUTATIONS_ARRAY, order='freq')       
SUGGESTED_MUTATIONS=np.zeros([1],dtype=object)
SUGGESTED_MUTATIONS[0]="These mutations may stabilize your protein since they differ from a consensus of over " + str(100*THRESHOLD) + "%"
if REDUCED_MAX_SEQS == '1':
	SUGGESTED_MUTATIONS = np.vstack((SUGGESTED_MUTATIONS,("Results based upon up to 200 BLAST results. BLAST timed out when trying to retrieve the requested number of hits (default request is 2000)."))) #add note to let user know that results are based on fewer sequences
for index in range(len(MUTATIONS_ARRAY[:,])): #for each suggested mutation
	SUGGESTED_MUTATIONS = np.vstack((SUGGESTED_MUTATIONS,("Change " + MUTATIONS_ARRAY[index,]['wt'] + " " + str(MUTATIONS_ARRAY[index,]['res']) + " to " + MUTATIONS_ARRAY[index,]['sug'] + " (found in " + str(int(-100*MUTATIONS_ARRAY[index,]['freq'])) + "% of similar proteins)" ))) #add new suggestion on to any existing "SUGGESTED_MUTATIONS"
print SUGGESTED_MUTATIONS

CONSENSUS=""
for index in range(len(CONSENSUS_SEQ[0,:])):
	CONSENSUS=CONSENSUS+str(CONSENSUS_SEQ[0,index])
CONSENSUS=">consensus_sequence",CONSENSUS # add header for FASTA format
np.savetxt(('./completed/'+SOURCE+'_consensus.fst'),CONSENSUS,delimiter="",fmt="%s") #save file with AA sequence of consensus sequence
np.savetxt(('./completed/'+SOURCE+'_mutations.txt'),SUGGESTED_MUTATIONS,delimiter=",",fmt="%s") #save file with suggested stabilizing mutations
sys.exit(0)
END_OF_PYTHON
mv ./processing/$SOURCE ./completed/$SOURCE
if ! [ "$KEEPTEMP" = "1" ]; then #if keeptemp is on, keep temp files
	rm ./processing/1BLAST.out ./processing/2BLAST-FASTA.out ./processing/3CDHIT.out ./processing/3CDHIT.out.clstr ./processing/4QUERY_CDHIT.out ./processing/5CLUSTAL.out ./processing/6SEQSONLY.out ./processing/7VERSION.out #cleaning up the temporary files
	rm ./processing/SEQ* #remove temp files
fi
echo "Run complete." 1>&2

