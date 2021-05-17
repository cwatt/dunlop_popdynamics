#!/usr/bin/env python3

# Read each line of a representative sequences fasta file and filter sequences of interested based on ASV label.
# Originally designed for use with PAPRICA
# Authored by Cassandra Wattenburger 10/05/20

import sys

Usage = """Subset Representative Sequences v1.0
By Cassi Wattenburger 10/06/20
Usage: subset_repsequences.py sequences.fasta listasvs.txt header
Arguements: 
Arg 1 is your fasta file (dereplicated or representative sequences only)
Arg 2 is the list of ASVs of interest in .txt format, each ASV on a separate line (no > at start of name)
Arg 3 is the number of header lines in the list (ie 1, enter 0 if none)
Arg 4 is the output file name (ie output.fasta)
Output: Creates a fasta file with only sequences from the ASV list"""

# Check for arguments
if len(sys.argv) < 5:
   print(Usage)

else:
   # Read in files/arguments
   InfileSeqsname = sys.argv[1] # FASTA file with representative sequences
   InfileListname = sys.argv[2] # List of ASVs
   Header = sys.argv[3] # Lines to skip in ASV list for header
   Outfilename = sys.argv[4] # Name of output file
   if Debug:
      print("Sequences file: " + InfileSeqsname)
      print("ASVs file: " + InfileListname)
      print("Skip " + Header + " line(s) in ASV list.")
      print("Output file: "+Outfilename)

   # Open and read in ASV list and sequence files
   with open(InfileListname) as n:
      InfileList = n.readlines()
   with open(InfileSeqsname) as f:
      InfileSeqs = f.readlines()

   # Read each line of ASV list
   LineNum = 0
   for LineList in InfileList:
      if LineNum < int(Header): # Skip designated number of lines
         LineNum += 1
      else:
         LineList = LineList.strip()
         LineNum += 1
         print("ASV to search: "+LineList)

         # Read fasta file line by line and search for ASV match
         LineSeqASV = ""
         Match = False
         for LineSeq in InfileSeqs:
            LineSeq = LineSeq.strip()
            if LineSeq[0] == '>': # > denotes sequence name
               LineSeqASV = LineSeq[1:]
               if LineList in LineSeqASV: # Check if ASV in list matches sequence name (> removed)
                  with open(Outfilename, "a") as o:
                     o.write(LineSeq + "\n")
                  print("Match found!")
                  Match = True
               else:
                  Match = False
            elif Match: # Save the next line if names matched (this is the associated sequence)
               with open(Outfilename, "a") as o:
                  o.write(LineSeq + "\n")
               break # We've found and saved the sequence, we can stop searching
            else: pass # Ignore unmatched sequences
# The End.
