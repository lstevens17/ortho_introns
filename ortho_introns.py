#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from os.path import basename, isfile, abspath, join
import sys

#########################
###### FUNCTIONS ########
#########################

## PARSE CDS ALIGNMENT
# We expect alignment in FASTA format
def parse_alignment ( CDS_alignment ):
    line_number = 0
    fasta_dict = {} # this will store our headers as keys and the sequence (with dashes) as values
    sequence_list = [] # this will store the nucleotide sequence as we parse the alignment
    for line in CDS_alignment:
        line_number += 1
        if line_number == 1:
            header = line.rstrip("\n").replace(">", "")
        else:
            if line.startswith(">"):
                fasta_dict[header] = "".join(sequence_list)
                sequence_list = []
                header = line.rstrip("\n").replace(">", "")
            else:
                sequence_list.append(line.rstrip("\n"))
    fasta_dict[header] = "".join(sequence_list)
    print fasta_dict
    return fasta_dict

## PARSE SPECIES PREFIX TO DBS
# We will read in a simple TSV with SPECIESPREFIX\tDATABASENAME and parse it into a dict
def parse_species_DB ( species_to_DB ):
    DB_dict = {}
    for line in species_to_DB:
        prefix = line.strip("\n").split("\t")[0]
        DB = line.rstrip("\n").split("\t")[1]
        DB_dict[prefix] = DB
        print DB_dict
        return DB_dict


# Run the code
def run_software():
   with open(CDS_alignment_file, 'r') as CDS_alignment:
       with open(species_to_DB_file, 'r') as species_to_DB:
           parse_alignment(CDS_alignment)
           parse_species_DB(species_to_DB)

if __name__ == "__main__":

   SCRIPT = basename(__file__)

   infile, value, outfile = '', 0,''
   try:
       CDS_alignment_file = sys.argv[1]
       species_to_DB_file = sys.argv[2]
   except:
       sys.exit("USAGE: ./%s %s %s" % (SCRIPT, "cds_alignment", "prefix_to_DB_file"))

   run_software()
