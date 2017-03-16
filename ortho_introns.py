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

## INSERT INTRONS INTO ALIGNMENT
def insert_introns(fasta_dict, DB_dict):
    intron_alignment_dict = {}
    for key in fasta_dict: # each key represents one species
        transcript_id = key
        species_prefix = key.split(".")[0] # all sequence headers are prefixed with 5-letter species prefix
        species_database = DB_dict[species_prefix]
        with open(species_database + ".CDS_positions.txt", 'r') as CDS_positions: # parse the CDS positions
            for line in CDS_positions:
                if line.startswith(transcript_id):
                    coord_list = line.rstrip("\n").split("\t")[1:] # we end up with a list of CDS lengths
        aligned_sequence = fasta_dict[transcript_id]
        base_number = 0
        intron_number = 0
        intron_containing_sequence = '' # this will store our new sequence
        for base in aligned_sequence:
            if not base == "-": # don't count dashes
                base_number += 1
                if base_number == int(coord_list[intron_number]):
                    intron_containing_sequence += base
                    intron_containing_sequence += "~"
                    intron_number += 1
                    base_number = 0 # this is because we are basing it off of CDS lengths, not intron positions
                else:
                    intron_containing_sequence += base
            else:
                intron_containing_sequence += base #this will be a -. but we want this in final dict
        intron_alignment_dict[key] = intron_containing_sequence
    return intron_alignment_dict

def print_intron_positions(intron_alignment_dict):
    for key in intron_alignment_dict:
        species_prefix = key.split(".")[0]
        aligned_sequence = intron_alignment_dict[key]
        alignment_position = 1
        intron_positions = []
        for position in aligned_sequence:
                if position == "~":
                        intron_positions.append(str(alignment_position))
                else:
                        alignment_position += 1
        print species_prefix + "\t" + "\t".join(intron_positions)

def print_intron_diagram(intron_alignment_dict):
    all_intron_positions = []
    for key in intron_alignment_dict:
        aligned_sequence = intron_alignment_dict[key]
        alignment_position = 1
        for position in aligned_sequence:
                if position == "~":
                        if not alignment_position in all_intron_positions:
                                all_intron_positions.append(str(alignment_position))
                else:
                        alignment_position += 1
    for key in intron_alignment_dict:
        species_prefix = key.split(".")[0]
        aligned_sequence = intron_alignment_dict[key]
        alignment_position = 1
        intron_positions = []
        for position in aligned_sequence:
                if position == "~":
                        intron_positions.append(str(alignment_position))
                else:
                        alignment_position += 1
        binary = ""
        for i in sorted(set(all_intron_positions)):
                if i in intron_positions:
                        binary += "|"
                else:
                        binary += "-"
        print species_prefix + "\t" + binary

# Run the code
def run_software():
   with open(CDS_alignment_file, 'r') as CDS_alignment:
       with open(species_to_DB_file, 'r') as species_to_DB:
           fasta_dict = parse_alignment(CDS_alignment)
           DB_dict = parse_species_DB(species_to_DB)
           intron_alignment_dict = insert_introns(fasta_dict, DB_dict)
           print_intron_positions(intron_alignment_dict)
           print_intron_diagram(intron_alignment_dict)

if __name__ == "__main__":

   SCRIPT = basename(__file__)

   infile, value, outfile = '', 0,''
   try:
       CDS_alignment_file = sys.argv[1]
       species_to_DB_file = sys.argv[2]
   except:
       sys.exit("USAGE: ./%s %s %s" % (SCRIPT, "cds_alignment", "prefix_to_DB_file"))

   run_software()
