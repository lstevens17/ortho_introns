#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Usage: ortho_intron.py --ortho <ortho_file> --cds <cds_dir> --proteins <protein_dir> --introns <intron_dir> --tree <newick_file> --config <config_file> --trim <trim_threshold> [-h|--help]

	Options:
		-h, --help	show this help message and exit

	Input files:
		--ortho		specify Orthogroup file eg. from OrthoFinder
		--cds		specify directory containing CDS FASTA files
		--protein	specify directory containing protein FASTA files
		--intron	specify directory containing intron position files
		--tree		specify newick tree
		--config	specify config file
		--trim		specify trim threshold (int)
"""
# remote test
#############
## MODULES ##
#############

from __future__ import division
from os.path import basename, isfile, abspath, join
from datetime import datetime
import sys
import os
import subprocess

try:
        import ete3
except ImportError:
        sys.exit("[ERROR] Module \'ete3\' is not available.")

try:
        from ete3 import Tree
except ImportError:
        sys.exit("[ERROR] Module \'ete3\' is not available.")

try:
	from docopt import docopt
except ImportError:
	sys.exit("[ERROR] Module \'docopt\' not available. ")

DEVNULL = open(os.devnull, 'w')

try:
	subprocess.check_call(['mafft-linsi', '-h'], stdout=DEVNULL, stderr=subprocess.STDOUT)
except subprocess.CalledProcessError:
	pass
except OSError:
	sys.exit("[ERROR] Please add \'MAFFT\' to path.")

try:
        subprocess.check_call(['tranalign', '-h'], stdout=DEVNULL, stderr=subprocess.STDOUT)
except subprocess.CalledProcessError:
        pass
except OSError:
        sys.exit("[ERROR] Please add \'tranalign\' to path.")

###############
## FUNCTIONS ##
###############

def parse_config(config_file):
	with open(config_file, 'r') as config:
		config_dict = {}
		for line in config:
			prefix = line.rstrip("\n").split("\t")[0]
			cdsfile = line.rstrip("\n").split("\t")[1]
			proteinfile = line.rstrip("\n").split("\t")[2]
			intronfile = line.rstrip("\n").split("\t")[3]
			config_dict[prefix] = [cdsfile, proteinfile, intronfile]
		print "\t " + str(len(config_dict)) + " species found in config file."
		return config_dict

def collect_single_copy_orthogroups(orthofile, config_dict):
	taxa_list = []
	for key in config_dict:
		taxa_list.append(key)
        with open(orthofile, 'r') as orthogroups:
                print "\t Retreiving single-copy orthogroups from " + orthofile + "..."
                single_copy_dict = {}
                for line in orthogroups:
                        present_taxa_list = []
                        cluster_id = line.rstrip("\n").split(": ")[0]
                        seq_list = line.rstrip("\n").split(": ")[1].split(" ")
                        for seq in seq_list:
                                present_taxa_list.append(seq.split(".")[0])
                                if not seq.split(".")[0] in taxa_list:
                                        sys.exit("[ERROR] Found species prefix in orthogroups file which is not in config file.")
                        if len(present_taxa_list) == len(taxa_list) and len(set(present_taxa_list)) == len(taxa_list):
                                single_copy_dict[cluster_id] = seq_list
                print "\t Indentified " + str(len(single_copy_dict)) + " single copy orthogroups."
                return single_copy_dict

def collect_single_copy_CDS_sequences(single_copy_dict, config_dict, cds_dir):
	print "\t Retrieving single-copy CDS sequences..."
	cds_sequence_dict = {}
	for species in config_dict:
		single_copy_seqs = []
		for key in single_copy_dict:
			for seq in single_copy_dict[key]:
				if seq.startswith(species):
					single_copy_seqs.append(seq)
		cds_file_name = config_dict[species][0]
		with open(cds_dir + cds_file_name) as cds_fasta:
			flag = 0
                	sequence = ''
			for line in cds_fasta:
                		if line.startswith(">") and line.rstrip("\n").replace(">", "") in single_copy_seqs:
					if flag == 1:
						cds_sequence_dict[header] = sequence
						sequence = ''
						header = line.rstrip("\n").replace(">", "")
					else:
                        			flag = 1
						header = line.rstrip("\n").replace(">", "")
                		elif line.startswith(">"):
                        		if flag == 1:
                                		cds_sequence_dict[header] = sequence
						sequence = ''
                                		flag = 0
                		else:
                        		if flag == 1:
                                		sequence += line.rstrip("\n")
			if flag == 1:
				cds_sequence_dict[header] = sequence
	if len(cds_sequence_dict) == (len(single_copy_dict) * len(config_dict)):
		print "\t All single-copy CDS sequences successfully retrieved."
	else:
		sys.exit("[ERROR] Did not find all single-copy CDS sequences.")
	return cds_sequence_dict

def collect_single_copy_protein_sequences(single_copy_dict, config_dict, protein_dir):
        print "\t Retrieving single-copy protein sequences..."
        protein_sequence_dict = {}
        for species in config_dict:
                single_copy_seqs = []
                for key in single_copy_dict:
                        for seq in single_copy_dict[key]:
                                if seq.startswith(species):
                                        single_copy_seqs.append(seq)
                protein_file_name = config_dict[species][1]
                with open(protein_dir + protein_file_name) as protein_fasta:
                        flag = 0
                        sequence = ''
                        for line in protein_fasta:
                                if line.startswith(">") and line.rstrip("\n").replace(">", "") in single_copy_seqs:
                                        if flag == 1:
                                                protein_sequence_dict[header] = sequence
						sequence = ''
                                                header = line.rstrip("\n").replace(">", "")
                                        else:
                                                flag = 1
                                                header = line.rstrip("\n").replace(">", "")
                                elif line.startswith(">"):
                                        if flag == 1:
                                                protein_sequence_dict[header] = sequence
						sequence = ''
                                                flag = 0
                                else:
                                        if flag == 1:
                                                sequence += line.rstrip("\n")
                        if flag == 1:
                                protein_sequence_dict[header] = sequence
        if len(protein_sequence_dict) == (len(single_copy_dict) * len(config_dict)):
                print "\t All single-copy protein sequences successfully retrieved."
	else:
		sys.exit("[ERROR] Did not find all single-copy protein sequences.")
        return protein_sequence_dict

def collect_single_copy_intron_information(single_copy_dict, config_dict, intron_dir):
	print "\t Retrieving intron position information for all single-copy sequences..."
        intron_position_dict = {}
        for species in config_dict:
                single_copy_seqs = []
                for key in single_copy_dict:
                        for seq in single_copy_dict[key]:
                                if seq.startswith(species):
                                        single_copy_seqs.append(seq)
                intron_file_name = config_dict[species][2]
                with open(intron_dir + intron_file_name) as intron_info:
			for line in intron_info:
				seq_id = line.rstrip("\n").split("\t")[0]
				intron_positions = line.rstrip("\n").split("\t")[1:]
				if seq_id in single_copy_seqs:
					intron_position_dict[seq_id] = intron_positions
	if len(intron_position_dict) == (len(single_copy_dict) * len(config_dict)):
		print "\t Intron positions successfully retrieved."
	else:
		sys.exit("[ERROR] Did not find intron positions for all single-copy genes.")
	return intron_position_dict

def write_protein_sequences(protein_sequence_dict, single_copy_dict):
	print "\t Writing protein sequences out to \'alignments\' dir..."
	aln_dir = "alignments"
	if not os.path.exists(aln_dir):
    		os.makedirs(aln_dir)
	for orthogroup in single_copy_dict:
		if not os.path.exists(aln_dir + "/" + orthogroup + ".proteins.fa"):
			with open(aln_dir + "/" + orthogroup + ".proteins.fa", 'w') as outfile:
				seq_list = single_copy_dict[orthogroup]
				for seq in seq_list:
					outfile.write(">" + seq + "\n" + protein_sequence_dict[seq] + "\n")

def align_protein_sequences(single_copy_dict):
	print "\t Aligning proteins with MAFFT..."
	aln_dir = "alignments"
	for orthogroup in single_copy_dict:
		if not os.path.exists(aln_dir + "/" + orthogroup + ".proteins.fa.aln"):
			subprocess.call("mafft-linsi " + aln_dir + "/" + orthogroup + ".proteins.fa > " + aln_dir + "/" + orthogroup + ".proteins.fa.aln", shell=True, stderr=DEVNULL)

def write_cds_sequences(cds_sequence_dict, single_copy_dict):
        print "\t Writing cds sequences out to \'alignments\' dir..."
        aln_dir = "alignments"
        if not os.path.exists(aln_dir):
                os.makedirs(aln_dir)
        for orthogroup in single_copy_dict:
                if not os.path.exists(aln_dir + "/" + orthogroup + ".cds.fa"):
                        with open(aln_dir + "/" + orthogroup + ".cds.fa", 'w') as outfile:
                                seq_list = single_copy_dict[orthogroup]
                                for seq in seq_list:
                                        outfile.write(">" + seq + "\n" + cds_sequence_dict[seq] + "\n")

def convert_alignments(single_copy_dict):
	print "\t Converting alignments with \'tranalign\'..."
	aln_dir = "alignments"
	for orthogroup in single_copy_dict:
                if not os.path.exists(aln_dir + "/" + orthogroup + ".cds.fa.aln"):
			tranalign_output = ''
                        subprocess.call("tranalign -asequence " + aln_dir + "/" + orthogroup + ".cds.fa -bsequence " + aln_dir + "/" + orthogroup + ".proteins.fa.aln -outseq " + aln_dir + "/" + orthogroup + ".cds.fa.aln", shell=True, stderr=DEVNULL)

def check_alignments(single_copy_dict):
	print "\t Checking alignments contain expected number of sequences..."
	aln_dir = "alignments"
	bad_orthogroups = []
	for orthogroup in single_copy_dict:
		alignseq_count = 0
		orthoseq_count = len(single_copy_dict[orthogroup])
		with open(aln_dir + "/" + orthogroup + ".cds.fa.aln", 'r') as cds_alignment:
			for line in cds_alignment:
				if line.startswith(">"):
					alignseq_count += 1
		if not alignseq_count == orthoseq_count:
			bad_orthogroups.append(orthogroup)
	for bad_orthogroup in bad_orthogroups:
		del single_copy_dict[bad_orthogroup]
	print "\t " + str(len(bad_orthogroups)) + " were excluded from downstream analysis..."
	return single_copy_dict

def parse_cds_alignments(single_copy_dict):
	print "\t Reading CDS alignments into dict..."
	aln_dir = "alignments"
	cds_alignment_dict = {}
	for orthogroup in single_copy_dict:
		with open(aln_dir + "/" + orthogroup + ".cds.fa.aln", 'r') as cds_alignment:
			line_number = 0
			sequence_list = []
			for line in cds_alignment:
				line_number += 1
                       		if line_number == 1:
                                	header = line.rstrip("\n").replace(">", "")
                        	else:
                                	if line.startswith(">"):
                                        	cds_alignment_dict[header] = "".join(sequence_list)
                                        	sequence_list = []
                                        	header = line.rstrip("\n").replace(">", "")
                                	else:
                                        	sequence_list.append(line.rstrip("\n"))
                	cds_alignment_dict[header] = "".join(sequence_list)
	return cds_alignment_dict

def determine_alignment_position_coverage(single_copy_dict, cds_alignment_dict, config_dict):
	alignment_coverage_dict = {}
	alignment_length_dict = {}
	taxa_count = len(config_dict)
	for orthogroup in single_copy_dict:
		orthoseq_list = single_copy_dict[orthogroup]
		total_length = 0
		for orthoseq in orthoseq_list:
			aligned_sequence = cds_alignment_dict[orthoseq]
			total_length += len(aligned_sequence)
		if total_length % taxa_count == 0:
			alignment_length = int(total_length / taxa_count)
			alignment_length_dict[orthogroup] = alignment_length
		else:
			sys.exit("[ERROR] Sequences in alignment are not all the same length")
		for position in range(0, alignment_length):
			species_count = 0
			for orthoseq in orthoseq_list:
				aligned_sequence = cds_alignment_dict[orthoseq]
				if not aligned_sequence[position] == "-":
					species_count += 1
			alignment_coverage_dict[orthogroup + "_" + str(position)] = species_count
	return alignment_coverage_dict, alignment_length_dict

def insert_introns(single_copy_dict, cds_alignment_dict, intron_position_dict):
        print "\t Inserting introns in alignments..."
	intron_alignment_dict = {}
	for cds_id in cds_alignment_dict:
		aligned_cds = cds_alignment_dict[cds_id]
		base_number = 0
		intron_number = 0
		intron_containing_sequence = ''
		coordinate_list = intron_position_dict[cds_id]
		intron_count = len(coordinate_list)
		for position in aligned_cds:
			if not position == "-":
				base_number += 1
				if intron_number + 1 <= intron_count: # +1 as intron number starts at 0, not 1
                                        if base_number == int(coordinate_list[intron_number]): # if current base number is equal to an intron position, then insert intron
                                                intron_containing_sequence += position # add the base to the output sequence
                                                intron_containing_sequence += "~" # add tilda to represent the intron
                                                intron_number += 1 # look for the next intron
                                        else:
                                                intron_containing_sequence += position # continue as normal
				else:
					intron_containing_sequence += position
			else:
				intron_containing_sequence += position
                intron_alignment_dict[cds_id] = intron_containing_sequence # store the sequence with inserted introns in dict
	return intron_alignment_dict

def check_intron_alignments(single_copy_dict, intron_alignment_dict):
	for orthogroup in single_copy_dict:
		orthoseq_list = single_copy_dict[orthogroup]
		total_length = 0
		for orthoseq in orthoseq_list:
			intron_alignment = intron_alignment_dict[orthoseq]
			alignment_length = 0
			for base in intron_alignment:
				if not base == "~":
					alignment_length += 1
			total_length += alignment_length
		taxa_count = len(orthoseq_list)
		if not total_length % taxa_count == 0:
			sys.exit("[ERROR] Alignments are not all the same length")

def convert_original_introns_to_binary(single_copy_dict, intron_alignment_dict):
        original_binary_intron_dict = {}
        for orthogroup in single_copy_dict:
                orthoseq_list = single_copy_dict[orthogroup]
                all_intron_sites = []
                for orthoseq in orthoseq_list:
                        aligned_sequence = intron_alignment_dict[orthoseq]
                        alignment_position = 0
                        for position in aligned_sequence:
                                if position == "~":
                                        all_intron_sites.append(alignment_position)
                                else:
                                        alignment_position += 1
                for orthoseq in orthoseq_list:
                        aligned_sequence = intron_alignment_dict[orthoseq]
                        alignment_position = 0
                        present_introns = []
                        for position in aligned_sequence:
                                if position == "~":
                                        present_introns.append(alignment_position)
                                else:
                                        alignment_position += 1
                        binary = ''
                        for site in sorted(set(all_intron_sites)):
                                if site in present_introns:
                                        binary += "1"
                                else:
                                        binary += "0"
                        original_binary_intron_dict[orthoseq] = binary
        return original_binary_intron_dict


def trim_intron_alignments(single_copy_dict, intron_alignment_dict, alignment_coverage_dict, alignment_length_dict, trim_threshold):
	print "\t Trimming start and end of alignment where there are < " + str(trim_threshold) + " aligned species..."
	trimmed_intron_alignment_dict = {} # output dict
	bad_orthogroup_list = []
	count_of_trimmed_introns = 0
	trimmed_intron_site_dict = {}
	for orthogroup in single_copy_dict:
		orthoseq_list = single_copy_dict[orthogroup]
		alignment_length = alignment_length_dict[orthogroup] # this is alignment length without introns!!!!
		# calculate number of bases to trim from start
		stop_flag = 0
		position = 0
		trim_start = 0
		complete_trim = 0
		while stop_flag == 0:
			if not position == alignment_length - 1: # but we are looping through a dict which has no idea of introns
				if alignment_coverage_dict[orthogroup + "_" + str(position)] < trim_threshold: 
					trim_start = position
				else:
					stop_flag = 1
			else:
				complete_trim = 1
				stop_flag = 1
			position += 1
		# calculate number of bases to trim from end
		stop_flag = 0
		position = alignment_length - 1 # again, fine, because all we're doing is collecting how many bases to trim
		trim_end = position
		if not complete_trim == 1:
			while stop_flag == 0:
				if alignment_coverage_dict[orthogroup + "_" + str(position)] < trim_threshold:
					trim_end = position
				else:
					stop_flag = 1
				position = position - 1
		if complete_trim == 1:
	                bad_orthogroup_list.append(orthogroup)
		# start trimming
		else:
			trimmed_intron_list = []
			for orthoseq in orthoseq_list: # for every sequence in orthogroup
				aligned_sequence = intron_alignment_dict[orthoseq] # get its sequence
				if not trim_start == 0: # assuming there are bases to be trimmed
					bases_to_remove = trim_start + 1 # calculate number of bases to be trimmed
					bases_trimmed = 0 # initiate trimmed base count
					position = 0
					while bases_trimmed < bases_to_remove: # until you've trimmed all necessary sequence
						if aligned_sequence[0] == "~": # if the first base in alignment is not an intron
							aligned_sequence = aligned_sequence[1:] # get rid of it
							trimmed_intron_list.append(position)
						else: # otherwise
							aligned_sequence = aligned_sequence[1:] # get rid of the intron, but don't consider it as a removed base
							bases_trimmed += 1
							position += 1
					start_trimmed_sequence = aligned_sequence # output trimmed sequence
				else: # otherwise
					start_trimmed_sequence = aligned_sequence # output aligned sequence as if it had been trimmed
				if not trim_end == alignment_length - 1: # if we also need to remove bases from the end
					bases_to_remove = alignment_length - trim_end # calculate number of bases to be trimmed
					bases_trimmed = 0
					while bases_trimmed < bases_to_remove:
						if start_trimmed_sequence[-1:] == "~": # if the last base in start_trimmed_sequence is not an intron
							start_trimmed_sequence = start_trimmed_sequence[:-1]
						else:
							start_trimmed_sequence = start_trimmed_sequence[:-1]
							bases_trimmed += 1
					end_trimmed_sequence = start_trimmed_sequence
				else:
					end_trimmed_sequence = start_trimmed_sequence
				trimmed_intron_alignment_dict[orthoseq] = end_trimmed_sequence
			if len(trimmed_intron_list) - len(sorted(set(trimmed_intron_list))) > 0:
				count_of_trimmed_introns += 1
			trimmed_intron_site_count = len(sorted(set(trimmed_intron_list)))
			trimmed_intron_site_dict[orthogroup] = trimmed_intron_site_count
	for orthogroup in bad_orthogroup_list:
		del single_copy_dict[orthogroup]
	print "\t Trimmed orthologous introns from " + str(count_of_trimmed_introns) + " orthogroups."
	print "\t Excluded " + str(len(bad_orthogroup_list)) + " orthogroups which had no aligned sequence which exceed threshold."
	return trimmed_intron_alignment_dict, trimmed_intron_site_dict

def calculate_trimmed_alignment_length(single_copy_dict, trimmed_intron_alignment_dict):
	trimmed_intron_alignment_length_dict = {}
	for orthogroup in single_copy_dict:
		aln_len_list = []
		orthoseq_list = single_copy_dict[orthogroup]
		taxa_count = len(orthoseq_list)
		total_length = 0
		for orthoseq in orthoseq_list:
			trimmed_intron_alignment = trimmed_intron_alignment_dict[orthoseq]
			alignment_length = 0
			for position in trimmed_intron_alignment:
				if not position == "~":
					alignment_length += 1
			total_length += alignment_length
			aln_len_list.append(alignment_length)
		if total_length % taxa_count == 0:
			trimmed_intron_alignment_length_dict[orthogroup] = int(total_length / taxa_count)
		else:
			print "[ERROR] Trimmed intron alignments are not all same length..."
	return trimmed_intron_alignment_length_dict

def convert_introns_to_binary(single_copy_dict, trimmed_intron_alignment_dict):
        print "\t Converting intron positions to binary format.."
        binary_intron_dict = {}
	orthogroup_intron_positions = {}
        for orthogroup in single_copy_dict:
                orthoseq_list = single_copy_dict[orthogroup]
		all_intron_sites = []
                for orthoseq in orthoseq_list:
			aligned_sequence = trimmed_intron_alignment_dict[orthoseq]
			alignment_position = 0
			for position in aligned_sequence:
				if position == "~":
					all_intron_sites.append(alignment_position)
				else:
					alignment_position += 1
		csv = ''
		for site in sorted(set(all_intron_sites)):
			csv += str(site) + ","
		orthogroup_intron_positions[orthogroup] = csv
		for orthoseq in orthoseq_list:
			aligned_sequence = trimmed_intron_alignment_dict[orthoseq]
			alignment_position = 0
			present_introns = []
			for position in aligned_sequence:
				if position == "~":
					present_introns.append(alignment_position)
				else:
					alignment_position += 1
			binary = ''
			for site in sorted(set(all_intron_sites)):
				if site in present_introns:
					binary += "1"
				else:
					binary += "0"
			binary_intron_dict[orthoseq] = binary
        return binary_intron_dict,orthogroup_intron_positions

def parse_tree(tree_file):
        print "\t Parsing newick tree file..."
	with open(tree_file, 'r') as newick_tree:
                newick = ''
                for line in newick_tree:
                        newick += line.rstrip("\n")
                tree = Tree(newick)
		node_number = 1
		for node in tree.traverse("preorder"):
			name = node.name
			if len(name) == 0:
				node.add_features(name=str(node_number))
				node_number += 1
                return tree

def define_gain_losses(binary_intron_dict, tree, orthogroup_intron_positions, trimmed_intron_alignment_length_dict, single_copy_dict, trimmed_intron_site_dict, original_binary_intron_dict):
	out_dir = "ortho_intron_output"
    	if not os.path.exists(out_dir):
    		os.makedirs(out_dir)
	print "\t Defining gain and loss events on phylogeny..."
	#fix this
	outfile = open("relative_intron_positions.txt", 'w')
	outfile2 = open("conserved_intron_positions.txt", 'w')
	outfile3 = open("unique_gain_positions.txt", 'w')
	outfile4 = open("lost_intron_ids.txt", "w")
	with open(out_dir + "/" + "lost_event_details.txt", 'w') as lost_event_details:
		for orthogroup in single_copy_dict:
			orthoseq_list = single_copy_dict[orthogroup]
			total_length = 0
			for orthoseq in orthoseq_list:
				binary_length = len(binary_intron_dict[orthoseq])
				total_length += binary_length
			if not binary_length == 0:
				if total_length % binary_length == 0:
					intron_site_count = binary_length
				else:
					sys.exit("[ERROR] Binary sequences are not all equal in length.")
				intron_site_count = intron_site_count - 1
				intron_site = 0
				while intron_site <= intron_site_count:
					species_with_intron = []
					species_without_intron = []
					for orthoseq in orthoseq_list:
						species_prefix = orthoseq.split(".")[0]
						if int(binary_intron_dict[orthoseq][intron_site]) == 1:
							species_with_intron.append(species_prefix)
						elif int(binary_intron_dict[orthoseq][intron_site]) == 0: # and those with an intron absent
	                                		species_without_intron.append(species_prefix)
						else:
							sys.exit("[ERROR] Something other than a 1 and a 0 in the binary - WHATTTT?")
					birth_node = tree.get_common_ancestor(species_with_intron)
					if birth_node.is_root() and len(species_with_intron) > 1:
						loss_count = 0
						try:
							existing_gain_count = birth_node.gain
							updated_gain_count = existing_gain_count + 1
							birth_node.add_features(gain=updated_gain_count)
						except AttributeError:
							birth_node.add_features(gain=1)
						for node in birth_node.iter_descendants("preorder"):
							if set(node.get_leaf_names()).issubset(set(species_without_intron)):
								try:
									existing_loss_count = node.loss
									updated_loss_count = existing_loss_count + 1
									node.add_features(loss=updated_loss_count)
								except AttributeError:
									node.add_features(loss=1)
								### HERE orthogroup_intron_positions, alignment_length_dict
								intron_position = int(orthogroup_intron_positions[orthogroup].split(",")[intron_site])
								alignment_length = trimmed_intron_alignment_length_dict[orthogroup]
								outfile.write(orthogroup + "\t" + str(intron_position/alignment_length) + "\n")
								# HERE
								loss_count += 1
								for leaf in node.get_leaf_names():
									species_without_intron.remove(leaf)
						if loss_count == 0:
							intron_position = int(orthogroup_intron_positions[orthogroup].split(",")[intron_site])
							alignment_length = trimmed_intron_alignment_length_dict[orthogroup]
							outfile2.write(orthogroup + "\t" + str(intron_position/alignment_length) + "\n")
					elif birth_node.is_root() and len(species_with_intron) == 1:
						for node in tree.traverse("preorder"):
							descendants = node.get_leaf_names()
							if len(descendants) == 1 and descendants[0] == species_with_intron[0]:
								try:
									existing_gain_count = node.gain
	                                        			updated_gain_count = existing_gain_count + 1
	                                        			node.add_features(gain=updated_gain_count)
	                                			except AttributeError:
	                                        			node.add_features(gain=1)
						intron_position = int(orthogroup_intron_positions[orthogroup].split(",")[intron_site])
						alignment_length = trimmed_intron_alignment_length_dict[orthogroup]
						outfile3.write(orthogroup + "\t" + str(intron_position/alignment_length) + "\n")
					elif not birth_node.is_root() and len(species_with_intron) > 1:
						loss_count = 0
						try:
							existing_gain_count = birth_node.gain
							updated_gain_count = existing_gain_count + 1
							birth_node.add_features(gain=updated_gain_count)
						except AttributeError:
							birth_node.add_features(gain=1)
						for node in birth_node.iter_descendants("preorder"):
							if set(node.get_leaf_names()).issubset(set(species_without_intron)):
								try:
	                                                        	existing_loss_count = node.loss
	                                                        	updated_loss_count = existing_loss_count + 1
	                                                        	node.add_features(loss=updated_loss_count)
	                                                	except AttributeError:
	                                                        	node.add_features(loss=1)
	                                                        intron_position = int(orthogroup_intron_positions[orthogroup].split(",")[intron_site])
	                                                        alignment_length = alignment_length_dict[orthogroup]
	                                                        outfile.write(orthogroup + "\t" + str(intron_position/alignment_length) + "\n")
								loss_event_species_list = []
	                                                	for leaf in node.get_leaf_names():
									loss_event_species_list.append(leaf)
	                                                        	species_without_intron.remove(leaf)
								loss_count += 1
								offset = trimmed_intron_site_dict[orthogroup]
								actual_intron_site = intron_site + offset
								species_with_lost_intron = []
								for orthoseq in orthoseq_list:
									original_binary = original_binary_intron_dict[orthoseq]
									intron_number = 1
									binary_site = 0
									for site in original_binary:
										if site == "1":
											if binary_site == actual_intron_site:
												species_with_lost_intron.append(orthoseq + ".i" + str(intron_number))
											intron_number += 1
										binary_site += 1
								loss_event_details.write(orthogroup + "\t" + str(intron_site) + "\t" + str(intron_position/alignment_length) + "\t" + ",".join(loss_event_species_list) + "\t" + ",".join(species_with_lost_intron) + "\n")
						if loss_count == 0:
							intron_position = int(orthogroup_intron_positions[orthogroup].split(",")[intron_site])
	                                                alignment_length = trimmed_intron_alignment_length_dict[orthogroup]
	                                                outfile2.write(orthogroup + "\t" + str(intron_position/alignment_length) + "\n")
					else:
						print "YOU BROKE ME!"
					intron_site += 1
			else:
				print orthogroup
		outfile.close()
		outfile2.close()
		print tree.get_ascii(attributes=["name", "gain", "loss"], show_internal=True)
##########
## MAIN ##
##########
args = docopt(__doc__)
orthofile = args['<ortho_file>']
cds_dir = args['<cds_dir>']
protein_dir = args['<protein_dir>']
intron_dir = args['<intron_dir>']
tree_file = args['<newick_file>']
config_file = args['<config_file>']
trim_threshold = int(args['<trim_threshold>'])

print
print "OrthoIntron.py"
print
print "[STEP 1] COLLECT SEQUENCE AND INTRON INFORMATION"
config_dict = parse_config(config_file)
single_copy_dict = collect_single_copy_orthogroups(orthofile, config_dict)
cds_sequence_dict = collect_single_copy_CDS_sequences(single_copy_dict, config_dict, cds_dir)
protein_sequence_dict = collect_single_copy_protein_sequences(single_copy_dict, config_dict, protein_dir)
intron_position_dict = collect_single_copy_intron_information(single_copy_dict, config_dict, intron_dir)
print "[STEP 2] DEFINE ORTHOLOGOUS INTRONS"
write_protein_sequences(protein_sequence_dict, single_copy_dict)
align_protein_sequences(single_copy_dict)
write_cds_sequences(cds_sequence_dict, single_copy_dict)
convert_alignments(single_copy_dict)
single_copy_dict = check_alignments(single_copy_dict)
cds_alignment_dict = parse_cds_alignments(single_copy_dict)
tuple = determine_alignment_position_coverage(single_copy_dict, cds_alignment_dict, config_dict)
alignment_coverage_dict = tuple[0]
alignment_length_dict = tuple[1]
intron_alignment_dict = insert_introns(single_copy_dict, cds_alignment_dict, intron_position_dict)
check_intron_alignments(single_copy_dict, intron_alignment_dict)

original_binary_intron_dict = convert_original_introns_to_binary(single_copy_dict, intron_alignment_dict)

trimming_output = trim_intron_alignments(single_copy_dict, intron_alignment_dict, alignment_coverage_dict, alignment_length_dict, trim_threshold)
trimmed_intron_alignment_dict = trimming_output[0]
trimmed_intron_site_dict = trimming_output[1]

trimmed_intron_alignment_length_dict = calculate_trimmed_alignment_length(single_copy_dict, trimmed_intron_alignment_dict)
binary_output = convert_introns_to_binary(single_copy_dict, trimmed_intron_alignment_dict)
binary_intron_dict = binary_output[0]
orthogroup_intron_positions = binary_output[1]
print "[STEP 3] DEFINE GAIN AND LOSS EVENTS ON TREE"
tree = parse_tree(tree_file)
define_gain_losses(binary_intron_dict, tree, orthogroup_intron_positions, trimmed_intron_alignment_length_dict, single_copy_dict, trimmed_intron_site_dict, original_binary_intron_dict)
print "[STATUS] OrthoIntron.py complete."
