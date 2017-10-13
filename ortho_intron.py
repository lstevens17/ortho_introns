#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from Bio import SeqIO, AlignIO
from os.path import basename, isfile, abspath, join
from os import listdir
import sys

try:
        from ete3 import Tree
except ImportError:
        sys.exit("[ERROR] Module \'ete3\' is not available.")

# cds_alignment_dir = directory of CDS alignments
# intron_info_dir = directory of intron information in tsv format
# alignment_dict = dictionary where keys are filenames (OGXXXX) and value are AlignIO aligment objects
# intron_length_dict = dictionary where keys are seqIDs and values are lists of intron lengths, in order
# intron_position_dict = dictionary where keys are seqIDs and values are lists of intron positions, in order
# intron_alignment_position_dict = dictonary where keys are seqIDs are values are lists of RELEVATIVE intron positions, in order
# presence_absence_dict = dictonary where keys are seqIDs and values are strings of BINARY, in order
# orthologous_length_dict = dictionary where keys are seqIDs and values are lists of lengths and dashes (showing absence), in order

######### MISSING DATA FILTER ############
min_species_with_data = int(sys.argv[4])
##########################################



def parse_cds_alignments(cds_alignment_dir):
    print "[+] Parsing CDS alignments..."
    alignment_dict = {}
    ## loop through and parse every file in cds_alignment_dir
    for file in listdir(cds_alignment_dir):
        with open(cds_alignment_dir + "/" + file, 'r') as alignment_file:
            alignment_dict[file.split(".")[0]] = AlignIO.read(alignment_file, "fasta") # ASSUMPTION: that the alignment file name can be split on the dot. Should be fine
    ## loop through every alignment now in dict and store orthogroup ID and seq IDs as a list
    orthogroup_dict = {}
    for orthogroup, alignment in alignment_dict.iteritems():
        seqid_list = []
        for seq in alignment:
            seqid_list.append(seq.id)
        orthogroup_dict[orthogroup] = seqid_list
    ## return dicts
    return alignment_dict, orthogroup_dict

def parse_intron_information(intron_info_dir):
    print "[+] Parsing intron information..."
    intron_position_dict = {}
    intron_length_dict = {}
    ## loop through every intron info file and store values in dict
    for file in listdir(intron_info_dir):
        with open(intron_info_dir + "/" + file, 'r') as intron_info_file:
                for line in intron_info_file:
                    seq_id = line.rstrip("\n").split("\t")[0]
                    int_list = line.rstrip("\n").split("\t")[1:]
                    if not seq_id in intron_position_dict:
                        intron_position_dict[seq_id] = int_list[:-1] # ASSUMPTION: the last intron position is actually the CDS length, and is therefore discarded
                    else:
                        intron_length_dict[seq_id] = int_list # ASSUMPTION: that the second line contains intron lengths
    ## return dicts
    return intron_length_dict, intron_position_dict

def determine_relevative_intron_positions(alignment_dict, intron_length_dict, intron_position_dict):
    print "[+] Determining relative introns positions..."
    intron_alignment_position_dict = {}
    orthogroup_count, orthogroup_total = 0, len(alignment_dict)
    ## loop through every alignment in dict and record relative intron position
    for OG, alignment in alignment_dict.iteritems():
        for seq in alignment:
            base_count, alignment_position, intron_count = 0, 0, 0
            for base in seq.seq:
                if not base == '-':
                    base_count += 1
                    alignment_position += 1
                    if intron_count < len(intron_position_dict[seq.id]): # to stop intron count exceeding length of intron list, replace with a while loop?
                        if base_count == int(intron_position_dict[seq.id][intron_count]):
                            try:
                                intron_alignment_position_dict[seq.id].append(str(alignment_position))
                            except KeyError:
                                intron_alignment_position_dict[seq.id] = [str(alignment_position)]
                            intron_count += 1
                else:
                    alignment_position += 1
        # progress
        orthogroup_count += 1
        sys.stdout.write('\r'+"\t"+str(orthogroup_count)+"/"+str(orthogroup_total)+" completed.")
        sys.stdout.flush()
    print
    ## return dicts
    ## plotting
    intron_count_dict = {}
    for seq, list in intron_alignment_position_dict.iteritems():
        try:
            intron_count_dict[seq.split(".")[0]] += len(list)
        except KeyError:
            intron_count_dict[seq.split(".")[0]] = len(list)
    for item, value in intron_count_dict.iteritems():
        print item + "\t" + str(value)
    return intron_alignment_position_dict

def defining_orthologous_introns(intron_alignment_position_dict, alignment_dict, intron_length_dict):
    print "[+] Defining orthologous introns..."
    presence_absence_dict = {}
    master_intron_position_dict, cds_len = {}, {}
    orthogroup_count, orthogroup_total = 0, len(alignment_dict)
    ## for each alignment, get a list of intron positions from all sequences and make a unique list
    for orthogroup, alignment in alignment_dict.iteritems():
        total_intron_list = []
        for seq in alignment:
            try:
                total_intron_list += intron_alignment_position_dict[seq.id]
            except KeyError:
                pass
        unique_intron_list = sorted(set(total_intron_list), key=int)
        ## for every sequence in the alignment
        for seq in alignment:
            presence_absence = ''
            alignment_position, seq_position = 0, 0
            ##  collect the intron positions for the current sequence
            try:
                intron_positions = intron_alignment_position_dict[seq.id]
            except KeyError:
                intron_positions = []
            ## for every position in sequence
            for position in seq.seq:
                alignment_position += 1
                if not position == '-':
                    seq_position += 1
                ## if current position is intron site
                if str(alignment_position) in unique_intron_list:
                    ## store presence, absence or missing in binary string
                    if not str(alignment_position) in intron_positions:
                        if position == '-':
                            presence_absence += 'X'
                        else:
                            presence_absence += '0'
                    else:
                        presence_absence += '1'
                    ## store position of intron (whether present, absent or missing) in dict
                    try:
                        master_intron_position_dict[seq.id].append(seq_position)
                    except KeyError:
                        master_intron_position_dict[seq.id] = [seq_position]
            ## for every sequence, store the CDS length (fix this) and binary in dicts
            cds_len[seq.id] = seq_position + 3
            #print seq.id.split(".")[0] + "\t" + presence_absence
            presence_absence_dict[seq.id] = presence_absence
        ## progress
        orthogroup_count += 1
        sys.stdout.write('\r'+"\t"+str(orthogroup_count)+"/"+str(orthogroup_total)+" completed.")
        sys.stdout.flush()
    print
    return presence_absence_dict, master_intron_position_dict, cds_len

def parse_tree(tree_file):
    print "[+] Parsing tree..."
    tree = Tree(tree_file)
    ## for every node, initiate gain and loss counts (and go back and do the same to the root)
    for node in tree.get_tree_root().iter_descendants("preorder"):
        node.add_features(gain_count=0, loss_count=0)
    tree.get_tree_root().add_features(gain_count=0, loss_count=0)
    return tree

def define_gain_and_loss_events(presence_absence_dict, orthogroup_dict, tree):
    root_of_caeno = tree.get_common_ancestor('CELEG', 'CMONO')
    print "[+] Inferring gain and loss events..."
    all_intron_list = []
    gain_record, loss_record = {}, {}
    gain_count, loss_count, specific_gain_count = 0, 0, 0
    orthogroup_count, orthogroup_total = 0, len(alignment_dict)
    ## for every orthogroup, get a count of the intron in the alignment
    for orthogroup, seqid_list in orthogroup_dict.iteritems():
        intron_count = len(presence_absence_dict[seqid_list[0]]) # ASSUMPTION: all binary strings in dict are the same - they should be!!
        ## for every intron in orthogroup
        for i in range(0, intron_count): # note: range stops 1 early, ie range(0, 21) generates numbers 0..20
            present_list, absent_list, missing_data_list = [], [], []
            ## generate lists of taxa which have an intron present, absent or missing data
            for seqid in seqid_list:
                if presence_absence_dict[seqid][i] == '1':
                    present_list.append(seqid.split(".")[0])
                elif presence_absence_dict[seqid][i] == '0':
                    absent_list.append(seqid.split(".")[0])
                else:
                    missing_data_list.append(seqid.split(".")[0])
            ## assuming the number of species without missing data is above threshold
            if len(present_list + absent_list) >= min_species_with_data:
                #if len(present_list) > 1:
                all_intron_list.append([orthogroup, i, present_list])
                ## if the intron is not unique to a single species
                if len(present_list) > 1:
                    ## record a gain
                    if tree.get_common_ancestor(present_list) == root_of_caeno:
                        gain_record[orthogroup + "_" + str(gain_count)] = [tree.get_common_ancestor(present_list).get_leaf_names(), i, "caeno_root"]
                    elif tree.get_common_ancestor(present_list).is_root() and present_list > 1:
                        gain_record[orthogroup + "_" + str(gain_count)] = [tree.get_common_ancestor(present_list).get_leaf_names(), i, "root"]
                    else:
                        gain_record[orthogroup + "_" + str(gain_count)] = [tree.get_common_ancestor(present_list).get_leaf_names(), i, "elsewhere"]
                    gain_count += 1
                    tree.get_common_ancestor(present_list).add_features(gain_count=tree.get_common_ancestor(present_list).gain_count + 1)
                    ## from gain/birth node, traverse tree hunting for loss events
                    for node in tree.get_common_ancestor(present_list).iter_descendants("preorder"):
                        if set(node.get_leaf_names()).issubset(set(absent_list + missing_data_list)) and not set(node.get_leaf_names()).issubset(set(missing_data_list)): # (if there is some missing data in there, then cool; if its all missing data, do nothing)
                            ## record loss
                            node.add_features(loss_count=node.loss_count + 1)
                            loss_record[orthogroup + "_" + str(loss_count)] = [node.get_leaf_names(), i]
                            loss_count += 1
                            ## remove species from being considered again
                            for leaf in node.get_leaf_names():
                                if leaf in absent_list:
                                    absent_list.remove(leaf)
                                else:
                                    missing_data_list.remove(leaf)
                ## else, if intron is unique to a species, record it as a gain at the leaf node
                else:
                    tree.get_leaves_by_name(present_list[0])[0].add_features(gain_count=tree.get_leaves_by_name(present_list[0])[0].gain_count + 1)
                    gain_record[orthogroup + "_" + str(gain_count)] = [tree.get_common_ancestor(present_list).get_leaf_names(), i, "tip"]
                    specific_gain_count += 1
     #               gain_record[orthogroup + "_" + str(gain_count)] = [tree.get_common_ancestor(present_list), i]
                    gain_count += 1
        ## progress
        orthogroup_count += 1
        sys.stdout.write('\r'+"\t"+str(orthogroup_count)+"/"+str(orthogroup_total)+" completed.")
        sys.stdout.flush()
    print
    print tree.get_ascii(attributes=["name", "gain_count", "loss_count"], show_internal=True)
    print "gains: ", gain_count, "species-specific gains:", specific_gain_count
    print "losses: ", loss_count
    return loss_record, gain_record, all_intron_list

def play_with_losses(loss_record, orthogroup_dict, alignment_dict):
    with open("relative_loss_positions.txt", "w") as outfile:
        for loss, loss_list in loss_record.iteritems():
            orthogroup, taxa_list, position = loss.split("_")[0], loss_list[0], loss_list[1]
            seq_list = orthogroup_dict[orthogroup]
            alignment = alignment_dict[orthogroup]
            sum_intron_position, alignment_len = 0, 0
            rel_pos_list = []
            for taxa in taxa_list:
                for seq in seq_list:
                    if seq.startswith(taxa):
                        sum_intron_position += (int(master_intron_position_dict[seq][position]) / cds_len[seq])
                        rel_pos_list.append(str(int(master_intron_position_dict[seq][position]) / cds_len[seq]))
            outfile.write(str(sum_intron_position / len(taxa_list)) + "\n")

def play_with_alls(all_intron_list, orthogroup_dict, alignment_dict):
    with open("all_intron_relative_positions.txt", "w") as outfile:
        for item in all_intron_list:
            orthogroup, position, taxa_list = item[0], item[1], item[2]
            seq_list = orthogroup_dict[orthogroup]
            alignment = alignment_dict[orthogroup]
            sum_intron_position, alignment_len = 0, 0
            rel_pos_list = []
            for taxa in taxa_list:
                for seq in seq_list:
                    if seq.startswith(taxa):
                        sum_intron_position += (int(master_intron_position_dict[seq][position]) / cds_len[seq])
                        rel_pos_list.append(str(int(master_intron_position_dict[seq][position]) / cds_len[seq]))
            outfile.write(str(sum_intron_position / len(taxa_list)) + "\n")

def play_with_gains(gain_record, orthogroup_dict, alignment_dict):
    with open("relative_gain_positions.txt", "w") as outfile:
        for gain, gain_list in gain_record.iteritems():
            orthogroup, taxa_list, position, place = gain.split("_")[0], gain_list[0], gain_list[1], gain_list[2]
            seq_list = orthogroup_dict[orthogroup]
            alignment = alignment_dict[orthogroup]
            sum_intron_position, alignment_len = 0, 0
            rel_pos_list = []
            for taxa in taxa_list:
                for seq in seq_list:
                    if seq.startswith(taxa):
                        sum_intron_position += (int(master_intron_position_dict[seq][position]) / cds_len[seq])
                        rel_pos_list.append(str(int(master_intron_position_dict[seq][position]) / cds_len[seq]))
            if place == 'elsewhere':
                outfile.write(str(sum_intron_position / len(taxa_list)) + "\n")
            elif place == 'tip':
                outfile.write(str(sum_intron_position / len(taxa_list)) + "\n")

def order_of_losses(loss_record, presence_absence_dict, orthogroup_dict):
    print "[+] Determining order of losses..."
    with open("losses_by_column.txt", 'w') as outfile:
        orthogroup_count, orthogroup_total = 0, len(alignment_dict)
        for orthogroup, seq_list in orthogroup_dict.iteritems():
            #print orthogroup
            ## first lets identify which columns we are removing (could this be done above)?
            orthogroup_binary_list = []
            for seq in seq_list:
                orthogroup_binary_list.append(presence_absence_dict[seq])
                intron_count = len(presence_absence_dict[seq])
                #print presence_absence_dict[seq]
            column_list = []
            for i in range(0, intron_count):
                existing_data_count = 0
                for binary in orthogroup_binary_list:
                    if not binary[i] == "X":
                        existing_data_count += 1
                if existing_data_count >= min_species_with_data:
                    column_list.append(i)
            actual_column = 0
            mapping_dict = {}
            for column in column_list:
                mapping_dict[column] = actual_column
                actual_column += 1
            #print mapping_dict
            losses_by_column = {}
            for loss, loss_list in loss_record.iteritems(): # SLOW!!!
                if loss.split("_")[0] == orthogroup:
                    taxa_loss, position = loss_list[0], loss_list[1]
                    actual_position = mapping_dict[position]
                    try:
                        losses_by_column[actual_position] += 1
                    except KeyError:
                        losses_by_column[actual_position] = 1
            out_list = []
            for y in range(0,len(column_list)):
                try:
                    count = losses_by_column[y]
                except KeyError:
                    count = 0
                out_list.append(str(count))
            outfile.write(orthogroup + ": " + " ".join(out_list) + "\n")
            ## progress
            orthogroup_count += 1
            sys.stdout.write('\r'+"\t"+str(orthogroup_count)+"/"+str(orthogroup_total)+" completed.")
            sys.stdout.flush()
    print

def print_supermatrix():
    with open("intron_supermatrix.fas", "w") as supermatrix_file:
        out_dict = {}
        for orthogroup, seq_list in orthogroup_dict.iteritems():
            for seq in seq_list:
                binary = presence_absence_dict[seq]
                taxa = seq.split(".")[0]
                try:
                    out_dict[taxa] += binary
                except KeyError:
                    out_dict[taxa] = binary
        for taxa in out_dict:
            supermatrix_file.write(">" + taxa + "\n")
            supermatrix_file.write(out_dict[taxa].replace("X", "?") + "\n")


if __name__ == "__main__":

    SCRIPT = basename(__file__)

    infile, value, outfile = '', 0,''
    try:
        cds_alignment_dir = sys.argv[1].rstrip("/")
        intron_info_dir = sys.argv[2].rstrip("/")
        tree_file = sys.argv[3]

    except:
        sys.exit("USAGE: ./%s %s %s %s" % (SCRIPT, "cds_alignment_dir", "intron_info_dir", "tree.nwk"))

    alignment_dict, orthogroup_dict = parse_cds_alignments(cds_alignment_dir)
    intron_length_dict, intron_position_dict = parse_intron_information(intron_info_dir)
    intron_alignment_position_dict = determine_relevative_intron_positions(alignment_dict, intron_length_dict, intron_position_dict)
    presence_absence_dict, master_intron_position_dict, cds_len = defining_orthologous_introns(intron_alignment_position_dict, alignment_dict, intron_length_dict)
    tree = parse_tree(tree_file)
    loss_record, gain_record, all_intron_list = define_gain_and_loss_events(presence_absence_dict, orthogroup_dict, tree)
    play_with_losses(loss_record, orthogroup_dict, alignment_dict)
    play_with_gains(gain_record, orthogroup_dict, alignment_dict)
    play_with_alls(all_intron_list, orthogroup_dict, alignment_dict)
    order_of_losses(loss_record, presence_absence_dict, orthogroup_dict)
    print_supermatrix()
