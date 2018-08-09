
import Bio
from Bio import Phylo
from Bio.Phylo import BaseTree
import numpy as np
import argparse
import sys
import itertools

MAX_ITERATIONS = 10
NEIGHBOURS_THRESHOLD = 100

def _load_ids(filename):
    with open(filename, 'r') as fh:
        IDs = [line.strip().split("\t") for line in fh.readlines()]
    if len(IDs[0]) == 2:
        IDs = dict(IDs)
    elif len(IDs[0]) == 1:
        IDs = list(itertools.chain(*IDs))
        IDs = dict(zip(IDs,IDs))
    return IDs

def _get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]

def _lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names

def _process_args(my_parser, my_args):
    mytree = Phylo.read(my_args.tree, "newick")
    myids = _load_ids(my_args.id)
    terminals=mytree.get_terminals()
    mytree_name_map = _lookup_by_names(mytree)
    my_query_nodes = [mytree_name_map[id] for id in myids]
    n_ids = len(myids)
    n_terminals = len(terminals)
    id_sample = list(myids.items())[0]
    if id_sample[0] == id_sample[1]:
        print("Query\tNeighbour\tDistance\tRank")
    else:
        print("Query\tOther_name\tNeighbour\tDistance\tRank")
    for q_node in my_query_nodes:
        iterations = 0
        # For each node, keep going up until there are neighbours + 1 required terminals

        # Parent of the query
        cur_parent = _get_parent(mytree, q_node)
        parent_terminals = cur_parent.get_terminals()
        n_parent_terminals = len(parent_terminals)

        # Find a parent where there are enough neighbours
        while (n_parent_terminals < my_args.neighbours + 1):

            cur_parent = _get_parent(mytree, cur_parent)
            parent_terminals = cur_parent.get_terminals()
            n_parent_terminals = len(parent_terminals)
            # if the current number of parent terminals is more than the threshold,
            # return the previous terminals
            if n_parent_terminals > NEIGHBOURS_THRESHOLD:
                parent_terminals = prev_parent_terminals
                break
            prev_parent_terminals = parent_terminals
            # If too many rounds have occurred, stop the search
            if iterations > MAX_ITERATIONS:
                break
            iterations += 1

        # Now go through the terminals and get the distance to the query node. Ignore the query node itself.
        node_distances = []
        for cur_terminal in parent_terminals:
            if cur_terminal.name != q_node.name:
                distance = mytree.distance(q_node, cur_terminal)
                node_distances.append((cur_terminal.name, distance))
        closest_n_nodes = sorted(node_distances, key = lambda x : x[1])[:my_args.neighbours]
        for rank, entry in enumerate(closest_n_nodes):
            # Get other name if any
            other_name = myids[q_node.name]
            if other_name == q_node.name:
                out_string = "{}\t{}\t{:0.3}\t{}".format(q_node.name, entry[0], entry[1], rank + 1)
            else:
                out_string = "{}\t{}\t{}\t{:0.3}\t{}".format(q_node.name,other_name, entry[0], entry[1], rank + 1)
            print(out_string)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find closest N neighbouring terminals from query nodes.')
    parser.add_argument('-t', '--tree', help='Input tree in newick format.', required=True)
    parser.add_argument('-id', help='Leaf IDs, one-per-line, to get neighbours of.', required=True)
    parser.add_argument('-neighbours', help='Number of neighbours to get for each leaf', required=False, default=1, type=int)
    # mytree="/Users/julianzaugg/Desktop/ACE/major_projects/skin_microbiome/21_gtdb/gtdb_aureus_080818/gtdb_phylogeny.wag_gamma.tree"
    # myids = "/Users/julianzaugg/Desktop/ACE/major_projects/skin_microbiome/21_gtdb/batch_gtdb_ids.txt"
    # myids = "/Users/julianzaugg/Desktop/ACE/major_projects/skin_microbiome/21_gtdb/gtdb_and_sample_ids.tsv"
    #
    # myargs = ["-t", mytree,
    #         "-id", myids,
    #           "-n", "20"]
    #args = parser.parse_args(myargs)

    args = parser.parse_args()

    _process_args(parser, args)