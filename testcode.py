from bio import Phylo
import pylab
import pandas as pd
from math import log10
from math import exp

class Likelihood:
    def __init__(self, dated_reconciled_tree_path, alignment_path):
        self.initial_freq = {'Gly': 0.089, 'Ala': 0.087, 'Leu': 0.085, 'Lys': 0.081, 'Ser': 0.070,
                             'Val': 0.065, 'Thr': 0.058, 'Pro': 0.051, 'Glu': 0.050, 'Asp': 0.047,
                             'Arg': 0.041, 'Asn': 0.040, 'Phe': 0.040, 'Gln': 0.038, 'Ile': 0.037,
                             'His': 0.034, 'Cys': 0.033, 'Tyr': 0.030, 'Met': 0.015, 'Trp': 0.010}
        self.amino_acid_name = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'Q': 'Gln', 'E': 'Glu',
                                'G': 'Gly', 'H': 'His', 'I': 'Ile',
                                'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr',
                                'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'}
        self.amino_acid_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T',
                                'W', 'Y', 'V']
        self.pam_matrix = {}
        self.tree = Phylo.read(dated_reconciled_tree_path, "newick")
        self.alignment = {}
        alignment_file = open(alignment_path, 'r')
        alignment_lines = alignment_file.readlines()
        i = 0
        while i < len(alignment_lines):
            self.alignment[alignment_lines[i].rstrip("\n").lstrip(">")] = alignment_lines[i + 1].strip("\n")
            i += 2
        self.alignment_length = len(alignment_lines[i - 1])
        alignment_file.close()

    def get_base_sub_prob(self, base_i, base_j, time):  # equation 7 in Felsenstein paper
        k = 0
        if base_i == base_j:
            k = 1
        return exp(-time)*k + (1-exp(-time))*self.pam_matrix[base_i][base_j]

    def get_all_possible_between_two(self, s_a, s_b, clade_a, clade_b):
        # get all possible pairs of amino acid between two sequences
        # s_a, s_b are two alignments
        # clade_a, clade_b are two sequences in the tree
        possible = []
        parent = s_a
        child = s_b
        if self.tree.distance(clade_a) > self.tree.distance(clade_b):
            parent = s_b
            child = s_a
        self.helper(parent, child, possible, 0, [])
        return possible

    def helper(self, parent, child, possible, i, current):  # helper function for the above function
        if i == self.alignment_length:
            possible.append(current)
        else:
            if parent[i] != '-' and child[i] != '-':
                current.append((parent[i], child[i]))
                self.helper(parent, child, possible, i+1, current)
                current.pop()
            elif parent[i] == '-' and child[i] != '-':
                for amino_acid in self.amino_acid_list:
                    current.append((amino_acid, child[i]))
                    self.helper(parent, child, possible, i + 1, current)
                    current.pop()
            elif parent[i] != '-' and child[i] == '-':
                for amino_acid in self.amino_acid_list:
                    current.append((parent[i], amino_acid))
                    self.helper(parent, child, possible, i + 1, current)
                    current.pop()
            else:
                for amino_acid in self.amino_acid_list:
                    for another_amino_acid in self.amino_acid_list:
                        current.append((amino_acid, another_amino_acid))
                        self.helper(parent, child, possible, i + 1, current)
                        current.pop()

    def calculate_substitution_likelihood(self):  # apply equation 1 of Felsenstein paper
        non_terminals = self.tree.get_nonterminals()
        terminals = self.tree.get_terminals()
        result = 1  # the result we want to return (the likelihood of the tree)

        for i in range(len(non_terminals)):
            for j in range(i+1, len(non_terminals)):
                # nested for loop to find all pairs of non terminal - non terminal sequence in the tree
                if self.tree.trace(non_terminals[i], non_terminals[j]) == 1:
                    time = self.tree.distance(non_terminals[i], non_terminals[j])

                    parent = child = "-" * self.alignment_length
                    possible = self.get_all_possible_between_two(parent, child, non_terminals[i], non_terminals[j])
                    total_prob = 0
                    for alignment_pair in possible:
                        alignment_prob = 1
                        for pair in alignment_pair:
                            alignment_prob *= self.get_base_sub_prob(pair[0], pair[1], time)
                        total_prob += alignment_prob
                    result *= total_prob

        for non_terminal in non_terminals:
            for terminal in terminals:
                # nested for loop to find all pairs of non terminal - terminal sequence in the tree
                if self.tree.trace(non_terminal, terminal) == 1:
                    time = self.tree.distance(non_terminal, terminal)
                    parent = "-" * self.alignment_length
                    child = self.alignment[str(non_terminal)]
                    possible = self.get_all_possible_between_two(parent, child, non_terminal, terminal)
                    total_prob = 0
                    for alignment_pair in possible:
                        alignment_prob = 1
                        for pair in alignment_pair:
                            alignment_prob *= self.get_base_sub_prob(pair[0], pair[1], time)
                        total_prob += alignment_prob
                    result *= total_prob
        return result  # this result haven't considered the value of prior probability yet. I will update later.
