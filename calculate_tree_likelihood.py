from bio import Phylo # this is a library that I told before
import pandas as pd
from math import exp


class Likelihood:
    def __init__(self, dated_reconciled_tree_path, alignment_path):
        # initial frequency of amino acid
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

        # used to establish an example pam matrix
        self.pam_matrix = {'A': {}, 'R': {}, 'N': {}, 'D': {}, 'C': {}, 'Q': {}, 'E': {}, 'G': {}, 'H': {}, 'I': {},
                  'L': {}, 'K': {}, 'M': {}, 'F': {}, 'P': {}, 'S': {}, 'T': {}, 'W': {}, 'Y': {}, 'V': {}}
        df = pd.read_csv('pam_matrix.csv')
        for index, row in df.iterrows():
            d = row.to_dict()
            a_acid = d['Unnamed: 0']
            for amino_acid in self.amino_acid_list:
                self.pam_matrix[amino_acid][a_acid] = d[amino_acid]

        self.tree = Phylo.read(dated_reconciled_tree_path, "newick")
        self.alignment = {} # contain all alignments of terminal species.
        alignment_file = open(alignment_path, 'r')
        alignment_lines = alignment_file.readlines()
        i = 0
        while i < len(alignment_lines):
            self.alignment[alignment_lines[i].rstrip("\n").lstrip(">")] = alignment_lines[i + 1].strip("\n")
            i += 2
        self.alignment_length = len(alignment_lines[i - 1])
        alignment_file.close()
        self.non_terminals = self.tree.get_nonterminals()
        self.terminals = self.tree.get_terminals()

    def get_base_sub_prob(self, base_i, base_j, time):  # equation 7 in Felsenstein paper
        k = 0
        if base_i == base_j:
            k = 1
        return exp(-time)*k + (1-exp(-time))*self.pam_matrix[base_i][base_j]

    def get_all_possible_sequences(self):
        # since some species have unknown alignment, or some amino acid in the species alignment is missing, this
        # function generates all possible alignments for each species.
        map_sequence = {}
        for non_terminal in self.non_terminals:
            all_sequences = []
            self.helper_get_all_non_terminal([], all_sequences)
            map_sequence[non_terminal] = all_sequences
        for terminal in self.terminals:
            alignment = self.alignment[str(terminal)]
            all_sequences = []
            self.helper_get_all_terminal([], all_sequences, alignment, 0)
            map_sequence[terminal] = all_sequences
        return map_sequence

    def helper_get_all_terminal(self, current, all_alignments, alignment, index):
        if len(current) == self.alignment_length:
            all_alignments.append(current)
        elif alignment[index] != '-':
            current.append(alignment[index])
            self.helper_get_all_terminal(current, all_alignments, alignment, index+1)
        else:
            for amino_acid in self.amino_acid_list:
                current.append(amino_acid)
                self.helper_get_all_terminal(current, all_alignments, alignment, index+1)
                current.pop()

    def helper_get_all_non_terminal(self, current, all_alignments):
        if len(current) == self.alignment_length:
            all_alignments.append(current)
        else:
            for amino_acid in self.amino_acid_list:
                current.append(amino_acid)
                self.helper_get_all_non_terminal(current, all_alignments)
                current.pop()

    def get_all_path(self): # get all possible pair of species to calculate
        # In the example of Felsenstein paper, all possible pair is (s_0, s_6), (s_6, s_1), (s_6, s_2), (s_0, s_8),
        # (s_8, s_3), (s_8, s_7), (s_7, s_4), (s_7, s_5)
        path = []
        for i in range(len(self.non_terminals)):
            for j in range(i+1, self.non_terminals):
                if self.tree.trace(self.non_terminals[i], self.non_terminals[j]) == 2:
                    if self.tree.distance(self.non_terminals[i]) < self.tree.distance(self.non_terminals[j]):
                        path.append((self.non_terminals[i], self.non_terminals[j]))
                    else:
                        path.append((self.non_terminals[j], self.non_terminals[i]))

        for non_terminal in self.non_terminals:
            for terminal in self.terminals:
                if self.tree.trace(non_terminal, terminal) == 2:
                    path.append((non_terminal, terminal))
        return path

    def get_prior_prob(self, alignment): # calculate pi s_O given the alignment of s_0
        # I just simply multiply all of the initial frequency of all amino acid in the alignment s_0 together
        result = 1
        for amino_acid in alignment:
            result *= self.initial_freq[self.amino_acid_name[amino_acid]]
        return result

    def calculate_likelihood(self): # calculate the likelihood of a tree
        path = self.get_all_path()
        result = 0
        # I am currently stuch here. Equation (1) turns out to be not easy to be implemented as I expected. I think
        # I will consider equation (3).
