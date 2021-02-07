from bio import Phylo
import pylab
import pandas as pd
from math import log10


class JTTMatrix:
    def __init__ (self, dated_reconciled_tree_path, alignment_path):
        self.amino_acid_name = {'A':'Ala', 'R':'Arg', 'N':'Asn', 'D':'Asp', 'C':'Cys', 'Q':'Gln', 'E':'Glu', 'G':'Gly', 'H':'His', 'I':'Ile',
                                 'L':'Leu', 'K':'Lys', 'M':'Met', 'F':'Phe', 'P':'Pro', 'S':'Ser', 'T':'Thr', 'W':'Trp', 'Y':'Tyr', 'V':'Val'}
        self.amino_acid_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        self.tree = Phylo.read(dated_reconciled_tree_path, "newick")
        self.alignment = {}
        alignment_file = open(alignment_path, 'r')
        alignment_lines = alignment_file.readlines()
        i = 0
        while i < len(alignment_lines):
            self.alignment[alignment_lines[i].rstrip("\n").lstrip(">")] = alignment_lines[i+1].strip("\n")
            i += 2
        self.alignment_length = len(alignment_lines[i-1])
        alignment_file.close()

    def get_species_list(self):
        species_list = []
        for key in self.alignment:
            species_list.append(key)
        return species_list

    def draw_tree(self):
        Phylo.draw(self.tree)
        pylab.show()

    def calculate_frequency(self):
        count_amino_acid = {'A':0, 'R':0, 'N':0, 'D':0, 'C':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0,
                            'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0}
        total = 0
        for alignment in self.alignment.values():
            for amino_acid in alignment:
                if amino_acid != "-":
                    count_amino_acid[amino_acid] += 1
                    total += 1
        for amino_acid in count_amino_acid:
            count_amino_acid[amino_acid] /= total
        return count_amino_acid

    def calculate_local_frequency(self, sequence):
        amino_acid_frequency = {'A':0, 'R':0, 'N':0, 'D':0, 'C':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0,
                                'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0}
        # count the number of each amino acid in the sequence
        for amino_acid in sequence:
            if amino_acid != '-':
                amino_acid_frequency[amino_acid] += 1
        # calculate the frequency of each amino acid in the sequence
        for amino_acid in amino_acid_frequency:
            amino_acid_frequency[amino_acid] /= self.alignment_length
        return amino_acid_frequency

    def construct_raw_PAM_matrix_and_get_alignment_mutation(self):
        # initialize the raw PAM matrix
        full_mutated_amino_acid = {'A': {}, 'R': {}, 'N': {}, 'D': {}, 'C': {}, 'Q': {}, 'E': {}, 'G': {}, 'H': {}, 'I': {},
                    'L': {}, 'K': {}, 'M': {}, 'F': {}, 'P': {}, 'S': {}, 'T': {}, 'W': {}, 'Y': {}, 'V': {}}
        for amino_acid in full_mutated_amino_acid:
            full_mutated_amino_acid[amino_acid] = {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0,
                                    'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0}
            full_mutated_amino_acid[amino_acid][amino_acid] = -1

        # initialize the counting of the number of mutations in each sequence
        alignment_mutation = {}
        for species_key in self.get_species_list():
            alignment_mutation[species_key] = 0
        species_list = self.get_species_list()

        # construct the PAM matrix and the counting of the number of mutations in each sequence
        for species in species_list:
            for another_species in species_list:
                if len(self.tree.trace(species, another_species)) == 2:  # check if two sequences are related
                    for i in range(len(self.alignment[species])):
                        if self.alignment[species][i] != "-" and self.alignment[another_species][i] != "-" \
                                and self.alignment[species][i] != self.alignment[another_species][i]:
                            # update the PAM matrix
                            full_mutated_amino_acid[self.alignment[species][i]][self.alignment[another_species][i]] += 1
                            full_mutated_amino_acid[self.alignment[another_species][i]][self.alignment[species][i]] += 1
                            # update the counting of the number of mutations in each sequence
                            alignment_mutation[species] += 1
                            alignment_mutation[another_species] += 1
        return full_mutated_amino_acid, alignment_mutation

    def get_total_exposure_of_amino_acid_to_mutation(self):
        result = {}
        alignment_mutation = self.construct_raw_PAM_matrix_and_get_alignment_mutation()[1]
        species_list = self.get_species_list()
        local_frequency = {}
        for species in species_list:
            local_frequency[species] = self.calculate_local_frequency(self.alignment[species])
        for amino_acid in self.amino_acid_list:
            total = 0
            for species in species_list:
                total += (local_frequency[species][amino_acid]*alignment_mutation[species]
                          / self.alignment_length*100)
            result[amino_acid] = total
        return result

    def get_normalized_frequency(self, input_set):
        result_set = input_set.copy()
        total = 0
        for amino_acid in result_set:
            total += result_set[amino_acid]
        for amino_acid in result_set:
            result_set[amino_acid] /= total
        return result_set


    def get_relative_mutability_and_sum_amino_acid(self):
        raw_pam_matrix = self.construct_raw_PAM_matrix_and_get_alignment_mutation()[0]
        exposure_of_amino_acid_to_mutation = self.get_total_exposure_of_amino_acid_to_mutation()
        # initialize the relative mutability of each amino acid
        relative_mutability = {}
        # initialize the counting of the number of mutations of each amino acid
        sum_amino_acid = {}
        for amino_acid in self.amino_acid_list:
            count = 0
            for another_amino_acid in self.amino_acid_list:
                if amino_acid != another_amino_acid:
                    count += raw_pam_matrix[amino_acid][another_amino_acid]
            if exposure_of_amino_acid_to_mutation[amino_acid] == 0:
                relative_mutability[amino_acid] = 0
            else:
                relative_mutability[amino_acid] = count/exposure_of_amino_acid_to_mutation[amino_acid]
            sum_amino_acid[amino_acid] = count
        return relative_mutability, sum_amino_acid

    def normalize_relative_mutability(self, relative_mutability):
        new_relative_mutability = relative_mutability.copy()
        multiplier = 100/new_relative_mutability['A']
        for amino_acid in new_relative_mutability:
            new_relative_mutability[amino_acid] *= multiplier
        return new_relative_mutability

    def construct_JTT_matrix(self):
        # initialize the raw JTT matrix
        jtt_matrix = {'A': {}, 'R': {}, 'N': {}, 'D': {}, 'C': {}, 'Q': {}, 'E': {}, 'G': {}, 'H': {}, 'I': {},
                    'L': {}, 'K': {}, 'M': {}, 'F': {}, 'P': {}, 'S': {}, 'T': {}, 'W': {}, 'Y': {}, 'V': {}}
        for amino_acid in jtt_matrix:
            jtt_matrix[amino_acid] = {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0,
                                    'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0}

        # construct the raw JTT matrix
        raw_pam_matrix = self.construct_raw_PAM_matrix_and_get_alignment_mutation()[0]
        old_relative_mutability, sum_amino_acid = self.get_relative_mutability_and_sum_amino_acid()
        relative_mutability = self.normalize_relative_mutability(old_relative_mutability)
        for amino_acid in self.amino_acid_list:
            for another_amino_acid in self.amino_acid_list:
                if amino_acid != another_amino_acid:
                    if sum_amino_acid[amino_acid] == 0:
                        jtt_matrix[amino_acid][another_amino_acid] = 0
                    else:
                        jtt_matrix[amino_acid][another_amino_acid] = (relative_mutability[amino_acid] *
                                    raw_pam_matrix[amino_acid][another_amino_acid])/sum_amino_acid[amino_acid]

        # construct the complete JTT matrix
        for amino_acid in self.amino_acid_list:
            # calculate the proportionality constant
            count = 0
            for another_amino_acid in self.amino_acid_list:
                count += jtt_matrix[amino_acid][another_amino_acid]
            count += relative_mutability[amino_acid]  # p_constant is the the proportionality constant
            p_constant = 0
            if count != 0:
                p_constant = 1/count
            for another_amino_acid in self.amino_acid_list:
                if amino_acid == another_amino_acid:  # calculate the value of diagonal element
                    jtt_matrix[amino_acid][another_amino_acid] = 1 - p_constant*relative_mutability[amino_acid]
                else:  # calculate the value of non-diagonal element
                    jtt_matrix[amino_acid][another_amino_acid] *= p_constant
        return jtt_matrix

    def scale_matrix(self, matrix, scale_factor):
        '''
        Scale the matrix by the desired scale factor
        Arguments:  matrix: A 2D matrix that needs to be scaled
                    scale_factor: A numeric value representing the scale factor
        '''
        new_matrix = matrix.copy()
        for amino_acid in new_matrix:
            for another_amino_acid in new_matrix:
                new_matrix[amino_acid][another_amino_acid] *= scale_factor
        return new_matrix

    def calculate_percent_difference(self, frequency, matrix):
        total = 0
        for amino_acid in self.amino_acid_list:
            total += (frequency[amino_acid]*matrix[amino_acid][amino_acid])
        return 100*(1-total)

    def construct_relatedness_odd_matrix(self):
        new_matrix = self.construct_JTT_matrix()
        normalized_frequency = self.get_normalized_frequency(self.get_total_exposure_of_amino_acid_to_mutation())
        for amino_acid in self.amino_acid_list:
            for another_amino_acid in self.amino_acid_list:
                new_matrix[amino_acid][another_amino_acid] /= normalized_frequency[another_amino_acid]
        return new_matrix

    def construct_MDM_matrix(self):
        new_matrix = self.construct_relatedness_odd_matrix()
        for amino_acid in self.amino_acid_list:
            for another_amino_acid in self.amino_acid_list:
                if new_matrix[amino_acid][another_amino_acid] != 0:
                    new_matrix[amino_acid][another_amino_acid] = 10 * log10(new_matrix[amino_acid][another_amino_acid])
        return new_matrix

    def display_result(self, matrix):
        result = {}
        for amino_acid in self.amino_acid_list:
            result[amino_acid] = []
            for another_amino_acid in self.amino_acid_list:
                result[amino_acid].append(matrix[amino_acid][another_amino_acid])
        return pd.DataFrame(result, index=self.amino_acid_list)



a = JTTMatrix("PF12720.dated.nwk", "PF12720.fasta")
a.draw_tree()
lst = a.tree.get_terminals()
lst1 = a.tree.get_nonterminals()
print(lst)
print(lst1)
print(a.tree.distance(lst[0], lst[1]))
print(a.tree.trace(lst[0], lst[0]))
print(a.tree.trace(lst[0], lst1[0]))
print(lst1[0] == lst1[3])
print(str(lst[1]))
'''
b = a.get_normalized_frequency(a.get_total_exposure_of_amino_acid_to_mutation())
c = a.construct_JTT_matrix()
print(a.calculate_percent_difference(b, c))
'''
'''
b = a.construct_JTT_matrix()
c = a.display_result(b)
c.to_csv('result_compare.csv')
print(a.calculate_frequency())
'''

'''
print(a.tree)
print(len(a.tree.trace('Aspergillus_oryzae|11096019-PF10056', 'Aspergillus_terreus|11096021-PF10056')) == 2)
'''
