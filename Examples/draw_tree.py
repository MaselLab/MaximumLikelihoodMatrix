from bio import Phylo
import pylab

tree_path = input("Enter the tree path:\n")
tree = Phylo.read(tree_path, "newick")
Phylo.draw(tree)
pylab.show()