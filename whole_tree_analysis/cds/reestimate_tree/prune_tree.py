from ete3 import Tree
from sys import argv

tree = Tree(argv[1], format=1)
species = argv[2].split(",")
new_tree_path = argv[3]

tree.prune(species)

with open(new_tree_path, 'w') as f:
    f.write(tree.write(format=1))