from ete3 import Tree
from sys import argv
# Загружаем деревья
t1 = Tree(argv[1], format = 1)
t2 = Tree(argv[2], format = 1)

# Сравнение топологий (игнорируя длины ветвей)
rf = t1.robinson_foulds(t2, unrooted_trees=True)[0]

if rf == 0:
    print("Topologies are identical")
else:
    print(f"Robinson-Foulds distance: {rf}")