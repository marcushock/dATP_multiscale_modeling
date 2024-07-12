# This script is used to combine the reaction pairs for different structures 
# within a metastable state. This ultimately collects all the reaction pair 
# information from a metatstbale state and then pools it into a set (unique only)
# of reaction pairs and writes a new file in each directory. Meant to be ran from the 
# metastable_X directory (where X could be 0,1,2) above the cluster_N direcoties. 

# Last edited 07/11/2024 by MTH


all_pairs = []
for i in range(1,15+1):
    with open(('cluster_{i}/reaction_pairs.xml').format(i=i)) as f:
        lines = f.read().splitlines()[1:-1]
    f.close()
    for line in lines:
        all_pairs.append(line)


with open(('cluster_{i}/reaction_pairs.xml').format(i=i)) as f:
        lines = f.read().splitlines()


for i in range(1,15+1):
    try:
        f = open(('cluster_{i}/new_reaction_pairs.xml').format(i=i), 'x')
    except:
        f = open(('cluster_{i}/new_reaction_pairs.xml').format(i=i), 'w')
    f.write('<rxn_pairs>\n')
    for line in set(all_pairs):
        f.write(line+'\n')
    f.write('</rxn_pairs>')
    f.close()



