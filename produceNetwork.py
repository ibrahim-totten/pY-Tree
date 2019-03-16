import networkx as nx
import matplotlib.pyplot as plt 
import random 
import csv

CHROMOSOME = 2
START_POINT = 3 #0
END_POINT = 4 #1
NAME = 0 #0
CM = 4 #5
PERCENT = 14

exclude = []
percent_upper_threshold = 25

def elim_done_segs(start_point, cur_list):
    new_list = []
    for i in range(len(cur_list)):
        #print(i)
        if int(cur_list[i][END_POINT]) > start_point:
            new_list.append(cur_list[i])
    return new_list


G = nx.Graph()

cur_list=[]
node_labels = {}

f = open('relDT190312.csv', 'r')
reader = csv.reader(f)
headerCheck = 0

for row in reader:
    row = tuple(row)
    if headerCheck == 0  or row[CHROMOSOME] == '' or row[0] in exclude or float(row[PERCENT][:-1]) >= percent_upper_threshold:
        headerCheck = 1
        continue
    #convert chromosome X
    if row[CHROMOSOME] == 'X':
        chrom = 23
    else:
        chrom = int(row[CHROMOSOME])
      
    cur_seg = (int(row[START_POINT]), int(row[END_POINT]))

    cur_list = elim_done_segs(int(row[START_POINT]), cur_list)

    G.add_node(row)
    node_labels[row] = row[NAME]

    for neighbour in cur_list:
        if row[CHROMOSOME] == neighbour[CHROMOSOME]:
            G.add_edge(row, neighbour)
        
    cur_list.append(row)



f.close()

print(nx.info(G))

nx.draw(G, labels = node_labels, with_labels = True, edge_color='g')

plt.show()