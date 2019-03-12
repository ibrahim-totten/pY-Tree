import matplotlib.pyplot as plt
import numpy
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

chromosome_lengths_GRCh37 = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566]
use_lengths_GRCh37 = [sum(chromosome_lengths_GRCh37)] + chromosome_lengths_GRCh37

centromere_positions_GRCh37 = [(121100000,124300000),]

#Assuming it is already sorted by starting point
def findUnbrokenSegments(chromosome):
    segments = []
    f = open('rel.csv', 'r')
    reader = csv.reader(f)
    headerCheck = 0

    for row in reader:
        if headerCheck == 0  or row[CHROMOSOME] == '' or row[0] in exclude or float(row[PERCENT][:-1]) >= percent_upper_threshold:
            headerCheck = 1
            continue
        #convert chromosome X
        if row[CHROMOSOME] == 'X':
            chrom = 23
        else:
            chrom = int(row[CHROMOSOME])

        cur_seg = (int(row[START_POINT]), int(row[END_POINT]))
        if chrom == chromosome:
            if len(segments) == 0:
                segments.append((cur_seg))
            else:
                if cur_seg[0] > segments[-1][1]:    #If the new segment begins after the last aggregate segment
                    segments.append((cur_seg))           #Create a new segment and add it to the end
                else:                               #If the new segment overlaps the last aggregate segment ends
                    segments[-1] = (segments[-1][0], cur_seg[1])    # Expand the last segment with this one

    f.close()
    return segments

#print(findUnbrokenSegments(13))

def plotSegments(chromosomes):
    plt.axes()
    for chromosome in chromosomes:
        segments = findUnbrokenSegments(chromosome)
        print(segments)
        rectangle = plt.Rectangle(((15000000*chromosome), 0), 10000000, use_lengths_GRCh37[chromosome], ec= 'grey',fc='grey')
        plt.gca().add_patch(rectangle)

        for seg in segments:
            rectangle = plt.Rectangle(((15000000*chromosome), seg[0]), 10000000, seg[1]-seg[0], ec='c',fc='b')
            plt.gca().add_patch(rectangle)

    plt.axis('scaled')
    plt.show()

#print All Chromosomes
plotSegments(range(1,24))
