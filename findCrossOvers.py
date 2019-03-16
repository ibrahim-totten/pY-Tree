import matplotlib.pyplot as plt
import numpy
import random
import csv

#curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip -c | grep acen > centromereLocs.tsv

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

#linear search
def search(suspect, lst):
    #print('===================')
    #print("Suspect: ", suspect)
    for i in range(len(lst)):
        #print(lst[i][0], lst[i][1])
        #print(suspect >= lst[i][0], suspect <= lst[i][1])
        #print("---------")
        if suspect >= lst[i][0] and suspect <= lst[i][1]: # and lst[i][1] - lst[i][0] != 0:
            return i
    return -1       

#INEFFICIENT, DO AGAIN WITH SEGMENT LAYERING AND DIVIDING
def buildHeatMap(chromosome):
    max_freq = 0
    segments = []   #(Starting point, ending point, frequency)
    freq = [0]*use_lengths_GRCh37[chromosome]
    f = open('relDT190312.csv', 'r')
    reader = csv.reader(f)
    headerCheck = 0


    segments.append((1, use_lengths_GRCh37[chromosome], 0))
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
            #make variable search alg later
            start_index_splice = search(cur_seg[0], segments)

            #middle of segment case: (Currently will make alligning SP and EPs be 'empty segments')
            #This must be split: cut off ending point (This remains the same)
            #First temporarily save the end of the segment before modifying it
            temp_start_index_end = segments[start_index_splice][1]
            segments[start_index_splice] = (segments[start_index_splice][0], cur_seg[0],segments[start_index_splice][2])

            #and inserting second half of split segment (This gets incremented)
            segments.insert(start_index_splice + 1, (cur_seg[0],  temp_start_index_end , segments[start_index_splice][2]))
   
            end_index_splice = search(cur_seg[1], segments)

            #This must be split: cut    (This gets incremented)
            #temporarily store the end
            temp_end_index_end = segments[end_index_splice][1]
            segments[end_index_splice] = (segments[end_index_splice][0], cur_seg[1],segments[end_index_splice][2])

            #and inserting second half of split segment (This remains the same)
            segments.insert(end_index_splice + 1, (cur_seg[1], temp_end_index_end, segments[end_index_splice][2]))

            #Increment all contained segments
            #print ("SI: " , start_index_splice , ", EI:" , end_index_splice)
            for incr in range(start_index_splice + 1, end_index_splice + 1):
                segments[incr] = (segments[incr][0], segments[incr][1], segments[incr][2] + 1)
                #update max frequency
                if segments[incr][2] > max_freq: 
                    max_freq = segments[incr][2]

    


    print (segments)
    f.close()
    return segments, max_freq

def heat_color(freq):
    #print("freq: ", freq)
    if freq == 0:
        return (0.0,0.0,0.0)
    if freq <0.2:
        return (1.0 - 5.0 * freq, 0.0 , 1.0)
    elif freq < 0.4:
        return (0.0, 5.0 * (freq - 0.2), 1.0)
    elif freq < 0.6:
        return (0, 1.0,1.0 - 5.0 * (freq - 0.4))
    elif freq < 0.8:
        return (5.0 * (freq - 0.6), 1.0, 0.0)
    else:
        return (1.0, 1.0 - 5.0 * (freq - 0.8),0.0)


def renderHeatMap(chromosomes):
    plt.axes()
    for chromosome in chromosomes:
        segments, max_freq = buildHeatMap(chromosome)
        print ("CHROMOSOME:",chromosome)
        print(segments)
        print("------------------------------------------------------------------")
        rectangle = plt.Rectangle(((15000000*chromosome), 0), 10000000, use_lengths_GRCh37[chromosome], ec= 'grey',fc='grey', joinstyle = 'round')
        plt.gca().add_patch(rectangle)

        for seg in segments:
            if seg[1]-seg[0] != 0 and max_freq != 0:
                #print("MAX_FREQ: ", max_freq)
                rectangle = plt.Rectangle(((15000000*chromosome), seg[0]), 10000000, seg[1]-seg[0], fc= heat_color(float(seg[2])/float(max_freq)))
                plt.gca().add_patch(rectangle)

    plt.axis('scaled')
    plt.show()

#Assuming it is already sorted by starting point
def findUnbrokenSegments(chromosome):
    segments = []
    f = open('relIT190312.csv', 'r')
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
        print ("CHROMOSOME:",chromosome)
        print(segments)
        print("------------------------------------------------------------------")
        rectangle = plt.Rectangle(((15000000*chromosome), 0), 10000000, use_lengths_GRCh37[chromosome], ec= 'grey',fc='grey', joinstyle = 'round')
        plt.gca().add_patch(rectangle)

        for seg in segments:
            rectangle = plt.Rectangle(((15000000*chromosome), seg[0]), 10000000, seg[1]-seg[0], ec='c',fc='b')
            plt.gca().add_patch(rectangle)

    plt.axis('scaled')
    plt.show()

#print All Chromosomes
#plotSegments(range(1,24))

#buildHeatMap(1)
renderHeatMap(range(1,24))

#buildHeatMap(22)