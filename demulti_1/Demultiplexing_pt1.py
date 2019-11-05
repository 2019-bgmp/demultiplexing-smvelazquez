#!/usr/bin/env python

### R1 = biological 1
### R2 = index1
### R3 = index 2
### R4 = biological 2

import numpy as np
import matplotlib.pyplot as plt
import gzip
import argparse


def get_args():
     parser = argparse.ArgumentParser()
     parser.add_argument("-f", "--FILEPATH", help="the absolute FILEPATH you want to analyze", required = True)
     parser.add_argument("-r", "--read_length", help = "the length of your read for your file", required=True, type = int)
     #parser.add_argument("-n", "--NUMBER_OF_RECORDS", help = "the number of records in your file", required = True, type = int)
     parser.add_argument('-H', "--HISTOGRAM_NAME", help = "the name of your inputted file will be generated into a pdf histogram", required = True)
     return parser.parse_args()

args = get_args()

def convert_phred(letter):
    ''' Converts a single character into a value (phred -33)'''
    x = ord(letter) -33
    return(x)

#READ_LENGTH = 101
#NUMBER_OF_RECORDS = ##
i = 0
record_count = 0

all_qscores = np.zeros(args.read_length)

with gzip.open(args.FILEPATH, 'rt') as fh:
    for line in fh:
        i +=1
        line = line.strip('\n')#line.decode("utf8").strip('\n')
        if i%4 == 0:
            record_count +=1
            for character in range(len(line)):
                letter = line[character] #takes the index value of the number and pulls out the actual letter associated with it
                y = convert_phred(letter)
                all_qscores[character] += y
        average_qscores = []
for item in range(len(all_qscores)):
    mean = all_qscores[item]/record_count
    average_qscores.append(mean)
    print(item, mean)

# for item in range(len(average_qscores)):
#     plt.bar(item, average_qscores[item])
#     plt.title('Average Quality Score Per Base Pair')
#     plt.xlabel('Base Pair')
#     plt.ylabel('Average Quality Score')
#     plt.savefig(HISTOGRAM_NAME + ".pdf")

plt.bar(range(len(average_qscores)), average_qscores)
plt.savefig(args.HISTOGRAM_NAME + ".pdf")
