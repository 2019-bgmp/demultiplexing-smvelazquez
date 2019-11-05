import argparse

def get_args():
     parser = argparse.ArgumentParser()
     parser.add_argument("-r1", "--FILEPATH_R1", help="the absolute FILEPATH for your first biological read", required = True)
     parser.add_argument("-r2", "--FILEPATH_R2", help="the absolute FILEPATH with your second biological read", required = True)
     parser.add_argument("-I1", "--FILEPATH_I1", help="the absolute FILEPATH for your first index", required = True)
     parser.add_argument("-I2", "--FILEPATH_I2", help="the absolute FILEPATH with your second index", required = True)
     parser.add_argument("-index_file", "--FILEPATH_INDEXES", help="the absolute FILEPATH a list of all of your indexes (should be indexes.txt)", required = True)
     return parser.parse_args()
args = get_args()


import gzip as gz


def convert_phred(letter):
    ''' Converts a single character into a value (phred -33)'''
    x = ord(letter) -33
    return(x)

def rev_comp(seq):
    """Takes a sequence and returns the reverse the complement"""
    valid_bases = ['A', 'T', 'G', 'C']
    comp_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    bases=[]
    for base in seq:
        if base in valid_bases:
            bases.append(comp_dict[base])
        else: bases.append(base)
    rev_seq=bases[::-1]
   #print(bases)
   #print(rev_seq)
   #print("".join(rev_seq))
    return "".join(rev_seq)

def avg_index_qual_score(seq):
    '''finds the average quality score for each index quality score line you pass it'''
    mean_scores = []
    x = 0
    for letter in seq:
        y = convert_phred(letter)
        mean_scores.append(y)

    return(sum(mean_scores)/len(mean_scores))

#UNIT TESTS:

#avg_index_qual_score(I1_l4) #31.125
#avg_index_qual_score(I2_l4) #26.375

file = open(args.FILEPATH_INDEXES, "r")



inner_list = []
barcode_reverse_list = []
barcodes_dic = {}
barcodes_counter = {}
for line in file:
    inner_list.append(line.strip().split()[4])

for barcode in inner_list[1:]:
    #rev_barcode = rev_comp(barcode)
    index_pair = barcode + "-" + barcode
    barcodes_dic[index_pair] = [gz.open(index_pair + "_R1.fastq.gz", "at"), gz.open(index_pair + "_R2.fastq.gz", "at")]
    barcodes_counter.setdefault(index_pair, 0)
    #print(index_pair)
#print(barcodes_dic)
#print(barcodes_counter)
barcodes_list = inner_list[1:]

for item in barcodes_list:
    rev = rev_comp(item)
    barcode_reverse_list.append(rev)
#print(barcode_reverse_list)

#print(barcodes_list)
#print(barcode_reverse_list)



R1_file = gz.open(args.FILEPATH_R1, 'rt')
R2_file = gz.open(args.FILEPATH_R2, 'rt')
I1_file = gz.open(args.FILEPATH_I1, 'rt')
I2_file = gz.open(args.FILEPATH_I2, 'rt')
undetermined_R1 = gz.open("Undetermined_reads_R1.fastq.gz", "wt")
undetermined_R2 = gz.open("Undetermined_reads_R2.fastq.gz", "wt")
index_R1 = gz.open("Index_hopped_reads_R1.fastq.gz", "wt")
index_R2 = gz.open("Index_hopped_reads_R2.fastq.gz", "wt")

count = 0
undetermine = 0
index_hopped = 0
index_paired = 0
for file in R1_file, R2_file, I1_file, I2_file:
    while True:
        R1_l1 = R1_file.readline().strip() #this pulls out the header line
        if R1_l1 == "":
            break
        R1_l2 = R1_file.readline().strip() #index line
        R1_l3 = R1_file.readline().strip() # +
        R1_l4 = R1_file.readline().strip()# quality score line

        R2_l1 = R2_file.readline().strip() #this pulls out the header line
        if R2_l1 == "":
            break
        R2_l2 = R2_file.readline().strip() #index line
        R2_l3 = R2_file.readline().strip() # +
        R2_l4 = R2_file.readline().strip()# quality score line

        I1_l1 = I1_file.readline().strip()#this pulls out the header line
        if I1_l1 == "":
            break
        I1_l2 = I1_file.readline().strip() #index line
        I1_l3 = I1_file.readline().strip() # +
        I1_l4 = I1_file.readline().strip()# quality score line

        I2_l1 = I2_file.readline().strip() #this pulls out the header line
        if I2_l1 == "":
            break
        I2_l2 = I2_file.readline().strip() #index line
        I2_l3 = I2_file.readline().strip() # +
        I2_l4 = I2_file.readline().strip()# quality score line

        #print(R1_l1, R1_l2, R1_l3, R1_l4)


        #appending the index-pair to the header line of forward/reverse biological reads:
        R1_l1 = R1_l1 + ":" + I1_l2 + "-" + I2_l2
        R2_l1 = R2_l1 + ":" + I1_l2 + "-" + I2_l2
        #print(R1_l1,"\n", R2_l1)

        #actually demultiplexing reads
        if I1_l2 in barcodes_list and I1_l2 == rev_comp(I2_l2):
            matched_pair = I1_l2 + "-" + rev_comp(I2_l2)
            #print(matched_pair)
            if matched_pair in barcodes_dic:
                #print('hi')
                #print(barcodes_dic[matched_pair][0])
                barcodes_dic[matched_pair][0].write(R1_l1 + "\n" + R1_l2 + "\n" + R1_l3 + "\n" + R1_l4 + "\n" )
                barcodes_dic[matched_pair][1].write(R2_l1 + "\n" + R2_l2 + "\n" + R2_l3 + "\n" + R2_l4  + "\n")
                barcodes_counter[matched_pair] += 1
                count +=1
                index_paired +=1
            #print("paired:", matched_pair)
        elif I1_l2 not in barcodes_list and rev_comp(I2_l2) not in barcode_reverse_list:
            undetermined_pair = I1_l2 + "-" + rev_comp(I2_l2)
            undetermined_R1.write(R1_l1 + "\n" + R1_l2 + "\n" + R1_l3 + "\n" + R1_l4 + "\n")
            undetermined_R2.write(R2_l1 + "\n" + R2_l2 + "\n" + R2_l3 + "\n" + R2_l4 + "\n")
            count +=1
            undetermine +=1
            #print("undetermined:", undetermined_pair)
        elif (avg_index_qual_score(I1_l4) < 30) or (avg_index_qual_score(I2_l4) < 30):
            undetermined_R1.write(R1_l1 + "\n" + R1_l2 + "\n" + R1_l3 + "\n" + R1_l4 + "\n")
            undetermined_R2.write(R2_l1 + "\n" + R2_l2 + "\n" + R2_l3 + "\n" + R2_l4 + "\n")
            undetermine +=1
            count +=1
        elif rev_comp(I2_l2) != I1_l2:
            index_hopped_pair = I1_l2 + "-" + rev_comp(I2_l2)
            #print("index_hopped:", index_hopped_pair)
            index_R1.write(R1_l1 + "\n" + R1_l2 + "\n" + R1_l3 + "\n" + R1_l4 + "\n")
            index_R2.write(R2_l1 + "\n" + R2_l2 + "\n" + R2_l3 + "\n" + R2_l4 + "\n")
            count +=1
            index_hopped +=1
        #print(R1_l1,'\n', R2_l1)

#print(barcodes_counter, sum_barcodes)

for matched_pair in barcodes_dic:
    barcodes_dic[matched_pair][0].close()
    barcodes_dic[matched_pair][1].close()
undetermined_R1.close()
undetermined_R2.close()
index_R1.close()
index_R2.close()
print(barcodes_counter)
print("done")


for key, value in barcodes_counter.items():
   print(key, round((value/count)*100,1),"%")
print("Reads"+":",count)
print("undetermine"+":",round((undetermine/count)*100,1), "%")
print("index_hopped"+":",round((index_hopped/count)*100,1), "%")
print("index_paired"+":",round((index_paired/count)*100,1), "%")
