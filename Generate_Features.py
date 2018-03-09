
#!/usr/bin/python
#!/usr/bin/python
# coding=utf-8

#Created by Parth Patel, DBI @ University of Delaware, Newark, Delaware 19717
#Date created: 03/9/2018

##This script generates Sequence and Positional features: 

import os,sys,subprocess,itertools,re,csv,math
import sys
import numpy as np
if sys.version < '2.6':
    print ('You are using a version of Python that this program does not support. Please update to the latest version!')
    sys.exit(1)
from repDNA.nac import Kmer


########################## GLOBAL VARIABLES ###############################

positiveSet="Pos.txt"
negativeSet="Neg.txt"
Size=33 #33 #max size of sRNA seen in the postive and negative sets

###########################################################################



###############################################################################################################
################################# Sequence Features Start #####################################################
###############################################################################################################

def make_kmer_feature_vector(seq,n,outputFile):

    """Initialize Kmer and make k-mer vector for each sequence

    Args:
        Sequence: Name of the file to be read
        n: Number of kmer to be generated (n=5 generates kmer features from 1 to 5.)
        Outputfile: File to write k-mers

    """

    output=outputFile
    seq=[seq.strip()]

    if bNorm=="T":

        kmer=Kmer(k=n,normalize=True)
    else:
        kmer=Kmer(k=n)
    vector=kmer.make_kmer_vec(seq)
    for e in vector[0]:
        # if bNorm=="T":
        #     norm_e=float(e)/len(seq)
        #     norm_e=str(norm_e)
        #     output.write("%s\t" % (norm_e))
        # else:
        output.write("%s\t" % (e))

       
       
def make_kmer_list (k, alphabet):

    """Make k-mer labels or headers for each sequence

    Args:
        k: Number of kmer to be generated. (k=3 generates kmer labels for 3 mers.)
        alphabets: [A,C,G,T]

    Returns:
        The kmer labels  k-mers.

    """

    # Base case.
    if (k == 1):
        return(alphabet)

    # Handle k=0 from user.
    if (k == 0):
        return([])

    # Error case.
    if (k < 1):
        sys.stderr.write("Invalid k=%d" % k)
        sys.exit(1)

    # Precompute alphabet length for speed.
    alphabet_length = len(alphabet)

    # Recursive call.
    return_value = []
    for kmer in make_kmer_list(k-1, alphabet):
        for i_letter in range(0, alphabet_length):
            return_value.append(kmer + alphabet[i_letter])
              
    return(return_value)


def make_upto_kmer_list (k_values,alphabet):
  

    """Compute the k-mer labels or headers for each value of k

    Args:
        k_values: list of k-mer values. [i.e., k_values =[1,2,3,4,5])
        alphabets: [A,C,G,T]

    Returns:
        The kmer labels from 1 to k kmers.

    """    

    return_value = []
    for k in k_values:
        return_value.extend(make_kmer_list(k, alphabet))
    return(return_value)
  

def GC(seq):
    
    
     """Calculates G+C content, returns the percentage (float between 0 and 100). 
     Args:
        seq: sequence to generate the GC content
        
     Returns:
        The GC content for a seqeunce.
   
     Note that this will return zero for an empty sequence. 

     """ 
     
     try: 
         gc = sum(seq.count(x) for x in ['G', 'C', 'g', 'c']) 
         return gc*100.0/len(seq) 
     except ZeroDivisionError: 
         return 0.0     


##############################################################################################################
################################# Shanon Entropy #############################################################
##############################################################################################################

def ShanonEntropy(sequence): 

    """Calculates length Normalized Shanon Entropy, returns it. 
     Args:
        seq: sequence to calculate the shannon entropy
        
     Returns:
        The length normalized shannon entropy for a seqeunce.

     """ 

    seqList = list(sequence)
    nucleotides = list(set(seqList)) # list of symbols in the string

    
    # Calculate the frequency of each nucleotide in the sequence
    frequencyList = []
    for symbol in nucleotides:
        i = 0
        for sym in seqList:
            if sym == symbol:
                i += 1
        # print sequence,frequencyList
        frequencyList.append(float(i) / len(seqList))
    # print 'Frequencies of alphabet symbols:'
    
    # Shannon entropy
    ent = 0.0
    for freq in frequencyList:
        # print ent,freq
        ent = ent + freq * math.log(freq, 2)
    ent = -ent/len(sequence)
    # print 'Shannon entropy:'

    return ent

##############################################################################################################
################################# Positional Features ########################################################
##############################################################################################################

def generatePositionalLabels(size):

    """Generates postional labels for a sequence, returns list of labels 
     Args:
        Size: max size of sRNA seen in the postive and negative sets
        
     Returns:
        reutrns list of positional labels upto the size of the longest sequence

     """ 


    alphabets=['A','C','G','T']
    n=0
    i=1
    position_labels=[]
    while (n < len(alphabets)*size): #4*33=132

        if(n%4==0):
    #        print  i,alphabets[0],"->",master_positional_vector_phasiRNA[0][n]
            label=str(i)+alphabets[0]
            position_labels.append(label)
        elif(n%4==1):
    #        print  i,alphabets[1],"->",master_positional_vector_phasiRNA[0][n]
            label=str(i)+alphabets[1]
            position_labels.append(label)
        elif(n%4==2):
    #        print  i,alphabets[2],"->",master_positional_vector_phasiRNA[0][n]
            label=str(i)+alphabets[2]
            position_labels.append(label)
        elif(n%4==3):
    #        print  i,alphabets[3],"->",master_positional_vector_phasiRNA[0][n]
            label=str(i)+alphabets[3]
            position_labels.append(label)
        else:
            continue
        
        if (n%4==3):
                i+=1
    
        n=n+1
        
    return position_labels


def GeneratePositionalFeatures(Sequence,Size,output1):

    """Generates postional features for a sequence, returns list of positional features
     Args:        
        Sequence: Name of the file to be read
        Size: Max size of sRNA seen in the postive and negative sets
        Outputfile: File to write k-mers

     """ 

    sequence = Sequence
    n=Size 
    output=output1

    # vector storing the information of frequencies of a1,c1,g1,t1,a2,c2,g2,t2..a21,c21,g21,t21 for a sequence.    
    master_positional_vector= np.zeros((1,4*n),int)
    #print len(master_positional_vector_phasiRNA)


    seq_len=len(sequence)
    # print seq_len
    i=0       
    for y in range (0,n):
        
        if(y<seq_len):
    #        print y,"->",sequence[y]
            if (sequence [y] == 'A'):
                master_positional_vector[i][4*y]=1            
            elif (sequence [y] == 'C'):
                master_positional_vector[i][4*y+1]=1
            elif (sequence [y] == 'G'):
                master_positional_vector[i][4*y+2]=1
            elif (sequence [y] == 'T'):
                master_positional_vector[i][4*y+3]=1
            else:
                continue    
                              

    row,columns=master_positional_vector.shape


    alphabets=['A','C','G','T']

    L=0
    i=1
    position_labels=[]
    while (L < columns):

        if(L%4==0): #A
           output.write("%s\t" % (master_positional_vector[0][L]))
            
        elif(L%4==1): #C
           output.write("%s\t" % (master_positional_vector[0][L]))
            
        elif(L%4==2):#G
           output.write("%s\t" % (master_positional_vector[0][L]))
            
        elif(L%4==3):#T
            output.write("%s\t" % (master_positional_vector[0][L]))
            
        else:
            continue
        
        if (L%4==3): #Jump to
                i+=1
        
        L=L+1
            


def main(input_File,output_File,classValue,n,Size):
    
    """ Main function generates sequence plus strucutral features. 
     Args:
        input_File: feature file Name
        output_File: output file Name
        classValue: Class label "yes" or "no"
        n: Number of k-mers to generate
        Size: Max size of sRNA seen in the postive and negative sets  
    """

    INPUT= open(input_File).readlines()
    
    alphabet = "ACGT"
    kmer_list= make_upto_kmer_list(n,alphabet)
   
    output= open(output_File,'w')
    for k in kmer_list:
        output.write('%s\t'% (k))

    print ("Sequence Features Generated!")
    
    
    output.write('gc_content\tShannonEntropy\t')

    position_labels= generatePositionalLabels(Size)

    for ent in position_labels:
        output.write('%s\t'% (ent))
    
    output.write('class\n')
    
    
    # count=0 #loop counter
    for lines in INPUT: #reading each line from the file 
        line = lines.split('\n')

        # Compute Sequence features from 1 to n k-mers.
        for i in n:# for each value of k count k-mers
            make_kmer_feature_vector(line[0],i,output)            
        seq=line[0].strip()

        # Compute GC content
        gc_per=GC(seq)

        #Compute Shannon Entropy
        ShanonEnt=ShanonEntropy(seq)      
        

        output.write("%s\t%s\t" % (gc_per,ShanonEnt))     

        #Compute Positional Features
        GeneratePositionalFeatures(seq,Size,output)

        output.write("%s\n" % (classValue))
        
        
    output.close()

    print ("Positional Features Generated!")


def concateFiles(POS_out_file_name,NEG_out_file_name,FINALFILE_NAME):

    """Parser for concatinating postive and negative feature files. 
     Args:
        POS_out_file_name: positive feature file Name
        NEG_out_file_name: positive feature file Name
        FINALFILE_NAME: Final file Name after merging both files.
        
    """

    filenames = [POS_out_file_name,NEG_out_file_name]
    i=1
    with open(FINALFILE_NAME, 'w') as outfile:
        for fname in filenames:                            
            with open(fname) as infile:
                if(i==len(filenames)): 
                    next(infile)   #Skip the duplicate headers while merging two feature files.
                    for line in infile:
                        outfile.write(line)
                else:
                    for line in infile:
                        outfile.write(line)
            i+=1 
    
    txt_file = FINALFILE_NAME
    csv_file = FINALFILE_NAME+".csv"
    in_txt = csv.reader(open(txt_file, "rb"), delimiter = '\t')
    out_csv = csv.writer(open(csv_file, 'wb'))

    out_csv.writerows(in_txt)


if __name__ == '__main__':

    """ 
    Calls the main function

    """

    bNorm="T" # if you want to Nomrlaize frequency of k_mers by frequency of all kmers look at repDNA manual for better explanantion.
    k_values=[1,2,3,4,5]

    #bBothFeatures="F" # generating and concatinating features for both datasets.
    POS_classvalue="yes"
    NEG_classvalue="no"

    positiveSet_file=positiveSet
    negativeSet_file=negativeSet
  
    
    print ("Generating POS - FEATURE SET")
    POS_out_file_name="Positive_feat.txt"
    main(positiveSet_file,POS_out_file_name,POS_classvalue,k_values,Size)
    
    print ("Generating NEG - FEATURE SET")  
    NEG_out_file_name="Negative_feat.txt"
    main(negativeSet_file,NEG_out_file_name,NEG_classvalue,k_values,Size)

    print ("Merging Features Files")
    FINALFILE_NAME= "Features"
    concateFiles(POS_out_file_name,NEG_out_file_name,FINALFILE_NAME)
    
    # removing files
    os.remove(POS_out_file_name)
    os.remove(NEG_out_file_name)
    os.remove(FINALFILE_NAME)
    
print ("Done..Enjoy Some Machine Learning bozo!")   
