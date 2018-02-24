
#!/usr/bin/python
#!/usr/bin/python
# coding=utf-8

#Created by Parth Patel, DBI @ University of Delaware, Newark, Delaware 19717
#Date created: 12/3/2017

##This script generates Sequence and Structural Features: 2nd-ary structure from miRNA precursor using RNAFold and parse the 2nd strcutre and generate structural features (dP,dG,MFEI1,MFEI2,MFEI3,MFEI4,NEFE) using RNAplot from Vienna RNA package.

import os,sys,subprocess,itertools,re,csv
import sys
if sys.version < '2.6':
    print ('You are using a version of Python that this program does not support. Please update to the latest version!')
    sys.exit(1)
from repDNA.nac import Kmer


########################## GLOBAL VARIABLES ###############################

positiveSet="Pos.txt"
negativeSet="Neg.txt"

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

###############################################################################################################
################################# Strucutral Features Start ###################################################
###############################################################################################################

#List of total basepairs in a predicted secondary strucutre of a sequence.
total_bp=[]

#List of number of stems in a predicted secondary strucutre of a sequence.
num_stems=[]

#List of number of loops in a predicted secondary strucutre of a sequence.
num_loops=[]

    
# def split_at(s, c, n):

#     words = s.split(c)
#     return c.join(words[:n]), c.join(words[n:])  
    
    
def RNAFOLD(Seq):

    """Calculates secondary strucuture of a seqeunce by calling the RNAfold program. 
     Args:
        seq: sequence to generate the GC content
        
     Returns:
        RNAfold.out file
    """
    os.environ['Seq']=Seq
    os.system("echo $Seq| RNAfold -d2 -p > RNAfold.out")

def RNAfoldParser(File):

    """Parser for RNAfold output file to compute the structural features. 
     Args:
        File: RNAfold.out file
        
     Returns: 
        Computed structural features: dP,dG,MFEI1,MFEI2,MFEI3,MFEI4,NEFE
    """

    fl_in = open(File,'r')
    lines = fl_in .readlines() #read SHORT SEQ. FILE INTO ARRAY

    #Compute total number of base pairs link - Dot bracket notation (http://rna.tbi.univie.ac.at/help.html)
    dot_notation=lines[1].split(" ", 1)[0]
    total_bases= dot_notation.count('(')
    total_bp.append(total_bases) # To compute average total basepairs, we need this array.
    # print "Total_bases=",total_bases

    #Compute total number of stems: n_stems is the number of stems in the secondary structure. A stem is a structural motif containing more than three contiguous base pairs
    dot_notation_parse=[(k, sum(1 for _ in vs)) for k, vs in itertools.groupby(dot_notation)] # e.g. [('(', 9), ('.', 3), ('(', 6), ('.', 9), (')', 6), ('.', 8), ('(', 6), ('.', 7), (')', 6), ('.', 2), (')', 9)]
    n_stems=0
    for i in range (0,len(dot_notation_parse)):
        if dot_notation_parse[i][0]=='(': # find all list with "(" which means there is a paired base
            if dot_notation_parse[i][1]>=3: #look for more than three contiguous base pairs
                n_stems+=1
    num_stems.append(n_stems) # To compute average stems, we need this array.

    #Compute total number of loops
    loops = []
    n_loops=0 # number of loops
    p = re.compile('[(][.]{1,100}[)]')
    loopsIter = p.finditer(dot_notation)
    # print loopsIter
    for m in loopsIter:    
        loops.append(m.span())
    n_loops=len(loops)
    num_loops.append(n_loops)
    # print "loops",len(loops)

    # Compute dP = 'Normalized base-pairing propensity'
    dP=float(total_bases)/len(lines[0].strip())
    # print "dP = ", dP
    #parse MFE
    lines[1]=lines[1].split(" ", 1)
    mfe= lines[1][1].split("(",1)[1].split(")",1)[0]
    # print "mfe = ",mfe

    #compute dG=norm_mfe
    dG=float(mfe)/len(lines[0].strip())
    # print len(lines[0].strip())
    # print "dG = ", dG

    #compute MFEI1 - 'Index 1 based on the minimum free energy'
    seq=lines[0].strip()
    g = seq.count("G")
    c = seq.count("C")
    t = seq.count("U")
    a = seq.count("A")
    gcCount = g+c
    totalBaseCount = g+c+t+a    
    gcFraction = float(gcCount) / totalBaseCount

    ##Compute MEFE1 - 'Index 1 based on the minimum free energy'
    if(gcFraction>0):
        MFEI1= dG/(gcFraction*100)
        # print "MFEI1 = ",MFEI1
    else:
        MFEI1 = dG

    ##Compute MEFE2 - 'Index 2 based on the minimum free energy'
    if(n_stems>0):
        MFEI2= float(dG)/n_stems
        # print "MFEI2 = ",MFEI2,n_stems
    else:
        MFEI2 = 0 

    ##Compute MEFE2 - 'Index 2 based on the minimum free energy'
    if(n_loops>0):
        MFEI3= float(dG)/n_loops
        # print "MFEI3 = ",MFEI3,n_loops
    else:
        MFEI3 = 0 
    
    #Compute MEFE4 - 'Index 4 based on the minimum free energy'
    if(total_bases>0):
        MFEI4= float(mfe)/total_bases
        # print "MFEI4 = ",MFEI4
    else:
        MFEI4 = 0 

    #parse EFE
    lines[2]=lines[2].split(" ", 1)
    efe= lines[2][1].split("[",1)[1].split("]",1)[0]
    # print "efe = ",efe

    #compute NEFE
    NEFE=float(efe)/len(lines[0].strip())
    # print "NEFE = ",NEFE

      
    fl_in.close()
    return dP,dG,MFEI1,MFEI2,MFEI3,MFEI4,NEFE


##############################################################################################################
################################# Strucutral Features End ####################################################
##############################################################################################################


def main(input_File,output_File,classValue,n):
    
    """ Main function generates sequence plus strucutral features. 
     Args:
        input_File: feature file Name
        output_File: output file Name
        classValue: Class label "yes" or "no"
        n: Number of k-mers to generate     
    """

    INPUT= open(input_File).readlines()
    
    alphabet = "ACGT"
    kmer_list= make_upto_kmer_list(n,alphabet)
   
    output= open(output_File,'w')
    for k in kmer_list:
        output.write('%s\t'% (k))

    print ("Sequence Features Generated!")
    
    output.write('gc_content\tdp\tdG\tMFEI1\tMFEI2\tMFEI3\tMFEI4\tNEFE\tclass\n')

    
    # count=0 #loop counter
    for lines in INPUT: #reading each line from the file 
        line = lines.split('\n')

        # Compute Sequence features from 1 to n k-mers.
        for i in n:# for each value of k count k-mers
            make_kmer_feature_vector(line[0],i,output)            
        seq=line[0].strip()

        # Compute GC content
        gc_per=GC(seq)
        
        # Compute secondary strucutre using RNAfold program.
        RNAFOLD(seq)

        # Compute secondary strucutre  features using a RNAfold parser.
        dP,dG,MFEI1,MFEI2,MFEI3,MFEI4,NEFE=RNAfoldParser("RNAfold.out")
        # count+=1

        output.write("%s\t" % (gc_per))
        output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (dP,dG,MFEI1,MFEI2,MFEI3,MFEI4,NEFE,classValue))
        
    output.close()

    print ("Strucutal_features Generated!")


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
    main(positiveSet_file,POS_out_file_name,POS_classvalue,k_values)
    
    print ("Generating NEG - FEATURE SET")  
    NEG_out_file_name="Negative_feat.txt"
    main(negativeSet_file,NEG_out_file_name,NEG_classvalue,k_values)

    print ("Merging Features Files")
    FINALFILE_NAME= "Features.txt"
    concateFiles(POS_out_file_name,NEG_out_file_name,FINALFILE_NAME)
    
    # removing files
    os.remove(POS_out_file_name)
    os.remove(NEG_out_file_name)
    os.remove(FINALFILE_NAME)
    
print ("Done..Enjoy Some Machine Learning bozo!")   
