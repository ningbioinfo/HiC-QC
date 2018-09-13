# HiC-QC

import os
from glob import glob
import re
import locale
locale.setlocale(locale.LC_ALL, 'en_US')
import subprocess
import pandas as pd
import argparse

'''
output format
colnames:
"Sequenced_Read_Pairs", "Ligations", "Unmapped", "Low_Mapping_Qual", "Unique_Aligned_Pairs", "Valid_Contacts",
"Duplicate_Contacts", "Intra_Fragment", "Inter_Chromosomal", "Intra_Chromosomal", "Intra_Short_Range",
"Intra_Long_Range", "Read_Pair_Type"
'''

def formatting(value):
    x = locale.format("%d", int(value), grouping=True)
    return x

def percent_formatting1(value,percent):
    x = formatting(value)+" ("+percent+")"
    return x

def percent_formatting2(value,percent1,percent2):
    x = formatting(value)+" ("+percent1+" / "+percent2+")"
    return x



parser = argparse.ArgumentParser(description='HiC QC script')
parser.add_argument('--data', nargs='?', help='path to the output directory of HiC-Pro pipeline')
parser.add_argument('--out', default='HiC_QC.csv', help='The QC output forms, by default its named HiC_QC.csv')


## read arguments
args = vars(parser.parse_args())

## get input
start_dir = args["data"]

## get hic experiment names
names = os.listdir(start_dir+"/bowtie_results/bwt2/")


## set output layout
c0 = ["Sequenced_Read_Pairs", "Normal_Reads", "Chimeric_Unambiguous","Ligations", "Unmapped", "Low_Mapping_Qual", "Unique_Aligned_Pairs", "Valid_Contacts",
      "Duplicate_Contacts", "Intra_Fragment", "Inter_Chromosomal", "Intra_Chromosomal", "Intra_Short_Range (< 20kb)",
      "Intra_Long_Range (> 20kb)", "Read_Pair_Type (L-I-O-R)"]

c_last = ["-", "-", "-", "30% - 40%", "less than 10%", "less than 10%", "-", "-", "less than 10%", "1% - 5%", "aroung or less than 20%", "aroung 60 - 70%",
        "around 20%", "at least 15%, good if more than 40%", "roughly 25% each"]

output = {'Metrics':c0, 'Recommand':c_last}

## get info from the results of HiC-Pro
for name in names:
    # Sequenced_Read_Pairs, Unmapped, Unique_aligned_pairs, Low_Mapping_Qual
    # these can be find in *.mpairstat file in bwt2


    #def get_first_4_row(start_dir, sample):
    mpairstat_files = []
    pattern   = name+"*.mpairstat"

    for dir,_,_ in os.walk(start_dir):
        mpairstat_files.extend(glob(os.path.join(dir,pattern)))

    #get numbers

    for line in open(mpairstat_files[0], 'r'):
        if re.search('^Total_pairs_processed', line):
            sequenced_read_pairs = line.strip().split('\t')[1]
        if re.search('^Unmapped_pairs', line):
            unmapped = line.strip().split('\t')[1]
            unmapped_percent = line.strip().split('\t')[2]+"%"
        if re.search('^Unique_paired_alignments', line):
            unique_aligned_pairs = line.strip().split('\t')[1]
            unique_aligned_pairs_percent = line.strip().split('\t')[2]+"%"
            unique_aligned_pairs_percent_of_bam = "100%"
        if re.search('^Low_qual_pairs', line):
            low_qual_pair = line.strip().split('\t')[1]
            low_qual_pair_percent = line.strip().split('\t')[2]+"%"

    # Valid_contact, Duplicate_contact, Intra-fragment, Inter_chromosomal,
    # Intra_chromosomal, Intra_short_range, Intra_long_range, Read_pair_type

    # most of them can be find in *.mergestat and *.mRSstat

    mergestat_files = []
    mRSstat_files = []
    pattern1  = name+"*.mergestat"
    pattern2  = name+"*.mRSstat"

    for dir,_,_ in os.walk(start_dir):
        mergestat_files.extend(glob(os.path.join(dir,pattern1)))
    for dir,_,_ in os.walk(start_dir):
        mRSstat_files.extend(glob(os.path.join(dir,pattern2)))


    for line in open(mergestat_files[0], 'r'):
        if re.search('^valid_interaction\t', line):
            valid_contact = line.strip().split('\t')[1]
            valid_contact_percent_of_all = str(round(int(valid_contact)/float(sequenced_read_pairs)*100,2))+"%"
            valid_contact_percent_of_bam = str(round(int(valid_contact)/float(unique_aligned_pairs)*100,2))+"%"
        if re.search('^valid_interaction_rmdup', line):
            duplicate_contact = int(valid_contact)-int(line.strip().split('\t')[1])
            duplicate_contact_percent_of_all = str(round(int(duplicate_contact)/float(sequenced_read_pairs)*100,2))+"%"
            duplicate_contact_percent_of_bam = str(round(int(duplicate_contact)/float(unique_aligned_pairs)*100,2))+"%"
        if re.search('^trans_interaction', line):
            inter_chrom = line.strip().split('\t')[1]
            inter_chrom_percent_of_all = str(round(int(inter_chrom)/float(sequenced_read_pairs)*100,2))+"%"
            inter_chrom_percent_of_bam = str(round(int(inter_chrom)/float(unique_aligned_pairs)*100,2))+"%"
        if re.search('^cis_interaction', line):
            intra_chrom = line.strip().split('\t')[1]
            intra_chrom_percent_of_all = str(round(int(intra_chrom)/float(sequenced_read_pairs)*100,2))+"%"
            intra_chrom_percent_of_bam = str(round(int(intra_chrom)/float(unique_aligned_pairs)*100,2))+"%"
        if re.search('^cis_shortRange', line):
            intra_short = line.strip().split('\t')[1]
            intra_short_percent_of_all = str(round(int(intra_short)/float(sequenced_read_pairs)*100,2))+"%"
            intra_short_percent_of_bam = str(round(int(intra_short)/float(unique_aligned_pairs)*100,2))+"%"
        if re.search('^cis_longRange', line):
            intra_long = line.strip().split('\t')[1]
            intra_long_percent_of_all = str(round(int(intra_long)/float(sequenced_read_pairs)*100,2))+"%"
            intra_long_percent_of_bam = str(round(int(intra_long)/float(unique_aligned_pairs)*100,2))+"%"

    for line in open(mRSstat_files[0], 'r'):
        if re.search('^Dangling_end_pairs', line):
            intra_fragment = line.strip().split('\t')[1]
            intra_fragment_percent_of_all = str(round(int(intra_fragment)/float(sequenced_read_pairs)*100,2))+"%"
            intra_fragment_percent_of_bam = str(round(int(intra_fragment)/float(unique_aligned_pairs)*100,2))+"%"
        if re.search('_RR', line):
            L_type = line.strip().split('\t')[1]
            L_type_percent = str(round(int(L_type)/float(valid_contact)*100,2))+"%"
        if re.search('_FF', line):
            R_type = line.strip().split('\t')[1]
            R_type_percent = str(round(int(R_type)/float(valid_contact)*100,2))+"%"
        if re.search('_RF', line):
            O_type = line.strip().split('\t')[1]
            O_type_percent = str(round(int(O_type)/float(valid_contact)*100,2))+"%"
        if re.search('_FR', line):
            I_type = line.strip().split('\t')[1]
            I_type_percent = str(round(int(I_type)/float(valid_contact)*100,2))+"%"

    read_pair_type = L_type_percent+"-"+I_type_percent+"-"+O_type_percent+"-"+R_type_percent

    # Ligation

    config_file = []
    pattern = "config*.txt"

    for dir,_,_ in os.walk(start_dir):
        config_file.extend(glob(os.path.join(dir,pattern)))

    for line in open(config_file[0], 'r'):
        if re.search('^LIGATION_SITE', line):
            ligation_site = line.strip().split(' = ')[1]

    path_to_raw = start_dir+"rawdata/"+name+"/"

    #get ligation site from read

    cmd = 'zgrep -c ' + ligation_site + ' ' + path_to_raw + '*fastq.gz'
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process_out = process.communicate()[0].decode('ascii')

    ligation_seq = []
    for i in process_out.strip().split('\n'):
        ligation_seq.append(i.split(':')[1])


    ligation_percent1 = str(round(int(ligation_seq[0])/float(sequenced_read_pairs)*100,2))+"%"
    ligation_percent2 = str(round(int(ligation_seq[1])/float(sequenced_read_pairs)*100,2))+"%"
    ligation_percent = ligation_percent1+"(R1) - "+ligation_percent2+"(R2)"


    # normal reads & chimeric reads
    # get pair name
    config_file = []
    pattern = "config*.txt"

    for dir,_,_ in os.walk(start_dir):
        config_file.extend(glob(os.path.join(dir,pattern)))

    for line in open(config_file[0], 'r'):
        if re.search('^PAIR1_EXT', line):
            pair1_name = line.strip().split(' = ')[1]
        if re.search('^PAIR2_EXT', line):
            pair2_name = line.strip().split(' = ')[1]

    # get file 1.global mapped 2. global unmapped trimmed
    #1.
    global_mapped_bam = []
    global_bam_dir = start_dir+"bowtie_results/bwt2_global/"+name+"/"
    pattern1 = "*"+pair1_name+"*.bwt2glob.bam"
    pattern2 = "*"+pair2_name+"*.bwt2glob.bam"

    for dir,_,_ in os.walk(global_bam_dir):
        global_mapped_bam.extend(glob(os.path.join(dir,pattern1)))
        global_mapped_bam.extend(glob(os.path.join(dir,pattern2)))

    normal_read = []
    for bam in global_mapped_bam:
        cmd = 'samtools view '+bam+' | wc -l'
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process_out = process.communicate()[0].decode('ascii')
        normal_read.append(process_out.strip())

    normal_read = formatting(normal_read[0])+" ("+str(round(int(normal_read[0])/float(sequenced_read_pairs)*100,2))+"%,R1) - "+formatting(normal_read[1])+" ("+str(round(int(normal_read[1])/float(sequenced_read_pairs)*100,2))+"%,R2)"

    # 2.
    chimeric_reads = []
    chimeric_trim_reads = []
    pattern1 = "*"+pair1_name+"*.unmap.fastq"
    pattern2 = "*"+pair1_name+"*.unmap_trimmed.fastq"
    pattern3 = "*"+pair2_name+"*.unmap.fastq"
    pattern4 = "*"+pair2_name+"*.unmap_trimmed.fastq"

    for dir,_,_ in os.walk(global_bam_dir):
        chimeric_reads.extend(glob(os.path.join(dir,pattern1)))
        chimeric_reads.extend(glob(os.path.join(dir,pattern3)))
        chimeric_trim_reads.extend(glob(os.path.join(dir,pattern2)))
        chimeric_trim_reads.extend(glob(os.path.join(dir,pattern4)))

    chimeric_unambiguous = []
    for i in [0,1]:
        cmd = 'diff '+chimeric_trim_reads[i]+' '+chimeric_reads[i]+' | grep -c "^<"'
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process_out = process.communicate()[0].decode('ascii')
        chimeric_unambiguous.append(int(int(process_out.strip())/2))

    chimeric_unambiguous = formatting(chimeric_unambiguous[0])+" ("+str(round(int(chimeric_unambiguous[0])/float(sequenced_read_pairs)*100,2))+"%,R1) - "+formatting(chimeric_unambiguous[1])+" ("+str(round(int(chimeric_unambiguous[1])/float(sequenced_read_pairs)*100,2))+"%,R2)"


    c1 = [formatting(sequenced_read_pairs), normal_read, chimeric_unambiguous, ligation_percent,
                  percent_formatting1(unmapped,unmapped_percent),
                  percent_formatting1(low_qual_pair, low_qual_pair_percent),
                  percent_formatting2(unique_aligned_pairs, unique_aligned_pairs_percent, unique_aligned_pairs_percent_of_bam),
                  percent_formatting2(valid_contact, valid_contact_percent_of_all, valid_contact_percent_of_bam),
                  percent_formatting2(duplicate_contact, duplicate_contact_percent_of_all, duplicate_contact_percent_of_bam),
                  percent_formatting2(intra_fragment, intra_fragment_percent_of_all, intra_fragment_percent_of_bam),
                  percent_formatting2(inter_chrom, inter_chrom_percent_of_all, inter_chrom_percent_of_bam),
                  percent_formatting2(intra_chrom, intra_chrom_percent_of_all, intra_chrom_percent_of_bam),
                  percent_formatting2(intra_short, intra_short_percent_of_all, intra_short_percent_of_bam),
                  percent_formatting2(intra_long, intra_long_percent_of_all, intra_long_percent_of_bam),
                  read_pair_type]

    output.update({name: c1})


df = pd.DataFrame(output)
names.insert(0, 'Metrics')
names.insert(1, 'Recommand')
output_df = df[names]
output_df.to_csv(args["out"], sep = "\t", index = False)
