# HiC-QC

import os
from os import listdir
from glob import glob
import re
import locale
locale.setlocale(locale.LC_ALL, 'en_US')
import subprocess
import pandas as pd
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')


'''
output format
colnames:
"Sequenced_Read_Pairs", "Ligations", "Unmapped", "Low_Mapping_Qual", "Unique_Aligned_Pairs", "Valid_Contacts",
"Duplicate_Contacts", "Intra_Fragment", "Inter_Chromosomal", "Intra_Chromosomal", "Intra_Short_Range",
"Intra_Long_Range", "Read_Pair_Type"
'''

'''
9 Oct update: QC the coverage of data, measure and plot the distribution of average raw read count against distance.
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
parser.add_argument('--excel', default='0', help='Set this to 1 if you want an excel output')
parser.add_argument('--excel_name', default='HiC_QC.xlsx', help='The output excel file name of the QC result, by default it is HiC_QC.xlsx')
parser.add_argument('--tad', default='0', help='Set this to 1 if you want to plot TADs')
parser.add_argument('--tad_chr', default='chr1', help='by default it will call tads for chr1, if you want the whole genome, set this to "all", but it will take a while ')



## read arguments
args = vars(parser.parse_args())

## get input

start_dir = args["data"] + "/"

print("Your intput directory is ", start_dir)

## get hic experiment names
sample_names = os.listdir(start_dir+"/bowtie_results/bwt2/")

print("You have ", len(sample_names), " samples, and they are: ", sample_names)

## set output layout
c0 = ["Sequenced_Read_Pairs", "Normal_Reads", "Chimeric_Unambiguous","Ligations", "Unmapped", "Low_Mapping_Qual", "Unique_Aligned_Pairs", "Valid_Contacts",
      "Duplicate_Contacts", "Intra_Fragment", "Inter_Chromosomal", "Intra_Chromosomal", "Intra_Short_Range (< 20kb)",
      "Intra_Long_Range (> 20kb)", "Read_Pair_Type (L-I-O-R)"]

c_last = ["-", "-", "-", "30% - 40%", "less than 10%", "less than 10%", "-", "-", "less than 10%", "less than 20%", "aroung or less than 20%", "aroung 60 - 70%",
        "around 20%", "at least 15%, good if more than 40%", "roughly 25% each"]

output = {'Metrics':c0, 'Recommand':c_last}


## get info from the results of HiC-Pro
print("Extracting info from HiC-Pro output.")
for name in sample_names:
    # Sequenced_Read_Pairs, Unmapped, Unique_aligned_pairs, Low_Mapping_Qual
    # these can be find in *.mpairstat file in bwt2


    #def get_first_4_row(start_dir, sample):
    mpairstat_files = []
    pattern = name+"*.mpairstat"

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
    pattern = "*.config"

    for dir,_,_ in os.walk(start_dir):
        config_file.extend(glob(os.path.join(dir,pattern)))

    for line in open(config_file[0], 'r'):
        if re.search('^LIGATION_SITE', line):
            ligation_line = line.strip().split(' = ')


    if len(ligation_line) == 1:
        ligation_seq = "NA"
        ligation_percent = "NA"
    else:
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
    pattern = "*.config"

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
    if len(ligation_line) == 1:
        chimeric_reads = "NA"
        chimeric_trim_reads = "NA"
        chimeric_unambiguous = "NA"
    else:
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

    # distance vs. average read count check
    print("plotting...")

    samplemxdir = start_dir + '/' + 'hic_results/matrix/' + name
    resolutions = listdir(samplemxdir + "/iced")
    for res in resolutions:
        index = samplemxdir + "/raw/" + res + "/" + name + "_" + res + "_abs.bed"
        matrix = samplemxdir + "/iced/" + res + "/" + name + "_" + res + "_iced.matrix"
    # build the index
        binlib = {}
        with open(index, 'r') as bed:
            for line in bed:
                cor = [line.strip().split("\t")[0], int(line.strip().split("\t")[1]), int(line.strip().split("\t")[2])]
                binname = line.strip().split("\t")[3]
                binlib[binname]=cor


    # subsitute by index
        interaction = pd.read_csv(matrix, header=None, sep= "\t")
        interaction.columns = ['bin1', 'bin2', 'read_count']

        chr1 = []
        start1 = []
        end1 = []
        for i in interaction['bin1']:
            key = str(int(i))
            chr1.append(binlib[key][0])
            start1.append(binlib[key][1])
            end1.append(binlib[key][2])

        chr2 = []
        start2 = []
        end2 = []
        for i in interaction['bin2']:
            key = str(int(i))
            chr2.append(binlib[key][0])
            start2.append(binlib[key][1])
            end2.append(binlib[key][2])

        interaction['chr1']=chr1
        interaction['start1']=start1
        interaction['end1']=end1
        interaction['chr2']=chr2
        interaction['start2']=start2
        interaction['end2']=end2

        interaction = interaction[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'read_count']]

    # get only cis interactions
        cis_int = interaction[interaction.chr1 == interaction.chr2]

    # get a distance column (kb)
        distance = (cis_int.start2 - cis_int.start1)/1000
        distance = [int(x) for x in distance.tolist()]
        cis_int = cis_int.assign(dis = distance)

    # distance vs. average interaction count
        dva = {}
        for i in list(set(distance)):
            df = cis_int.loc[cis_int['dis'] == i]
            b = len(df)
            a = sum(df.read_count)
            c = a/b

            dva[i] = c

    # plot
        plt.scatter(list(dva.keys()),list(dva.values()))
        plt.yscale('log')
        plt.xlabel("Distance")
        plt.ylabel("average raw read count")
        plt.title(name + " at " + res + " resolution")

    # save plot
        plt.savefig(name + "_at_" + res + "_resolution.png")

    # clear plot
        plt.clf()


    if args['tad'] == '1':
        src_dir = os.path.dirname(os.path.abspath(__file__))
        try:
            open(src_dir+'/callTADs.R','r')
        except:
            print('can not find the callTADs.R to call TADs')
            exit()

    #### Call TADs
    # call tad and plot.
    # need : 1. matrix file 2. bin file
        print("calling TADs....")

    # matrix file
        samplemxdir = start_dir + '/' + 'hic_results/matrix/' + name
        resolutions = listdir(samplemxdir + "/iced")
    ## highest resolution
        res = str(min(list(map(int,resolutions))))
        mfile = samplemxdir + '/iced/' + res + '/' + name + '_' + res + '_iced.matrix'
    # bin file
        bfile = samplemxdir + '/raw/' + res + '/' + name + '_' + res + '_abs.bed'
    # chr
        tad_chr = "chr1"
    # tad.bed
        tad_bed = name + '_' + res + '_TADs.bed'
    # tad.plot
        tad_plot = name + '_' + res + '_TADs.pdf'

    # run R script
        cmd = " ".join(['Rscript', src_dir+'/callTADs.R', mfile, bfile, res, tad_chr, tad_bed, tad_plot])
        process = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'wb'))
        process.communicate()



# write output table to txt
df = pd.DataFrame(output)
df.to_csv(args["out"], sep = "\t", index = False)

# excel output
if args['excel'] == '1':
    print("Writing QC results to excel file...")
    writer = pd.ExcelWriter(args['excel_name'], engine="xlsxwriter")
    df.to_excel(writer, sheet_name="Sheet1", index=False)

    workbook = writer.book
    worksheet = writer.sheets['Sheet1']

    # Green fill with dark green text.
    greenformat = workbook.add_format({'bg_color':   '#C6EFCE',
                               'font_color': '#006100'})
    # Light red fill with dark red text.
    redformat = workbook.add_format({'bg_color':   '#FFC7CE',
                               'font_color': '#9C0006'})


    for name in sample_names:
        print(name)
        col_idx = sample_names.index(name) + 2
        # conditions
        # Ligation
        if 30 <= float(df.iloc[3][name].split('%')[0]) <= 40 and 30 <= float(df.iloc[3][name].split('%')[1].split(' - ')[1]) <= 40:
            worksheet.conditional_format(3+1,col_idx,3+1,col_idx,{'type':'no_blanks', 'format':greenformat})
        else:
            worksheet.conditional_format(3+1,col_idx,3+1,col_idx,{'type':'no_blanks', 'format':redformat})

        # unmapped
        if float(df.iloc[4][name].split('(')[1].rstrip('%)')) <= 10:
            worksheet.conditional_format(4+1,col_idx,4+1,col_idx,{'type':'no_blanks', 'format':greenformat})
        else:
            worksheet.conditional_format(4+1,col_idx,4+1,col_idx,{'type':'no_blanks', 'format':redformat})

        # mapping qual
        if float(df.iloc[5][name].split('(')[1].rstrip('%)')) <= 0:
            worksheet.conditional_format(5+1,col_idx,5+1,col_idx,{'type':'no_blanks', 'format':greenformat})
        else:
            worksheet.conditional_format(5+1,col_idx,5+1,col_idx,{'type':'no_blanks', 'format':redformat})

        # duplicate
        if float(df.iloc[8][name].split(' / ')[1].rstrip('%)')) <= 10:
            worksheet.conditional_format(8+1,col_idx,8+1,col_idx,{'type':'no_blanks', 'format':greenformat})
        else:
            worksheet.conditional_format(8+1,col_idx,8+1,col_idx,{'type':'no_blanks', 'format':redformat})

        # intra fragment
        if float(df.iloc[9][name].split('%')[0].split('(')[1]) <= 20:
            worksheet.conditional_format(9+1,col_idx,9+1,col_idx,{'type':'no_blanks', 'format':greenformat})
        else:
            worksheet.conditional_format(9+1,col_idx,9+1,col_idx,{'type':'no_blanks', 'format':redformat})

        # inter chrom
        if 15 <=float(df.iloc[10][name].split(' / ')[1].rstrip('%)')) <= 25:
            worksheet.conditional_format(10+1,col_idx,10+1,col_idx,{'type':'no_blanks', 'format':greenformat})
        else:
            worksheet.conditional_format(10+1,col_idx,10+1,col_idx,{'type':'no_blanks', 'format':redformat})

        # intra chrom
        if 55 <= float(df.iloc[11][name].split(' / ')[1].rstrip('%)')) <= 75:
            worksheet.conditional_format(11+1,col_idx,11+1,col_idx,{'type':'no_blanks', 'format':greenformat})
        else:
            worksheet.conditional_format(11+1,col_idx,11+1,col_idx,{'type':'no_blanks', 'format':redformat})

        # intra short
        if 15 <= float(df.iloc[12][name].split(' / ')[1].rstrip('%)')) <= 25:
            worksheet.conditional_format(12+1,col_idx,12+1,col_idx,{'type':'no_blanks', 'format':greenformat})
        else:
            worksheet.conditional_format(12+1,col_idx,12+1,col_idx,{'type':'no_blanks', 'format':redformat})

        # intra long
        if float(df.iloc[13][name].split(' / ')[1].rstrip('%)')) >= 15:
            worksheet.conditional_format(13+1,col_idx,13+1,col_idx,{'type':'no_blanks', 'format':greenformat})
        else:
            worksheet.conditional_format(13+1,col_idx,13+1,col_idx,{'type':'no_blanks', 'format':redformat})

        # rpt
        rpt = df.iloc[14][name].rstrip('%').split('%-')
        if 24 <= float(rpt[0]) <= 26 and 24 <= float(rpt[1]) <= 26 and 24 <= float(rpt[2]) <= 26 and 24 <= float(rpt[3]) <= 26:
            worksheet.conditional_format(14+1,col_idx,14+1,col_idx,{'type':'no_blanks', 'format':greenformat})
        else:
            worksheet.conditional_format(14+1,col_idx,14+1,col_idx,{'type':'no_blanks', 'format':redformat})

    writer.save()

print("Finished!")
