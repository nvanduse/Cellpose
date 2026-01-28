#!/bin/bash
#SBATCH --mail-type=ALL                 #directs Slurm to send job-related email when an event of the specified type(s) occurs; valid type values include all, begin, end, and fail.
#SBATCH --mail-user=nvanduse@iu.edu     #indicates the email address to which Slurm will send job-related mail      
#SBATCH --time=12:00:00                 #requested compute time in hrs:min:seconds
#SBATCH --mem=12G                       #memory requested
#SBATCH -A r00226                      #specifies that the job should run in the general partition
#SBATCH --cpus-per-task=12           #number of nodes requested. Its confusing but this is cores, not CPUs. Quartz nodes have 2 CPUs with 64 cores each. So the max number you can request here is 128. Parralelized parts of the script will use one core per input file.

#uncommment if desired:
###SBATCH --nodes=1                     #requests that a minimum of one node be allocated to this job.
###SBATCH --ntasks-per-node=12          #specifies the number of tasks that should be launched per node.

#General purpose compute nodes have 256gb of ram, two 64-core CPUs, and four 480GB solid state drives.
#module avail                               #shows what modules are available.
#squeue -u nvanduse -p dl -t PENDING        #put in your username and it will show all the jobs you are running 
#scancel -n my_job                          #will cancel your job it it's called "my_job"
#scancel 8990                               #will cancel job #8990

#####################################################################################################
#BEFORE STARTING

#ADD A TAB DELIMITED TWO COLUMN TEXT FILE ("region_annotation_and_sequence.txt") THAT HAS REGION ANNOTATION AND REGION SEQUENCE TO THE WORKING DIRECTORY. MAKE SURE IT HAS UNIX LINE ENDINGS.
    #this file has some annotation like region coordinates and region type in column A, and region barcode in column B, with no headers.  Example:
    #Myh7_M_Chr14:54968085-54968271_wt1_1   CTTCAGCTCCCACCCTATCTACTG
    #in this case the barcode is a 24bp sequence composed of the 4-bp interval barcode (differentiates wildtype and mutant tiles) + the first 20bp of each interval

#adjust the merged sequence length (PEAR), and the file name regular expressions (two lines in PEAR block, one line in Extract gRNA block)
#verify that reads have the expected structure by completing the block below.
#####################################################################################################
#VERIFY READS HAVE THE EXPECTED SEQUENCE

#MedGenome sequencing is 2x150bp paired end sequencing. 
#The reference sequence for merged R1/R2 reads is 260 bp (PEAR auto-trims off the ends of reads that extend into the opposite adapter).

#Reference amplicon (without adapters) --> 130bp.  Reads will extend 20bp past the end, into the adapter, but PEAR will auto-trim this off.
#cttgtggaaaggacgaaacaccgNNNNNNNNNNNNNNNNNNNNgttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaNNNNNNNNNNGGAGTTATGTGGGTCCCTAG
#Actual assembled amplicon:
#CTTGTGGAAAGGACGAAACACCGTACAAGCAGGACACCGATTAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAATATTAAAAGCGGAGTTATGTGGGTCCCTAG

#check and confirm that the reads have the anticipated sequence:
    ### srun -p general -A r00226 --pty bash        #move to an interactive session on a compute node
    ### gunzip -c Fyco_P4_m1_R1.fastq.gz | head -n 200 > Fyco_P4_m1_R1_first200lines.txt          #save the first 200 lines from each file to a text file. Note the lines at the very beginning are often worse quality than the majority of reads, so don't worry too much if they aren't all exactly as expected.
    ### gunzip -c Fyco_P4_m1_R2.fastq.gz | head -n 200 > Fyco_P4_m1_R2_first200lines.txt

#Example R1 reads:
#CTTGTGGAAAGGACGAAACACCGTATATGTTGGATCAGTCGAGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAACTTCGAGGGCGGAGTTATGTGGGTCCCTAGAGATCGGAAGAGCACACGTC
#CTTGTGGAAAGGACGAAACACCGGTGCTGGACGTGAATGTCCGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAACTTATTGGGGGAGTTATGTGGGTCCCTAGAGATCGGAAGAGCACACGTC
#CTTGTGGAAAGGACGAAACACCGAGGAGTTATGTGGGTCCCTAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATTCCGTTTTATTATTGAAAATTTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

#Example R2 reads:
#CTAGGGACCCACATAACTCCGCCCTCGAAGTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAACCTCGACTGATCCAACATATACGGTGTTTCGTCCTTTCCACAAGAGATCGGAAGAGCGTCGGGT
#CTAGGGACCCACATAACTCCCCCAATAAGTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAACCGGACATTCACGTCCAGCACCGGTGTTTCGTCCTTTCCACAAGAGATCGGAAGAGCGTCGTGT
#CTAGGGACCCACATAACTCCTCGGTGTTTCGTCCTTTCCACAAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAATTTGTTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG


###########################################################################################################################################################
#set up environment 
module load gnu-parallel/20210222       #need to load to enable parallel processing 
module load python/3.11.4                #needs to be loaded for pip to work
#if necessary install packages by running the commmands below (first time only; I think most of these are already included in Python 3.11.4 module):
#pip install PyPDF2
#pip install seaborn
#pip install pandas
#pip install matplotlib
#pip install numpy

######################################################################################################
#ANALYSIS PARAMETERS (UPDATE AS NEEDED; 2x150bp ONLY)

forward_pattern="*_R1.fastq"
reverse_pattern="*_R2.fastq"

# Splicing/cryptic-acceptor motifs (all searched within forward reads)
motif_correct_splicing="gcggccgc"
motif_failed_splicing="gcgatcgcagcccatatatg"
motif_cryptic_1="agggactttccattgacgtc"
motif_cryptic_2="gtatcatatgccaagtacgc"
motif_cryptic_3="tattagtcatcgctattacc"
motif_cryptic_4="gtgccgccatggtgagcaag"
motif_cryptic_5="ctacaagacccgcgccgagg"
motif_cryptic_6="gtatcaaggttacaagacag"
motif_cryptic_7="tctctctgcctattggtcta"

# Barcode extraction target (2x150bp only)
barcode_anchor_150="gcggccgc"                   # 2x150bp anchor (NotI)

# UMI anchor in reverse read (reverse-complement orientation in R2)
umi_left_rc="acataactcc"
umi_right_rc="acatgaagtt"

# Export motif/anchor variables so GNU parallel subshells see them
export motif_correct_splicing barcode_anchor_150 umi_left_rc umi_right_rc

######################################################################################################
#1) COUNT RAW READS CONTAINING SPECIFIC SEQUENCES (FORWARD READ ONLY)

splicing_stats_file="splicing_site_counts.tsv"
echo -e "sample\ttotal_reads\tcorrect_splicing\tfailed_splicing\tcryptic_1\tcryptic_2\tcryptic_3\tcryptic_4\tcryptic_5\tcryptic_6\tcryptic_7\tother\tcorrect_pct\tfailed_pct\tcryptic1_pct\tcryptic2_pct\tcryptic3_pct\tcryptic4_pct\tcryptic5_pct\tcryptic6_pct\tcryptic7_pct\tother_pct" > "$splicing_stats_file"

for r1 in $forward_pattern; do
    sample_name="${r1%.fastq}"
    sample_name="${sample_name%_R1}"
    awk -v sample="$sample_name" \
        -v m1="$motif_correct_splicing" \
        -v m2="$motif_failed_splicing" \
        -v m3="$motif_cryptic_1" \
        -v m4="$motif_cryptic_2" \
        -v m5="$motif_cryptic_3" \
        -v m6="$motif_cryptic_4" \
        -v m7="$motif_cryptic_5" \
        -v m8="$motif_cryptic_6" \
        -v m9="$motif_cryptic_7" \
        'BEGIN{FS="\t"; total=0; c1=0; c2=0; c3=0; c4=0; c5=0; c6=0; c7=0; c8=0; c9=0; other=0;}
        NR%4==2{
            seq=tolower($0);
            total++;
            if(index(seq,m1)) {c1++; next}
            if(index(seq,m2)) {c2++; next}
            if(index(seq,m3)) {c3++; next}
            if(index(seq,m4)) {c4++; next}
            if(index(seq,m5)) {c5++; next}
            if(index(seq,m6)) {c6++; next}
            if(index(seq,m7)) {c7++; next}
            if(index(seq,m8)) {c8++; next}
            if(index(seq,m9)) {c9++; next}
            other++;
        }
        END{
            if(total==0){total=1}
            printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", \
            sample,total,c1,c2,c3,c4,c5,c6,c7,c8,c9,other, \
            (c1/total)*100,(c2/total)*100,(c3/total)*100,(c4/total)*100,(c5/total)*100,(c6/total)*100,(c7/total)*100,(c8/total)*100,(c9/total)*100,(other/total)*100
        }' "$r1" >> "$splicing_stats_file"
done

######################################################################################################
#2) FILTER READS WITH NotI SITE IN FORWARD READ; COPY MATCHING REVERSE READS

filter_fastq_pairs() {
    r1="$1"
    r2="${r1/_R1/_R2}"
    if [ ! -f "$r2" ]; then
        echo "No matching reverse file found for $r1"
        return
    fi

    # Write uncompressed filtered FASTQ to avoid gz write complexity in awk
    out_r1="filtered_${r1}"
    out_r2="filtered_${r2}"

    # Stream paired FASTQ files with paste so each line has R1<TAB>R2
    paste "$r1" "$r2" | \
    awk -v motif="$motif_correct_splicing" -v out1="$out_r1" -v out2="$out_r2" '
        BEGIN{OFS="\n"}
        NR%4==1{r1_h=$1; r2_h=$2}
        NR%4==2{r1_s=$1; r2_s=$2}
        NR%4==3{r1_p=$1; r2_p=$2}
        NR%4==0{
            r1_q=$1; r2_q=$2;
            if(index(tolower(r1_s), motif)){
                print r1_h, r1_s, r1_p, r1_q >> out1;
                print r2_h, r2_s, r2_p, r2_q >> out2;
            }
        }'
}

export -f filter_fastq_pairs
parallel -j "$SLURM_CPUS_PER_TASK" filter_fastq_pairs ::: $forward_pattern

######################################################################################################
#3) EXTRACT BARCODES AND UMIs FROM FILTERED FILES
#   Output: "barcode-UMI_<filtered_R1_file_base>.txt"

process_filtered_pair() {
    r1="$1"
    r2="${r1/_R1/_R2}"
    if [ ! -f "$r2" ]; then
        echo "No matching reverse file found for $r1"
        return
    fi

    base_name="${r1%.fastq}"
    clean_name="${base_name#filtered_}"
    clean_name="${clean_name%_R1}"
    out_file="barcode-UMI_${clean_name}.txt"
    tmp_file="tmp_${clean_name}.txt"

    # 2x150bp only
    barcode_anchor="$barcode_anchor_150"
    barcode_len=20

    # Extract barcode and UMI per read-pair; only output when both are present
    paste "$r1" "$r2" | \
    awk -v anchor="$barcode_anchor" -v blen="$barcode_len" \
        -v left_rc="$umi_left_rc" -v right_rc="$umi_right_rc" '
        BEGIN{OFS="\t"}
        NR%4==1{r1_h=$1; r2_h=$2}
        NR%4==2{r1_s=$1; r2_s=$2}
        NR%4==3{r1_p=$1; r2_p=$2}
        NR%4==0{
            r1_q=$1; r2_q=$2;
            r1_seq=tolower(r1_s);
            r2_seq=tolower(r2_s);

            bpos=index(r1_seq, anchor);
            if(bpos){
                barcode=substr(r1_seq, bpos+length(anchor), blen);
            } else {
                next;
            }

            # UMI is between left_rc and right_rc in the reverse read
            if(match(r2_seq, left_rc "[acgtn]{10}" right_rc)){
                umi=substr(r2_seq, RSTART+length(left_rc), 10);
            } else {
                next;
            }

            if(length(barcode)==blen && length(umi)==10){
                print barcode, umi;
            }
        }' > "$tmp_file"

    # Filter for correct format and write final output
    awk -F"\t" -v blen="$barcode_len" '
        length($1)==blen && length($2)==10 {print $0}
    ' "$tmp_file" > "$out_file"

    rm "$tmp_file"
}

export -f process_filtered_pair
parallel -j "$SLURM_CPUS_PER_TASK" process_filtered_pair ::: filtered_*_R1.fastq

#####################################################################################################
#3b) LOG LINE COUNTS FOR BARCODE-UMI FILES

log_file="processing_log.txt"
echo -e "file\tlines" > "$log_file"
for file in barcode-UMI_*.txt; do
    line_count=$(wc -l < "$file")
    echo -e "$file\t$line_count" >> "$log_file"
done

#####################################################################################################
#4) DEDUPLICATE BARCODE/UMIs
#Only PCR duplicates should have the same barcode and the same UMI.

for file in barcode-UMI_*.txt; do
    sort -u "$file" -o "deduplicated_$file"
done

# Record deduplicated line counts in the log
for file in deduplicated_barcode-UMI_*.txt; do
    line_count=$(wc -l < "$file")
    echo -e "$file\t$line_count" >> "$log_file"
done

######################################################################################################
#5) SENSOR COUNTING FROM DEDUPLICATED FILES
#Uses sensor_pool_barcodes_first20.csv with headers: "sensor" and "barcode (first 20bp)"

python3 - <<'PYTHON_SCRIPT'
import pandas as pd
import glob
import os

barcode_table = pd.read_csv("sensor_pool_barcodes_first20.csv")
barcode_table = barcode_table.rename(columns=lambda c: c.strip())

barcode_col = "barcode (first 20bp)"
if barcode_col not in barcode_table.columns:
    raise ValueError(f"Expected column '{barcode_col}' in sensor_pool_barcodes_first20.csv")

# Normalize barcode strings to lowercase for robust matching
barcode_table[barcode_col] = barcode_table[barcode_col].astype(str).str.lower().str.strip()
barcode_to_sensor = barcode_table.set_index(barcode_col)

count_frames = []
for f in sorted(glob.glob("deduplicated_barcode-UMI_*.txt")):
    sample_name = os.path.basename(f).replace("deduplicated_barcode-UMI_", "").replace(".txt", "")
    df = pd.read_csv(f, sep="\t", header=None, names=["barcode", "umi"])
    df["barcode"] = df["barcode"].astype(str).str.lower().str.strip()
    counts = df["barcode"].value_counts()
    merged = barcode_to_sensor.copy()
    merged[sample_name] = merged.index.map(counts).fillna(0).astype(int)
    count_frames.append(merged[[sample_name]])

raw_counts = barcode_table.copy()
for frame in count_frames:
    raw_counts = raw_counts.join(frame, on=barcode_col)

raw_counts.to_csv("raw_counts.csv", index=False)

# 6) RPM normalization
rpm = raw_counts.copy()
count_cols = [c for c in rpm.columns if c not in ["sensor", barcode_col]]
for col in count_cols:
    total = rpm[col].sum()
    rpm[col] = (rpm[col] / total * 1_000_000) if total else 0

rpm.to_csv("RPM_counts.csv", index=False)
PYTHON_SCRIPT

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
#Switch to Python for easier manipulation of the data table

python3 - <<'PYTHON_SCRIPT'
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import PyPDF2
import os
import numpy as np


#############
# 1) Read the readcount_table.txt into pandas and remove duplicate columns --> reformat the header columns and save as tab delimitted txt. Will be input for MAGeCK.

data = pd.read_csv("readcount_table.txt", delimiter="\t")                # Read in the table
data.columns = data.columns.str.split('.').str[0]                        # Remove text from column headers that comes after the "."
data = data.loc[:, ~data.columns.duplicated()]                           # Remove duplicate columns if any
data['annotation'] = data.iloc[:, 0] + "-" + data.iloc[:, 1]             # Create the 'annotation' column by concatenating the first two columns with a "-"
data['gene'] = data.iloc[:, 0].str.split('_').str[0]                     # Create the 'gene' column by extracting the text before the underscore in the first column
data = data.drop(data.columns[[0, 1]], axis=1)                           # Now remove the original first two columns (columns 0 and 1)
columns = ['annotation', 'gene'] + [col for col in data.columns if col not in ['annotation', 'gene']]        # Reorder the columns to put 'annotation' and 'gene' at the beginning
data = data[columns]
data.to_csv("readcount_table_for_MAGeCK.txt", sep='\t', index=False)     # Save the DataFrame to a tab-delimited text file

#############
# 2) Read the readcount_table.txt into pandas and remove duplicate columns --> save as csv. Will be input for subsequent steps.

data = pd.read_csv("readcount_table.txt", delimiter="\t")                         # Read in the table
data.columns = data.columns.str.split('.').str[0]                                 # remove text from column headers that comes after the "."
data = data.loc[:, ~data.columns.duplicated()]                                    # deduplicate (will remove duplicated annotation columns)
data.to_csv("readcount_table_deduplicatedColumns.csv", index=False)               # Save the DataFrame to a CSV file

#############
# 3) Normalize the counts to counts-per-million (CPM)

data = pd.read_csv("readcount_table_deduplicatedColumns.csv")           # Read in the de-duplicated table
data_columns = data.columns[2:]                                         # Exclude the first two columns (Annotation and gRNA-seq)
cpm_data = data.copy()

for col in data_columns:
    cpm_data[col] = (data[col] / data[col].sum()) * 1_000_000

cpm_data.to_csv("readcount_table_CPM.csv", index=False)               # Save the CPM table to a CSV file

##############
# 4) Plot a histogram for each column using the log2(counts+1) transformed counts-per-million data

pdf_files = []
with pd.option_context("mode.use_inf_as_na", True):  # Handle infinite values as NaNs
    for col in data_columns:
        fig, ax = plt.subplots()
        log_data = np.log2(cpm_data[col] + 1)                                    # Apply the log2(counts+1) transformation to the column data
        sns.histplot(data=log_data, ax=ax, bins=40, kde=False)                   # plot with 40 bins
        ax.set_xlabel("log2(Counts per Million + 1)")
        ax.set_ylabel("# of tiles")
        ax.set_title(col)

        # Save histogram to a PDF file
        pdf_filename = f"{col}_histogram.pdf"                                   # pdf file name
        with plt.rc_context(rc={"figure.subplot.bottom": 0.2}):
            fig.savefig(pdf_filename, bbox_inches="tight")
        
        pdf_files.append(pdf_filename)
        plt.close(fig)

# Merge individual PDFs into a single PDF document
output_pdf = PyPDF2.PdfMerger()

for pdf_file in pdf_files:
    output_pdf.append(pdf_file)

# Save the merged PDF document
with open("histograms_merged.pdf", "wb") as f:
    output_pdf.write(f)

# Clean up individual PDF files
for pdf_file in pdf_files:
    os.remove(pdf_file)
#############################################

PYTHON_SCRIPT
