#!/gpfs/projects/bsc40/mschikora/anaconda3/envs/VarCall_CNV_env/bin/python

# This script is intended to run SIFT (https://sift.bii.a-star.edu.sg/sift4g/SIFT4G_codes.html#Database) given a vcf, reference genome and a gtf file

import argparse, os
import pandas as pd
import numpy as np
from argparse import RawTextHelpFormatter
import copy as cp
import pickle
import string
import shutil 
from Bio import SeqIO
from Bio import SeqRecord
import random
import sys
from collections import ChainMap

# define paths, depending if you are in the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    is_cluster = False    
else:
    is_cluster = True
    ParentDir = "/gpfs/projects/bsc40/mschikora"

# import functions
FunDir = "%s/scripts/VarCall_CNV_ReadProcessing"%ParentDir; sys.path.append(FunDir)
import functions as fun


# define paths
EnvDir = "/gpfs/projects/bsc40/mschikora/anaconda3/envs/VarCall_CNV_env"
sift4g = "%s/bin/sift4g"%EnvDir
nr_protein_database = "%s/databases/uniref90/uniref90.fasta"%ParentDir
perl = "%s/bin/perl"%EnvDir
scripts_to_build_SIFT_db_dir = "%s/software/scripts_to_build_SIFT_db"%ParentDir
gffread = "%s/software/gffread/gffread"%ParentDir

# parse cmd arguments
description = """
This script is intended to run SIFT (https://sift.bii.a-star.edu.sg/sift4g/SIFT4G_codes.html#Database) given a vcf, reference genome and a gtf file. It takes into consideration the possibility of concurrent runs. If another run is working it will not run.

It should be run from the VarCall_CNV_env, in /gpfs/projects/bsc40/mschikora/anaconda3/envs/VarCall_CNV_env

This pipeline can be tested running:

/home/mschikora/samba/scripts/VarCall_CNV_ReadProcessing/run_sift.py -r %s/Cglabrata_antifungals/data/Cglabrata_genomes_and_annotations/C_glabrata_CBS138_current_chromosomes.fasta -g %s/Cglabrata_antifungals/data/Cglabrata_genomes_and_annotations/C_glabrata_CBS138_current_features.gff -v %s/scripts/VarCall_CNV_ReadProcessing/test_dir/RUN1_CST34_2G_FLZ_VarCallresults/freebayes_ploidy1_out/output.filt.norm_gatk.vcf --replace

This script is based on the tutorial described in  https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB/blob/master/README.md#configFile (section Making a SIFT database from local genomic and gene annotation file (.gtf))

"""%(ParentDir, ParentDir, ParentDir)

parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

parser.add_argument("-r", "--ref", dest="ref", required=True, help="Reference genome. Has to end with .fasta")
#parser.add_argument("-g", "--gtf", dest="gtf", required=True, help="gtf annotation file")
parser.add_argument("-g", "--gff", dest="gff", required=False, default=None, help="gff annotation file")
parser.add_argument("-v", "--vcf", dest="vcf", required=True, help="vcf variants file")
parser.add_argument("--replace", dest="replace", action="store_true", help="Replace existing files")

parser.add_argument("-mchr", "--mitochondrial_chromosome", dest="mitochondrial_chromosome", default="mito_C_glabrata_CBS138", type=str, help="The name of the mitochondrial chromosome. This is important if you have mitochondrial proteins for which to annotate the impact of nonsynonymous variants, as the mitochondrial genetic code is different. This should be the same as in the gff")
parser.add_argument("-mcode", "--mitochondrial_code", dest="mitochondrial_code", default=3, type=int, help="The code of the NCBI mitochondrial genetic code. For yeasts it is 3. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi")


opt = parser.parse_args()


###############
# GENERAL THINGS
###############

# map each chromosome to a a numeric (or Mito) name. This is to avoid problems with the pipelines

#gtf_df = fun.load_gtf(opt.gtf)
gff_df =  fun.load_gtf(opt.gff)
all_chromosomes = sorted(set(gff_df.chromosome).difference({opt.mitochondrial_chromosome}))
originalChrName_to_changedChrName = {c : str(I+1) for I,c in enumerate(all_chromosomes)}
originalChrName_to_changedChrName[opt.mitochondrial_chromosome] = "Mito" # Mito is identified as mitochondrial by the make-SIFT-db-all.pl file

################################
###### MAKE SIFT DATABASE ######
################################

# make the directory where all the contents will be saved
siftDB_outdir = "%s/%s_SIFTdb"%("/".join(opt.gff.split("/")[0:-1]), opt.gff.split("/")[-1].split(".")[0])
if opt.replace is True: shutil.rmtree(siftDB_outdir)
if not os.path.isdir(siftDB_outdir): os.mkdir(siftDB_outdir)

# first create the SIFT database in the directory where the opt.gff is provided
print("\nCreating SIFT database. This is created under the opt.gff file")

# first obtain the translation of all the CDSs in the gff

# obtain the cds
cds_merged_file = "%s/all_cds_merged.fasta"%siftDB_outdir
fun.run_cmd("%s %s -g %s -x %s"%(gffread, opt.gff, opt.ref, cds_merged_file))

# map each transcript ID to a chromosome and a translation code
chromosome_to_transcriptIDs = {chrom : set([[x.split('"')[1] for x in a.split(";") if "transcript_id" in x][0] for a in set(gff_df[gff_df.chromosome==chrom]["attribute"])]) for chrom in set(gff_df.chromosome)}
chromosome_to_translationCode = {c:1 for c in set(gff_df.chromosome)}
chromosome_to_translationCode[opt.mitochondrial_chromosome] = opt.mitochondrial_code
transcriptID_to_chromosome = ChainMap(*[{t:chrom for t in transcripts} for chrom, transcripts in chromosome_to_transcriptIDs.items()])

# translate them according to the mitocode
aa_merged_file = "%s/all_cds_translated.fasta"%siftDB_outdir
SeqIO.write([SeqRecord.SeqRecord(seq=s.seq.translate(table=chromosome_to_translationCode[transcriptID_to_chromosome[s.id]]), id=s.id, description="") for s in SeqIO.parse(cds_merged_file, "fasta")], aa_merged_file, "fasta")

# build the sift4g database
database_outdir = "%s/database"%siftDB_outdir
if not os.path.isdir(database_outdir): os.mkdir(database_outdir)

# check if all transcripts have file in sift db
transcripts_with_db = set([x.split(".")[0] for x in os.listdir(database_outdir) if not fun.file_is_empty("%s/%s"%(database_outdir, x))])
all_transcripts = set([s.id for s in SeqIO.parse(aa_merged_file, "fasta")])
remaining_transcripts = all_transcripts.difference(transcripts_with_db)

if len(remaining_transcripts)>0: 
    print("Running sift database")
    fun.run_cmd("%s -q %s -d %s --threads 4 --out %s"%(sift4g, aa_merged_file, nr_protein_database, database_outdir))


##########
# RUN SIFT FOR THE GENES IN THE VCF. THIS HAS TO BE IMPLEMENTED, HOWEVER, SIFT ONLY CONSIDERS SNPs. USE PROVEAN FOR INDELS
###########################

