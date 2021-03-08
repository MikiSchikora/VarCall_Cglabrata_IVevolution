
        # generate the raw bedpe with only PASS rearrangements. 

        bedpe_onlyPASS = "%s.PASS"%(bedpe_file)
        fun.run_cmd("grep 'PASS||PASS' %s > %s"%(bedpe_file, bedpe_onlyPASS))
        nlines_onlyPASS = len(open(bedpe_onlyPASS, "r").readlines())
        bedpe_raw_onlyPASS = "%s.raw"%bedpe_onlyPASS
        fun.run_cmd("cut -f1-10 %s | tail -n %i > %s"%(bedpe_onlyPASS, nlines_onlyPASS-1, bedpe_raw_onlyPASS))
        os.unlink(bedpe_onlyPASS)

        # generate the raw BEDPE with all rearrangements
        nlines = len(open(bedpe_file, "r").readlines())
        bedpe_raw = "%s.raw"%bedpe_file
        fun.run_cmd("cut -f1-10 %s | tail -n %i > %s"%(bedpe_file, nlines-1, bedpe_raw))

        print("running clove")

        # calculate the median coverage of non mitochondrial genes
        median_cov = np.median(df_cov.loc[set(df_cov.index).difference({opt.mitochondrial_chromosome}), "coverage"])

        # run clove putting the median coverage and a window of +-0.9 of the median as thresholds for del and dup
        fun.run_cmd("%s -jar %s -i %s BEDPE -b %s -o %s -c %i %i"%(java, clove, bedpe_raw, sorted_bam, output_vcf_clove_tmp, median_cov, 1))
        #fun.run_cmd("%s -jar %s -i %s BEDPE -b %s -o %s -c %i %i"%(java, clove, bedpe_raw_onlyPASS, sorted_bam, output_vcf_clove_tmp, median_cov, 1)) #### DEBUG
        os.rename(output_vcf_clove_tmp, output_vcf_clove)

### QUALITY CONTROL OF THE BAM FILE ### THIS IS NOT CRITICAL

"""
bamqc_outdir = "%s/bamqc_out"%opt.outdir
if fun.file_is_empty("%s/qualimapReport.html"%bamqc_outdir) or opt.replace is True:
    print("Running bamqc to analyze the bam alignment")
    bamqc_cmd = "%s bamqc -bam %s -outdir %s -nt %i"%(qualimap, sorted_bam, bamqc_outdir, opt.threads); fun.run_cmd(bamqc_cmd)
"""



######### GATK



    # run GATK
    gatk_out = "%s/output.raw.vcf"%outdir_gatk; gatk_out_tmp = "%s.tmp"%gatk_out
    if fun.file_is_empty(gatk_out) or opt.replace is True:

        print("Running GATK HaplotypeCaller...")
        gatk_cmd = "%s HaplotypeCaller -R %s -I %s -O %s -ploidy %i --genotyping-mode DISCOVERY --emit-ref-confidence NONE --stand-call-conf 30 --native-pair-hmm-threads %i > %s.log"%(gatk, opt.ref, sorted_bam, gatk_out_tmp, opt.ploidy, opt.threads, gatk_out); fun.run_cmd(gatk_cmd)
        os.rename(gatk_out_tmp, gatk_out)

        # rename the index as well
        os.rename("%s.tmp.idx"%gatk_out, "%s.idx"%gatk_out)

    # variant filtration. There's a field called filter that has the FILTER argument
    gatk_out_filtered = "%s/output.filt.vcf"%outdir_gatk; gatk_out_filtered_tmp = "%s.tmp"%gatk_out_filtered
    if fun.file_is_empty(gatk_out_filtered) or opt.replace is True:

        print("Running GATK HaplotypeCaller Variant filtration...")

        # this depends on the ploidy. If ploidy is 2 you don't want to filter out heterozygous positions
        if opt.ploidy==1: filterHeterozygous = '-G-filter-name "heterozygous" -G-filter "isHet == 1"'
        else: filterHeterozygous = ''

        gatk_filt_cmd = '%s VariantFiltration -V %s -O %s -cluster 5 -window 20 %s --filter-name "BadDepthofQualityFilter" -filter "DP <= %i || QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" > %s.log'%(gatk, gatk_out, gatk_out_filtered_tmp, filterHeterozygous ,opt.coverage, gatk_out_filtered); fun.run_cmd(gatk_filt_cmd)
        os.rename(gatk_out_filtered_tmp, gatk_out_filtered)

        # rename the index as well
        os.rename("%s.tmp.idx"%gatk_out_filtered, "%s.idx"%gatk_out_filtered)



###########
# RUN gridss_purple_linx PIPELINE
##########

if opt.run_sv_gridss_purple_linx is True:

    # debug the fact that you don't have any callers
    if len(filtered_vcf_results)==0: raise ValueError("You have to run a variant caller with -sv_gridss_purple_linx")

    # debug the fact that you did not provide the reference fastq files
    if opt.fastq1_ref is None or opt.fastq1_ref is None: raise ValueError("You have to provide fastq1_ref and fastq2_ref to run -sv_gridss_purple_linx")

    # define a filtered vcf that has the variants of the sample
    if opt.caller=="all": sample_variants_vcf = [x for x in filtered_vcf_results if "HaplotypeCaller" in x][0]
    else: sample_variants_vcf = [x for x in filtered_vcf_results if opt.caller in x][0]

    # define the outdir where everything will be stored
    sv_gridss_purple_linx_outdir = "%s/gridss_purple_linx"%opt.outdir

    # run the pipeline
    df = fun.run_gridss_purple_linx(sorted_bam, sample_variants_vcf, opt.fastq1_ref, opt.fastq2_ref, sv_gridss_purple_linx_outdir, opt.ref, name_sample, opt.coverage, opt.gff, opt.mitochondrial_chromosome, threads=opt.threads, replace=False)





def run_freebayes_parallel(outdir_freebayes, ref, sorted_bam, ploidy, coverage, threads=4, replace=False):

    """IT GIVES SOME ERRORS"""

    # make the dir if not already done
    if not os.path.isdir(outdir_freebayes): os.mkdir(outdir_freebayes)

    #run freebayes
    freebayes_output ="%s/output.raw.vcf"%outdir_freebayes; freebayes_output_tmp = "%s.tmp"%freebayes_output
    if file_is_empty(freebayes_output) or replace is True:

        print("running freebayes in parallel")

        # define the regions file with the location of each chromosome. This is to parallelize on them
        regions_content = ["%s:0-%i"%(chrom.id, len(chrom.seq)) for chrom in SeqIO.parse(ref, "fasta")]
        regions_file = "%s.regions_for_freebayes"%freebayes_output
        open(regions_file, "w").write("\n".join(regions_content)+"\n")

        # remove the previous
        if not file_is_empty(freebayes_output_tmp): os.unlink(freebayes_output_tmp)

        # run the parallel freebayes
        run_cmd("%s %s %i -f %s -p %i --min-coverage %i --haplotype-length -1 -b %s -v %s"%(freebayes_parallel, regions_file, threads, ref, ploidy, coverage, sorted_bam, freebayes_output_tmp))
        os.rename(freebayes_output_tmp, freebayes_output)

    # filter the freebayes by quality
    freebayes_filtered = "%s/output.filt.vcf"%outdir_freebayes; freebayes_filtered_tmp = "%s.tmp"%freebayes_filtered
    if file_is_empty(freebayes_filtered) or replace is True:
        print("filtering freebayes")
        cmd_filter_fb = '%s -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" --tag-pass PASS %s > %s'%(vcffilter, freebayes_output, freebayes_filtered_tmp); run_cmd(cmd_filter_fb)
        os.rename(freebayes_filtered_tmp, freebayes_filtered)

    return freebayes_filtered



def run_cloveOnBedpe_and_process_output(bedpe_file, name, sorted_bam, replace=False, median_coverage=10, median_coverage_dev=1, check_coverage=True):

    """This function runs clove on a bedpe file writing files under name. . Then it writes, under outfile_clove, a table with insertions, deletions, inversions, tandemDuplications, translocations and other events. This is useful to compare against the output of simulateSV."""

    # first generate the clove outfile
    outfile_clove = "%s.%s.clove.vcf"%(bedpe_file, name)
    run_clove_filtered_bedpe(bedpe_file, outfile_clove, sorted_bam, replace=replace, median_coverage=median_coverage, median_coverage_dev=median_coverage_dev, check_coverage=check_coverage)

    # now get a dataframe with the clove output
    df_clove = get_clove_output(outfile_clove)

    # add the filter of the coverage. When the check_coverage is False, all TAN and DEL are set to PASS
    if check_coverage is True: df_clove["FILTERconsideringCoverage"] = df_clove["FILTER"]
    else: 

        def get_filter_considering_TANDEL(r):

            if r["SVTYPE"] in {"TAN", "DEL"}: return "PASS"
            else: return r["FILTER"]

        df_clove["FILTERconsideringCoverage"] = df_clove.apply(get_filter_considering_TANDEL, axis=1)


    print(df_clove[["FILTER", "FILTERconsideringCoverage", "SVTYPE"]], df_clove.keys())



def get_bedpe_for_clovebalTRA_5with3_sorted(r, chr_to_len):

    """Takes a row of the df_balTRA_5with3 df and returns a bedpe row, sorted"""

    # define the chrA (the first alphabetically)
    chrA = min(r["#CHROM"], r["CHR2"])
    chrB = max(r["#CHROM"], r["CHR2"])


    # define the choordinates
    if chrA==r["#CHROM"]: 
        startA = 0; endA = r["POS"]
        startB = r["END"]; endB = chr_to_len[r["CHR2"]]-1

    else: 
        startA = 0; endA = r["END"]
        startB = r["POS"]; endB = chr_to_len[r["#CHROM"]]-1


    return pd.Series({"chrA":chrA, "startA":startA, "endA":endA, "chrB":chrB, "startB":startB, "endB":endB, "Balanced":True})

def benchmark_processedSVs_against_knownSVs_suitable_for_compareSV_of_RSVSim(svtype_to_predsvfile, know_SV_dict, fileprefix, replace=False, analysis_benchmarking=False, tolerance_bp=50):

    """Takes two dictionaries that map some SVfiles. It runs, for all the types in svtype_to_predsvfile, a benchmarking against the known ones, writing a file under fileprefix. It returns a df of this benchmark, created with the benchmark_SVs_against_knownSVs_R script."""

    # define outfile, with the name of the benchmarked events
    outfile = "%s.benchmark%s.tab"%(fileprefix, "_".join(sorted(svtype_to_predsvfile)))

    if file_is_empty(outfile) or replace is True:
        print("benchmarking for %s"%outfile)

        # define the backbone cmd
        cmd = "%s --outfile_benchmark %s --tolerance_bp %i "%(benchmark_SVs_against_knownSVs_R, outfile, tolerance_bp)

        # now define extra arguments with the provided dict
        for svtype, predsvfile in svtype_to_predsvfile.items(): 
            if svtype in know_SV_dict: cmd += "--known_%s %s --predicted_%s %s "%(svtype, know_SV_dict[svtype], svtype, predsvfile)

        # run it
        Rstd = "%s.std"%outfile
        cmd = "%s > %s 2>&1"%(cmd, Rstd)
        run_cmd(cmd)

        remove_file(Rstd)

    # load into df
    df_benchmark = pd.read_csv(outfile, sep="\t")

    ### ANALYSIS OF THE BENCHMARK IN RELATIONSHIP TO UNCLASSIFIED BPs ####
    if analysis_benchmarking is True: analyse_benchmarking_processedSVs_against_knownSVs_suitable_for_compareSV_of_RSVSim(df_benchmark, svtype_to_predsvfile, know_SV_dict, fileprefix, unassigned_BPs_name="remaining_noTANDEL", tolerance_bp=tolerance_bp)

        
    return df_benchmark


def get_coverage_df_INVTX_ITX(df_bedpe, outdir, sorted_bam, reference_genome, window_l, replace=False, run_in_parallel=True):

    """Takes a dataframe bedpe, runs clove on it and calculates the coverage for all the ITX* and INVTX* region"""

    # define final output
    filename_df_INVTX_ITX = "%s/df_INVTX_ITX.py"%outdir

    if file_is_empty(filename_df_INVTX_ITX) or replace is True:
    #if True: # debug
    
        # write a bedpe that has the only necessary names
        raw_bedpe = "%s/raw_breakpoints.bedpe"%outdir
        important_fields = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2']
        df_bedpe[important_fields].to_csv(raw_bedpe, sep="\t", header=False, index=False)

        # run clove on this output
        outfile_clove = "%s.clove.vcf"%(raw_bedpe)
        run_clove_filtered_bedpe(raw_bedpe, outfile_clove, sorted_bam, replace=replace, median_coverage=10, median_coverage_dev=1, check_coverage=False)
        df_clove = get_clove_output(outfile_clove)

        # get the df that has ITX* or INVTX*
        df_INVTX_ITX = df_clove[df_clove.SVTYPE.isin({"ITX1", "IXT2", "INVTX1", "INVTX2"})]

        # get the coverage regions
        df_all = calculate_coverage_df_clove_regions_5end_3end(df_INVTX_ITX, "%s/regions_INVTX_ITX_calculate_coverage"%outdir, sorted_bam, reference_genome, window_l, replace=replace, run_in_parallel=run_in_parallel)

        # remove files
        for f in [raw_bedpe, outfile_clove]: remove_file(f)

        # save
        save_object(df_all, filename_df_INVTX_ITX)

    else: df_all = load_object(filename_df_INVTX_ITX)

    return df_all

def get_2bit_from_fasta(fasta_file):

    """Takes a fasta file and converts it to 2 bit"""

    # create a folder where the converted files will be
    fasta_dir = "/".join(fasta_file.split("/")[0:-1])
    fasta_name = fasta_file.split("/")[-1]
    twobit_conversion_dir = "%s/%s_twobit_conversion"%(fasta_dir, fasta_name.split(".")[0]); make_folder(twobit_conversion_dir)

    # split the fasta across chromosomes
    seqs_folder = "%s/chromosomes"%twobit_conversion_dir; make_folder(seqs_folder)
    chromosome_names = []
    for seq in SeqIO.parse(fasta_file, "fasta"):

        # change the description and ID
        seq.name = seq.id; seq.description = seq.id

        # change to upper
        seq.seq = seq.seq.upper()

        # write
        SeqIO.write([seq], "%s/%s.fa"%(seqs_folder, seq.id), "fasta")

        # keep
        chromosome_names.append(seq.id)

    chromosome_names_str = 'c("%s")'%('","'.join(chromosome_names))

    # generate the 2bit file with R
    r_script = "%s/convert_to_twoBit.R"%twobit_conversion_dir
    open(r_script, "w").write('library(BSgenome)\nforgeSeqFiles(%s, prefix="", suffix=".fa", seqs_srcdir="%s", seqs_destdir="%s", ondisk_seq_format="2bit")\n'%(chromosome_names_str, seqs_folder, twobit_conversion_dir))
    run_cmd("Rscript %s"%r_script)

    return "%s/single_sequences.2bit"%twobit_conversion_dir

def generate_BSgenome_package(genome, mitochondrial_chromosome, replace=False):

    """This function will generate a BSgenome package given a fasta file, generating the files under the same dir"""

    print("Making BS genome package for %s"%genome)

    # create a folder
    genome_dir = "/".join(genome.split("/")[0:-1])
    genome_name = genome.split("/")[-1]
    BSgenome_dir = "%s/%s_BSgenome"%(genome_dir, genome_name.split(".")[0]); make_folder(BSgenome_dir)

    # install the BS package --> this is only necessary the first time
    # run_cmd("R CMD INSTALL %s"%BSgenome_source_tarball)

    # split each chromosome into a file
    seqs_folder = "%s/chromosomes"%BSgenome_dir; make_folder(seqs_folder)
    chromosome_names = []
    all_chromosomes = []
    for seq in SeqIO.parse(genome, "fasta"):

        # change the description and ID
        seq.name = seq.id; seq.description = seq.id

        # change to upper
        seq.seq = seq.seq.upper()

        # write
        SeqIO.write([seq], "%s/%s.fa"%(seqs_folder, seq.id), "fasta")

        # keep
        chromosome_names.append(seq.id)
        all_chromosomes.append(seq)

    # write the corrected genome
    correct_genome = "%s/correct_genome.fasta"%BSgenome_dir
    SeqIO.write(all_chromosomes, correct_genome, "fasta")

    # generate a package name (only alphanumeric)
    bs_object_name = "".join([x for x in genome_name.split(".")[0] if x.isalnum()])
    package_name = "BSgenome.%s.GabaldonLab.Version1"%bs_object_name

    # checj if the package exists
    loading_file = "%s/checking_loading.R"%BSgenome_dir
    open(loading_file, "w").write("library(%s, quietly=TRUE)"%package_name + "\n")
    try: run_cmd("Rscript %s"%loading_file); existing_package = True
    except: existing_package = False

    if existing_package is False or replace is True:

        # make the seed file that will be used to construct the package
        seedfile_content = ["Package: %s"%package_name,
                            "Title: Full genome sequences for %s"%bs_object_name,
                            "Description: Full genome %s at %s"%(bs_object_name, get_date()),
                            "BSgenomeObjname: %s"%bs_object_name,
                            "seqs_srcdir: %s"%seqs_folder,
                            "Version: 1.0.0",
                            "License: GPL-2",
                            "Author: Miquel Angel Schikora",
                            "Maintainer: Miki <mikischikora@gmail.com>",
                            "organism: Genus species subspecies",
                            "common_name: organism",
                            "provider: GabaldonLab",
                            "organism_biocview: Genera_species",
                            "provider_version: Version1",
                            "release_date: Apr. 2019",
                            "release_name: release name",
                            "source_url: https://www.r-project.org",
                            "SrcDataFiles: random.2bit from https://www.r-project.org",
                            'seqnames: c("%s")'%('","'.join(chromosome_names)),
                            'mseqnames: c("%s")'%mitochondrial_chromosome]

        seed_file = "%s/seedfileBSgenome.txt"%BSgenome_dir
        open(seed_file, "w").write("\n".join(seedfile_content) + "\n")

        # remove previosuly generated files
        package_SourceTree_dir = "%s/%s"%(BSgenome_dir, package_name)
        if os.path.isdir(package_SourceTree_dir): shutil.rmtree(package_SourceTree_dir)

        # change the directory
        os.chdir(BSgenome_dir)

        # generate the forgeBSgenomeDataPkg with R
        forge_package_script = "%s/forge_package.R"%BSgenome_dir
        open(forge_package_script, "w").write('library(BSgenome)\nforgeBSgenomeDataPkg("%s", seqs_srcdir="%s", destdir="%s", verbose=TRUE)\n'%(seed_file, seqs_folder, BSgenome_dir))
        run_cmd("Rscript %s"%forge_package_script)

        # build the package
        run_cmd("R CMD build %s"%package_SourceTree_dir)
        tarball = "%s/%s_1.0.0.tar.gz"%(BSgenome_dir, package_name)

        # check it, without pdf generation
        run_cmd("R CMD check --no-manual --no-build-vignettes %s"%tarball)

        # install it
        run_cmd("R CMD INSTALL %s"%tarball)

    # return the name of the package, which will be useful for further usage
    return package_name

def generate_PON_folder_StructuralVariants(pon_outdir, sorted_bam, ref, mitochondrial_chromosome,  name_sample, threads = 4, replace=False):

    """This funtion creates a directory with breakpoints for a bam file. This directory is called PON in the gridss pipeline, and it is created in the way that the normal set of samples was created in the gridss-purple-linx pipeline """

    # make the output directory if it does not exist
    if not os.path.isdir(pon_outdir): os.mkdir(pon_outdir)

    # run gridss to get the structural variants
    gridss_outdir = "%s/gridss_outdir"%pon_outdir
    gridss_VCFoutput = run_gridss_raw(gridss_outdir, sorted_bam, ref, threads=threads, replace=replace)

    print(gridss_VCFoutput)

    # THIS DOES NOT WORK PROPERLY

def generate_PON_folder_StructuralVariants_using_scripts_of_Cameron2019_GridssPurpleLinx(pon_outdir, sorted_bam, ref, mitochondrial_chromosome,  name_sample, threads = 4, replace=False):

    """This funtion creates a directory with breakpoints for a bam file. This directory is called PON in the gridss pipeline, and it is created in the way that the normal set of samples was created in the gridss-purple-linx pipeline """

    # make the output directory if it does not exist
    if not os.path.isdir(pon_outdir): os.mkdir(pon_outdir)

    # run gridss to get the structural variants
    gridss_outdir = "%s/gridss_outdir"%pon_outdir
    gridss_VCFoutput = run_gridss_raw(gridss_outdir, sorted_bam, ref, threads=threads, replace=replace)

    # generate BS package
    bsgenome_package_name = generate_BSgenome_package(ref, mitochondrial_chromosome, replace=replace)

    # change the genomic package loaded by libgridss.R
    libgridss_file = "%s/libgridss.R"%libgridss_dir
    libgridss_content = []

    for line in open(libgridss_file, "r").readlines():

        if line.startswith("library(BSgenome."): libgridss_content.append("library(%s, quietly=TRUE)\n"%bsgenome_package_name)
        elif "vcf = readVcf(paste0(directory, filename)" in line: libgridss_content.append('\tvcf = readVcf(paste0(directory, filename), "%s")\n'%(bsgenome_package_name.split(".")[-1]))
        else: libgridss_content.append(line)

    open(libgridss_file, "w").write("".join(libgridss_content))

    # change the gridss_pon package
    create_gridss_pon_content = []
    for line in open(create_gridss_pon, "r").readlines():

        if "full_vcf = readVcf(vcf_file" in line: create_gridss_pon_content.append('  full_vcf = readVcf(vcf_file, "%s")\n'%(bsgenome_package_name.split(".")[-1]))
        else:  create_gridss_pon_content.append(line)

    open(create_gridss_pon, "w").write("".join(create_gridss_pon_content))

    # generate the gridss pon file
    print("Running pon formation")
    
    # remove previous files
    for file in os.listdir(pon_outdir):
        path = "%s/%s"%(pon_outdir, file)
        if os.path.isfile(path): os.unlink(path)

    # delet the cache
    cache_folder = "%s/Rcache"%pon_outdir
    if os.path.isdir(cache_folder): shutil.rmtree(cache_folder)

    # create the file that has the list of samples to use
    open("%s/gridss_pon_samples.txt"%pon_outdir, "w").write("%s\n"%name_sample)

    # run the pon formation
    # THIS FUNCTION HAS NEVER WORKED
    run_cmd("Rscript %s --pondir %s --scriptdir %s --batchsize 1 --normalordinal 1 --input %s"%(create_gridss_pon, pon_outdir, libgridss_dir, gridss_VCFoutput))

def run_gridss_purple_linx(sample_sorted_bam, sample_variants_vcf, fastq1_ref, fastq2_ref, outdir, ref_genome, name_sample, coverage, gff, mitochondrial_chromosome, threads=4, replace=False, name_ref="reference_strain"):

    """Takes several inputs and runs the gridss_purple_linx pipeline. This calls the variants that are in a sample (as indicated by sample_sorted_bam) against a reference sample
    - sample_sorted_bam is the sorted_bam of interest (it has to be indexed and sorted)
    - sample_variants_vcf is a vcf with the variants of the sample. These are important for a proper working of PURPLE

    In the default mode, this does not run any blacklisting nore considering that there are some repeats or viral sequences.

    Acc
    
    """

    # make the outdir if not already done
    if not os.path.isdir(outdir): os.mkdir(outdir)

    # for fastq1_ref, create a folder in the same dir
    ref_dir, ref_filename = ["/".join(fastq1_ref.split("/")[0:-1]), fastq1_ref.split("/")[-1].split(".")[0]]
    ref_dir = "%s/varcall_%s"%(ref_dir, ref_filename)
    if not os.path.isdir(ref_dir): os.mkdir(ref_dir)

    # create the sorted bamfile
    ref_bamfile = "%s/aligned_reads.bam"%ref_dir
    ref_sorted_bam = "%s.sorted"%ref_bamfile
    ref_index_bam = "%s.bai"%ref_sorted_bam

    # run bwa mem to get the bam for the reference bam
    run_bwa_mem(fastq1_ref, fastq2_ref, ref_genome, ref_dir, ref_bamfile, ref_sorted_bam, ref_index_bam, name_sample, threads=threads, replace=replace)

    # get only SNPs that PASS the filters in the sample
    sampleSNPs_PASS_vcf = get_onlySNPs_from_vcf_filtering(sample_variants_vcf)

    # link the reference genome to be under ref_dir
    reference_genome_under_ref_dir = "reference_genome.fasta"
    if file_is_empty("%s/%s"%(ref_dir, reference_genome_under_ref_dir)): run_cmd("ln -s %s %s/%s"%(ref_genome, ref_dir, reference_genome_under_ref_dir))

    # run GATK to identify heterozygous positions in the reference, which will be used by AMBER
    outdir_gatk = "%s/HaplotypeCaller_ploidy2"%ref_dir
    reference_het_vcf = run_gatk_HaplotypeCaller(outdir_gatk, ref_genome, ref_sorted_bam, ploidy=2, threads=threads, coverage=coverage, replace=replace)

    # get the bed file with the heterozygous positions
    bed_reference_heterozygous_positions = "heterozygous_position_reference.bed"
    reference_path_het_positions_bed = "%s/%s"%(ref_dir, bed_reference_heterozygous_positions)
    write_bed_heterozygous_positions(reference_het_vcf, reference_path_het_positions_bed, replace=replace)

    # get the annotation of the repeats' location and softlink under ref_dir
    repeat_masker_outfile = run_repeat_masker(ref_genome, threads=threads, replace=replace)
    repeat_masker_file = "repeat_masker.fasta.out"
    if file_is_empty("%s/%s"%(ref_dir, repeat_masker_file)): run_cmd("ln -s %s %s/%s"%(repeat_masker_outfile, ref_dir, repeat_masker_file))

    # create the gc profiles
    gcprofiles_file = "gcprofile.1000bp.cnp" # the human example can be found in https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd?path=%2FHMFTools-Resources%2FCobalt
    gcprofiles_path = "%s/%s"%(ref_dir, gcprofiles_file)
    generate_gcprofiles(ref_genome, gcprofiles_path, replace=replace, threads=threads, window_l=1000)

    # create a folder with the pon SV (these are the ones that are already in the reference)
    gridss_pon_dir = "%s/gridss_pon"%ref_dir
    generate_PON_folder_StructuralVariants(gridss_pon_dir, ref_sorted_bam, ref_genome, mitochondrial_chromosome, name_sample, threads=threads, replace=replace)

    """
    THIS IS HOW THE PON HAVE TO BE CREATED
    Variants called were annotated as “PON” based on overlap with a panel of normals (PON).
    A GRIDSS 2 PON was constructed from the 40x coverage WGS matched normals for 3,972
    patients using the gridss.GeneratePonBedpe utility. If multiple samples for a patient existed,
    only the normal for the first sample was included in the PON. Variants were aggregated
    across samples and the following filters were applied:
    - Existing GRIDSS 2 filters were cleared
    - All imprecise calls were filtered
    - Breakpoint variants with GRIDSS QUAL score less than 75 were filtered
    - Single breakend variants with a GRIDSS QUAL score of less least 428 were filtered
    A breakpoint BEDPE and single breakend BED file were generated requiring at least two
    supporting PON samples. A 2bp margin of error was allowed for breakpoint variants. If a
    breakpoint variant overlapped a PON breakpoint at one breakend and was within 2bp of
    overlapping on the other side, the bounds of the PON breakpoint was expanded to include
    the new variant. For breakpoint variants with microhomology, the variant was considered
    overlapping if any homologous position was overlapping. 
    """


    ######### EMPTY FILES ##########

    # create an empty list of blacklisted regions
    blacklist = "blacklisted_regions.bed"
    open("%s/%s"%(ref_dir, blacklist), "w").write("")

    ################################

    # create an empty structure of gridss pon variants 


    print("WARNING: running without blacklisting repeats nore ")


    # run the cmd of the pipeline
    run_cmd("%s --normal_bam %s -t %s --sample %s --snvvcf %s -o %s --threads %s --normal_sample %s --tumour_sample %s --ref_dir %s --install_dir %s --reference %s --repeatmasker %s --blacklist %s --bafsnps %s --gcprofile %s "%(gridss_purple_linx_exec, ref_sorted_bam, sample_sorted_bam, name_sample, sampleSNPs_PASS_vcf, outdir, threads, name_ref, name_sample, ref_dir, software_dir, reference_genome_under_ref_dir, repeat_masker_file, blacklist, bed_reference_heterozygous_positions, gcprofiles_file))

    #run_cmd("bash %s --normal_bam %s --sample %s --snvvcf %s --output_dir %s --threads %s --normal_sample %s --tumour_sample %s --ref_dir %s --install_dir %s --reference %s"%(gridss_purple_linx_exec, ref_sorted_bam, name_sample, sampleSNPs_PASS_vcf, outdir, threads, name_ref, name_sample, ref_dir, gridss_purple_linx_dir, "refgenomes/reference_genome"))

    print(sampleSNPs_PASS_vcf)


def get_bedpe_for_clove_unbalTRA_5with5_3with3(df_clove, chr_to_len, outdir, reference_genome, sorted_bam, replace=False, median_coverage=10, threshold_p=0.1):

    """Takes a df of CLOVE (that has ITX2) and returns the umbalanced translocations where the 5' is replacnig another 5' or a 3' is replacing another 3'. The df_clove should be filtered so that it does not include any ITX* that may contribute to balanced 5-to-5 or 3-to-3 translocations, or to complex interchromosomal duplications / translocations (which are insertions in simulateSV)

    There are 2 options: 

    - 5' of A replaces 5' of B: The coverage of B 5' should be low, while the coverage of A 5' should be high
    - 3' of A replaces 3' of B: The same as before, but inverted

    It returns a df with the bedpe of the unbalanced translocation, so that  the B region is removed
    """

    print("getting bedpe for 5'5' or 3'3' unbalanced translocations")
    make_folder(outdir)

    # define a df
    df_ITX = df_clove[df_clove.SVTYPE=="ITX2"]

    # chamge types
    df_ITX["POS"] = df_ITX.POS.apply(int); df_ITX["END"] = df_ITX.END.apply(int)

    # return empty dataframe if it is empty
    if len(df_ITX)==0: return pd.DataFrame()

    #print(df_ITX, df_ITX.keys()) # this is to debug

    # add the 3' and 5' regions of each breakpoint
    df_ITX["5end_#CHROM"] = df_ITX.apply(lambda r: (r["#CHROM"], 0, r["POS"]), axis=1)
    df_ITX["3end_#CHROM"] = df_ITX.apply(lambda r: (r["#CHROM"], r["POS"], chr_to_len[r["#CHROM"]]), axis=1)
    df_ITX["5end_CHR2"] = df_ITX.apply(lambda r: (r["CHR2"], 0, r["END"]), axis=1)
    df_ITX["3end_CHR2"] = df_ITX.apply(lambda r: (r["CHR2"], r["END"], chr_to_len[r["CHR2"]]), axis=1)

    # write a bed file with the region
    bed_dict = {}
    for region in ["5end_#CHROM", "3end_#CHROM", "5end_CHR2", "3end_CHR2"]:
        for IDX in df_ITX.index:
            bed_dict[(region, IDX)] = dict(zip(["chr", "start", "end"], df_ITX.loc[IDX, region]))

    df_bed = pd.DataFrame(bed_dict).transpose()[["chr", "start", "end"]]
    regions_bed = "%s/unbalanced_translocations_5with5_or_3with3.bed"%outdir
    df_bed.to_csv(regions_bed, sep="\t", header=False, index=False)

    # get the coverage_df of the interesting regions
    coverage_df = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir, sorted_bam, windows_file=regions_bed, replace=replace, window_l=0, run_in_parallel=False, delete_bams=False), sep="\t")
    coverage_df["relative_coverage"] = coverage_df.mediancov_1 / median_coverage
    coverage_df = coverage_df.merge(df_bed, left_on=["#chrom", "start", "end"], right_on=["chr", "start", "end"], validate="one_to_one")
    coverage_df.index = [(c, s, e) for c, s, e in coverage_df[["#chrom", "start", "end"]].values]

    # add columns of coverage of each region
    for region in ["5end_#CHROM", "3end_#CHROM", "5end_CHR2", "3end_CHR2"]: df_ITX["cov_%s"%region] = [float(coverage_df.loc[{coords}, "relative_coverage"]) for coords in df_ITX[region]]

    def define_unbalanced_TRA(r):

        """Takes a row of df_ITX and defines which is the unbalanced translocation"""

        # initialize a dict that will store probabilities of each theory being true
        p_to_hypothesis = {}
        hypothesis_to_bedpe = {}

        for hyp in {"5end", "3end"}:

            # calculate the probabilities as the difference in relative coverage to 

            # when the 5' replaces another 5', then CHR2 is the chrA and #CHROM is chrB
            if hyp=="5end":
                ChrA_id = "CHR2"
                ChrB_id = "#CHROM"

            # the other way arround
            else:
                ChrA_id = "#CHROM"
                ChrB_id = "CHR2"

            # define the duplicated and deleted_regions
            duplicated_region_cov = r["cov_%s_%s"%(hyp, ChrA_id)] 
            deleted_region_cov = r["cov_%s_%s"%(hyp, ChrB_id)]

            # define the probabilities of each event
            p_duplication = set_float_to_maximum(duplicated_region_cov - 1) # the higher the better, up to 1
            p_deletion = 1 - deleted_region_cov # the lower the better
            p_to_hypothesis[np.mean([p_duplication, p_deletion])] = hyp

            # define the region coords
            regionA = r["%s_%s"%(hyp, ChrA_id)]
            regionB = r["%s_%s"%(hyp, ChrB_id)]

            # keep the bedpe
            hypothesis_to_bedpe[hyp] = {"ChrA": regionA[0], "StartA": regionA[1],"EndA": regionA[2], "ChrB": regionB[0], "StartB": regionB[1], "EndB": regionB[2], "Balanced":False}

        # get the best bedpe
        best_p = max(p_to_hypothesis.keys())
        best_bedpe = hypothesis_to_bedpe[p_to_hypothesis[best_p]]
        best_bedpe["probability"] = best_p

        return pd.Series(best_bedpe)

    # get unbalanced TRA df and filter
    df_unbal_TRA = df_ITX.apply(define_unbalanced_TRA, axis=1)
    df_unbal_TRA = df_unbal_TRA[df_unbal_TRA.probability>=threshold_p]
    
    return df_unbal_TRA


def get_bedpe_for_cloveROW_unbalTRA_5with3_3with5(r, chr_to_len):

    """Takes a row of a clove DF that is unmatched INVTX1 or INVTX2 and returns a bedpe of the umblanaced translocation"""

    # define chromosomes, which are always the same ones (this is based on a couple of examples)
    ChrA = r["CHR2"]
    ChrB = r["#CHROM"]

    print("unbal_TRA_metadata", r["#CHROM"], r["POS"], r["CHR2"], r["END"], r["SVTYPE"])

    # 5' replaces 3'
    if r["SVTYPE"]=="INVTX1":

        StartA = 0; EndA = r["END"]
        StartB = r["POS"]; EndB = chr_to_len[ChrB]

    # 3' replaces 5'
    elif r["SVTYPE"]=="INVTX2":

        StartA = r["END"]; EndA = chr_to_len[ChrA]
        StartB = 0; EndB = r["POS"]
    
    else: raise ValueError("%s does not correspond to unbalanced translocations with 5' to 3' or 3' to 5'")

    # return the BEDPE
    return pd.Series({"ChrA": ChrA, "StartA": StartA,"EndA": EndA, "ChrB": ChrB, "StartB": StartB, "EndB": EndB, "Balanced":False})

def get_alignment_two_genomes_with_minimap2(query_genome, reference_genome, replace=False, threads=4):

    """Takes a query genome (or a multifasta) and aligns it agains the reference_genome, writing under query_genome. This is fine if the sequence divergence is of <15%"""

    # define outdir
    query_name = get_file(query_genome).split(".")[0]
    ref_name = get_file(reference_genome).split(".")[0]
    outdir = "%s/%s_alignedVS_%s"%(get_dir(query_genome), query_name, ref_name); make_folder(outdir)

    # define the outfiles
    samfile = "%s/aligned_genome.sam"%(outdir)
    bamfile = "%s/aligned_genome.bam"%(outdir)
    sorted_bam = "%s.sorted"%bamfile

    print("running minimap2 to align genomes")

    if file_is_empty(sorted_bam) or replace is True:

        # run minimap2 into samfile
        if file_is_empty(samfile) or replace is True:
            samfile_tmp = "%s.tmp"%samfile
            run_cmd("%s -a -o %s -t %i %s %s"%(minimap2, samfile_tmp, threads, reference_genome, query_genome))
            os.rename(samfile_tmp, samfile)

        # convert to bam
        if file_is_empty(bamfile) or replace is True:
            bamfile_tmp = "%s.tmp"%bamfile
            run_cmd("%s view -Sbu %s > %s"%(samtools, samfile, bamfile_tmp))
            os.rename(bamfile_tmp , bamfile)

        # sort bam
        if file_is_empty(sorted_bam) or replace is True:
            sorted_bam_tmp = "%s.tmp"%sorted_bam
            sort_bam(bamfile, sorted_bam_tmp, threads=threads)
            os.rename(sorted_bam_tmp , sorted_bam)

    # index bam
    if file_is_empty("%s.bai"%sorted_bam) or replace is True: index_bam(sorted_bam, threads=threads)

    # remove files
    for f in [bamfile, samfile]: remove_file(f)

    return sorted_bam



def simulate_longReads_per_chromosome(chr_obj, fraction_chr_per_read=0.5, fraction_chr_window_increment=0.01, n_extra_reads=30):

    """This function takes a chromosome object (SeqRecord) and it generates reads that are as long as fraction_chr_per_read, in a window that increases each time by fraction_chr_window. It returns a list of objects, each of which has a chromosomeID_readID"""

    # define metrics
    len_window_increment = int(len(chr_obj)*fraction_chr_window_increment)
    len_read = int(len(chr_obj)*fraction_chr_per_read)
    len_chromosome = len(chr_obj)

    # debug 
    if len_read<=len_window_increment and not (fraction_chr_per_read==1.0 and fraction_chr_window_increment==1.0): raise ValueError("the read length has to be longer that the window increment")

    # go through each window
    all_reads = []
    startRead = 0
    while startRead<len_chromosome:

        # define the end of the read
        endRead = startRead + len_read

        # keep the read
        read = chr_obj[startRead:endRead]
        ID = "%s_readStart%i"%(read.id, startRead)
        read.id = ID; read.name = ""; read.description = ""
        all_reads.append(read)

        # increment for the next window
        startRead += len_window_increment

    # add some reads of the begining
    for I in range(n_extra_reads):

         # keep the read
        read = chr_obj[0:len_read]
        ID = "%s_readStart0_id%i"%(read.id, I)
        read.id = ID; read.name = ""; read.description = ""
        all_reads.append(read)


    return all_reads





lakhflfahlahf
# first simulate reads in parallel for each chromosome
inputs_fn = [(chr_obj, fraction_chr_per_read, fraction_chr_window_increment, n_extra_reads) for chr_obj in SeqIO.parse(query_genome, "fasta")]
with multiproc.Pool(multiproc.cpu_count()) as pool:
    all_reads_list = pool.starmap(simulate_longReads_per_chromosome, inputs_fn)
    pool.close()

#all_reads_list = list(map(lambda x: simulate_longReads_per_chromosome(x[0], x[1], x[2]), inputs_fn))
simulated_reads = "%s/%s_simulated_long_reads_fractChr%.2f_fractWindow%.3f_nextraReads%i.fasta"%(working_dir, get_file(query_genome), fraction_chr_per_read, fraction_chr_window_increment, n_extra_reads)
if file_is_empty(simulated_reads) or replace is True: SeqIO.write(make_flat_listOflists(all_reads_list), simulated_reads, "fasta")

# define the sorted bam and the numbers of reads
sorted_bam = "%s.%s.coordsorted.bam"%(simulated_reads.rstrip(".fasta"), aligner)
index_bam = "%s.bai"%sorted_bam

# remove all intermediate files
for f in os.listdir(working_dir):
    filename = "%s/%s"%(working_dir, f)
    if filename not in {sorted_bam, index_bam, simulated_reads} and ".coordsorted." not in filename: remove_file(filename); delete_folder(filename)

#run_cmd("%s alignment %s %s %s --min_sv_size %i --skip_genotyping --minimum_score 0 --minimum_depth 0 --max_sv_size %i --trans_destination_partition_max_distance 1000000000 --trans_partition_max_distance 10000000 --trans_sv_max_distance 1000000 --min_mapq 0 --segment_gap_tolerance 10000000 --segment_overlap_tolerance 100000000 --partition_max_distance 10000000 --cluster_max_distance 1.0 --del_ins_dup_max_distance 10 --types BND"%(svim, working_dir, sorted_bam, reference_genome, min_sv_size, max_sv_size))

run_cmd("%s reads %s %s %s --min_sv_size %i --max_sv_size %i --cores %i --aligner %s --minimum_depth 0 --min_mapq 0 --skip_genotyping --nanopore"%(svim, working_dir, simulated_reads, reference_genome, min_sv_size, max_sv_size, threads, aligner))



    ########## SIMULATE LONG READS AS IF THEY WERE PAIRED ###############

    # simulate long reads, but paired, with wgsim
    readPairs = int(100/fraction_chr_per_read) # these are arround 100 reads per position
    median_insert_size = 100
    median_insert_size_sd = 30
    error_rate = 0.0

    # define the simulated reads output
    ID = "nreads%i_insertsize%i+-%i_error%.3f"%(readPairs, median_insert_size, median_insert_size_sd, error_rate)
    all_reads_1 = "%s/wgsim_%s_all_reads1.fq.gz"%(working_dir, ID)
    all_reads_2 = "%s/wgsim_%s_all_reads2.fq.gz"%(working_dir, ID)

    if any([file_is_empty(x) for x in [all_reads_1, all_reads_2]]) or replace is True:

        print("getting reads in parallel")      
        inputs_fn = [(query_genome, chr_obj.id, 0, len(chr_obj), readPairs, int(len(chr_obj))*fraction_chr_per_read, working_dir,  median_insert_size, median_insert_size_sd, replace, error_rate) for chr_obj in SeqIO.parse(query_genome, "fasta")]

        with multiproc.Pool(multiproc.cpu_count()) as pool:
            all_pairedReads_tupleFiles = pool.starmap(run_wgsim_pairedEnd_for_window, inputs_fn)
            pool.close()

        # concatenate all reads into one
        print("concatenating all reads in one")
        all_reads_1_tmp = "%s.tmp"%all_reads_1; all_reads_2_tmp = "%s.tmp"%all_reads_2;
        run_cmd("cat %s/*%s_read1.fq.gz > %s"%(working_dir, ID, all_reads_1_tmp)) 
        run_cmd("cat %s/*%s_read2.fq.gz > %s"%(working_dir, ID, all_reads_2_tmp)) 

        # remove all the reads
        print("Removing")
        for f1, f2 in all_pairedReads_tupleFiles: remove_file(f1); remove_file(f2) 

        # keep
        os.rename(all_reads_1_tmp, all_reads_1); os.rename(all_reads_2_tmp, all_reads_2)


    # align them with bwa mem, (or maybe minimap2)
    bamfile = "%s/wgsim_%s_aligned_reads.bam"%(working_dir, ID)
    sorted_bam = "%s.sorted"%bamfile
    index_bam = "%s.bai"%sorted_bam

    run_bwa_mem(all_reads_1, all_reads_2, reference_genome, working_dir, bamfile, sorted_bam, index_bam, name_sample="wgsim_%s"%ID, threads=threads, replace=replace)

    # plot the coverage across genome
    plot_coverage_across_genome_pairedEndReads(sorted_bam, reference_genome, replace=replace)

    ############################################################################

    ######### LONG READS SIMULATION ##########
    
    """
    # define the output
    ID = "fractChr%.5f_fractWind%.5f"%(fraction_chr_per_read, fraction_chr_window_increment)
    all_reads_1 = "%s/uniformSim_%s_all_reads1.fasta"%(working_dir, ID)
    all_reads_2 = "%s/uniformSim_%s_all_reads2.fasta"%(working_dir, ID)

    if any([file_is_empty(x) for x in [all_reads_1, all_reads_2]]) or replace is True:
        print("simulating long reads")

        # first simulate reads in parallel for each chromosome, as if they were paired
        inputs_fn = [(chr_obj, fraction_chr_per_read, fraction_chr_window_increment) for chr_obj in SeqIO.parse(query_genome, "fasta")]
        with multiproc.Pool(multiproc.cpu_count()) as pool:
            all_reads_list_tuples = pool.starmap(simulate_pairedEndReads_per_chromosome_uniform, inputs_fn)
            pool.close()

        #all_reads_list_tuples = list(map(lambda x: simulate_pairedEndReads_per_chromosome_uniform(x[0], x[1], x[2]), inputs_fn))

        # get the reads writen
        SeqIO.write(make_flat_listOflists([r[0] for r in all_reads_list_tuples]), all_reads_1, "fasta")
        SeqIO.write(make_flat_listOflists([r[1] for r in all_reads_list_tuples]), all_reads_2, "fasta")

    # run bwa mem
    bamfile = "%s/uniformSim_%s_aligned_reads.bam"%(working_dir, ID)
    sorted_bam = "%s.sorted"%bamfile
    index_bam = "%s.bai"%sorted_bam

    run_bwa_mem(all_reads_1, all_reads_2, reference_genome, working_dir, bamfile, sorted_bam, index_bam, name_sample="uniformSim_%s"%ID, threads=threads, replace=replace)

    # plot the coverage across genome
    plot_coverage_across_genome_pairedEndReads(sorted_bam, reference_genome, replace=replace)
    """

    ########################################



def get_SVtables_from_SVIM_to_RSVSim(svim_svtype_to_SVtable, sorted_bam, outdir, replace=False, threads=4):

    """Takes a dict of SV variation as output by SVIM and writes SVtables in a way that is equivalent as RSVSim, written under outdir. The sorted bam should span the deletions"""

    make_folder(outdir)

    # define all the svtypes as dfs
    svim_svtype_to_df = {svtype: pd.read_csv(SVtable, sep="\t") for svtype, SVtable in svim_svtype_to_SVtable.items()}
    all_svim_svtypes = set(svim_svtype_to_SVtable)

    # add the name of each event
    for svtype, df in svim_svtype_to_df.items(): 
        df["name_event"] = ["%s_%i"%(svtype, I+1) for I in range(len(df))]
        print(svtype)

    ldahflh

    # ren

    # calculate the coverage of all the regions
    all_regions_df

    ##### INSERTIONS ####
    if "int_duplications_source" in all_svim_svtypes:

        # define the dataframes
        df_from = svim_svtype_to_df["int_duplications_source"]
        df_to = svim_svtype_to_df["int_duplications_dest"]



        print(df_from, df_to)

        # defnie the end and the start and the end of each in


    # load the 


    ####################

    afkjfjkdakjdf

def run_svim_for_reads(reads, reference_genome, outdir,  threads=4, replace=False, min_sv_size=50, max_sv_size=1000000000000000000000, aligner="ngmlr", read_lengths=[kb*1000 for kb in [0.1, 0.3, 0.5, 0.7, 1, 2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 100000000]], coverage=50, max_fraction_chromosome=1.0, min_read_len=100, is_nanopore=True):

    """Takes some reads and a reference genome and runs svim"""

    # define the outfiles
    svType_to_file = {svtype : "%s/candidates/candidates_%s.corrected.bed"%(outdir, svtype) for svtype in {"breakends", "deletions", "int_duplications_dest", "int_duplications_source", "inversions", "novel_insertions", "tan_duplications_dest", "tan_duplications_source"}}

    # get the name of the sorted bam
    sorted_bam_long = "%s/%s.%s.coordsorted.bam"%(outdir, get_file(reads).rstrip(".fasta"), aligner)
    sorted_bam_long_idx = "%s.bai"%sorted_bam_long

    # change the name so that it is shorter, this is good for making further folders
    sorted_bam_short = "%s/aligned_reads.sorted.bam"%outdir
    sorted_bam_short_idx = "%s.bai"%sorted_bam_short

    if any([not os.path.isfile(x) for x in svType_to_file.values()]) or file_is_empty(sorted_bam_short) or file_is_empty(sorted_bam_short_idx) or replace is True:
    #if True:
     
        # make the folder
        delete_folder(outdir); make_folder(outdir)

        # run svim with few filters
        svim_cmd = "%s reads %s %s %s --min_sv_size %i --max_sv_size %i --cores %i --aligner %s --minimum_depth 1 --min_mapq 0 --skip_genotyping"%(svim, outdir, reads, reference_genome, min_sv_size, max_sv_size, threads, aligner)
        if is_nanopore is True: svim_cmd += " --nanopore"
        run_cmd(svim_cmd)
        
        os.rename(sorted_bam_long, sorted_bam_short)
        os.rename(sorted_bam_long_idx, sorted_bam_short_idx)

        # plotting the coverage
        plot_coverage_across_genome_pairedEndReads(sorted_bam_short, reference_genome, replace=replace, window_l=5000)

        #### ADD HEADER TO TABLES ####

        # define the column names
        col3_Bnd_IntDup_TanDup = "svtype;partner_dest;std_pos_across_cluster;std_span_across_cluster"
        col3_Del_Inv_Ins = "svtype;std_pos_across_cluster;std_span_across_cluster"

        svtype_to_col3_name = {"breakends":col3_Bnd_IntDup_TanDup, "deletions":col3_Del_Inv_Ins, "int_duplications_dest":col3_Bnd_IntDup_TanDup, "int_duplications_source":col3_Bnd_IntDup_TanDup, "inversions":col3_Del_Inv_Ins, "novel_insertions":col3_Del_Inv_Ins, "tan_duplications_dest":col3_Bnd_IntDup_TanDup, "tan_duplications_source":col3_Bnd_IntDup_TanDup}

        colnamesDict_InsDelTanInv = {0:"Chr", 1:"Start", 2:"End", 4:"score", 5:"evidence_deleted_origin", 6:"signatures_making_this_candidate"}
        colnamesDict_Bnd = {0:"Chr", 1:"Start", 2:"End", 4:"score", 5:"signatures_making_this_candidate"}

        # rewrite the candidates adding header
        candidates_dir = "%s/candidates"%outdir
        for file in os.listdir(candidates_dir):
            svtype = file.split(".")[0].split("candidates_")[1]

            # define the colnames for this svtype
            if svtype=="breakends": colnames_dict = colnamesDict_Bnd
            else: colnames_dict = colnamesDict_InsDelTanInv
            colnames_dict[3] = svtype_to_col3_name[svtype]

            # get the df
            filename = "%s/%s"%(candidates_dir, file)
            df = pd.read_csv(filename, sep="\t", header=-1).rename(columns=colnames_dict)

            # write
            df.to_csv(svType_to_file[svtype], sep="\t", index=False, header=True)

        ############################

    return svType_to_file, sorted_bam_short

def generate_tables_of_SV_between_genomes_svim(query_genome, reference_genome, replace=False, min_sv_size=50, max_sv_size=1000000000000000000000, threads=4, aligner="ngmlr", read_lengths=[kb*1000 for kb in [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1, 2, 3, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 100, 100000000]], coverage=100, max_fraction_chromosome=1.0, min_read_len=100, type_return="RSVSim"):

    """This function gets the structural variants with svim simulating reads of 10kb of the query genome, and the running svim to get structural variations """

    # initialize time
    start_time = time.time()

    # first run svim under the outdir of the aligned reads
    working_dir = "%s/svim_ouptut_against_%s"%(get_dir(query_genome), get_file(reference_genome))
    if replace is True: delete_folder(working_dir); 
    make_folder(working_dir)

    ############ READS SIMULATION ###########
    
    # define simulation parms
    def get_int_or_float_as_text(number):

        if int(number)==float(number): return "%i"%number
        else: return "%.2f"%number

    simulated_reads = "%s/%s_simulated_long_reads_readLengths%skb_coverage%ix_maxFractChrom%.3f_minReadLen%ibp.fasta"%(working_dir, get_file(query_genome), "-".join([get_int_or_float_as_text(x/1000) for x in read_lengths]), coverage, max_fraction_chromosome, min_read_len)

    if file_is_empty(simulated_reads) or replace is True:

        print("simulating reads")

        # first simulate reads in parallel for each chromosome
        inputs_fn = [(chr_obj, read_lengths, coverage, max_fraction_chromosome, min_read_len) for chr_obj in SeqIO.parse(query_genome, "fasta")]
        with multiproc.Pool(multiproc.cpu_count()) as pool:
            all_reads_list = pool.starmap(simulate_longReads_per_chromosome_uniform, inputs_fn)
            pool.close()

        print("writing reads")
        simulated_reads_tmp = "%s.tmp"%simulated_reads
        SeqIO.write(make_flat_listOflists(all_reads_list), simulated_reads_tmp, "fasta")
        os.rename(simulated_reads_tmp, simulated_reads)

        print("there are %i reads simulated"%(len(list(SeqIO.parse(simulated_reads, "fasta")))))
        
    ##########################################

    # RUN ALIGNMENT AND VAR CALLING
    outdir_svim = "%s.%s.svim_outdir"%(simulated_reads, aligner)
    svim_svType_to_file, sorted_bam_svim = run_svim_for_reads(simulated_reads, reference_genome, outdir_svim,  threads=threads, replace=replace, aligner=aligner, read_lengths=read_lengths, coverage=coverage, max_fraction_chromosome=max_fraction_chromosome, min_read_len=min_read_len, min_sv_size=min_sv_size, max_sv_size=max_sv_size, is_nanopore=True)

    print("---- It took %s seconds to generate the known vars for %s ----"%(time.time() - start_time, query_genome))

    # return the calculated elements if it is from svim
    if type_return=="svim": return svim_svType_to_file 

    # format as SVIM if indicated:
    elif type_return=="RSVSim":

        # get dict as in RSVim
        outdir_RSVSim = "%s/SVtables_RSVSim"%outdir_svim
        return get_SVtables_from_SVIM_to_RSVSim(svim_svType_to_file, sorted_bam_svim, outdir_RSVSim, replace=replace, threads=threads)   


########## SIMULATE LONG READS AS IF THEY WERE PAIRED ###############

"""
# define the simulated reads output
ID_all = "coverage%i_insertsize%i+-%i_error%.3f_readLengths%skb"%(coverage, insert_size, median_insert_size_sd, error_rate, "-".join([get_int_or_float_as_text(x/1000) for x in read_lengths]))
all_reads_1 = "%s/wgsim_%s_all_reads1.fq.gz"%(working_dir, ID_all)
all_reads_2 = "%s/wgsim_%s_all_reads2.fq.gz"%(working_dir, ID_all)

# define the genome length
genome_len = sum([len(c) for c in SeqIO.parse(query_genome, "fasta")])

# initialize the tuples of reads that will be removed
all_pairedReads_tupleFiles = []
all_reads_1_individualFiles = []
all_reads_2_individualFiles = []

# generate the reads
if any([file_is_empty(x) for x in [all_reads_1, all_reads_2]]) or replace is True:

    # go through each read length
    for read_length in read_lengths:

        # define the number of readPairs
        readPairs = int(((genome_len/read_length)*coverage)/2)
        print("simulating %i read pairs for coverage %i and read length %i"%(readPairs, coverage, read_length))
        
        # define the inputs of wgsim            
        inputs_fn = [(query_genome, chr_obj.id, 0, len(chr_obj), readPairs, read_length, working_dir,  median_insert_size, median_insert_size_sd, replace, error_rate) for chr_obj in SeqIO.parse(query_genome, "fasta")]

        with multiproc.Pool(multiproc.cpu_count()) as pool:
            all_pairedReads_tupleFiles += pool.starmap(run_wgsim_pairedEnd_for_window, inputs_fn)
            pool.close()

        # keep
        all_reads_1_individualFiles += [x[0] for x in all_pairedReads_tupleFiles]
        all_reads_2_individualFiles += [x[1] for x in all_pairedReads_tupleFiles]

    # concatenate all reads into one
    print("concatenating all reads in one")
    all_reads_1_tmp = "%s.tmp"%all_reads_1; all_reads_2_tmp = "%s.tmp"%all_reads_2;
    run_cmd("cat %s > %s"%(working_dir, " ".join(all_reads_1_individualFiles), all_reads_1_tmp)) 
    run_cmd("cat %s > %s"%(working_dir, " ".join(all_reads_2_individualFiles), all_reads_2_tmp)) 

    # remove all the reads
    print("Removing")
    for f1, f2 in all_pairedReads_tupleFiles: remove_file(f1); remove_file(f2) 

    # keep
    os.rename(all_reads_1_tmp, all_reads_1); os.rename(all_reads_2_tmp, all_reads_2)

# align them with bwa mem, (or maybe minimap2)
bamfile = "%s/wgsim_%s_aligned_reads.bam"%(working_dir, ID)
sorted_bam = "%s.sorted"%bamfile
index_bam = "%s.bai"%sorted_bam

run_bwa_mem(all_reads_1, all_reads_2, reference_genome, working_dir, bamfile, sorted_bam, index_bam, name_sample="wgsim_%s"%ID, threads=threads, replace=replace)

# plot the coverage across genome
plot_coverage_across_genome_pairedEndReads(sorted_bam, reference_genome, replace=replace)

"""

############################################################################

def simulate_longReads_per_chromosome_uniform(chr_obj, read_lengths = [1000, 5000, 10000, 15000, 20000], coverage=100, max_fraction_chromosome=0.2, min_read_len=500):

    """This function takes a chromosome object (SeqRecord) and it generates reads that are as long as read_length, covering the whole chromosome. This works in a way that the coverage is equally achieved by each of the read lengths"""

    # define metrics
    len_chromosome = len(chr_obj)
    print("working on chromosome %s"%chr_obj.id)

    # define the read_lenghts to investigate
    read_lengths = [int(x) for x in sorted(read_lengths + [int(len_chromosome*max_fraction_chromosome)]) if x<=int(len_chromosome*max_fraction_chromosome)]

    # define the coverage that has to be achieved by each length of reads
    per_readLenth_coverage = int(coverage/len(read_lengths))

    # define a function that takes a start and a length and returns the read
    def get_read(start, length, chr_obj, extraID):

        # define the end of the read
        end = start + length

        if (len_chromosome-start)<=0: raise ValueError("There was an error with the read generation")

        # keep the read 
        read = chr_obj[start:end]
        read.id = "%s_readStart%i_%s"%(chr_obj.id, start, extraID); read.name = ""; read.description = ""

        return read

    # initialize reads
    all_reads = []

    # go through each read_length
    for read_len in read_lengths:

        all_reads_readLen = []

        # define the size of the increasing window that has to be reached to generate the per_readLenth_coverage
        len_window_increment = int( len_chromosome / ((len_chromosome/read_len)*per_readLenth_coverage) )
        #print("The len of the window is %i for read length %i"%(len_window_increment, read_len))

        # define the start
        start = 0

        # go through each window
        while (len_chromosome - start) > min_read_len : 

            # get the read
            read = get_read(start, read_len, chr_obj, "readLen%i_intrachromosomal_read"%read_len)
            all_reads_readLen.append(read)

            # move the window
            start += len_window_increment

        # add the reads for the beginning of the chr, which would be otherwise uncovered
        for I in range(per_readLenth_coverage): all_reads_readLen.append(get_read(0, read_len, chr_obj, "readLen%i_fromStart%i"%(I, read_len)))

        # keep 
        all_reads += all_reads_readLen  

        #print("Simulating %i reads, which yields a coverage of %i. The expected coverage is %i"%(len(all_reads_readLen), (len(all_reads_readLen)*read_len/len_chromosome), per_readLenth_coverage))

    return all_reads

    def get_coverage_per_window_for_chromosomeDF(chromosome_id, destination_dir, windows_bed, sorted_bam, replace, window_l):

    """Takes a chromosome id, a destination dir where to write files, a windows file (provided by generate_coverage_per_window_file_parallel) and a sorted bam and it generates a dataframe with the coverage stats"""

    # define the output coverage file
    outfile = "%s/%s_coverage_windows%ibp.tab"%(destination_dir, chromosome_id, window_l); outfile_tmp = "%s.tmp"%outfile
    print(outfile)
    print("running bamstats for %s"%chromosome_id)
        
    # generate a randomID
    randID = id_generator(25)

    # define a file for the coverage
    windows_bed_chromsome = "%s.%s.%s.bed"%(windows_bed, chromosome_id, randID)
    run_cmd("grep $'%s\t' %s > %s"%(chromosome_id, windows_bed, windows_bed_chromsome))

    # if there is nothing, return an empty df
    bamstats_fields = ["#chrom", "start", "end", "length", "sample", "mincov", "maxcov", "avgcov_1", "mediancov_1", "nocoveragebp_1", "percentcovered_1"]
    if file_is_empty(windows_bed_chromsome): return pd.DataFrame(columns=bamstats_fields)

    # remove previously existing files
    #if file_is_empty(outfile) or replace is True:
    if True:

        # generate the bam file for this chromosome (and index)
        sorted_bam_chr = "%s.%s.bam"%(sorted_bam, chromosome_id)
        sorted_bam_chr_index = "%s.bai"%sorted_bam_chr
        if file_is_empty(sorted_bam_chr) or file_is_empty(sorted_bam_chr_index) or replace is True: 

            # get bam for chr
            sorted_bam_chr_tmp = "%s.%s"%(sorted_bam_chr, randID); remove_file(sorted_bam_chr_tmp)
            run_cmd("%s view -b %s %s > %s"%(samtools, sorted_bam, chromosome_id, sorted_bam_chr_tmp))

            # get the index for the tmp (generates a tmp.bai file)
            run_cmd("%s index -@ 1 %s"%(samtools, sorted_bam_chr_tmp))
            sorted_bam_chr_index_tmp = "%s.bai"%sorted_bam_chr_tmp; remove_file(sorted_bam_chr_index_tmp)

            os.rename(sorted_bam_chr_tmp, sorted_bam_chr)
            os.rename(sorted_bam_chr_index_tmp, sorted_bam_chr_index)

        # generate a df with the coverage statistics using bedtools genomecov
        genomecov_outdf = "%s.coverage_per_positionDF.py"%sorted_bam_chr
        if file_is_empty(genomecov_outdf) or replace is True: 

            # run genomecov to get df of coverage
            coverage_per_position_bed = "%s.coverage_per_positionDF.bed"%sorted_bam_chr; coverage_per_position_bed_tmp = "%s.tmp.%s"%(coverage_per_position_bed, randID)
            if file_is_empty(coverage_per_position_bed) or replace is True:
                run_cmd("%s genomecov -ibam %s -dz | grep '%s' | cut -f2,3 > %s "%(bedtools, sorted_bam_chr, chromosome_id, coverage_per_position_bed_tmp)) # it does not output for all the positions in the genome
                os.rename(coverage_per_position_bed_tmp, coverage_per_position_bed)

            # get the length of the chromosome
            chr_len_line = "%s.chromosome_length.tab"%sorted_bam_chr
            run_cmd("%s view -h %s | grep -m 1 $'@SQ\tSN:%s' > %s"%(samtools, sorted_bam_chr, chromosome_id, chr_len_line))
            chr_len = int(open(chr_len_line, "r").readlines()[0].split("\t")[2].split(":")[1])
            all_positions = set(list(range(chr_len)))

            # load into df 
            df_coverage_per_position = pd.read_csv(coverage_per_position_bed, sep="\t", header=-1, names=["position", "coverage"])

            # add the missing positions
            all_positions_covered = set(df_coverage_per_position.position)
            all_positions_uncovered = list(all_positions.difference(all_positions_covered))
            df_coverage_per_position_missing = pd.DataFrame({"position":all_positions_uncovered, "coverage":[0]*len(all_positions_uncovered)})
            df_coverage_per_position = df_coverage_per_position.append(df_coverage_per_position_missing).sort_values(by="position")

            # save
            genomecov_outdf_tmp = "%s.tmp.%s"%(genomecov_outdf, randID)
            save_object(df_coverage_per_position, genomecov_outdf_tmp)
            os.rename(genomecov_outdf_tmp, genomecov_outdf)

        else: df_coverage_per_position = load_object(genomecov_outdf)

        # generate a df that has all the bamstats_fields
        df_windows = pd.read_csv(windows_bed_chromsome, sep="\t", header=-1, names=["chromosome", "start", "end", "ID"])

        # add other trivial measurements
        df_windows["length"] = df_windows.end - df_windows.start
        df_windows["sample"] = ["coverage_calc"]*len(df_windows)
        df_windows = df_windows.rename(columns={"chromosome":"#chrom"})

        def get_statistics(r):

            """Takes a window with start and end and returns the statistics"""

            # define an array where each of the positions is a coverage for one position in the window
            cov_list = df_coverage_per_position[(df_coverage_per_position.position>=r["start"]) & (df_coverage_per_position.position<=r["end"])].coverage.values

            # define uncovered regions
            nocoveragebp_1 = int((cov_list<=1).sum()) # the number of uncovered bp by more than 1 read
            percentcovered_1 = np.divide(nocoveragebp_1, r["length"])*100

            # return as a series
            return pd.Series({"mincov":cov_list.min(), "maxcov":cov_list.max(), "avgcov_1": np.mean(cov_list), "mediancov_1": np.median(cov_list), "nocoveragebp_1": nocoveragebp_1, "percentcovered_1": percentcovered_1})


        # add coverage statistcis
        df_windows[["mincov", "maxcov", "avgcov_1", "mediancov_1", "nocoveragebp_1", "percentcovered_1"]] = df_windows.apply(get_statistics, axis=1)

        df_windows = df_windows[bamstats_fields]

        return df_windows
        print(df_windows)

        khagfagk




        # run the bam stats into outfile_tmp, try it 10 times or quit
        for I in range(10):

            remove_file(outfile_tmp)

            # run bamstats, and keep whether it worked
            try:
                bamstats_cmd = "%s -jar %s -B %s -cov 1 %s > %s"%(JAVA, bamstats04_jar, windows_bed_chromsome, sorted_bam_chr, outfile_tmp)
                run_cmd(bamstats_cmd)
                correct_bamstats_run = True

            except: correct_bamstats_run = False

            # check that the bamstats are as expected
            try: 
                check_df = pd.read_csv(outfile_tmp, sep="\t")[bamstats_fields] # tgese are some fields that should be in the outfile
                correct_bamstats_out = True
            
            except: correct_bamstats_out = False 

            # decide whether to try again
            if correct_bamstats_run and correct_bamstats_out: break
            else: continue

        correct_bamstats = (correct_bamstats_run and correct_bamstats_out)
        if correct_bamstats is False: raise ValueError("%s did not work"%bamstats_cmd)
            
        # remove files and rename
        #os.unlink(sorted_bam_chr); os.unlink("%s.bai"%sorted_bam_chr)
        os.rename(outfile_tmp, outfile)

    # remove
    os.unlink(windows_bed_chromsome)

    return pd.read_csv(outfile, sep="\t")




            print(sys.getsizeof(series))


            positions_array = np.array(series.index)

            print(series)

            

            df_coverage_per_position = pd.DataFrame()
            df_coverage_per_position["position"] = series.index
            df_coverage_per_position["coverage"] = series

            print(sys.getsizeof(df_coverage_per_position))

            lashasldhsad


            print(df_coverage_per_position)


            dfgfghfhd
            print("hu")
            df_coverage_per_position = pd.read_csv(coverage_per_position_bed, sep="\t", header=None, names=["position", "coverage"], dtype={"position":np.int16, "coverage":np.int16})

            # add the missing positions
            print("adding missing positions")
            all_positions_covered = set(df_coverage_per_position.position)
            all_positions_uncovered = list(all_positions.difference(all_positions_covered))

            print("geting the missing positions")
            df_coverage_per_position_missing = pd.DataFrame({"position":all_positions_uncovered, "coverage":[0]*len(all_positions_uncovered)})
            df_coverage_per_position = df_coverage_per_position.append(df_coverage_per_position_missing).sort_values(by="position")

            # remove the file
            remove_file(coverage_per_position_bed)

            # save
            print("saving df")
            genomecov_outdf_tmp = "%s.tmp.%s"%(genomecov_outdf, randID)
            save_object(df_coverage_per_position, genomecov_outdf_tmp)
            os.rename(genomecov_outdf_tmp, genomecov_outdf)



            lsfkhfldhk

            print(genomecov_out_regionsDF)

            ljfahkfjakfeag

            # get the length of the chromosome
            chr_len_line = "%s.chromosome_length.tab"%sorted_bam_chr
            run_cmd("%s view -H %s -@ 0 | grep -m 1 $'@SQ\tSN:%s' > %s"%(samtools, sorted_bam_chr, chromosome_id, chr_len_line))
            chr_len = int(open(chr_len_line, "r").readlines()[0].split("\t")[2].split(":")[1])
            all_positions = set(list(range(chr_len)))

            # load into dict 
            position_to_coverage = {p:c for p, c in map(lambda l: [int(x) for x in l.strip().split("\t")], open(coverage_per_position_bed, "r").readlines())}

            # add missing positions
            all_positions_covered = set(position_to_coverage)
            all_positions_uncovered = all_positions.difference(all_positions_covered)
            for p in all_positions_uncovered: position_to_coverage[p] = 0

            # sort dict by key
            position_to_coverage = collections.OrderedDict(sorted(position_to_coverage.items()))

            # get as series
            coverage_series = pd.Series(position_to_coverage)

            # remove the file
            remove_file(coverage_per_position_bed)

            # save
            genomecov_out_series_tmp = "%s.tmp.%s"%(genomecov_out_series, randID)

            save_object(coverage_series, genomecov_out_series_tmp)
            os.rename(genomecov_out_series_tmp, genomecov_out_series)

        else: coverage_series = load_object(genomecov_out_series)


        def get_statistics(r):

            """Takes a window with start and end and returns the statistics"""
            
            print("\n\n")
            print(r) ## debug here parallelization
            #print(coverage_series.index)
            #print((coverage_series.index>=r["start"]))
            #print((coverage_series.index<=r["end"]))


            

            # define an array where each of the positions is a coverage for one position in the window
            print(sys.getsizeof(coverage_series))
            positions_larger_than_start = coverage_series.index[(coverage_series.index>=r["start"])]
            print("positions larger than start", positions_larger_than_start)
            #print(coverage_series.index[((coverage_series.index>=r["start"]) & (coverage_series.index<=r["end"]))])

            cov_list = coverage_series.index[(coverage_series.index<=r["end"])].values

            #print(cov_list, coverage_series)

            # define uncovered regions
            nocoveragebp_1 = int((cov_list<=1).sum()) # the number of uncovered bp by more than 1 read
            percentcovered_1 = 100 - (np.divide(nocoveragebp_1, r["length"])*100)

            # return as a series
            return pd.Series({"mincov":cov_list.min(), "maxcov":cov_list.max(), "avgcov_1": np.mean(cov_list), "mediancov_1": np.median(cov_list), "nocoveragebp_1": nocoveragebp_1, "percentcovered_1": percentcovered_1})




            
def get_coverage_per_window_for_chromosomeDF(chromosome_id, destination_dir, windows_bed, sorted_bam, replace, window_l):

    """Takes a chromosome id, a destination dir where to write files, a windows file (provided by generate_coverage_per_window_file_parallel) and a sorted bam and it generates a dataframe with the coverage stats"""

    # define the output coverage file
    outfile = "%s/%s_coverage_windows%ibp.tab"%(destination_dir, chromosome_id, window_l); outfile_tmp = "%s.tmp"%outfile
    print("running bamstats for %s"%chromosome_id)
        
    # generate a randomID
    randID = id_generator(25)

    # define a file for the coverage
    windows_bed_chromsome = "%s.%s.%s.bed"%(windows_bed, chromosome_id, randID)
    run_cmd("grep $'%s\t' %s > %s"%(chromosome_id, windows_bed, windows_bed_chromsome))

    # if there is nothing, return an empty df
    bamstats_fields = ["#chrom", "start", "end", "length", "sample", "mincov", "maxcov", "avgcov_1", "mediancov_1", "nocoveragebp_1", "percentcovered_1"]
    if file_is_empty(windows_bed_chromsome): return pd.DataFrame(columns=bamstats_fields)

    # define the output bams
    sorted_bam_chr = "%s.%s.bam"%(sorted_bam, chromosome_id)
    sorted_bam_chr_index = "%s.bai"%sorted_bam_chr

    # define the filename of a line that will contain the line of the bam that has the length of the chromosome
    chr_len_line = "%s.chromosome_length.tab"%sorted_bam_chr

    # remove previously existing files
    if file_is_empty(outfile) or replace is True:

        # generate the bam file for this chromosome (and index)
        if file_is_empty(sorted_bam_chr) or file_is_empty(sorted_bam_chr_index) or replace is True: 

            # get bam for chr
            sorted_bam_chr_tmp = "%s.%s"%(sorted_bam_chr, randID); remove_file(sorted_bam_chr_tmp)
            run_cmd("%s view -b %s %s > %s"%(samtools, sorted_bam, chromosome_id, sorted_bam_chr_tmp))

            # get the index for the tmp (generates a tmp.bai file)
            sorted_bam_chr_index_tmp = "%s.bai"%sorted_bam_chr_tmp; remove_file(sorted_bam_chr_index_tmp)
            run_cmd("%s index -@ 1 %s"%(samtools, sorted_bam_chr_tmp))

            os.rename(sorted_bam_chr_tmp, sorted_bam_chr)
            os.rename(sorted_bam_chr_index_tmp, sorted_bam_chr_index)

        # generate a df with the coverage statistics using bedtools genomecov
        genomecov_out_series = "%s.coverage_per_positionSeries.py"%sorted_bam_chr
        if file_is_empty(genomecov_out_series) or replace is True: 

            # run genomecov to get df of coverage
            coverage_per_position_bed = "%s.coverage_per_position.bed"%sorted_bam_chr; coverage_per_position_bed_tmp = "%s.tmp.%s"%(coverage_per_position_bed, randID)
            if file_is_empty(coverage_per_position_bed) or replace is True:
                run_cmd("%s genomecov -ibam %s -dz | grep '%s' | cut -f2,3 > %s "%(bedtools, sorted_bam_chr, chromosome_id, coverage_per_position_bed_tmp)) # it does not output for all the positions in the genome
                os.rename(coverage_per_position_bed_tmp, coverage_per_position_bed)

            # get the length of the chromosome
            run_cmd("%s view -H %s -@ 0 | grep -m 1 $'@SQ\tSN:%s' > %s"%(samtools, sorted_bam_chr, chromosome_id, chr_len_line))
            chr_len = int(open(chr_len_line, "r").readlines()[0].split("\t")[2].split(":")[1])
            all_positions = set(list(range(chr_len)))

            # load into dict 
            position_to_coverage = {p:c for p, c in map(lambda l: [int(x) for x in l.strip().split("\t")], open(coverage_per_position_bed, "r").readlines())}

            # add missing positions
            all_positions_covered = set(position_to_coverage)
            all_positions_uncovered = all_positions.difference(all_positions_covered)
            for p in all_positions_uncovered: position_to_coverage[p] = 0

            # sort dict by key
            position_to_coverage = collections.OrderedDict(sorted(position_to_coverage.items()))

            # get as series
            coverage_series = pd.Series(position_to_coverage); del position_to_coverage

            # remove the file
            remove_file(coverage_per_position_bed)

            # save
            genomecov_out_series_tmp = "%s.tmp.%s"%(genomecov_out_series, randID)

            save_object(coverage_series, genomecov_out_series_tmp)
            os.rename(genomecov_out_series_tmp, genomecov_out_series)

        else: coverage_series = load_object(genomecov_out_series)

        # generate a df that has all the bamstats_fields
        df_windows = pd.read_csv(windows_bed_chromsome, sep="\t", header=None, names=["chromosome", "start", "end", "ID"])

        # add other trivial measurements
        df_windows["length"] = df_windows.end - df_windows.start
        df_windows["sample"] = ["coverage_calc"]*len(df_windows)
        df_windows = df_windows.rename(columns={"chromosome":"#chrom"})

        # add coverage statistcis
        df_windows[["mincov", "maxcov", "avgcov_1", "mediancov_1", "nocoveragebp_1", "percentcovered_1"]] = df_windows.apply(lambda r: get_statistics_of_coverage(r, coverage_series), axis=1)
        del coverage_series

        # write
        df_windows[bamstats_fields].to_csv(outfile_tmp, sep="\t", header=True, index=False)
        os.rename(outfile_tmp, outfile)

    else: df_windows = pd.read_csv(outfile, sep="\t")

    # remove unnecessary files
    print("removing files")
    for f in [sorted_bam_chr, sorted_bam_chr_index, windows_bed_chromsome, chr_len_line]: remove_file(f)

    return df_windows[bamstats_fields]



print(var_to_GT)
khghkdsfd




# add for multialleles
for var in df.index:
    chrom, pos, ref, alt = var
    filterTag = df.loc[{var}, "FILTER"][0]

    # get all alternative alleles
    alt_alleles = alt.split(",")

    # get all the alternative frequencies
    altAllele_to_freq = dict(zip(alt_alleles, df.loc[{var}, "alternative_allelle_frequencies"][0]))

    # get all the GT tags
    altAllele_to_GT = dict(zip(alt_alleles, df.loc[{var}, "GT"][0].split(",")))

    # map them to frequs and tags
    for real_alt, freq in altAllele_to_freq.items(): 
        var_to_frequency[(chrom, pos, ref, real_alt)] = freq
        var_to_filter[(chrom, pos, ref, real_alt)] = filterTag
        var_to_GT[(chrom, pos, ref, real_alt)] = altAllele_to_GT[real_alt]

# add alternative frequencies for different representations of the variants
print("adding for different representations of the vars")
for var, freq in cp.deepcopy(var_to_frequency).items():

    # change
    chrom, pos, ref, alt = var
    mod_pos, mod_ref, mod_alt = leftTrimVariant(pos, ref, alt)

    # keep if changed
    if (mod_pos, mod_ref, mod_alt)!=(pos, ref, alt): 
        var_to_frequency[(chrom, mod_pos, mod_ref, mod_alt)] = freq
        var_to_filter[(chrom, mod_pos, mod_ref, mod_alt)] = var_to_filter[var]
        var_to_GT[(chrom, mod_pos, mod_ref, mod_alt)] = var_to_GT[var]


#### chromosome calculation at 01 january 2020:

def get_statistics_of_coverage(r, coverage_series):

    """Takes a window with start and end and returns the statistics"""

    positions = coverage_series.index[(coverage_series.index>=r["start"]) & (coverage_series.index<=r["end"])]
    cov_list = coverage_series[positions].values

    # define uncovered regions
    nocoveragebp_1 = int((cov_list<=1).sum()) # the number of uncovered bp by more than 1 read
    percentcovered_1 = 100 - (np.divide(nocoveragebp_1, r["length"])*100)

    # define basic stats
    mincov = cov_list.min()
    maxcov = cov_list.max()
    avgcov_1 = np.mean(cov_list)
    mediancov_1 = np.median(cov_list)

    # delete objects
    del cov_list

    # return as a series
    return pd.Series({"mincov":mincov, "maxcov":maxcov, "avgcov_1":avgcov_1, "mediancov_1": mediancov_1, "nocoveragebp_1": nocoveragebp_1, "percentcovered_1": percentcovered_1})



def get_coverage_per_window_for_chromosomeDF(chromosome_id, destination_dir, windows_bed, sorted_bam, replace, window_l):

    """Takes a chromosome id, a destination dir where to write files, a windows file (provided by generate_coverage_per_window_file_parallel) and a sorted bam and it generates a dataframe with the coverage stats"""

    # define the output coverage file
    outfile = "%s/%s_coverage_windows%ibp.tab"%(destination_dir, chromosome_id, window_l); outfile_tmp = "%s.tmp"%outfile
    print("running bamstats for %s"%chromosome_id)
        
    # generate a randomID
    randID = id_generator(25)

    # define a file for the coverage
    windows_bed_chromsome = "%s.%s.%s.bed"%(windows_bed, chromosome_id, randID)
    run_cmd("grep $'%s\t' %s > %s"%(chromosome_id, windows_bed, windows_bed_chromsome))

    # if there is nothing, return an empty df
    bamstats_fields = ["#chrom", "start", "end", "length", "sample", "mincov", "maxcov", "avgcov_1", "mediancov_1", "nocoveragebp_1", "percentcovered_1"]
    if file_is_empty(windows_bed_chromsome): return pd.DataFrame(columns=bamstats_fields)

    # define the output bams
    sorted_bam_chr = "%s.%s.bam"%(sorted_bam, chromosome_id)
    sorted_bam_chr_index = "%s.bai"%sorted_bam_chr

    # define the filename of a line that will contain the line of the bam that has the length of the chromosome
    chr_len_line = "%s.chromosome_length.tab"%sorted_bam_chr

    # remove previously existing files
    if file_is_empty(outfile) or replace is True:

        # generate the bam file for this chromosome (and index)
        if file_is_empty(sorted_bam_chr) or file_is_empty(sorted_bam_chr_index) or replace is True: 

            # get bam for chr
            sorted_bam_chr_tmp = "%s.%s"%(sorted_bam_chr, randID); remove_file(sorted_bam_chr_tmp)
            run_cmd("%s view -b %s %s > %s"%(samtools, sorted_bam, chromosome_id, sorted_bam_chr_tmp))

            # get the index for the tmp (generates a tmp.bai file)
            sorted_bam_chr_index_tmp = "%s.bai"%sorted_bam_chr_tmp; remove_file(sorted_bam_chr_index_tmp)
            run_cmd("%s index -@ 1 %s"%(samtools, sorted_bam_chr_tmp))

            os.rename(sorted_bam_chr_tmp, sorted_bam_chr)
            os.rename(sorted_bam_chr_index_tmp, sorted_bam_chr_index)

        # generate a df with the coverage statistics using bedtools genomecov
        genomecov_out_series = "%s.coverage_per_positionSeries.py"%sorted_bam_chr
        if file_is_empty(genomecov_out_series) or replace is True: 
            print("calculating coverage per position")

            # run genomecov to get df of coverage
            coverage_per_position_bed = "%s.coverage_per_position.bed"%sorted_bam_chr; coverage_per_position_bed_tmp = "%s.tmp.%s"%(coverage_per_position_bed, randID)
            if file_is_empty(coverage_per_position_bed) or replace is True:
                run_cmd("%s genomecov -ibam %s -dz | grep '%s' | cut -f2,3 > %s "%(bedtools, sorted_bam_chr, chromosome_id, coverage_per_position_bed_tmp)) # it does not output for all the positions in the genome
                os.rename(coverage_per_position_bed_tmp, coverage_per_position_bed)

            # get the length of the chromosome
            run_cmd("%s view -H %s -@ 0 | grep -m 1 $'@SQ\tSN:%s' > %s"%(samtools, sorted_bam_chr, chromosome_id, chr_len_line))
            chr_len = int(open(chr_len_line, "r").readlines()[0].split("\t")[2].split(":")[1])
            all_positions = set(list(range(chr_len)))

            # load into dict 
            position_to_coverage = {p:c for p, c in map(lambda l: [int(x) for x in l.strip().split("\t")], open(coverage_per_position_bed, "r").readlines())}

            # add missing positions
            all_positions_covered = set(position_to_coverage)
            all_positions_uncovered = all_positions.difference(all_positions_covered)
            for p in all_positions_uncovered: position_to_coverage[p] = 0

            # sort dict by key
            position_to_coverage = collections.OrderedDict(sorted(position_to_coverage.items()))

            # get as series
            coverage_series = pd.Series(position_to_coverage); del position_to_coverage

            # remove the file
            remove_file(coverage_per_position_bed)

            # save
            genomecov_out_series_tmp = "%s.tmp.%s"%(genomecov_out_series, randID)

            save_object(coverage_series, genomecov_out_series_tmp)
            os.rename(genomecov_out_series_tmp, genomecov_out_series)

        else: coverage_series = load_object(genomecov_out_series)

        # generate a df that has all the bamstats_fields
        print("coverage for windows obtention")
        df_windows = pd.read_csv(windows_bed_chromsome, sep="\t", header=None, names=["chromosome", "start", "end", "ID"])


        # add other trivial measurements
        df_windows["length"] = df_windows.end - df_windows.start
        df_windows["sample"] = ["coverage_calc"]*len(df_windows)
        df_windows = df_windows.rename(columns={"chromosome":"#chrom"})

        print("getting windows") #print(df_windows, coverage_series)

        # add coverage statistcis
        df_windows[["mincov", "maxcov", "avgcov_1", "mediancov_1", "nocoveragebp_1", "percentcovered_1"]] = df_windows.apply(lambda r: get_statistics_of_coverage(r, coverage_series), axis=1)

        print("windows got")

        # write
        df_windows[bamstats_fields].to_csv(outfile_tmp, sep="\t", header=True, index=False)
        os.rename(outfile_tmp, outfile)

    else: df_windows = pd.read_csv(outfile, sep="\t")

    # remove unnecessary files
    print("removing files")
    for f in [sorted_bam_chr, sorted_bam_chr_index, windows_bed_chromsome, chr_len_line]: remove_file(f)

    return df_windows[bamstats_fields]



#### COLORS

def closest_colour(requested_colour):
    min_colours = {}
    for key, name in webcolors.css3_hex_to_names.items():
        r_c, g_c, b_c = webcolors.hex_to_rgb(key)
        rd = (r_c - requested_colour[0]) ** 2
        gd = (g_c - requested_colour[1]) ** 2
        bd = (b_c - requested_colour[2]) ** 2
        min_colours[(rd + gd + bd)] = name
    return min_colours[min(min_colours.keys())]

def get_colour_name(requested_colour):
    try:
        closest_name = actual_name = webcolors.rgb_to_name(requested_colour)
    except ValueError:
        closest_name = closest_colour(requested_colour)
        actual_name = None
    return actual_name, closest_name

def hex_to_rgb(value):

    """Gets a hex and returns an RGB tuple"""

    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def hex_to_color(value): 

    """takes a hex value and returns the color"""

    # get RGB
    rgb = hex_to_rgb(value)

    #get color
    return get_colour_name(rgb)[1]

def define_colorbar(min, max, n=4, rgb_pos=0, type="btw_two_colors", color_from="red", color_to="blue"):
    
    """ Returns a dictionary in which the keys are a linspace array between min and max (with n steps) and the values are HEX color codes, progressively 
    rgb_pos indicates which of the rgb codes will change"""
    
    color_dict = {}
    
    if type=="btw_two_colors":
        
        # define the colors
        colors = list(Color(color_from).range_to(Color(color_to),n))
        
        for I, value in enumerate(np.linspace(min, max, n)):
            color_dict[value] = colors[I].get_hex()
        
    else:
        
        # this is a white to gray progression
        
        rgb_array = np.linspace(0, 255, n) # generate a linspace of the rgb color that will change
        rgb = [0, 0, 0] # this is going to be used
        color_dict = {}
        
        for I, value in enumerate(np.linspace(min, max, n)):
            
            # define the color
            rgb[0] = int(rgb_array[I]); rgb[1] = int(rgb_array[I]); rgb[2] = int(rgb_array[I]); # from black to gray
            #rgb[1] = int(rgb_array[I]); # from black to 
            
            hex_color = '#%02x%02x%02x'%(rgb[0], rgb[1], rgb[2])
            
            # add to dictionary
            color_dict[value] = hex_color
            
    return color_dict, list(color_dict.keys())

def get_colorbar(iterable, color_from="red", color_to="blue"): 
    
    """ function to take list of elements and map them to different colors """

    return {k:hex_to_color(v) for k,v in dict(zip(iterable, define_colorbar(1,2, n=len(iterable), color_from=color_from, color_to=color_to)[0].values())).items()}


#####


def get_clusters_overlapping_vars(svtype_to_svDF):

    """This function takes a dictionary that maps svtype to the df and returns a list of sets. Each of them represents a cluster of variants that are overlapping. For variants where there is more than one position it is considered overlapping if at least one of the boxes overlaps with the others. """

    print("getting clusters of overlapping vars")

    # define a function that changes the START by the END if it is -1
    def arrange_minus1_START(r):

        if r["Start"] == -1: return r["End"]
        else: return r["Start"]

    # initialize a df that will have Chr, Start, End, ID for all vars
    all_regions_df = pd.DataFrame()

    # go through each svtype
    for svtype, svDF in svtype_to_svDF.items():

        # simple types, only one box
        if svtype in {"deletions", "tandemDuplications", "inversions"}: all_regions_df = all_regions_df.append(svDF[["Chr", "Start", "End", "ID"]])

        # two boxes
        elif svtype in {"insertions", "translocations"}:

            # get two dfs
            dfA = svDF.rename(columns={"ChrA":"Chr", "StartA":"Start", "EndA":"End"})
            dfB = svDF.rename(columns={"ChrB":"Chr", "StartB":"Start", "EndB":"End"})

            # keep them
            for df in [dfA, dfB]: all_regions_df = all_regions_df.append(df[["Chr", "Start", "End", "ID"]])

        elif svtype=="remaining":
            
            # get the CHROMdf, which has POS as the same
            dfCHROM = svDF.rename(columns={"#CHROM":"Chr", "POS":"Start"}); dfCHROM["End"] = dfCHROM["Start"]

            # get the other df
            dfCHR2 = svDF.rename(columns={"CHR2":"Chr", "START":"Start", "END":"End"}); dfCHR2["Start"] = dfCHR2.apply(arrange_minus1_START, axis=1)

            # keep them
            for df in [dfCHROM, dfCHR2]: all_regions_df = all_regions_df.append(df[["Chr", "Start", "End", "ID"]])
            
        else: raise ValueError("%s is not a valid svtype"%svtype)

    # now identify clusters of overlapping IDs, and save the IDs

    # map each chr to it's regions df
    chrom_to_regionsDF = {chrom : all_regions_df[all_regions_df.Chr==chrom] for chrom in set(all_regions_df["Chr"])}

    # initialize vars
    list_clusters = [] # a list of sets of IDs, each of the corresponding to a cluster of SVs

    # go through each region
    for qChr, qStart, qEnd, qID in all_regions_df[["Chr", "Start", "End", "ID"]].values:

        # first check if the ID is already in any of the previosuly run clusters
        if any([qID in cluster for cluster in list_clusters]): continue

        # get the df that is of this region
        df = chrom_to_regionsDF[qChr]

        # find overlapping SVs
        overlappingIDs = set(df[~((df.End<qStart) | (df.Start>qEnd))]["ID"])

        # if there are overlaps add them to the cluster_list with itself
        if overlappingIDs!={qID} and len(overlappingIDs)>0: list_clusters.append(overlappingIDs.union({qID}))

    all_SVs = set(all_regions_df.ID)
    if len(list_clusters)>0: all_SVs_in_cluster = set.union(*list_clusters)
    else: all_SVs_in_cluster = set()

    print("There are %i clusters of SVs involving %i of %i SVs"%(len(list_clusters), len(all_SVs_in_cluster), len(all_SVs)))

    return list_clusters

def get_clusters_overlapping_vars(svtype_to_svDF):

    """This function takes a dictionary that maps svtype to the df and returns a list of sets. Each of them represents a cluster of variants that are overlapping. For variants where there is more than one position it is considered overlapping if at least one of the boxes overlaps with the others. """

    print("getting clusters of overlapping vars")

    # define a function that changes the START by the END if it is -1
    def arrange_minus1_START(r):

        if r["Start"] == -1: return r["End"]
        else: return r["Start"]

    # initialize a df that will have Chr, Start, End, ID for all vars
    all_regions_df = pd.DataFrame()

    # go through each svtype
    for svtype, svDF in svtype_to_svDF.items():

        # simple types, only one box
        if svtype in {"deletions", "tandemDuplications", "inversions"}: all_regions_df = all_regions_df.append(svDF[["Chr", "Start", "End", "ID"]])

        # two boxes
        elif svtype in {"insertions", "translocations"}:

            # get two dfs
            dfA = svDF.rename(columns={"ChrA":"Chr", "StartA":"Start", "EndA":"End"})
            dfB = svDF.rename(columns={"ChrB":"Chr", "StartB":"Start", "EndB":"End"})

            # keep them
            for df in [dfA, dfB]: all_regions_df = all_regions_df.append(df[["Chr", "Start", "End", "ID"]])

        elif svtype=="remaining":
            
            # get the CHROMdf, which has POS as the same
            dfCHROM = svDF.rename(columns={"#CHROM":"Chr", "POS":"Start"}); dfCHROM["End"] = dfCHROM["Start"]

            # get the other df
            dfCHR2 = svDF.rename(columns={"CHR2":"Chr", "START":"Start", "END":"End"}); dfCHR2["Start"] = dfCHR2.apply(arrange_minus1_START, axis=1)

            # keep them
            for df in [dfCHROM, dfCHR2]: all_regions_df = all_regions_df.append(df[["Chr", "Start", "End", "ID"]])
            
        else: raise ValueError("%s is not a valid svtype"%svtype)

    # now identify clusters of overlapping IDs, and save the IDs

    # map each chr to it's regions df
    chrom_to_regionsDF = {chrom : all_regions_df[all_regions_df.Chr==chrom] for chrom in set(all_regions_df["Chr"])}

    # initialize vars
    list_clusters = [] # a list of sets of IDs, each of the corresponding to a cluster of SVs

    # go through each region
    for qChr, qStart, qEnd, qID in all_regions_df[["Chr", "Start", "End", "ID"]].values:

        # first check if the ID is already in any of the previosuly run clusters
        if any([qID in cluster for cluster in list_clusters]): continue

        # get the df that is of this region
        df = chrom_to_regionsDF[qChr]

        # find overlapping SVs
        overlappingIDs = set(df[~((df.End<qStart) | (df.Start>qEnd))]["ID"])

        # if there are overlaps add them to the cluster_list with itself
        if overlappingIDs!={qID} and len(overlappingIDs)>0: list_clusters.append(overlappingIDs.union({qID}))

    all_SVs = set(all_regions_df.ID)
    if len(list_clusters)>0: all_SVs_in_cluster = set.union(*list_clusters)
    else: all_SVs_in_cluster = set()

    print("There are %i clusters of SVs involving %i of %i SVs"%(len(list_clusters), len(all_SVs_in_cluster), len(all_SVs)))

    return list_clusters

   # run samtools mpileup, keeping only the first 4 columns (which are chromosome, position, base and coverage)
    mpileup_output = "%s/aligned_reads_bed_first_columns.mpileup"%cnv_outdir; mpileup_output_tmp = "%s.tmp"%mpileup_output
    if fun.file_is_empty(mpileup_output) or opt.replace is True or opt.replace_cnv is True:
    #if True:

        print("Running samtools mpileup")

        def get_files_mpileup_one_chr(cnv_outdir): return sorted(["%s/%s"%(cnv_outdir, file) for file in os.listdir(cnv_outdir) if file.endswith("mpileup_one_chr")])

        # remove previously generated files
        for x in get_files_mpileup_one_chr(cnv_outdir): os.unlink(x)

        # parallelized: generate mpileups for all chromosomes
        cmd_mpileup = '%s view -H %s | grep "\@SQ" | sed "s/^.*SN://g" | cut -f 1 | xargs -I {} -n 1 -P %i sh -c "%s mpileup -BQ0 -a -d 100000000 -f %s -r \"{}\" %s | cut -f1,2,3,4 > %s/\"{}\".mpileup_one_chr"'%(samtools, sorted_bam, opt.threads, samtools, opt.ref, sorted_bam, cnv_outdir); fun.run_cmd(cmd_mpileup)

        # merge all in one
        fun.run_cmd("cat %s > %s"%(" ".join(get_files_mpileup_one_chr(cnv_outdir)), mpileup_output_tmp))

        # unparallelized --> this is only computing for all regions, not the entire chromosome
        #cmd_mpileup = "%s mpileup -a -l %s -f %s %s | cut -f1,2,3,4 > %s"%(samtools, bed_file_regions, opt.ref, sorted_bam, mpileup_output_tmp); fun.run_cmd(cmd_mpileup)

        os.rename(mpileup_output_tmp, mpileup_output)

    # mpileup does not always finish correctly, something you should test
    if fun.file_is_empty(mpileup_output): raise ValueError("mpileup was incorrectly performed")

    # create a python dataframe that contains, for each chromosome, the position and the coverage
    dataframe_mpileup = "%s_df.py"%mpileup_output
    if fun.file_is_empty(dataframe_mpileup) or opt.replace is True or opt.replace_cnv is True:  
    #if True:
        print("Creating mpileup dataframe")
        fun.save_object(pd.read_csv(mpileup_output, sep="\t", header=-1, names=["chr", "position", "base", "coverage"]).set_index("chr"), dataframe_mpileup)

    # create a table in which each gene has the median per-position coverage
    gene_to_coverage_file = "%s/gene_to_coverage_genes.tab"%cnv_outdir
    if fun.file_is_empty(gene_to_coverage_file) or opt.replace is True or opt.replace_cnv is True: # or opt.replace_cnv is True:
    #if True:
        print("Getting gene to coverage file")
        fun.write_coverage_per_gene(mpileup_output, bed_file, gene_to_coverage_file)

    # same as before, now for regions that are +-10,000 away of the gene. It will be useful for further normalisation
    gene_to_coverage_file_regions = "%s/gene_to_coverage_regions.tab"%cnv_outdir
    if fun.file_is_empty(gene_to_coverage_file_regions) or opt.replace is True or opt.replace_cnv is True:
    #if True:
        print("Getting gene to coverage file for windows")
        fun.write_coverage_per_gene(mpileup_output, bed_file_regions, gene_to_coverage_file_regions)

# GATK LeftAlignAndTrimVariants normalisation (only do it if there are no ambiguous characters in your genome)

# define ambiguous nucleotides present in the reference
ambiguous_nts = set.union(*[set(str(chrom.seq).upper()) for chrom in SeqIO.parse(opt.ref, "fasta")]).difference({"A", "C", "G", "T", "N"})

if len(ambiguous_nts)==0:

    # normalize with the GATK
    for unnormalised_vcf in filtered_vcf_results:

        # define the normalised output
        folder = "/".join(unnormalised_vcf.split("/")[0:-1])
        normalised_vcf = "%s/output.filt.norm_gatk.vcf"%folder; normalised_vcf_tmp = "%s.tmp"%normalised_vcf

        # generate an unifyed representation of the vcfs
        if fun.file_is_empty(normalised_vcf) or opt.replace is True:
            print("Running GATK LeftAlignAndTrimVariants for vcf %s"%unnormalised_vcf)
            cmd_normalise = "%s LeftAlignAndTrimVariants -R %s --variant %s -O %s"%(gatk, opt.ref, unnormalised_vcf, normalised_vcf_tmp); fun.run_cmd(cmd_normalise)
            os.rename(normalised_vcf_tmp, normalised_vcf)

        # keep
        all_normalised_vcfs.add(normalised_vcf)

    print("GATK Normalisation is done")

else: print("Warning, GATK LeftAlignAndTrimVariants cannot handle the ambiguous nucleotides found in the provided genome: ", ambiguous_nts)



   ############## KNOWN VARS OPTIMISATION ##########
    if known_genomes_withSV_and_shortReads_table is not None:
    #if False is True: # debug

        # go through each genome, which will be regarded as genomeID in genomeID_to_knownSVdict
        for genomeCounter, (assembly, reads1, reads2, ID) in enumerate(pd.read_csv(known_genomes_withSV_and_shortReads_table, sep="\t")[["assembly", "R1", "R2", "ID"]].values):

            print(ID, assembly, reads1, reads2)
            
            # get a dict that maps structural variation type to a table that has it, in the same format as in RsvSim            
            know_SV_dict = generate_tables_of_SV_between_genomes_gridssClove(assembly, reference_genome, replace=replace, threads=threads, expected_ploidy=expected_ploidy)

                        
            # this is to generate the parameter optimisation from known data

            # define an outdir for the benchmark. Everything will be written under the reads of 
            benchmark_outdir = "%s_withReads_%s.benchmarking_SV"%(assembly, get_file(reads1)); make_folder(benchmark_outdir)
            
            #### ALIGN READS ######
            realData_bamfile = "%s/aligned_reads.bam"%(benchmark_outdir)
            realData_sorted_bam = "%s.sorted"%(realData_bamfile)
            realData_index_bam = "%s.bai"%(realData_sorted_bam)

            run_bwa_mem(reads1, reads2, reference_genome, benchmark_outdir, realData_bamfile, realData_sorted_bam, realData_index_bam, name_sample="real_data_withSV", threads=threads, replace=replace)
            #######################

            # define parms of of insert size
            realData_median_insert_size, realData_median_insert_size_sd  = get_insert_size_distribution(realData_sorted_bam, replace=replace, threads=threads)

            # get a df with a benchmark of many different parameters. This will also report some plots with the 
            benchmarking_df = benchmark_GridssClove_for_knownSV(realData_sorted_bam, reference_genome, know_SV_dict, benchmark_outdir, range_filtering=range_filtering_benchmark, expected_AF=1.0, replace=replace, threads=threads, median_insert_size=realData_median_insert_size, median_insert_size_sd=realData_median_insert_size_sd, window_l=window_l, mitochondrial_chromosome=mitochondrial_chromosome)

            # add some parms and keep
            ploidy = "consensus_ref"
            benchmarking_df["genomeID"] = [genomeID]*len(benchmarking_df)
            benchmarking_df["ploidy"] = [ploidy]*len(benchmarking_df)
            df_benchmark_all = df_benchmark_all.append(benchmarking_df, sort=True)

            genomeID = "realData_%s"%ID


            # keep the known SVs
            genomeID_to_knownSVdict[genomeID] = know_SV_dict

            #if genomeCounter==1: break # debug (only consider two genomes)
        
    #############################################



    all_df_benchmark_longReads = pd.DataFrame()



            ### CHECK HOW THE SV from ALIGNMENT WORKS ON THESE REARRANGED GENOMES ####
            if check_SVfromALNpipeline_simulated_genomes is True:
                print("checking how the alignment-based checking of SV works")

                # get the variants from simulating reads
                predicted_tablesSV = generate_tables_of_SV_between_genomes_gridssClove(rearranged_genome, reference_genome, replace=replace, threads=threads, expected_ploidy=expected_ploidy)

                # get a df of benchmarking
                fileprefix = "%s/rearranged_genome_benchmarking_SV"%simType_dir
                df_benchmark_longReads = benchmark_processedSVs_against_knownSVs_inHouse(predicted_tablesSV, know_SV_dict, fileprefix, replace=replace, analysis_benchmarking=True, tolerance_bp=50)

                # plot the results of the benchmark
                filename = "%s/rearranged_genome_benchmarking_SV_bars.pdf"%simType_dir
                plot_bars_single_df_benchmark(df_benchmark_longReads, filename)

                # keep
                df_benchmark_longReads["genomeID"] = [genomeID]*len(df_benchmark_longReads)
                df_benchmark_longReads["simType"] = [simType]*len(df_benchmark_longReads)
                all_df_benchmark_longReads = all_df_benchmark_longReads.append(df_benchmark_longReads)

            
            ###########################################################################



    # plot the results of all the benchmarking of long reads
    plots_dir = "%s/plots"%outdir; make_folder(plots_dir)
    for simType in sorted(set(all_df_benchmark_longReads.simType)):

        # get df
        simType_df = all_df_benchmark_longReads[all_df_benchmark_longReads.simType==simType]

        # get a file
        filename = "%s/rearranged_genome_benchmarking_SV_bars_for_longReads_%s.pdf"%(plots_dir, simType)

        # get plot
        plot_bars_single_df_benchmark(simType_df, filename)



def get_and_report_filtering_accuracy_across_genomes_and_ploidies(df_benchmark, genomeID_to_knownSVdict, outdir, known_genomes_withSV_and_shortReads_table, PlotsDir, reference_genome, replace=False, consider_integrated_filtering=True, threads=4, mitochondrial_chromosome="mito_C_glabrata_CBS138"):

    """This function takes a df that has the benchmarking info (each line is a set of filtering parameters) and a "genomeID" and "ploidy" fields, which indicate the unique genomes. The idea is to pick, for each genomeID and ploidy combination, the best filtering set (highest Fscore) for each svtype and test for each svtype of genomeID_to_knownSVdict. If the genomeID is not in df_benchmark.genomeID it will run the whole gridss pipeline for the desired filtering sets, considering known_genomes_withSV_and_shortReads_table """

    # general vars
    
    ##### add an overal accuracy measurement for df_benchmark  ########
    if consider_integrated_filtering is True:

        # define the expected df field
        df_benchmark_with_integratedInfo_file = "%s/df_benchmark_all_with_integratedInfo.py"%outdir
        #df_benchmark_with_integratedInfo_file = "%s/df_benchmark_all_with_integratedInfo_filters0.py"%outdir # DEBUG


        if file_is_empty(df_benchmark_with_integratedInfo_file) or replace is True:
        #if True: # debug

            # define the fileds related to a filtering ID
            fields_filtering_ID = ['bedpe', 'benchmarkID','genomeID', 'gridss_VCFoutput', 'gridss_maxcoverage', 'gridss_regionsToIgnoreBed', 'ploidy']

            # get as filtered df
            print("getting all the variants integrated for each set of filters")
            df_benchmark_allSVtypes = df_benchmark.groupby(fields_filtering_ID, as_index=True).apply(get_integrated_benchmarking_fields_series_for_setFilters_df)

            # add the indices as fields, to match those of df_benchmark
            for I, field in enumerate(fields_filtering_ID): df_benchmark_allSVtypes[field] =  df_benchmark_allSVtypes.index.get_level_values(I)
            df_benchmark_allSVtypes = df_benchmark_allSVtypes.set_index("svtype", drop=False)

            # append to df benchmark
            df_benchmark = df_benchmark.append(df_benchmark_allSVtypes[list(df_benchmark.keys())])

            # keep
            save_object(df_benchmark, df_benchmark_with_integratedInfo_file)

        else: df_benchmark = load_object(df_benchmark_with_integratedInfo_file)

    #####################################################################

    # get, for each combination of genomeID, ploidy, svtype the best and less conservative filterset into a list
    print("Getting list of best filters for each genome, ploidy and svtype")
    IDs_sepparate_measurements = ["genomeID", "ploidy", "svtype"]
    df_best_filters = df_benchmark.groupby(IDs_sepparate_measurements).apply(get_best_less_conservative_row_df_benchmark)
    for I, field in enumerate(IDs_sepparate_measurements): df_best_filters[field] =  df_best_filters.index.get_level_values(I)

    # add a field of rounded max coverage (some samples are less than and the others more than 50000)
    median_non50000 = int(np.median(df_best_filters[df_best_filters.gridss_maxcoverage!=50000].gridss_maxcoverage))
    def get_maxcoverage(x):

        if x==50000: return int(x) # the maximum
        elif abs(median_non50000-x)<1000: return median_non50000 # the median of the low-coverage regions
        else: 
            print("warning: There are some values very far away of the rounded median"); return int(x) # something else

    df_best_filters["rounded_gridss_maxcoverage"] = df_best_filters.gridss_maxcoverage.apply(get_maxcoverage)
    df_benchmark["rounded_gridss_maxcoverage"] = df_benchmark.gridss_maxcoverage.apply(get_maxcoverage)

    # define the combinations of regions_to_ignore and max_coverage
    df_regionsIgnore_maxCov = df_best_filters[["gridss_regionsToIgnoreBed", "rounded_gridss_maxcoverage"]].drop_duplicates()

    ####### GENERATE A DF WITH THE INFO OF EACH SV SET TO BE TESTED THROUGH df_best_filters ############

    # initialize dicts to run cross-benchmarking
    genomeIDandPlody_to_info = {}

    # initialize a folder were files of the benchmarking will be stored
    cross_benchmarking_files_dir = "%s/cross_benchmarking_files"%outdir; make_folder(cross_benchmarking_files_dir)


    ## add the info for real data
    for genomeCounter, (assembly, reads1, reads2, ID) in enumerate(pd.read_csv(known_genomes_withSV_and_shortReads_table, sep="\t")[["assembly", "R1", "R2", "ID"]].values):
        print("running bwa mem and gridss for %s"%ID)

        # ALIGN READS #
        gridss_outdir = "%s_gridss_runs_against_%s"%(reads1, get_file(reference_genome)); make_folder(gridss_outdir)
        print(gridss_outdir)
        realData_bamfile = "%s/aligned_reads.bam"%(gridss_outdir)
        realData_sorted_bam = "%s.sorted"%(realData_bamfile)
        realData_index_bam = "%s.bai"%(realData_sorted_bam)

        run_bwa_mem(reads1, reads2, reference_genome, gridss_outdir, realData_bamfile, realData_sorted_bam, realData_index_bam, name_sample="reads_aligned_for_GRIDSS", threads=threads, replace=replace)
        
        # initialize test dict
        test_gridss_info_dict = {}

        # run gridss for all possible filters
        for gridss_regionsToIgnoreBed, gridss_maxcoverage in df_regionsIgnore_maxCov.values:

            # define an outdir
            gridss_outdir_filterSet = "%s/gridss_maxcoverage%i_regionsToIgnore%s"%(gridss_outdir, gridss_maxcoverage, get_file(gridss_regionsToIgnoreBed)); make_folder(gridss_outdir_filterSet)

            # generate gridss vcf
            gridss_VCFoutput = run_gridss_and_annotateSimpleType(realData_sorted_bam, reference_genome, gridss_outdir_filterSet, replace=replace, threads=threads, blacklisted_regions=gridss_regionsToIgnoreBed, maxcoverage=gridss_maxcoverage)

            # keep into test_gridss_info_dict
            test_gridss_info_dict.setdefault(gridss_regionsToIgnoreBed, {}).setdefault(gridss_maxcoverage, gridss_VCFoutput)

        # keep in the dict
        genomeID = "realData_%s"%ID
        ploidy = "consensus_ref"
        knownSVdict = genomeID_to_knownSVdict[genomeID]
        processing_dir = "%s/%s_%s"%(cross_benchmarking_files_dir, genomeID, ploidy); make_folder(processing_dir)

        # calculate insert size metrics
        median_insert_size, median_insert_size_sd  = get_insert_size_distribution(realData_sorted_bam, replace=replace, threads=multiproc.cpu_count())

        # calculate median coverage per 10000 bp coverage
        outdir_coverage = "%s/%s_calculating_coverage_perwindow"%(outdir, genomeID); make_folder(outdir_coverage)
        coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir_coverage, realData_sorted_bam, windows_file="none", replace=replace, window_l=5000), sep="\t")
        median_coverage = np.median(coverage_df[coverage_df["#chrom"]!=mitochondrial_chromosome].mediancov_1); print("The median coverage is %i"%median_coverage)


        genomeIDandPlody_to_info[(genomeID, ploidy)] = {"test_SVdict":knownSVdict, "outdir":processing_dir, "df_filters_train":df_best_filters, "test_gridss_info_dict":test_gridss_info_dict, "sorted_bam":realData_sorted_bam, "median_coverage":median_coverage, "median_insert_size":median_insert_size, "median_insert_size_sd":median_insert_size_sd}

    ## add the info for all genomes and ploidies found in df_best_filters
    for genomeID, ploidy in df_best_filters[["genomeID", "ploidy"]].drop_duplicates().values:

        # initialize test dict
        test_gridss_info_dict = {}

        # go through the gridss filterings
        for gridss_regionsToIgnoreBed, gridss_maxcoverage in df_regionsIgnore_maxCov.values:

            # find in df_benchmark the corresponding value
            df_ben_int = df_benchmark[(df_benchmark.genomeID==genomeID) & (df_benchmark.ploidy==ploidy) & (df_benchmark.rounded_gridss_maxcoverage==gridss_maxcoverage) & (df_benchmark.gridss_regionsToIgnoreBed==gridss_regionsToIgnoreBed)][["gridss_VCFoutput", "sorted_bam", "median_coverage", "median_insert_size", "median_insert_size_sd"]].drop_duplicates()

            # debug
            if len(df_ben_int)!=1: raise ValueError("There are not only one gridss vcfs with the given genomeID and ploidy")
            gridss_VCFoutput = df_ben_int.gridss_VCFoutput.iloc[0]

            # get the known vars
            knownSVdict = genomeID_to_knownSVdict[genomeID]

            # keep into test_gridss_info_dict
            test_gridss_info_dict.setdefault(gridss_regionsToIgnoreBed, {}).setdefault(gridss_maxcoverage, gridss_VCFoutput)

        # keep in the dict
        knownSVdict = genomeID_to_knownSVdict[genomeID]
        processing_dir = "%s/%s_%s"%(cross_benchmarking_files_dir, genomeID, ploidy); make_folder(processing_dir)
        genomeIDandPlody_to_info[(genomeID, ploidy)] = {"test_SVdict":knownSVdict, "outdir":processing_dir, "df_filters_train":df_best_filters, "test_gridss_info_dict":test_gridss_info_dict, "sorted_bam":df_ben_int.sorted_bam.iloc[0], "median_coverage":df_ben_int.median_coverage.iloc[0], "median_insert_size":df_ben_int.median_insert_size.iloc[0], "median_insert_size_sd":df_ben_int.median_insert_size_sd.iloc[0]}

    ###########################################################################################

    # run a function in parallel that will take a genome and ploidy combination and evaluate the accuracy of all the filters in df_best_filters. first prepare input as list of tuples
    list_inputs = [(d["test_SVdict"], d["outdir"], d["df_filters_train"], d["test_gridss_info_dict"], genomeID, ploidy, d["sorted_bam"], reference_genome, d["median_coverage"], d["median_insert_size"], d["median_insert_size_sd"],  replace) for (genomeID, ploidy), d in genomeIDandPlody_to_info.items()] # 0 is for debug


    # test withot parallelization
    #list_cross_benchmarking_dfs = list(map(lambda x: get_benchmarking_df_for_testSVs_from_trainSV_filterSets(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11]), list_inputs))

    # run in parallel
    with multiproc.Pool(multiproc.cpu_count()) as pool:
        list_cross_benchmarking_dfs = pool.starmap(get_benchmarking_df_for_testSVs_from_trainSV_filterSets, list_inputs) # needs if __name__=="__main__" 
        #list_cross_benchmarking_dfs = pool.map(lambda x: get_benchmarking_df_for_testSVs_from_trainSV_filterSets(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11]), list_inputs)
        
        pool.close()
        pool.terminate()
    

    # concatenate all the dfs
    df_cross_benchmark = pd.concat(list_cross_benchmarking_dfs)

    # add the simulation type for train and test
    def add_simulation_name_and_type(genomeID, tag):

        if "biased_towards_repeats" in genomeID: 
            simName = genomeID.split("_simType_biased_towards_repeats")[0]
            simType = "biased_towards_repeats"

        elif "uniform" in genomeID: 
            simName = genomeID.split("_simType_uniform")[0]
            simType = "uniform"

        elif "realData" in genomeID: 
            simName = genomeID
            simType = "realData"

        return pd.Series({"%s_simName"%tag : simName, "%s_simType"%tag : simType})

    df_cross_benchmark[["train_simName", "train_simType"]] = df_cross_benchmark.train_genomeID.apply(lambda x: add_simulation_name_and_type(x, "train"))
    df_cross_benchmark[["test_simName", "test_simType"]] = df_cross_benchmark.test_genomeID.apply(lambda x: add_simulation_name_and_type(x, "test"))


    return df_cross_benchmark






    JLDKLKLSFDKLFD


# a function to copy some parts

def transform_cut_and_paste_to_copy_and_paste_insertions(reference_genome, rearranged_genome, insertions_file, svtype_to_svDF):

    """ This function takes a rearranged genome and reinserts the copy-and-paste insertions where they should be """

    print("reinserting-copy-and-paste insertions into %s"%insertions_file)

    # load df and keep the copy-and-paste insertions
    df = pd.read_csv(insertions_file, sep="\t")
    df = df[df.Copied]

    if len(df)>0:

        # define an unmodified genome
        rearranged_genome_unmodified = "%s.unmodified.fasta"%rearranged_genome
        rearranged_genome_unmodified_tmp = "%s.tmp"%rearranged_genome_unmodified

        if file_is_empty(rearranged_genome_unmodified):

            # if the unmodified tmps is writen, replace the rearranged_genome with it
            if not file_is_empty(rearranged_genome_unmodified_tmp): os.rename(rearranged_genome_unmodified_tmp, rearranged_genome)

            # get the rearranged genome seq
            chr_to_rearrangedSeq = {seq.id: str(seq.seq) for seq in SeqIO.parse(rearranged_genome, "fasta")}
            all_rearranged_chromosomes_together = "".join(chr_to_rearrangedSeq.values())

            # get the seq
            chr_to_refSeq = {seq.id: str(seq.seq) for seq in SeqIO.parse(reference_genome, "fasta")}

            # define the length of each chrom
            chr_to_lenSeq = {chrom : len(seq) for chrom, seq in chr_to_refSeq.items()}

            # define all the positions with breakpoints
            df_positions = pd.concat([get_breakpoint_positions_df_in_svDF(svDF) for svtype, svDF in svtype_to_svDF.items()])
            chr_to_bpPositions = dict(df_positions.groupby("Chr").apply(lambda df_c: set(df_c["Pos"])))

            # add the ends of the chromosome, and convert to np array
            for chrom, lenSeq in chr_to_lenSeq.items(): 

                chr_to_bpPositions[chrom].update({1, lenSeq})
                chr_to_bpPositions[chrom] = np.array(sorted(chr_to_bpPositions[chrom]))

            # add the closest breakpoint position of ChrA in the reference
            df["closest_5'breakpoint_position"] = df.apply(lambda r: find_nearest(chr_to_bpPositions[r["ChrA"]][chr_to_bpPositions[r["ChrA"]]<(r["StartA"])], r["StartA"]), axis=1)

            df["closest_3'breakpoint_position"] = df.apply(lambda r: find_nearest(chr_to_bpPositions[r["ChrA"]][chr_to_bpPositions[r["ChrA"]]>(r["EndA"])], r["EndA"]), axis=1)

            # get the 5' sequence (from one position after the closest breakpoint to the position before the breakpoint)
            df["5'sequence"] = df.apply(lambda r: chr_to_refSeq[r["ChrA"]][r["closest_5'breakpoint_position"]:r["StartA"]-1], axis=1)

            # get the 3' sequence (from the position after End to the position before the closest breakpoint)
            df["3'sequence"] = df.apply(lambda r: chr_to_refSeq[r["ChrA"]][r["EndA"]:r["closest_3'breakpoint_position"]-1], axis=1)

            # get the deleted sequence (from the start to the end)
            df["deleted_sequence"] = df.apply(lambda r: chr_to_refSeq[r["ChrA"]][r["StartA"]-1:r["EndA"]], axis=1)

            # change the chromosome seq in the sequence 
            for I, (chrA, seq5, seq3, del_seq) in enumerate(df[["ChrA", "5'sequence", "3'sequence", "deleted_sequence"]].values):
                print("copy-paste-insertion %i.."%I)

                # all seq
                ref_seq = seq5+del_seq+seq3

                # conformation in the rearranged chromosome
                rearranged_seq = seq5+seq3

                # check that the rearranged seq appears once in the genome and the ref seq in the ref genome. And they do not cross.
                chrA_refSeq = chr_to_refSeq[chrA]
                if not(chrA_refSeq.count(ref_seq)==1 and chrA_refSeq.count(rearranged_seq)==0 and all_rearranged_chromosomes_together.count(rearranged_seq)==1 and all_rearranged_chromosomes_together.count(ref_seq)==0): raise ValueError("The sequence is not unique")

                # go through each chrom of the rearranged seqs
                for chrom in chr_to_rearrangedSeq.keys():

                    # get the rearranged sequence
                    seq = cp.deepcopy(chr_to_rearrangedSeq[chrom])

                    # if the rearrangement sequence is in this chromosome, change it
                    if rearranged_seq in seq: 

                        # update the chr_to_rearrangedSeq so that it contains the reference sequence (copied)
                        chr_to_rearrangedSeq[chrom] = seq.replace(rearranged_seq, ref_seq)
                        break

            # get the rearranged genome into the file
            seq_records_list = [SeqRecord(Seq(seq), id=chrom, name=chrom, description=chrom) for chrom, seq in chr_to_rearrangedSeq.items()]

            # write the unmodified one
            print("writing")
            run_cmd("cp %s %s.tmp"%(rearranged_genome, rearranged_genome_unmodified_tmp))
            os.rename("%s.tmp"%rearranged_genome_unmodified_tmp, rearranged_genome_unmodified_tmp)

            # write the modified genome
            SeqIO.write(seq_records_list, rearranged_genome, "fasta")

            # write the modified genome
            os.rename(rearranged_genome_unmodified_tmp, rearranged_genome_unmodified)

        else: print("the insertions have already been modified")


