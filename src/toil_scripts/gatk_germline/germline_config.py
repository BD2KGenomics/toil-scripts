import textwrap


def generate_config():
    return textwrap.dedent("""
        ##############################################################################################################
        # GATK Germline Pipeline configuration file
        # Variant databases can be obtained through the GATK resource bundle:
        # https://software.broadinstitute.org/gatk/guide/article?id=1213
        # http://gatkforums.broadinstitute.org/gatk/discussion/1259/what-vqsr-training-sets-arguments-should-i-use-for-my-specific-project
        ##############################################################################################################
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline
        # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank.
        ##############################################################################################################
        # Required: Number of cores per job
        cores:

        # Required: Java heap size (human readable bytes format i.e. 10G)
        xmx:

        # Required: Approximate input file size (human readable bytes format)
        file-size:

        # Required: S3 URL or local path to output directory
        output-dir:

        # Required: Input BAM file is sorted (Default: False)
        sorted:

        # Required: URL or local path to reference genome FASTA file
        genome-fasta:

        # Optional: URL or local path to reference genome index (Default: None)
        genome-fai:

        # Optional: URL or local path to reference genome sequence dictionary (Default: None)
        genome-dict:

        # Required for VQSR: URL or local path to 1000G SNP resource file (Default: None)
        g1k_snp:

        # Required for preprocessing: URL or local path to 1000G INDEL resource file (Default: None)
        g1k_indel:

        # Required for VQSR: URL or local path HapMap resource file (Default: None)
        hapmap:

        # Required for VQSR: URL or local path Omni resource file (Default: None)
        omni:

        # Required for VQSR: URL or local path to Mills resource file (Default: None)
        mills:

        # Required for VQSR: URL or local path to dbSNP resource file (Default: None)
        dbsnp:

        # Required for FASTQ samples: Align FASTQs or Realign BAM file (Default: False)
        run-bwa:

        # Optional. Trim adapters (Default: False)
        trim:

        # Required for BWA alignment: URL or local path to BWA index file prefix.amb (Default: None)
        amb:

        # Required for BWA alignment: URL or local path to BWA index file prefix.ann (Default: None)
        ann:

        # Required for BWA alignment: URL or local path to BWA index file prefix.bwt (Default: None)
        bwt:

        # Required for BWA alignment: URL or local path to BWA index file prefix.pac (Default: None)
        pac:

        # Required for BWA alignment: URL or local path to BWA index file prefix.sa (Default: None)
        sa:

        # Required for ALT-aware alignment: URL or local path to alternate contigs (Default: None)
        alt:

        # Optional: Run GATK Preprocessing (Default: False)
        preprocess:

        # Optional: Stops after GATK Preprocessing (Default: False)
        preprocess-only:

        # Required: GATK annotations used to filter SNPs
        snp-filter-annotations:

        # Required: GATK annotations used to filter INDELS
        indel-filter-annotations:

        # Required for hard filtering: Name of SNP hard filter for VCF header
        snp_filter_name:

        # Required for hard filtering: SNP JEXL filter expression
        snp_filter_expression:

        # Required for hard filtering: Name of INDEL hard filter for VCF header
        indel_filter_name:

        # Required for hard filtering: INDEL JEXL filter expression
        indel_filter_expression:

        # Optional: Run GATK VQSR (Default: False)
        run-vqsr:

        # Optional: Merges all samples into a single GVCF for genotyping and filtering (Default: False)
        joint-genotype:

        # Optional: Run Oncotator (Default: False)
        run-oncotator:

        # Required for Oncotator: URL or local path to Oncotator database (Default: None)
        oncotator-db:

        # Optional: Suffix added to output filename (i.e. .toil)
        suffix:

        # Optional: Path to key file for SSE-C Encryption (Default: None)
        ssec:

        # Optional: Allow seq dict incompatibility (Default: False)
        unsafe-mode:
        """[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #   There are 2-4 tab-separated columns: UUID, URL, URL2, and Read Group
        #
        #   UUID        Unique sample identifier
        #   URL         URL (http://, ftp://, file://, s3://) pointing to the input FASTQ or BAM file
        #   URL2        (Optional) URL (http://, ftp://, file://, s3://) pointing to paired FASTQ file
        #   RG-line     (Optional) Read group information for FASTQ sample: @RG\tID:foo\tSM:bar
        #
        #   Example below. Lines beginning with # are ignored.
        #
        #   UUID_BAM      file:///path/to/sample.bam
        #                           OR
        #   UUID_FASTQ    file:///path/to/sample.1.fq   file:///path/to/sample.2.fq   @RG\tID:foo\tSM:bar
        #
        #   Place your samples below, one per line.
    """[1:])
