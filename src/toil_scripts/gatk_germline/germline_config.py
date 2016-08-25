import textwrap


def generate_config():
    return textwrap.dedent("""
        ##############################################################################################################
        # GATK Germline Pipeline configuration file
        # Variant databases can be obtained through the GATK resource bundle:
        # https://software.broadinstitute.org/gatk/guide/article?id=1213
        ##############################################################################################################
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline
        # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank.
        ##############################################################################################################
        # Required: Number of cores per job
        cores:

        # Required: Java heap size (human readable format)
        xmx:

        # Required: Approximate input file size (human readable format)
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

        # Optional: URL or local path to 1000 Genomes INDELs resource file (Default: None)
        phase:

        # Optional: URL or local path to Mills INDELs resource file (Default: None)
        mills:

        # Optional: URL or local path to dbSNP resource file (Default: None)
        dbsnp:

        # Optional: URL or local path HapMap resource file (Default: None)
        hapmap:

        # Optional: URL or local path Omni resource file (Default: None)
        omni:

        # Optional: Align FASTQs or Realign BAM file (Default: False)
        run-bwa:

        # Optional. Trim adapters (Default: False)
        trim:

        # Optional: URL or local path to BWA index file prefix.amb (Default: None)
        amb:

        # Optional: URL or local path to BWA index file prefix.ann (Default: None)
        ann:

        # Optional: URL or local path to BWA index file prefix.bwt (Default: None)
        bwt:

        # Optional: URL or local path to BWA index file prefix.pac (Default: None)
        pac:

        # Optional: URL or local path to BWA index file prefix.sa (Default: None)
        sa:

        # Optional: URL or local path to alternate contigs (Default: None)
        # Necessary for ALT-aware alignment
        alt:

        # Optional: Run GATK Preprocessing (Default: False)
        preprocess:

        # Optional: Only runs GATK Preprocessing steps (Default: False)
        preprocess-only:

        # Optional: Run GATK VQSR (Default: False)
        run-vqsr:

        # Optional: Joint Genotyping (Default: False)
        joint:

        # Optional: Run Oncotator (Default: False)
        run-oncotator:

        # Optional: URL or local path to Oncotator database (Default: None)
        oncotator-db:

        # Optional: Suffix added to output filename (i.e. .toil)
        suffix:

        # Optional: Path to CSV file containing Synapse username and password (Default: None)
        synapse_credentials:

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
