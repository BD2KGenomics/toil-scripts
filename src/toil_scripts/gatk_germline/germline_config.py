import textwrap


def generate_config():
    return textwrap.dedent("""
        # GATK Germline Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline
        # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank.
        ##############################################################################################################
        # Required: Reference Genome URL
        genome-fasta:

        # Required: Reference assembly (i.e. hg19)
        assembly:

        # Optional: Reference Genome Index URL
        genome-fai:

        # Optional: Reference Genome Sequence Dictionary URL
        genome-dict:

        # Optional: 1000G INDELs URL
        phase:

        # Optional: Mills INDELs URL
        mills:

        # Optional: dbSNP URL
        dbsnp:

        # Optional: HapMap URL
        hapmap:

        # Optional: Omni URL
        omni:

        # Optional: Align FASTQs or Realign BAM file (Boolean)
        run-bwa:

        # Optional. BWA trims adapters (Boolean)
        trim: True

        # Optional: BWA Index amb URL
        amb:

        # Optional: BWA Index ann URL
        ann:

        # Optional: BWA Index bwt URL
        bwt:

        # Optional: BWA Index pac URL
        pac:

        # Optional: BWA Index sa URL
        sa:

        # Optional: BWA Index alt URL
        alt:

        # Optional: Run GATK Preprocessing (Boolean)
        preprocess:

        # Optional: Run GATK VQSR (Boolean)
        run-vqsr:

        # Optional: Joint Calling (Boolean)
        joint:

        # Optional: Run Oncotator (Boolean)
        run-oncotator:

        # Optional: Oncotator Database URL
        oncotator-db:

        # Optional: Approximate input file size (Default: 10G)
        file-size:

        # BAM file is already sorted (Boolean)
        sorted:

        # Number of cores per job (Default: 8)
        cores:

        # Java heap size (Default: 10G)
        xmx:

        # Optional: Suffix to be added to output filename (i.e. .toil)
        suffix:

        # Optional: Path to output directory (PATH/URL)
        output-dir:

        # Optional: Path to CSV file containing Synapse username and password
        synapse_credentials:

        # Optional: Path to key file for SSE-C Encryption
        ssec:

        # Optional: Allow seq dict incompatibility (Boolean)
        unsafe-mode:
        """[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #   There are 2 tab-separated columns: UUID and URL
        #
        #   UUID        This should be a unique identifier for the sample to be processed
        #   URL         A URL (http://, ftp://, file://, s3://) pointing to the input SAM or BAM file
        #
        #   Example below. Lines beginning with # are ignored.
        #
        #   UUID_1    file:///path/to/sample.bam
        #                   OR
        #   UUID_1    file:///path/to/sample.1.fq @RG\tID:foo\tSM:bar
        #
        #   Place your samples below, one per line.
    """[1:])