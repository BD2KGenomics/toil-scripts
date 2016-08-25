import textwrap


def generate_config():
    return textwrap.dedent("""
        # GATK Germline Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline
        # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank.
        ##############################################################################################################
        # Required: Number of cores per job
        cores:

        # Required: Java heap size
        xmx:

        # Required: URL or local path to output directory
        output-dir:

        # Required: Approximate input file size
        file-size:

        # Required: Input BAM file is sorted (False)
        sorted:

        # Required: URL or local path to reference genome FASTA file
        genome-fasta:

        # Optional: URL or local path to reference genome index (False)
        # if not present will be generated
        genome-fai:

        # Optional: URL or local path to reference genome sequence dictionary
        # if not present will be generated
        genome-dict:

        # Optional: URL or local path to 1000G INDELs resource file
        phase:

        # Optional: URL or local path to Mills INDELs resource file
        mills:

        # Optional: URL or local path to dbSNP resource file
        dbsnp:

        # Optional: URL or local path HapMap resource file
        hapmap:

        # Optional: URL or local path Omni resource file
        omni:

        # Optional: Align FASTQs or Realign BAM file (False)
        run-bwa:

        # Optional. Trim adapters (False)
        trim: True

        # Optional: URL or local path to BWA index file prefix.amb
        amb:

        # Optional: URL or local path to BWA index file prefix.ann
        ann:

        # Optional: URL or local path to BWA index file prefix.bwt
        bwt:

        # Optional: URL or local path to BWA index file prefix.pac
        pac:

        # Optional: URL or local path to BWA index file prefix.sa
        sa:

        # Optional: URL or local path to ALT contigs. Necessary for ALT-aware alignment
        alt:

        # Optional: Run GATK Preprocessing (False)
        preprocess:

        # Optional: Run GATK VQSR (False)
        run-vqsr:

        # Optional: Joint Calling (False)
        joint:

        # Optional: Run Oncotator (False)
        run-oncotator:

        # Optional: URL or local path to Oncotator database
        oncotator-db:


        # Optional: Suffix added to output filename (i.e. .toil)
        suffix:

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
