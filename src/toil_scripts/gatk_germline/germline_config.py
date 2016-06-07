import textwrap


def generate_config():
    return textwrap.dedent("""
        # GATK Germline Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline
        # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank.
        ##############################################################################################################
        # Required: Reference Genome URL (hg19.fa)
        genome-fasta:

        # Optional: Reference Genome Index URL (hg19.fa.fai)
        genome-fai:

        # Optional: Reference Genome Sequence Dictionary (hg19.dict)
        genome-dict:

        # Optional: URL (1000G_phase1.indels.hg19.sites.fixed.vcf)
        phase:

        # Optional: URL (Mills_and_1000G_gold_standard.indels.hg19.sites.vcf)
        mills:

        # Optional: URL (dbsnp_138.hg19.excluding_sites_after_129.vcf)
        dbsnp:

        # Optional: URL (hapmap_3.3.hg19.sites.vcf)
        hapmap:

        # Optional: URL (1000G_omni2.5.hg19.sites.vcf)
        omni:

        # Optional: (boolean) Run BWA on fastqs
        run-bwa:

        # Optional. If true, BWA trims adapters
        trim: True

        # Optional: Reference fasta file (amb) -- if not present will be generated
        amb:

        # Optional: Reference fasta file (ann) -- If not present will be generated
        ann:

        # Optional: Reference fasta file (bwt) -- If not present will be generated
        bwt:

        # Optional: Reference fasta file (pac) -- If not present will be generated
        pac:

        # Optional: Reference fasta file (sa) -- If not present will be generated
        sa:

        # Optional: Alternate file for reference build (alt). Necessary for alt aware alignment
        alt:

        # Optional: (boolean) Run GATK Preprocessing
        preprocess:

        # Optional: (boolean) Run GATK VQSR
        run-vqsr:

        # Optional: (boolean) Joint Calling
        joint:

        # Optional: Approximate input file size
        file-size: 50G

        # Memory allocation for Java option Xmx
        xmx: 100G

        # Optional: (string) Suffix to be added to final output
        suffix:

        # Optional: (string) Path to output directory
        output-dir:

        # Optional: (string) Path to Key File for SSE-C Encryption
        ssec:

        # Optional: (boolean) Set to True to allow seq dict incompatibility
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