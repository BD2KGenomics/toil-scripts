import os
import subprocess

from toil_scripts.lib.programs import docker_call


def gatk_select_variants(job, mode, vcf_id, config):
    """
    Isolate variant types using GATK SelectVariants

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param mode str: variant type (i.e. SNP or INDEL)
    :param vcf_id str: VCF FileStoreID
    :param config Namespace: Configuration and shared FileStoreIDs
    :return: VCF FileStoreID
    :rtype: str
    """
    job.fileStore.logToMaster('Running GATK SelectVariants: %s' % mode)
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'genome.fa': config.genome_fasta,
              'genome.fa.fai': config.genome_fai,
              'genome.dict': config.genome_dict,
              'input.vcf': vcf_id}
    for name, file_store_id in inputs.iteritems():
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))

    command = ['-T', 'SelectVariants',
               '-R', 'genome.fa',
               '-V', 'input.vcf',
               '-o', 'output.vcf',
               '-selectType', mode]

    outputs = {'output.vcf': None}
    docker_call(work_dir=work_dir,
                env={'JAVA_OPTS': '-Djava.io.tmpdir=/data/ -Xmx{}'.format(job.memory)},
                parameters=command,
                tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'output.vcf'))


def split_vcf_by_name(job, sample_names, vcf_id, config):
    """
    Splits a VCF using a list of sample names.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param list sample_names: Sample names
    :param str vcf_id: VCF FileStoreID
    :param Namespace config: Configuration and shared FileStoreIDs
    :return: VCF FileStoreIDs
    :rtype: dictionary
    """
    job.fileStore.logToMaster('Splitting VCF with GATK SelectVariants')
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'genome.fa': config.genome_fasta,
              'genome.fa.fai': config.genome_fai,
              'genome.dict': config.genome_dict,
              'input.vcf': vcf_id}
    for name, file_store_id in inputs.iteritems():
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))

    vcfs = {}
    for sample in sample_names:
        output = '{}.vcf'.format(sample)
        command = ['-T', 'SelectVariants',
                   '-R', 'genome.fa',
                   '-V', 'input.vcf',
                   '-o', output,
                   '--sample_name', sample]
        outputs = {output: None}
        docker_call(work_dir=work_dir,
                    env={'JAVA_OPTS': '-Djava.io.tmpdir=/data/ -Xmx{}'.format(job.memory)},
                    parameters=command,
                    tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                    inputs=inputs.keys(),
                    outputs=outputs)
        vcfs[sample] = job.fileStore.writeGlobalFile(os.path.join(work_dir, output))
    return vcfs


def gatk_variant_filtration(job, mode, vcf_id, config):
    """
    Filters VCF using GATK recommended filters. VCF must contain a single variant
    type: SNPs or INDELs.

    SNP Filter:
    QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0

    INDEL Filter:
    QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0

    :param job: Toil Job instance
    :param mode str: variant type
    :param vcf_id str: VCF FileStoreID
    :param config Namespace: Configuration and shared FileStoreIDs
    :return: Filtered VCF FileStoreID
    :rtype: str
    """
    mode = mode.upper()
    job.fileStore.logToMaster('Apply %s Filter' % mode)
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'genome.fa': config.genome_fasta,
              'genome.fa.fai': config.genome_fai,
              'genome.dict': config.genome_dict,
              'input.vcf': vcf_id}
    for name, file_store_id in inputs.iteritems():
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))

    # Recommended GATK hard filters:
    # https://software.broadinstitute.org/gatk/documentation/article?id=2806
    # Filters are different for SNPs and INDELs. Please refer to GATK documentation.
    if mode == 'SNP':
        expression = '"QD < 2.0 || FS > 60.0 || MQ < 40.0 || ' \
                     'MQRankSum < -12.5 || ReadPosRankSum < -8.0"'

    elif mode == 'INDEL':
        expression = '"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"'

    else:
        raise ValueError('Variant filter modes can be SNP or INDEL, got %s' % mode)

    command = ['-T', 'VariantFiltration',
               '-R', 'genome.fa',
               '-V', 'input.vcf',
               '--filterExpression', expression,
               '--filterName', 'GATK_Germline_Hard_Filter_%s' % mode,   # Documents filter in header
               '-o', 'filtered_variants.vcf']

    outputs = {'filtered_variants.vcf': None}
    docker_call(work_dir=work_dir,
                env={'JAVA_OPTS': '-Djava.io.tmpdir=/data/ -Xmx{}'.format(job.memory)},
                parameters=command,
                tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs)

    outpath = os.path.join(work_dir, 'filtered_variants.vcf')

    # Fix extra quotation marks in FILTER line
    sed_cmd = 's/"{}"/{}/'.format(expression, expression)
    subprocess.call(['sed', '-i', sed_cmd, outpath])
    return job.fileStore.writeGlobalFile(outpath)


def gatk_variant_recalibrator(job, mode, vcf_id, config):
    """
    Variant quality score recalibration for SNP or INDEL variants. VQSR must be run twice because
    SNP and INDEL VQSR use distinct models.

    :param JobFunctionWrappingJob job: Job instance
    :param str mode: Determines variant recalibration mode (SNP or INDEL)
    :param str vcf_id: VCF FileStoreID
    :param Namespace config: Input parameters and shared FileStoreIDs
    :return: recalibration table, tranches, plot FileStoreIDs
    :rtype: tuple
    """
    mode = mode.upper()
    if mode not in {'INDEL', 'SNP'}:
        raise ValueError('Variant recalibration mode must be INDEL or SNP, got %s' % mode)

    job.fileStore.logToMaster('Running GATK VariantRecalibrator ({} Mode)'.format(mode.upper()))
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'genome.fa': config.genome_fasta,
              'genome.fa.fai': config.genome_fai,
              'genome.dict': config.genome_dict,
              'input.vcf': vcf_id}

    # Refer to GATK documentation for description of recommended parameters:
    # https://software.broadinstitute.org/gatk/documentation/article?id=1259
    # https://software.broadinstitute.org/gatk/documentation/article?id=2805

    # This base command includes parameters for both INDEL and SNP VQSR.
    # DP is a recommended annotation, but does not work well with exome data
    command = ['-T', 'VariantRecalibrator',
               '-R', 'genome.fa',
               '-input', 'input.vcf',
               '-nt', str(job.cores),
               '--maxGaussians', '4',
               '-an', 'QualByDepth',
               '-an', 'FisherStrand',
               '-an', 'StrandOddsRatio',
               '-an', 'ReadPosRankSum',
               '-an', 'MQRankSum',
               '-an', 'InbreedingCoeff',
               '-tranche', '100.0',
               '-tranche', '99.9',
               '-tranche', '99.0',
               '-tranche', '90.0',
               '-recalFile', 'output.recal',
               '-tranchesFile', 'output.tranches',
               '-rscriptFile', 'output.plots.R']

    # These parameters and resource files are specific for SNP VQSR.
    if mode == 'SNP':
        command += ['-resource:hapmap,known=false,training=true,truth=true,prior=15.0', 'hapmap.vcf',
                    '-resource:omni,known=false,training=true,truth=true,prior=12.0', 'omni.vcf',
                    '-resource:dbsnp,known=true,training=false,truth=false,prior=2.0', 'dbsnp.vcf',
                    '-resource:1000G,known=false,training=true,truth=false,prior=10.0', '1000G.vcf',
                    '-an', 'RMSMappingQuality',
                    '-mode', 'SNP']
        inputs['hapmap.vcf'] = config.hapmap
        inputs['omni.vcf'] = config.omni
        inputs['dbsnp.vcf'] = config.dbsnp
        inputs['1000G.vcf'] = config.phase

    # These parameters and resource files are specific for INDEL VQSR.
    if mode == 'INDEL':
        command += ['-resource:mills,known=false,training=true,truth=true,prior=12.0', 'mills.vcf',
                    '-resource:dbsnp,known=true,training=false,truth=false,prior=2.0', 'dbsnp.vcf',
                    '-mode', 'INDEL']
        inputs['mills.vcf'] = config.mills
        inputs['dbsnp.vcf'] = config.dbsnp

    if config.unsafe_mode:
        command.extend(['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'])

    # Delay reading in files until function is configured
    for name, file_store_id in inputs.iteritems():
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))

    outputs = {'output.recal': None, 'output.tranches': None, 'output.plots.R': None}
    docker_call(work_dir=work_dir,
                env={'JAVA_OPTS': '-Djava.io.tmpdir=/data/ -Xmx{}'.format(job.memory)},
                parameters=command,
                tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs)
    recal_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'output.recal'))
    tranches_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'output.tranches'))
    plots_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'output.plots.R'))
    return recal_id, tranches_id, plots_id


def gatk_apply_variant_recalibration(job, mode, vcf_id, recal_id, tranches_id, config):
    """
    Apply variant quality score recalibration.

    :param JobFunctionWrappingJob job: Toil Job instance
    :param str vcf_id: VCF FileStoreID
    :param str recal_id: Recalibration table FileStoreID
    :param str tranches_id: Tranches FileStoreID
    :param Namespace config: Pipeline configuration options and shared files
    :return: Recalibrated VCF file store ID
    :rtype: str
    """
    mode = mode.upper()
    job.fileStore.logToMaster(
        'Running GATK ApplyRecalibration ({} Mode): {}'.format(mode, config.uuid))
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'genome.fa': config.genome_fasta,
              'genome.fa.fai': config.genome_fai,
              'genome.dict': config.genome_dict,
              'input.vcf': vcf_id,
              'recal': recal_id,
              'tranches': tranches_id}
    for name, file_store_id in inputs.iteritems():
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))

    # GATK recommended parameters:
    # https://software.broadinstitute.org/gatk/documentation/article?id=2805
    command = ['-T', 'ApplyRecalibration',
               '-nt', str(job.cores),
               '-mode', mode,
               '-R', 'genome.fa',
               '-input', 'input.vcf',
               '-o', 'vqsr.vcf',
               '-ts_filter_level', '99.0',
               '-recalFile', 'recal',
               '-tranchesFile', 'tranches']

    if config.unsafe_mode:
        command.extend(['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'])

    outputs = {'vqsr.vcf': None}
    docker_call(work_dir=work_dir,
                env={'JAVA_OPTS': '-Djava.io.tmpdir=/data/ -Xmx{}'.format(job.memory)},
                parameters=command,
                tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'vqsr.vcf'))


def gatk_combine_variants(job, vcfs, config):
    """
    Combine a dictionary of VCF FileStoreIDs from the same sample

    :param JobFunctionWrappingJob job: Toil Job instance
    :param dict vcfs: Dictionary of VCF FileStoreIDs
    :param Namespace config: Pipeline configuration options and shared files
    :return: Merged VCF FileStoreID
    :rtype: str
    """
    job.fileStore.logToMaster('Running GATK CombineVariants')
    work_dir = job.fileStore.getLocalTempDir()

    inputs = {'genome.fa': config.genome_fasta,
              'genome.fa.fai': config.genome_fai,
              'genome.dict': config.genome_dict}
    inputs.update(vcfs)
    for name, file_store_id in inputs.iteritems():
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))

    command = ['-T', 'CombineVariants',
               '-R', '/data/genome.fa',
               '-o', '/data/merged.vcf',
               '--genotypemergeoption', 'UNSORTED']  # For merging VCFs from same sample

    for uuid, vcf_id in vcfs.iteritems():
        command.extend(['--variant', os.path.join('/data', uuid)])

    outputs = {'merged.vcf': None}
    docker_call(work_dir=work_dir,
                env={'JAVA_OPTS': '-Djava.io.tmpdir=/data/ -Xmx{}'.format(job.memory)},
                parameters=command,
                tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'merged.vcf'))
