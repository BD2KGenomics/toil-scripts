import os

from toil_scripts.lib import require
from toil_scripts.lib.programs import docker_call


def run_cutadapt(job, r1_id, r2_id, fwd_3pr_adapter, rev_3pr_adapter):
    """
    Adapter triming for RNA-seq data

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str r1_id: FileStoreID of fastq read 1
    :param str r2_id: FileStoreID of fastq read 2 (if paired data)
    :param str fwd_3pr_adapter: Adapter sequence for the forward 3' adapter
    :param str rev_3pr_adapter: Adapter sequence for the reverse 3' adapter (second fastq pair)
    :return: R1 and R2 FileStoreIDs
    :rtype: tuple
    """
    work_dir = job.fileStore.getLocalTempDir()
    if r2_id:
        require(rev_3pr_adapter, "Paired end data requires a reverse 3' adapter sequence.")
    # Retrieve files
    parameters = ['-a', fwd_3pr_adapter,
                  '-m', '35']
    if r1_id and r2_id:
        job.fileStore.readGlobalFile(r1_id, os.path.join(work_dir, 'R1.fastq'))
        job.fileStore.readGlobalFile(r2_id, os.path.join(work_dir, 'R2.fastq'))
        parameters.extend(['-A', rev_3pr_adapter,
                           '-o', '/data/R1_cutadapt.fastq',
                           '-p', '/data/R2_cutadapt.fastq',
                           '/data/R1.fastq', '/data/R2.fastq'])
    else:
        job.fileStore.readGlobalFile(r1_id, os.path.join(work_dir, 'R1.fastq'))
        parameters.extend(['-o', '/data/R1.fastq', '/data/R1.fastq'])
    # Call: CutAdapt
    docker_call(tool='quay.io/ucsc_cgl/cutadapt:1.9--6bd44edd2b8f8f17e25c5a268fedaab65fa851d2',
                work_dir=work_dir, parameters=parameters)
    # Write to fileStore
    if r1_id and r2_id:
        r1_cut_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1_cutadapt.fastq'))
        r2_cut_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R2_cutadapt.fastq'))
    else:
        r1_cut_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1_cutadapt.fastq'))
        r2_cut_id = None
    return r1_cut_id, r2_cut_id


def run_samtools_faidx(job, ref_id):
    """
    Use Samtools to create reference index file

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str ref_id: FileStoreID for the reference genome
    :return: FileStoreID for reference index
    :rtype: str
    """
    job.fileStore.logToMaster('Created reference index')
    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.readGlobalFile(ref_id, os.path.join(work_dir, 'ref.fasta'))
    command = ['faidx', 'ref.fasta']
    docker_call(work_dir=work_dir, parameters=command,
                tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e')
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'ref.fasta.fai'))


def run_samtools_index(job, bam_id):
    """
    Runs samtools index to create (.bai) files

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str bam_id: FileStoreID of the bam file
    :return: BAM index FileStoreID
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.readGlobalFile(bam_id, os.path.join(work_dir, 'sample.bam'))
    # Call: index the bam
    parameters = ['index', '/data/sample.bam']
    docker_call(work_dir=work_dir, parameters=parameters,
                tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e')
    # Write to fileStore
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sample.bam.bai'))


def run_picard_create_sequence_dictionary(job, ref_id):
    """
    Use Picard-tools to create reference dictionary

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str ref_id: FileStoreID for the reference genome
    :return: FileStoreID for reference dictionary
    :rtype: str
    """
    job.fileStore.logToMaster('Created reference dictionary')
    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.readGlobalFile(ref_id, os.path.join(work_dir, 'ref.fasta'))
    command = ['CreateSequenceDictionary', 'R=ref.fasta', 'O=ref.dict']
    docker_call(work_dir=work_dir, parameters=command,
                tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e')
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'ref.dict'))


def run_gatk_preprocessing(job, cores, bam, bai, ref, ref_dict, fai, phase, mills, dbsnp, mem='10G', unsafe=False):
    """
    Convenience method for grouping together GATK preprocessing

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param int cores: Maximum number of cores on a worker node
    :param str bam: Sample BAM FileStoreID
    :param str bai: Bam Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str ref_dict: Reference dictionary FileStoreID
    :param str fai: Reference index FileStoreID
    :param str phase: Phase VCF FileStoreID
    :param str mills: Mills VCF FileStoreID
    :param str dbsnp: DBSNP VCF FileStoreID
    :param str mem: Memory value to be passed to children. Needed for CI tests
    :param bool unsafe: If True, runs gatk UNSAFE mode: "-U ALLOW_SEQ_DICT_INCOMPATIBILITY"
    :return: BAM and BAI FileStoreIDs from Print Reads
    :rtype: tuple(str, str)
    """
    rtc = job.wrapJobFn(run_realigner_target_creator, cores, bam, bai, ref, ref_dict, fai, phase, mills, mem, unsafe)
    ir = job.wrapJobFn(run_indel_realignment, rtc.rv(), bam, bai, ref, ref_dict, fai, phase, mills, mem, unsafe)
    br = job.wrapJobFn(run_base_recalibration, cores, ir.rv(0), ir.rv(1), ref, ref_dict, fai, dbsnp, mem, unsafe)
    pr = job.wrapJobFn(run_print_reads, cores, br.rv(), ir.rv(0), ir.rv(1), ref, ref_dict, fai, mem, unsafe)
    # Wiring
    job.addChild(rtc)
    rtc.addChild(ir)
    ir.addChild(br)
    br.addChild(pr)
    return pr.rv(0), pr.rv(1)


def run_realigner_target_creator(job, cores, bam, bai, ref, ref_dict, fai, phase, mills, mem, unsafe=False):
    """
    Creates intervals file needed for indel realignment

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param int cores: Maximum number of cores on a worker node
    :param str bam: Sample BAM FileStoreID
    :param str bai: Bam Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str ref_dict: Reference dictionary FileStoreID
    :param str fai: Reference index FileStoreID
    :param str phase: Phase VCF FileStoreID
    :param str mills: Mills VCF FileStoreID
    :param str mem: Memory value to be passed to children. Needed for CI tests
    :param bool unsafe: If True, runs gatk UNSAFE mode: "-U ALLOW_SEQ_DICT_INCOMPATIBILITY"
    :return: FileStoreID for the processed bam
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_ids = [ref, fai, ref_dict, bam, bai, phase, mills]
    inputs = ['ref.fasta', 'ref.fasta.fai', 'ref.dict', 'sample.bam', 'sample.bam.bai', 'phase.vcf', 'mills.vcf']
    for file_store_id, name in zip(file_ids, inputs):
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))
    # Call: GATK -- RealignerTargetCreator
    parameters = ['-T', 'RealignerTargetCreator',
                  '-nt', str(cores),
                  '-R', '/data/ref.fasta',
                  '-I', '/data/sample.bam',
                  '-known', '/data/phase.vcf',
                  '-known', '/data/mills.vcf',
                  '--downsampling_type', 'NONE',
                  '-o', '/data/sample.intervals']
    if unsafe:
        parameters.extend(['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'])
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs,
                outputs={'sample.intervals': None},
                work_dir=work_dir, parameters=parameters, env=dict(JAVA_OPTS='-Xmx{}'.format(mem)))
    # Write to fileStore
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sample.intervals'))


def run_indel_realignment(job, intervals, bam, bai, ref, ref_dict, fai, phase, mills, mem, unsafe=False):
    """
    Creates realigned bams using the intervals file from previous step

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str intervals: Indel interval FileStoreID
    :param str bam: Sample BAM FileStoreID
    :param str bai: Bam Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str ref_dict: Reference dictionary FileStoreID
    :param str fai: Reference index FileStoreID
    :param str phase: Phase VCF FileStoreID
    :param str mills: Mills VCF FileStoreID
    :param str mem: Memory value to be passed to children. Needed for CI tests
    :param bool unsafe: If True, runs gatk UNSAFE mode: "-U ALLOW_SEQ_DICT_INCOMPATIBILITY"
    :return: FileStoreID for the processed bam
    :rtype: tuple(str, str)
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_ids = [ref, fai, ref_dict, intervals, bam, bai, phase, mills]
    inputs = ['ref.fasta', 'ref.fasta.fai', 'ref.dict', 'sample.intervals',
              'sample.bam', 'sample.bam.bai', 'phase.vcf', 'mills.vcf']
    for file_store_id, name in zip(file_ids, inputs):
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))
    # Call: GATK -- IndelRealigner
    parameters = ['-T', 'IndelRealigner',
                  '-R', '/data/ref.fasta',
                  '-I', '/data/sample.bam',
                  '-known', '/data/phase.vcf',
                  '-known', '/data/mills.vcf',
                  '-targetIntervals', '/data/sample.intervals',
                  '--downsampling_type', 'NONE',
                  '-maxReads', str(720000), # Taken from MC3 pipeline
                  '-maxInMemory', str(5400000), # Taken from MC3 pipeline
                  '-o', '/data/sample.indel.bam']
    if unsafe:
        parameters.extend(['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'])
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs,
                outputs={'sample.indel.bam': None, 'sample.indel.bai': None},
                work_dir=work_dir, parameters=parameters, env=dict(JAVA_OPTS='-Xmx{}'.format(mem)))
    # Write to fileStore
    indel_bam = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sample.indel.bam'))
    indel_bai = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sample.indel.bai'))
    return indel_bam, indel_bai


def run_base_recalibration(job, cores, indel_bam, indel_bai, ref, ref_dict, fai, dbsnp, mem, unsafe=False):
    """
    Creates recal table used in Base Quality Score Recalibration

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param int cores: Maximum number of cores on a worker node
    :param str indel_bam: Indel interval FileStoreID
    :param str indel_bai: Bam Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str ref_dict: Reference dictionary FileStoreID
    :param str fai: Reference index FileStoreID
    :param str dbsnp: DBSNP VCF FileStoreID
    :param str mem: Memory value to be passed to children. Needed for CI tests
    :param bool unsafe: If True, runs gatk UNSAFE mode: "-U ALLOW_SEQ_DICT_INCOMPATIBILITY"
    :return: FileStoreID for the processed bam
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_ids = [ref, fai, ref_dict, indel_bam, indel_bai, dbsnp]
    inputs = ['ref.fasta', 'ref.fasta.fai', 'ref.dict', 'sample.indel.bam', 'sample.indel.bai', 'dbsnp.vcf']
    for file_store_id, name in zip(file_ids, inputs):
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))
    # Call: GATK -- IndelRealigner
    parameters = ['-T', 'BaseRecalibrator',
                  '-nct', str(cores),
                  '-R', '/data/ref.fasta',
                  '-I', '/data/sample.indel.bam',
                  '-knownSites', '/data/dbsnp.vcf',
                  '-o', '/data/sample.recal.table']
    if unsafe:
        parameters.extend(['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'])
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs,
                outputs={'sample.recal.table': None},
                work_dir=work_dir, parameters=parameters, env=dict(JAVA_OPTS='-Xmx{}'.format(mem)))
    # Write output to file store
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sample.recal.table'))


def run_print_reads(job, cores, table, indel_bam, indel_bai, ref, ref_dict, fai, mem, unsafe=False):
    """
    Creates BAM that has had the base quality scores recalibrated

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param int cores: Maximum number of cores on host node
    :param str table: Recalibration table FileStoreID
    :param str indel_bam: Indel interval FileStoreID
    :param str indel_bai: Bam Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str ref_dict: Reference dictionary FileStoreID
    :param str fai: Reference index FileStoreID
    :param str mem: Memory value to be passed to children. Needed for CI tests
    :param bool unsafe: If True, runs gatk UNSAFE mode: "-U ALLOW_SEQ_DICT_INCOMPATIBILITY"
    :return: FileStoreID for the processed bam
    :rtype: tuple(str, str)
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_ids = [ref, fai, ref_dict, table, indel_bam, indel_bai]
    inputs = ['ref.fasta', 'ref.fasta.fai', 'ref.dict', 'sample.recal.table',
              'sample.indel.bam', 'sample.indel.bai']
    for file_store_id, name in zip(file_ids, inputs):
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))
    # Call: GATK -- PrintReads
    parameters = ['-T', 'PrintReads',
                  '-nct', str(cores),
                  '-R', '/data/ref.fasta',
                  '--emit_original_quals',
                  '-I', '/data/sample.indel.bam',
                  '-BQSR', '/data/sample.recal.table',
                  '-o', '/data/sample.bqsr.bam']
    if unsafe:
        parameters.extend(['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'])
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs,
                outputs={'sample.bqsr.bam': None, 'sample.bqsr.bai': None},
                work_dir=work_dir, parameters=parameters, env=dict(JAVA_OPTS='-Xmx{}'.format(mem)))
    # Write ouptut to file store
    bam_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sample.bqsr.bam'))
    bai_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sample.bqsr.bai'))
    return bam_id, bai_id
