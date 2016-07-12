import os
from glob import glob

from toil_scripts.tools import get_mean_insert_size
from toil_scripts.lib.files import tarball_files
from toil_scripts.lib.programs import docker_call


def run_mutect(job, normal_bam, normal_bai, tumor_bam, tumor_bai, ref, ref_dict, fai, cosmic, dbsnp):
    """
    Calls MuTect to perform variant analysis

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str normal_bam: Normal BAM FileStoreID
    :param str normal_bai: Normal BAM index FileStoreID
    :param str tumor_bam: Tumor BAM FileStoreID
    :param str tumor_bai: Tumor BAM Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str ref_dict: Reference dictionary FileStoreID
    :param str fai: Reference index FileStoreID
    :param str cosmic: Cosmic VCF FileStoreID
    :param str dbsnp: DBSNP VCF FileStoreID
    :return: MuTect output (tarball) FileStoreID
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_ids = [normal_bam, normal_bai, tumor_bam, tumor_bai, ref, fai, ref_dict, cosmic, dbsnp]
    file_names = ['normal.bam', 'normal.bai', 'tumor.bam', 'tumor.bai', 'ref.fasta',
                  'ref.fasta.fai', 'ref.dict', 'cosmic.vcf', 'dbsnp.vcf']
    for file_store_id, name in zip(file_ids, file_names):
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))
    # Call: MuTect
    parameters = ['--analysis_type', 'MuTect',
                  '--reference_sequence', 'ref.fasta',
                  '--cosmic', '/data/cosmic.vcf',
                  '--dbsnp', '/data/dbsnp.vcf',
                  '--input_file:normal', '/data/normal.bam',
                  '--input_file:tumor', '/data/tumor.bam',
                  '--tumor_lod', str(10), # Taken from MC3 pipeline
                  '--initial_tumor_lod', str(4.0), # Taken from MC3 pipeline
                  '--out', 'mutect.out',
                  '--coverage_file', 'mutect.cov',
                  '--vcf', 'mutect.vcf']
    docker_call(work_dir=work_dir, parameters=parameters,
                tool='quay.io/ucsc_cgl/mutect:1.1.7--e8bf09459cf0aecb9f55ee689c2b2d194754cbd3')
    # Write output to file store
    output_file_names = ['mutect.vcf', 'mutect.cov', 'mutect.out']
    output_file_paths = [os.path.join(work_dir, x) for x in output_file_names]
    tarball_files('mutect.tar.gz', file_paths=output_file_paths, output_dir=work_dir)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'mutect.tar.gz'))


def run_muse(job, cores, normal_bam, normal_bai, tumor_bam, tumor_bai, ref, ref_dict, fai, dbsnp):
    """
    Calls MuSe to find variants

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param int cores: Maximum number of cores on host node
    :param str normal_bam: Normal BAM FileStoreID
    :param str normal_bai: Normal BAM index FileStoreID
    :param str tumor_bam: Tumor BAM FileStoreID
    :param str tumor_bai: Tumor BAM Index FileStoreID
    :param str tumor_bai: Tumor BAM Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str ref_dict: Reference genome dictionary FileStoreID
    :param str fai: Reference index FileStoreID
    :param str dbsnp: DBSNP VCF FileStoreID
    :return: MuSe output (tarball) FileStoreID
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_ids = [normal_bam, normal_bai, tumor_bam, tumor_bai, ref, ref_dict, fai, dbsnp]
    file_names = ['normal.bam', 'normal.bai', 'tumor.bam', 'tumor.bai',
                  'ref.fasta', 'ref.dict', 'ref.fasta.fai', 'dbsnp.vcf']
    for file_store_id, name in zip(file_ids, file_names):
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))
    # Call: MuSE
    parameters = ['--mode', 'wxs',
                  '--dbsnp', '/data/dbsnp.vcf',
                  '--fafile', '/data/ref.fasta',
                  '--tumor-bam', '/data/tumor.bam',
                  '--tumor-bam-index', '/data/tumor.bai',
                  '--normal-bam', '/data/normal.bam',
                  '--normal-bam-index', '/data/normal.bai',
                  '--outfile', '/data/muse.vcf',
                  '--cpus', str(cores)]
    docker_call(tool='quay.io/ucsc_cgl/muse:1.0--6add9b0a1662d44fd13bbc1f32eac49326e48562',
                work_dir=work_dir, parameters=parameters)
    # Return fileStore ID
    tarball_files('muse.tar.gz', file_paths=[os.path.join(work_dir, 'muse.vcf')], output_dir=work_dir)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'muse.tar.gz'))


def run_pindel(job, cores, normal_bam, normal_bai, tumor_bam, tumor_bai, ref, fai):
    """
    Calls Pindel to compute indels / deletions

    :param JobFunctionWrappingJob job: Passed automatically by Toil
    :param int cores: Maximum number of cores on host node
    :param str normal_bam: Normal BAM FileStoreID
    :param str normal_bai: Normal BAM index FileStoreID
    :param str tumor_bam: Tumor BAM FileStoreID
    :param str tumor_bai: Tumor BAM Index FileStoreID
    :param str ref: Reference genome FileStoreID
    :param str fai: Reference index FileStoreID
    :return: Pindel output (tarball) FileStoreID
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_ids = [normal_bam, normal_bai, tumor_bam, tumor_bai, ref, fai]
    file_names = ['normal.bam', 'normal.bai', 'tumor.bam', 'tumor.bai', 'ref.fasta', 'ref.fasta.fai']
    for file_store_id, name in zip(file_ids, file_names):
        job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, name))
    # Create Pindel config
    with open(os.path.join(work_dir, 'pindel-config.txt'), 'w') as f:
        for bam in ['normal', 'tumor']:
            f.write('/data/{} {} {}\n'.format(bam + '.bam', get_mean_insert_size(work_dir, bam + '.bam'), bam))
    # Call: Pindel
    parameters = ['-f', '/data/ref.fasta',
                  '-i', '/data/pindel-config.txt',
                  '--number_of_threads', str(cores),
                  '--minimum_support_for_event', '3',
                  '--report_long_insertions', 'true',
                  '--report_breakpoints', 'true',
                  '-o', 'pindel']
    docker_call(tool='quay.io/ucsc_cgl/pindel:0.2.5b6--4e8d1b31d4028f464b3409c6558fb9dfcad73f88',
                work_dir=work_dir, parameters=parameters)
    # Collect output files and write to file store
    output_files = glob(os.path.join(work_dir, 'pindel*'))
    tarball_files('pindel.tar.gz', file_paths=output_files, output_dir=work_dir)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'pindel.tar.gz'))
