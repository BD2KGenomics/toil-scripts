import os

import subprocess

from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.urls import download_url


def run_star(job, cores, r1_id, r2_id, star_index_url):
    """
    Performs alignment of fastqs to bam via STAR

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param int cores: Number of cores to run star with
    :param str r1_id: FileStoreID of fastq (pair 1)
    :param str r2_id: FileStoreID of fastq (pair 2 if applicable, else pass None)
    :param str star_index_url: STAR index tarball
    :param bool paired: If True, treats the sample as
    :return: FileStoreID from RSEM
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    download_url(url=star_index_url, name='starIndex.tar.gz', work_dir=work_dir)
    subprocess.check_call(['tar', '-xvf', os.path.join(work_dir, 'starIndex.tar.gz'), '-C', work_dir])
    os.remove(os.path.join(work_dir, 'starIndex.tar.gz'))
    # Determine tarball structure - star index contains are either in a subdir or in the tarball itself
    star_index = os.path.join('/data', os.listdir(work_dir)[0]) if len(os.listdir(work_dir)) == 1 else '/data'
    # Parameter handling for paired / single-end data
    parameters = ['--runThreadN', str(cores),
                  '--genomeDir', star_index,
                  '--outFileNamePrefix', 'rna',
                  '--outSAMtype', 'BAM', 'SortedByCoordinate',
                  '--outSAMunmapped', 'Within',
                  '--quantMode', 'TranscriptomeSAM',
                  '--outSAMattributes', 'NH', 'HI', 'AS', 'NM', 'MD',
                  '--outFilterType', 'BySJout',
                  '--outFilterMultimapNmax', '20',
                  '--outFilterMismatchNmax', '999',
                  '--outFilterMismatchNoverReadLmax', '0.04',
                  '--alignIntronMin', '20',
                  '--alignIntronMax', '1000000',
                  '--alignMatesGapMax', '1000000',
                  '--alignSJoverhangMin', '8',
                  '--alignSJDBoverhangMin', '1',
                  '--sjdbScore', '1']
    if r1_id and r2_id:
        job.fileStore.readGlobalFile(r1_id, os.path.join(work_dir, 'R1.fastq'))
        job.fileStore.readGlobalFile(r2_id, os.path.join(work_dir, 'R2.fastq'))
        parameters.extend(['--readFilesIn', '/data/R1.fastq', '/data/R2.fastq'])
    else:
        job.fileStore.readGlobalFile(r1_id, os.path.join(work_dir, 'R1_cutadapt.fastq'))
        parameters.extend(['--readFilesIn', '/data/R1.fastq'])
    # Call: STAR Mapping
    docker_call(tool='quay.io/ucsc_cgl/star:2.4.2a--bcbd5122b69ff6ac4ef61958e47bde94001cfe80',
                work_dir=work_dir, parameters=parameters)
    # Write to fileStore
    transcriptome_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rnaAligned.toTranscriptome.out.bam'))
    sorted_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rnaAligned.sortedByCoord.out.bam'))
    return transcriptome_id, sorted_id