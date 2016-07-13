import os

import subprocess

from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.urls import download_url


def run_star(job, cores, r1_id, r2_id, star_index_url, wiggle=False):
    """
    Performs alignment of fastqs to bam via STAR

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param int cores: Number of cores to run star with
    :param str r1_id: FileStoreID of fastq (pair 1)
    :param str r2_id: FileStoreID of fastq (pair 2 if applicable, else pass None)
    :param str star_index_url: STAR index tarball
    :param bool wiggle: If True, will output a wiggle file and return it
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
    if wiggle:
        parameters.extend(['--outWigType', 'bedGraph',
                           '--outWigStrand', 'Unstranded',
                           '--outWigReferencesPrefix', 'chr'])
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
    if wiggle:
        wiggle_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rnaSignal.UniqueMultiple.str1.out.bg'))
        return transcriptome_id, sorted_id, wiggle_id
    else:
        return transcriptome_id, sorted_id


def run_bwakit(job, config, threads, sort=True, trim=False):
    """
    Runs BWA-Kit to align a fastq file or fastq pair into a BAM file.

    :param JobFunctionWrappingJob job: Passed by Toil automatically
    :param Namespace config: A configuration object that holds strings as attributes.
        The attributes must be accessible via the dot operator.
        The config must have:
        config.r1               FileStoreID for 1st fastq file
        config.r2               FileStoreID for 2nd fastq file (or None if single-ended)
        config.ref              FileStoreID for the reference genome
        config.fai              FileStoreID for the reference index file
        config.amb              FileStoreID for the reference amb file
        config.ann              FileStoreID for the reference ann file
        config.bwt              FileStoreID for the reference bwt file
        config.pac              FileStoreID for the reference pac file
        config.sa               FileStoreID for the reference sa file
        config.alt              FileStoreID for the reference alt (or None)
        config.rg_line          The read group value to use (or None -- see below)
        config.library          Read group attribute: library
        config.platform         Read group attribute: platform
        config.program_unit     Read group attribute: program unit
        config.uuid             Read group attribute: sample ID

        If specifying config.rg_line, use the following format:
            BAM read group header line (@RG), as defined on page 3 of the SAM spec.
            Tabs should be escaped, e.g., @RG\\tID:foo\\tLB:bar...
            for the read group "foo" from sequencing library "bar".
            Multiple @RG lines can be defined, but should be split by an escaped newline \\n,
            e.g., @RG\\tID:foo\\t:LB:bar\\n@RG\\tID:santa\\tLB:cruz.

    :param int threads: Number of threads to use
    :param bool sort: If True, sorts the BAM
    :param bool trim: If True, performs adapter trimming
    :return: FileStoreID of BAM
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_names = ['r1.fq.gz', 'ref.fa.fai', 'ref.fa', 'ref.fa.amb', 'ref.fa.ann',
                  'ref.fa.bwt', 'ref.fa.pac', 'ref.fa.sa']
    ids = [config.r1, config.ref, config.fai, config.amb, config.ann, config.bwt, config.pac, config.sa]
    # If a fastq pair was provided
    if getattr(config, 'r2', None):
        file_names.insert(1, 'r2.fq.gz')
        ids.insert(1, config.r2)
    # If an alt file was provided
    if getattr(config, 'alt', None):
        file_names.append('ref.fa.alt')
        ids.append(config.alt)
    for fileStoreID, name in zip(ids, file_names):
        job.fileStore.readGlobalFile(fileStoreID, os.path.join(work_dir, name))
    # If a read group line was provided
    if getattr(config, 'rg_line', None):
        rg = config.rg_line
    # Otherwise, generate a read group line to place in the BAM.
    else:
        rg = "@RG\\tID:{0}".format(config.uuid)  # '\' character is escaped so bwakit gets passed '\t' properly
        rg_attributes = [config.library, config.platform, config.program_unit, config.uuid]
        for tag, info in zip(['LB', 'PL', 'PU', 'SM'], rg_attributes):
            rg += '\\t{0}:{1}'.format(tag, info)
    # BWA Options
    opt_args = []
    if sort:
        opt_args.append('-s')
    if trim:
        opt_args.append('-a')
    # Call: bwakit
    parameters = (['-t', str(threads),
                   '-R', rg] +
                  opt_args +
                  ['-o', '/data/aligned',
                   '/data/ref.fa',
                   '/data/r1.fq.gz'])
    if getattr(config, 'r2', None):  # If a fastq pair was provided
        parameters.append('/data/r2.fq.gz')
    mock_bam = config.uuid + '.bam'
    outputs = {'aligned.aln.bam': mock_bam}
    docker_call(tool='quay.io/ucsc_cgl/bwakit:0.7.12--528bb9bf73099a31e74a7f5e6e3f2e0a41da486e',
                parameters=parameters, inputs=file_names, outputs=outputs, work_dir=work_dir)

    # Either write file to local output directory or upload to S3 cloud storage
    job.fileStore.logToMaster('Aligned sample: {}'.format(config.uuid))
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'aligned.aln.bam'))
