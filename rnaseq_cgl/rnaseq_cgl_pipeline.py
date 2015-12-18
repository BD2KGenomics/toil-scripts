#!/usr/bin/env python2.7
"""
UCSC Computational Genomics Lab RNA-seq Pipeline
Author: John Vivian
Affiliation: UC Santa Cruz Genomics Institute

Please see the README.md in the same directory

Structure of RNA-Seq Pipeline (per sample)

                 0
                 |
                 1
                 |
                 2 - - - - > 7 -> 8
                / \
               3   6
               |
               4
               |
               5

0 = Download sample
1 = Merge fastqs
2 = CutAdapt
3 = STAR Alignment
4 = RSEM Quantification
5 = RSEM Post-processing
6 = Kallisto
7 = Consoliate output (into a tarball)
8 = upload results to S3 (optional)
=======================================
Dependencies
Curl:       apt-get install curl
Docker:     apt-get install docker.io # docker.io if using linux, o.w. just docker
Toil:       pip install git+https://github.com/BD2KGenomics/toil.git

Optional (if uploading results to S3)
Boto:       pip install boto
"""
import os
import argparse
import base64
from contextlib import closing
import glob
import hashlib
import errno
import subprocess
import shutil
import tarfile
import multiprocessing
from urlparse import urlparse
from toil.job import Job


def build_parser():
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-c', '--config', default=None,
                       help='Path to CSV. One sample per line with the format: uuid,url (see example config). '
                            'URLs follow the format of: https://www.address.com/file.tar or file:///full/path/file.tar.'
                            ' Samples must be tarfiles that contain fastq files.')
    group.add_argument('-d', '--dir', default=None,
                       help='Path to directory of samples. Samples must be tarfiles that contain fastq files. '
                            'The UUID for the sample will be derived from the file basename.')
    group.add_argument('-s', '--sample_urls', nargs='?', default=None, type=str,
                       help='Sample URLs (any number). Samples must be tarfiles that contain fastq files. '
                            'URLs follow the format of: https://www.address.com/file.tar or file:///full/path/file.tar.'
                            'The UUID for the sample will be derived from the file basename.')
    parser.add_argument('--starIndex', help='URL to download STAR Index built from HG38/gencodev23 annotation.',
                        default='https://s3-us-west-2.amazonaws.com/'
                                'cgl-pipeline-inputs/rnaseq_cgl/starIndex_hg38_no_alt.tar.gz')
    parser.add_argument('--kallistoIndex', help='URL to download Kallisto Index built from HG38 transcriptome.',
                        default='https://s3-us-west-2.amazonaws.com/'
                                'cgl-pipeline-inputs/rnaseq_cgl/kallisto_hg38.idx')
    parser.add_argument('--rsemRef', help='URL to download RSEM reference built from HG38/gencodev23.',
                        default='https://s3-us-west-2.amazonaws.com/'
                                'cgl-pipeline-inputs/rnaseq_cgl/rsem_ref_hg38_no_alt.tar.gz')
    parser.add_argument('--fwd_3pr_adapter', help="Sequence for the FWD 3' Read Adapter.", default='AGATCGGAAGAG')
    parser.add_argument('--rev_3pr_adapter', help="Sequence for the REV 3' Read Adapter.", default='AGATCGGAAGAG')
    parser.add_argument('--ssec', default=None, help='Path to Key File for SSE-C Encryption')
    parser.add_argument('--output_dir', default=None, help='full path where final results will be output')
    parser.add_argument('--s3_dir', default=None, help='S3 Directory, starting with bucket name. e.g.: '
                                                       'cgl-driver-projects/ckcc/rna-seq-samples/')
    parser.add_argument('--sudo', dest='sudo', action='store_true',
                        help='Docker usually needs sudo to execute locally, but not when running Mesos or when '
                             'the user is a member of a Docker group.')
    parser.set_defaults(sudo=False)
    return parser


# Convenience functions used in the pipeline
def mkdir_p(path):
    """
    It is Easier to Ask for Forgiveness than Permission
    """
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def flatten(x):
    """
    Flattens a nested array into a single list

    x: list/tuple       The nested list/tuple to be flattened.
    """
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result


def humansize(nbytes):
    """
    Returns human readable format
    """
    suffixes = ['B', 'K', 'M', 'G', 'T', 'P']
    if nbytes == 0:
        return '0 B'
    i = 0
    while nbytes >= 1024 and i < len(suffixes)-1:
        nbytes /= 1024.
        i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])


def generate_unique_key(master_key_path, url):
    """
    master_key_path: str    Path to the BD2K Master Key (for S3 Encryption)
    url: str                S3 URL (e.g. https://s3-us-west-2.amazonaws.com/bucket/file.txt)

    Returns: str            32-byte unique key generated for that URL
    """
    with open(master_key_path, 'r') as f:
        master_key = f.read()
    assert len(master_key) == 32, 'Invalid Key! Must be 32 characters. ' \
                                  'Key: {}, Length: {}'.format(master_key, len(master_key))
    new_key = hashlib.sha256(master_key + url).digest()
    assert len(new_key) == 32, 'New key is invalid and is not 32 characters: {}'.format(new_key)
    return new_key


def download_encrypted_file(job, input_args, name):
    """
    Downloads encrypted files from S3 via header injection

    input_args: dict    Input dictionary defined in main()
    name: str           Symbolic name associated with file
    """
    work_dir = job.fileStore.getLocalTempDir()
    key_path = input_args['ssec']
    file_path = os.path.join(work_dir, name)
    url = input_args[name]

    with open(key_path, 'r') as f:
        key = f.read()
    if len(key) != 32:
        raise RuntimeError('Invalid Key! Must be 32 bytes: {}'.format(key))

    key = generate_unique_key(key_path, url)

    encoded_key = base64.b64encode(key)
    encoded_key_md5 = base64.b64encode(hashlib.md5(key).digest())
    h1 = 'x-amz-server-side-encryption-customer-algorithm:AES256'
    h2 = 'x-amz-server-side-encryption-customer-key:{}'.format(encoded_key)
    h3 = 'x-amz-server-side-encryption-customer-key-md5:{}'.format(encoded_key_md5)
    try:
        subprocess.check_call(['curl', '-fs', '--retry', '5', '-H', h1, '-H', h2, '-H', h3, url, '-o', file_path])
    except OSError:
        raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)
    return job.fileStore.writeGlobalFile(file_path)


def download_from_url(job, url):
    """
    Simple curl request made for a given url

    url: str    URL to download
    """
    work_dir = job.fileStore.getLocalTempDir()
    file_path = os.path.join(work_dir, os.path.basename(url))
    if not os.path.exists(file_path):
        try:
            subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', url, '-o', file_path])
        except OSError:
            raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)
    return job.fileStore.writeGlobalFile(file_path)


def read_from_filestore(job, work_dir, ids, *filenames):
    """
    Reads file from fileStore and writes it to working directory.


    work_dir: str       working directory
    ids: dict           dict of fileStore IDs
    *filenames: str     filenames to be read from jobStore
    """
    for filename in filenames:
        if not os.path.exists(os.path.join(work_dir, filename)):
            job.fileStore.readGlobalFile(ids[filename], os.path.join(work_dir, filename))


def docker_call(work_dir, tool_parameters, tool, java_opts=None, sudo=False, outfile=None):
    """
    Makes subprocess call of a command to a docker container.


    tool_parameters: list   An array of the parameters to be passed to the tool
    tool: str               Name of the Docker image to be used (e.g. quay.io/ucsc_cgl/samtools)
    java_opts: str          Optional commands to pass to a java jar execution. (e.g. '-Xmx15G')
    outfile: file           Filehandle that stderr will be passed to
    sudo: bool              If the user wants the docker command executed as sudo
    """
    base_docker_call = 'docker run --log-driver=none --rm -v {}:/data'.format(work_dir).split()
    if sudo:
        base_docker_call = ['sudo'] + base_docker_call
    if java_opts:
        base_docker_call = base_docker_call + ['-e', 'JAVA_OPTS={}'.format(java_opts)]
    try:
        if outfile:
            subprocess.check_call(base_docker_call + [tool] + tool_parameters, stdout=outfile)
        else:
            subprocess.check_call(base_docker_call + [tool] + tool_parameters)
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status. Check error logs.')
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


def copy_to_output_dir(work_dir, output_dir, uuid=None, files=list()):
    """
    A list of files to move from work_dir to output_dir.

    work_dir: str       Current working directory
    output_dir: str     Output directory for files to go
    uuid: str           UUID to "stamp" onto output files
    files: list         List of files to iterate through
    """
    for fname in files:
        if uuid is None:
            shutil.copy(os.path.join(work_dir, fname), os.path.join(output_dir, fname))
        else:
            shutil.copy(os.path.join(work_dir, fname), os.path.join(output_dir, '{}.{}'.format(uuid, fname)))


def tarball_files(work_dir, tar_name, uuid=None, files=None):
    """
    Tars a group of files together into a tarball

    work_dir: str       Current Working Directory
    tar_name: str       Name of tarball
    uuid: str           UUID to stamp files with
    files: str(s)       List of filenames to place in the tarball from working directory
    """
    with tarfile.open(os.path.join(work_dir, tar_name), 'w:gz') as f_out:
        for fname in files:
            if uuid:
                f_out.add(os.path.join(work_dir, fname), arcname=uuid + '.' + fname)
            else:
                f_out.add(os.path.join(work_dir, fname), arcname=fname)


def parse_config(path_to_config):
    """
    Parses config file. Returns list of samples: [ [uuid1, url1], [uuid2, url2] ... ]
    """
    samples = []
    with open(path_to_config, 'r') as f:
        for line in f.readlines():
            if not line.isspace():
                sample = line.strip().split(',')
                assert len(sample) == 2, 'Error: Config file is inappropriately formatted. Read documentation.'
                assert sample[1].endswith('.tar'), 'Error: Samples must be tarfiles: {}'.format(sample[1])
                samples.append(sample)
    return samples


# Job Functions
def download_shared_files(job, input_args):
    """
    Downloads and stores shared input files in the FileStore

    input_args: dict        Dictionary of input arguments (from main())
    """
    shared_files = ['starIndex.tar.gz', 'rsem_ref_hg38.tar.gz', 'kallisto_hg38.idx']
    shared_ids = {}
    for f in shared_files:
        shared_ids[f] = job.addChildJobFn(download_from_url, input_args[f], disk='25G').rv()
    job.addFollowOnJobFn(parse_input_samples, shared_ids, input_args)


def parse_input_samples(job, shared_ids, input_args):
    """
    Launches pipeline for each sample.

    shared_ids: dict        Dictionary of fileStore IDs
    input_args: dict        Dictionary of input arguments
    """
    config = input_args['config']
    sample_dir = input_args['dir']
    sample_urls = input_args['sample_urls']
    samples = None
    if config:
        samples = parse_config(config)
    if sample_dir:
        files = os.listdir(sample_dir)
        samples = [[os.path.splitext(os.path.basename(f))[0], 'file://' + os.path.abspath(f)] for f in files]
    if sample_urls:
        samples = [[os.path.splitext(os.path.basename(f))[0], f] for f in sample_urls]
    for sample in samples:
        job.addChildJobFn(download_sample, shared_ids, input_args, sample)


def download_sample(job, ids, input_args, sample):
    """
    Defines variables unique to a sample that are used in the rest of the pipelines

    ids: dict           Dictionary of fileStore IDS
    input_args: dict    Dictionary of input arguments
    sample: tuple       Contains uuid and sample_url
    """
    uuid, url = sample
    # Update values unique to sample
    sample_input = dict(input_args)
    sample_input['uuid'] = uuid
    sample_input['sample.tar'] = url
    if sample_input['output_dir']:
        sample_input['output_dir'] = os.path.join(input_args['output_dir'], uuid)
    sample_input['cpu_count'] = multiprocessing.cpu_count()
    job_vars = (sample_input, ids)
    # Download or locate local file and place in the jobStore
    if urlparse(url).scheme == 'file':
        ids['sample.tar'] = job.fileStore.writeGlobalFile(urlparse(url).path)
    else:
        if sample_input['ssec']:
            ids['sample.tar'] = job.addChildJobFn(download_encrypted_file, sample_input, 'sample.tar', disk='25G').rv()
        else:
            ids['sample.tar'] = job.addChildJobFn(download_from_url, sample_input['sample.tar'], disk='25G').rv()
    job.addFollowOnJobFn(static_dag_launchpoint, job_vars)


def static_dag_launchpoint(job, job_vars):
    """
    Statically define pipeline

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    a = job.wrapJobFn(merge_fastqs, job_vars, disk='70G').encapsulate()
    b = job.wrapJobFn(consolidate_output, job_vars, a.rv(), disk='2G')
    # Take advantage of "encapsulate" to simplify pipeline wiring
    job.addChild(a)
    a.addChild(b)


def merge_fastqs(job, job_vars):
    """
    Untars input sample and concats the Read1 and Read2 groups together.

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    # Unpack variables
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    # I/O
    read_from_filestore(job, work_dir, ids, 'sample.tar')
    sample_tar = os.path.join(work_dir, 'sample.tar')
    # Untar File and concat
    subprocess.check_call(['tar', '-xvf', sample_tar, '-C', work_dir])
    os.remove(os.path.join(work_dir, 'sample.tar'))
    # TODO: Change for TCGA data _1 and _2
    r1_files = sorted(glob.glob(os.path.join(work_dir, '*R1*')))
    r2_files = sorted(glob.glob(os.path.join(work_dir, '*R2*')))
    with open(os.path.join(work_dir, 'R1.fastq'), 'w') as f1:
        p1 = subprocess.Popen(['zcat'] + r1_files, stdout=f1)
    with open(os.path.join(work_dir, 'R2.fastq'), 'w') as f2:
        p2 = subprocess.Popen(['zcat'] + r2_files, stdout=f2)
    p1.wait()
    p2.wait()
    # Write to fileStore
    ids['R1.fastq'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1.fastq'))
    ids['R2.fastq'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R2.fastq'))
    job.fileStore.deleteGlobalFile(ids['sample.tar'])
    # Start cutadapt step
    return job.addChildJobFn(cutadapt, job_vars, disk='150G').rv()


def cutadapt(job, job_vars):
    """
    Filters out adapters that may be left in the RNA-seq files

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    # Unpack variables
    input_args, ids = job_vars
    sudo = input_args['sudo']
    cores = input_args['cpu_count']
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve files
    read_from_filestore(job, work_dir, ids, 'R1.fastq', 'R2.fastq')
    # Call: CutAdapt
    parameters = ['-a', input_args['fwd_3pr_adapter'],
                  '-A', input_args['rev_3pr_adapter'],
                  '-m', '35',
                  '-o', '/data/R1_cutadapt.fastq',
                  '-p', '/data/R2_cutadapt.fastq',
                  '/data/R1.fastq', '/data/R2.fastq']
    docker_call(tool='quay.io/ucsc_cgl/cutadapt:1.9--6bd44edd2b8f8f17e25c5a268fedaab65fa851d2',
                work_dir=work_dir, tool_parameters=parameters, sudo=sudo)
    # Write to fileStore
    ids['R1_cutadapt.fastq'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1_cutadapt.fastq'))
    ids['R2_cutadapt.fastq'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R2_cutadapt.fastq'))
    # start STAR and Kallisto steps
    rsem_output = job.addChildJobFn(star, job_vars, cores=cores, disk='150G', memory='40G').rv()
    kallisto_output = job.addChildJobFn(kallisto, job_vars, cores=cores, disk='100G').rv()
    return rsem_output, kallisto_output


def kallisto(job, job_vars):
    """
    Performs RNA quantification via Kallisto

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    # Unpack variables
    input_args, ids = job_vars
    sudo = input_args['sudo']
    cores = input_args['cpu_count']
    uuid = input_args['uuid']
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve files
    read_from_filestore(job, work_dir, ids, 'R1_cutadapt.fastq', 'R2_cutadapt.fastq', 'kallisto_hg38.idx')
    # Call: Kallisto
    parameters = ['quant',
                  '-i', '/data/kallisto_hg38.idx',
                  '-t', str(cores),
                  '-o', '/data/',
                  '-b', '100',
                  '/data/R1_cutadapt.fastq',
                  '/data/R2_cutadapt.fastq']
    docker_call(tool='quay.io/ucsc_cgl/kallisto:0.42.4--35ac87df5b21a8e8e8d159f26864ac1e1db8cf86',
                work_dir=work_dir, tool_parameters=parameters, sudo=sudo)
    # Tar output files together and store in fileStore
    output_files = ['run_info.json', 'abundance.tsv', 'abundance.h5']
    tarball_files(work_dir, tar_name='kallisto.tar.gz', uuid=uuid, files=output_files)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'kallisto.tar.gz'))


def star(job, job_vars):
    """
    Performs alignment of fastqs to BAM via STAR

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    # Unpack variables
    input_args, ids = job_vars
    sudo = input_args['sudo']
    cores = input_args['cpu_count']
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve files
    read_from_filestore(job, work_dir, ids, 'R1_cutadapt.fastq', 'R2_cutadapt.fastq', 'starIndex.tar.gz')
    subprocess.check_call(['tar', '-zxvf', os.path.join(work_dir, 'starIndex.tar.gz'), '-C', work_dir])
    # Call: STAR Map
    parameters = ['--runThreadN', str(cores),
                  '--genomeDir', '/data/starIndex',
                  '--readFilesIn', '/data/R1_cutadapt.fastq', '/data/R2_cutadapt.fastq',
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
                  '--sjdbScore', '1',
                  '--outBAMcompression', '-1']
    docker_call(tool='quay.io/ucsc_cgl/star:2.4.2a--bcbd5122b69ff6ac4ef61958e47bde94001cfe80',
                work_dir=work_dir, tool_parameters=parameters, sudo=sudo)
    # Write to fileStore
    ids['transcriptome.bam'] = job.fileStore.writeGlobalFile(os.path.join(work_dir,
                                                                          'rnaAligned.toTranscriptome.out.bam'))
    # RSEM doesn't use more than 16 cores, so be efficient
    cores = 16 if cores >= 16 else cores
    return job.addChildJobFn(rsem, job_vars, cores=cores, disk='50G').rv()


def rsem(job, job_vars):
    """
    Runs RSEM to produce counts

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    cores = input_args['cpu_count']
    cores = 16 if cores >= 16 else cores
    sudo = input_args['sudo']
    # I/O
    read_from_filestore(job, work_dir, ids, 'transcriptome.bam', 'rsem_ref_hg38.tar.gz')
    subprocess.check_call(['tar', '-zxvf', os.path.join(work_dir, 'rsem_ref_hg38.tar.gz'), '-C', work_dir])
    output_prefix = 'rsem'
    # Call: RSEM
    parameters = ['--quiet',
                  '--no-qualities',
                  '--paired-end',
                  '-p', str(cores),
                  '--forward-prob', '0.5',
                  '--seed-length', '25',
                  '--fragment-length-mean', '-1.0',
                  '--bam', '/data/transcriptome.bam',
                  '/data/rsem_ref_hg38/hg38',
                  output_prefix]

    docker_call(tool='quay.io/ucsc_cgl/rsem:1.2.25--d4275175cc8df36967db460b06337a14f40d2f21',
                tool_parameters=parameters, work_dir=work_dir, sudo=sudo)
    os.rename(os.path.join(work_dir, output_prefix+'.genes.results'), os.path.join(work_dir, 'rsem_gene.tab'))
    os.rename(os.path.join(work_dir, output_prefix+'.isoforms.results'), os.path.join(work_dir, 'rsem_isoform.tab'))
    # Write to FileStore
    ids['rsem_gene.tab'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem_gene.tab'))
    ids['rsem_isoform.tab'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem_isoform.tab'))
    # Run child jobs
    return job.addChildJobFn(rsem_postprocess, job_vars).rv()


def rsem_postprocess(job, job_vars):
    """
    Parses RSEMs output to produce the separate .tab files (TPM, FPKM, counts) for both gene and isoform

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    uuid = input_args['uuid']
    sudo = input_args['sudo']
    # I/O
    read_from_filestore(job, work_dir, ids, 'rsem_gene.tab', 'rsem_isoform.tab')
    # Command
    sample = input_args['uuid']
    docker_call(tool='jvivian/rsem_postprocess', tool_parameters=[sample], work_dir=work_dir, sudo=sudo)
    # Tar output files together and store in fileStore
    output_files = ['rsem.genes.norm_counts.tab', 'rsem.genes.raw_counts.tab', 'rsem.genes.norm_fpkm.tab',
                    'rsem.genes.norm_tpm.tab', 'rsem.isoform.norm_counts.tab', 'rsem.isoform.raw_counts.tab',
                    'rsem.isoform.norm_fpkm.tab', 'rsem.isoform.norm_tpm.tab']
    tarball_files(work_dir, tar_name='rsem.tar.gz', uuid=uuid, files=output_files)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem.tar.gz'))


def consolidate_output(job, job_vars, output_ids):
    """
    Combine the contents of separate zipped outputs into one via streaming

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    output_ids: tuple   Nested tuple of all the output fileStore IDs
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    uuid = input_args['uuid']
    # Retrieve IDs
    rsem_id, kallisto_id = flatten(output_ids)
    # Retrieve output file paths to consolidate
    rsem_tar = job.fileStore.readGlobalFile(rsem_id, os.path.join(work_dir, 'rsem.tar.gz'))
    kallisto_tar = job.fileStore.readGlobalFile(kallisto_id, os.path.join(work_dir, 'kallisto.tar.gz'))
    # I/O
    out_tar = os.path.join(work_dir, uuid + '.tar.gz')
    # Consolidate separate tarballs into one as streams (avoids unnecessary untaring)
    with tarfile.open(os.path.join(work_dir, out_tar), 'w:gz') as f_out:
        for tar in [rsem_tar, kallisto_tar]:
            with tarfile.open(tar, 'r') as f_in:
                for tarinfo in f_in:
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        if tar == rsem_tar:
                            tarinfo.name = os.path.join(uuid, 'RSEM', os.path.basename(tarinfo.name))
                        else:
                            tarinfo.name = os.path.join(uuid, 'Kallisto', os.path.basename(tarinfo.name))
                        f_out.addfile(tarinfo, fileobj=f_in_file)
    # Move to output directory of selected
    if input_args['output_dir']:
        output_dir = input_args['output_dir']
        mkdir_p(output_dir)
        copy_to_output_dir(work_dir, output_dir, uuid=None, files=[uuid + '.tar.gz'])
    # Write output file to fileStore
    ids['uuid.tar.gz'] = job.fileStore.writeGlobalFile(out_tar)
    # If S3 bucket argument specified, upload to S3
    if input_args['s3_dir']:
        job.addChildJobFn(upload_to_s3, job_vars)


def upload_to_s3(job, job_vars):
    """
    If s3_dir is specified in arguments, file will be uploaded to S3 using boto.
    WARNING: ~/.boto credentials are necessary for this to succeed!

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    import boto
    from boto.s3.key import Key

    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    uuid = input_args['uuid']
    # Parse s3_dir
    s3_dir = input_args['s3_dir']
    bucket_name = s3_dir.split('/')[0]
    bucket_dir = '/'.join(s3_dir.split('/')[1:])
    # I/O
    read_from_filestore(job, work_dir, ids, 'uuid.tar.gz')
    uuid_tar = os.path.join(work_dir, 'uuid.tar.gz')
    # Upload to S3 via boto
    conn = boto.connect_s3()
    bucket = conn.get_bucket(bucket_name)
    k = Key(bucket)
    k.key = os.path.join(bucket_dir, uuid + '.tar.gz')
    k.set_contents_from_filename(uuid_tar)


def main():
    """
    This is a Toil pipeline for the UCSC CGL's RNA-seq pipeline.
    RNA-seq fastqs are combined, aligned, and quantified with 2 different methods.

    Please read the README.md located in the same directory for run instructions.
    """
    # Define Parser object and add to toil
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    # Store inputs from argparse
    inputs = {'config': args.config,
              'dir': args.dir,
              'sample_urls': args.sample_urls,
              'starIndex.tar.gz': args.starIndex,
              'rsem_ref_hg38.tar.gz': args.rsemRef,
              'kallisto_hg38.idx': args.kallistoIndex,
              'fwd_3pr_adapter': args.fwd_3pr_adapter,
              'rev_3pr_adapter': args.rev_3pr_adapter,
              'output_dir': args.output_dir,
              'ssec': args.ssec,
              's3_dir': args.s3_dir,
              'sudo': args.sudo,
              'uuid': None,
              'sample.tar': None,
              'cpu_count': None}

    # Sanity checks
    if args.ssec:
        assert os.path.isfile(args.ssec)
    if args.config:
        assert os.path.isfile(args.config)
    if args.dir:
        assert os.path.isdir(args.dir)
        assert os.listdir(args.dir) != []
    if args.sample_urls:
        assert args.sample_urls != []

    # Start Pipeline
    Job.Runner.startToil(Job.wrapJobFn(download_shared_files, inputs), args)


if __name__ == '__main__':
    main()