#!/usr/bin/env python2.7
"""
UNC Best Practice RNA-Seq Pipeline
Author: John Vivian
Affiliation: UC Santa Cruz Genomics Institute

Please see the README.md in the same directory

Tree Structure of RNA-Seq Pipeline (per sample)

    0---> 2
    |     |
    1     3 - - - - -> Consolidate Output -> Upload_to_S3
         / \
       *4   5
            |
            6
           / \
         *7   8
             / \
           *9   10
                |
                11
                |
                12
                |
               *13

0 = Start Node
1 = Download Sample
2 = Unzip
3 = Mapsplice
4 = Mapping Stats (not currently included)
5 = Add Read Groups
6 = Bamsort and Index
7 = Rseq-QC
8 = Sort Bam by Reference
9 = Exon Quantification
10 = Transcriptome
11 = Filter
12 = RSEM
13 = RSEM Post-Process

7,9,13 contribute to producing the final output

Dependencies
Curl:       apt-get install curl
Docker:     apt-get install docker.io # docker.io if using linux, o.w. just docker
Samtools:   apt-get install samtools
Unzip:      apt-get install unzip
Toil:       pip install git+https://github.com/BD2KGenomics/toil.git

Optional
Boto:       pip install boto
"""
import argparse
import base64
from collections import OrderedDict
from contextlib import closing
import glob
import hashlib
import os
import subprocess
import errno
import multiprocessing
import shutil
import tarfile
from urlparse import urlparse
from toil.job import Job

from toil_scripts import download_from_s3_url


def build_parser():
    parser = argparse.ArgumentParser(description=main.__doc__, add_help=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--config', default=None, help='Path to config. One sample per line, with the format: '
                                                      'uuid,url_to_sample.tar')
    group.add_argument('--input', default=None, help='Accepts a local sample: /path/to/sample.tar. Take note! The'
                                                     'UUID for this sample is derived from the name. So samples'
                                                     'should be in the form of uuid.tar.')
    group.add_argument('-f', '--config_fastq', default=None,
                       help='Path to CSV. One sample per line with the format: '
                            'uuid,file:///path/to/R_1.fastq,file:///path/to/R_2.fastq')
    parser.add_argument('--single_end_reads', default=False, action='store_true',
                        help='Set this flag if input data is non-paired (single end reads).')
    parser.add_argument('--unc', help='URL to unc_hg19.bed',
                        default='https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rna-seq/unc_hg19.bed')
    parser.add_argument('--fasta', help='URL to hg19_M_rCRS_ref.transcripts.fa',
                        default='https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rna-seq/hg19_M_rCRS_ref.transcripts.fa')
    parser.add_argument('--composite_exons', help='URL to composite_exons.bed',
                        default='https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rna-seq/composite_exons.bed')
    parser.add_argument('--normalize', help='URL to normalizeBedToolsExonQuant.pl',
                        default='https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rna-seq/normalizeBedToolsExonQuant.pl')
    parser.add_argument('--rsem_ref', help='RSEM_REF URL',
                        default='https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rna-seq/rsem_ref.zip')
    parser.add_argument('--chromosomes', help='Chromosomes Directory',
                        default='https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rna-seq/chromosomes.zip')
    parser.add_argument('--ebwt', help='EBWT Directory',
                        default='https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/rna-seq/ebwt.zip')
    parser.add_argument('--ssec', help='Path to Key File for SSE-C Encryption')
    parser.add_argument('--output_dir', default=None, help='full path where final results will be output')
    parser.add_argument('--upload_bam_to_s3', default=False, action='store_true',
                        help='uploads alignment bam to S3 directory specified.')
    parser.add_argument('--s3_dir', default=None, help='S3 Directory, starting with bucket name. e.g.: '
                                                       'cgl-driver-projects/ckcc/rna-seq-samples/')
    parser.add_argument('--sudo', dest='sudo', action='store_true', default=False,
                        help='Docker usually needs sudo to execute locally, but not when running Mesos or when '
                             'the user is a member of a Docker group.')
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


def which(program):
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


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
        if url.startswith('s3:'):
            download_from_s3_url(file_path, url)
        else:
            try:
                subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', url, '-o', file_path])
            except OSError:
                raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')
    assert os.path.exists(file_path)
    return job.fileStore.writeGlobalFile(file_path)


def return_input_paths(job, work_dir, ids, *args):
    """
    Given one or more strings representing file_names, return the paths to those files. Each item must be unpacked!

    work_dir: str       Current working directory
    ids: dict           Dictionary of fileStore IDs
    *args: str(s)       for every file in *args, place file in work_dir via FileStore
    """
    paths = OrderedDict()
    for name in args:
        if not os.path.exists(os.path.join(work_dir, name)):
            file_path = job.fileStore.readGlobalFile(ids[name], os.path.join(work_dir, name))
        else:
            file_path = os.path.join(work_dir, name)
        paths[name] = file_path
        if len(args) == 1:
            return file_path

    return paths.values()


def docker_path(filepath):
    """
    Given a path, returns that files path inside the docker mount directory (/data).
    """
    return os.path.join('/data', os.path.basename(filepath))


def docker_call(work_dir, tool_parameters, tool, java_opts=None, outfile=None, sudo=False):
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


# Job Functions
def program_checks(job, input_args):
    """
    Checks that dependency programs are installed.

    input_args: dict        Dictionary of input arguments (from main())
    """
    # Program checks
    for program in ['curl', 'docker', 'unzip', 'samtools']:
        assert which(program), 'Program "{}" must be installed on every node.'.format(program)
    job.addChildJobFn(download_shared_files, input_args)


def download_shared_files(job, input_args):
    """
    Downloads and stores shared inputs files in the FileStore

    input_args: dict        Dictionary of input arguments (from main())
    """
    shared_files = ['unc.bed', 'hg19.transcripts.fa', 'composite_exons.bed', 'normalize.pl', 'rsem_ref.zip',
                    'ebwt.zip', 'chromosomes.zip']
    shared_ids = {}
    for f in shared_files:
        shared_ids[f] = job.addChildJobFn(download_from_url, input_args[f]).rv()
    if input_args['config'] or input_args['config_fastq']:
        job.addFollowOnJobFn(parse_config_file, shared_ids, input_args)
    else:
        sample_path = input_args['input']
        uuid = os.path.splitext(os.path.basename(sample_path))[0]
        sample = (uuid, sample_path)
        job.addFollowOnJobFn(download_sample, shared_ids, input_args, sample)


def parse_config_file(job, ids, input_args):
    """
    Launches pipeline for each sample.

    shared_ids: dict        Dictionary of fileStore IDs
    input_args: dict        Dictionary of input arguments
    """
    samples = []
    config = input_args['config']
    with open(config, 'r') as f:
        for line in f.readlines():
            if not line.isspace():
                sample = line.strip().split(',')
                samples.append(sample)
    for sample in samples:
        job.addChildJobFn(download_sample, ids, input_args, sample)


def download_sample(job, ids, input_args, sample):
    """
    Defines variables unique to a sample that are used in the rest of the pipelines

    ids: dict           Dictionary of fileStore IDS
    input_args: dict    Dictionary of input arguments
    sample: tuple       Contains uuid and sample_url
    """
    if len(sample) == 2:
        uuid, sample_location = sample
        url1, url2 = None, None
    else:
        uuid, url1, url2 = sample
        sample_location = None
    # Update values unique to sample
    sample_input = dict(input_args)
    sample_input['uuid'] = uuid
    sample_input['sample.tar'] = sample_location
    if sample_input['output_dir']:
        sample_input['output_dir'] = os.path.join(input_args['output_dir'], uuid)
    sample_input['cpu_count'] = multiprocessing.cpu_count()
    job_vars = (sample_input, ids)
    # Download or locate local file and place in the jobStore
    if sample_input['input']:
        ids['sample.tar'] = job.fileStore.writeGlobalFile(os.path.abspath(sample_location))
    elif sample_input['config_fastq']:
        ids['R1.fastq'] = job.fileStore.writeGlobalFile(urlparse(url1).path)
        ids['R2.fastq'] = job.fileStore.writeGlobalFile(urlparse(url2).path)
    else:
        if sample_input['ssec']:
            ids['sample.tar'] = job.addChildJobFn(download_encrypted_file, sample_input, 'sample.tar', disk='25G').rv()
        else:
            ids['sample.tar'] = job.addChildJobFn(download_from_url, sample_input['sample.tar'], disk='25G').rv()
    job.addFollowOnJobFn(static_dag_launchpoint, job_vars)


def static_dag_launchpoint(job, job_vars):
    """
    Statically define jobs in the pipeline

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    input_args, ids = job_vars
    if input_args['config_fastq']:
        cores = input_args['cpu_count']
        a = job.wrapJobFn(mapsplice, job_vars, cores=cores, disk='130G').encapsulate()
    else:
        a = job.wrapJobFn(merge_fastqs, job_vars, disk='70 G').encapsulate()
    b = job.wrapJobFn(consolidate_output, job_vars, a.rv())
    # Take advantage of "encapsulate" to simplify pipeline wiring
    job.addChild(a)
    a.addChild(b)


def merge_fastqs(job, job_vars):
    """
    Unzips input sample and concats the Read1 and Read2 groups together.

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    cores = input_args['cpu_count']
    single_end_reads = input_args['single_end_reads']
    # I/O
    sample = return_input_paths(job, work_dir, ids, 'sample.tar')
    # Untar File
    # subprocess.check_call(['unzip', sample, '-d', work_dir])
    subprocess.check_call(['tar', '-xvf', sample, '-C', work_dir])
    # Remove large files before creating concat versions.
    os.remove(os.path.join(work_dir, 'sample.tar'))
    # Zcat files in parallel
    if single_end_reads:
        files = sorted(glob.glob(os.path.join(work_dir, '*')))
        with open(os.path.join(work_dir, 'R1.fastq'), 'w') as f1:
            subprocess.check_call(['zcat'] + files, stdout=f1)
        # FileStore
        ids['R1.fastq'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1.fastq'))
    else:
        r1_files = sorted(glob.glob(os.path.join(work_dir, '*R1*')))
        r2_files = sorted(glob.glob(os.path.join(work_dir, '*R2*')))
        with open(os.path.join(work_dir, 'R1.fastq'), 'w') as f1:
            p1 = subprocess.Popen(['zcat'] + r1_files, stdout=f1)
        with open(os.path.join(work_dir, 'R2.fastq'), 'w') as f2:
            p2 = subprocess.Popen(['zcat'] + r2_files, stdout=f2)
        p1.wait()
        p2.wait()
        # FileStore
        ids['R1.fastq'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1.fastq'))
        ids['R2.fastq'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R2.fastq'))
    job.fileStore.deleteGlobalFile(ids['sample.tar'])
    # Spawn child job
    return job.addChildJobFn(mapsplice, job_vars, cores=cores, disk='130 G').rv()


def mapsplice(job, job_vars):
    """
    Maps RNA-Seq reads to a reference genome.

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    # Unpack variables
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    cores = input_args['cpu_count']
    sudo = input_args['sudo']
    single_end_reads = input_args['single_end_reads']
    files_to_delete = ['R1.fastq']
    # I/O
    return_input_paths(job, work_dir, ids, 'ebwt.zip', 'chromosomes.zip')
    if single_end_reads:
        return_input_paths(job, work_dir, ids, 'R1.fastq')
    else:
        return_input_paths(job, work_dir, ids, 'R1.fastq', 'R2.fastq')
        files_to_delete.extend(['R2.fastq'])
    for fname in ['chromosomes.zip', 'ebwt.zip']:
        subprocess.check_call(['unzip', '-o', os.path.join(work_dir, fname), '-d', work_dir])
    # Command and call
    parameters = ['-p', str(cores),
                  '-s', '25',
                  '--bam',
                  '--min-map-len', '50',
                  '-x', '/data/ebwt',
                  '-c', '/data/chromosomes',
                  '-1', '/data/R1.fastq',
                  '-o', '/data']
    if not single_end_reads:
        parameters.extend(['-2', '/data/R2.fastq'])
    docker_call(tool='quay.io/ucsc_cgl/mapsplice:2.0.1.9--2296da365ead6b12ded9d9c7b7798fbc927cd66b',
                tool_parameters=parameters, work_dir=work_dir, sudo=sudo)
    # Write to FileStore
    for fname in ['alignments.bam', 'stats.txt']:
        ids[fname] = job.fileStore.writeGlobalFile(os.path.join(work_dir, fname))
    for fname in files_to_delete:
        job.fileStore.deleteGlobalFile(ids[fname])
    # Run child job
    # map_id = job.addChildJobFn(mapping_stats, job_vars).rv()
    if input_args['upload_bam_to_s3'] and input_args['s3_dir']:
        job.addChildJobFn(upload_bam_to_s3, job_vars)
    output_ids = job.addChildJobFn(add_read_groups, job_vars, disk='30 G').rv()
    return output_ids


def mapping_stats(job, job_vars):
    """
    This function is not currently in use.

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    # Unpack variables
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    # I/O
    return_input_paths(job, work_dir, ids, 'stats.txt')
    uuid = input_args['uuid']
    # Command
    docker_call(tool='jvivian/mapping_stats', tool_parameters=[uuid], work_dir=work_dir, sudo=sudo)
    # Zip output files and store
    output_files = ['{}_stats2.txt'.format(uuid), '{}_stats_all.txt'.format(uuid), '{}_mapping.tab'.format(uuid)]
    tarball_files(work_dir, tar_name='map.tar.gz', files=output_files)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'map.tar.gz'))


def add_read_groups(job, job_vars):
    """
    This function adds read groups to the headers

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    # I/O
    alignments = return_input_paths(job, work_dir, ids, 'alignments.bam')
    output = os.path.join(work_dir, 'rg_alignments.bam')
    # Command and callg
    parameter = ['AddOrReplaceReadGroups',
                 'INPUT={}'.format(docker_path(alignments)),
                 'OUTPUT={}'.format(docker_path(output)),
                 'RGSM={}'.format(input_args['uuid']),
                 'RGID={}'.format(input_args['uuid']),
                 'RGLB=TruSeq',
                 'RGPL=illumina',
                 'RGPU=barcode',
                 'VALIDATION_STRINGENCY=SILENT']
    docker_call(tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                tool_parameters=parameter, work_dir=work_dir, sudo=sudo)
    # Write to FileStore
    ids['rg_alignments.bam'] = job.fileStore.writeGlobalFile(output)
    # Run child job
    return job.addChildJobFn(bamsort_and_index, job_vars, disk='30 G').rv()


def bamsort_and_index(job, job_vars):
    """
    Sorts bam file and produces index file

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    # Unpack variables
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    # I/O
    rg_alignments = return_input_paths(job, work_dir, ids, 'rg_alignments.bam')
    output = os.path.join(work_dir, 'sorted.bam')
    # Command -- second argument is "Output Prefix"
    cmd1 = ['sort', docker_path(rg_alignments), docker_path('sorted')]
    cmd2 = ['index', docker_path(output)]
    docker_call(tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                tool_parameters=cmd1, work_dir=work_dir, sudo=sudo)
    docker_call(tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                tool_parameters=cmd2, work_dir=work_dir, sudo=sudo)
    # Write to FileStore
    ids['sorted.bam'] = job.fileStore.writeGlobalFile(output)
    ids['sorted.bam.bai'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'sorted.bam.bai'))
    # Run child job
    output_ids = job.addChildJobFn(sort_bam_by_reference, job_vars, disk='50 G').rv()
    rseq_id = job.addChildJobFn(rseq_qc, job_vars, disk='20 G').rv()
    return rseq_id, output_ids


def rseq_qc(job, job_vars):
    """
    QC module: contains QC metrics and information about the BAM post alignment

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    uuid = input_args['uuid']
    sudo = input_args['sudo']
    # I/O
    return_input_paths(job, work_dir, ids, 'sorted.bam', 'sorted.bam.bai')
    # Command
    docker_call(tool='jvivian/qc', tool_parameters=['/opt/cgl-docker-lib/RseqQC_v2.sh', '/data/sorted.bam', uuid],
                work_dir=work_dir, sudo=sudo)
    # Write to FileStore
    output_files = [f for f in glob.glob(os.path.join(work_dir, '*')) if 'sorted.bam' not in f]
    tarball_files(work_dir, tar_name='qc.tar.gz', uuid=None, files=output_files)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'qc.tar.gz'))


def sort_bam_by_reference(job, job_vars):
    """
    Sorts the bam by reference

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    # Unpack variables
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    # I/O
    sorted_bam, sorted_bai = return_input_paths(job, work_dir, ids, 'sorted.bam', 'sorted.bam.bai')
    output = os.path.join(work_dir, 'sort_by_ref.bam')
    # Call: Samtools
    ref_seqs = []
    handle = subprocess.Popen(["samtools", "view", "-H", sorted_bam], stdout=subprocess.PIPE).stdout
    for line in handle:
        if line.startswith("@SQ"):
            tmp = line.split("\t")
            chrom = tmp[1].split(":")[1]
            ref_seqs.append(chrom)
    handle.close()
    # Iterate through chromosomes to create mini-bams
    for chrom in ref_seqs:
        # job.addChildJobFn(sbbr_child, chrom, os.path.join(work_dir, chrom), sorted_bam)
        cmd_view = ["samtools", "view", "-b", sorted_bam, chrom]
        cmd_sort = ["samtools", "sort", "-m", "3000000000", "-n", "-", os.path.join(work_dir, chrom)]
        p1 = subprocess.Popen(cmd_view, stdout=subprocess.PIPE)
        subprocess.check_call(cmd_sort, stdin=p1.stdout)
    sorted_files = [os.path.join(work_dir, chrom) + '.bam' for chrom in ref_seqs]
    cmd = ["samtools", "cat", "-o", output] + sorted_files
    subprocess.check_call(cmd)
    # Write to FileStore
    ids['sort_by_ref.bam'] = job.fileStore.writeGlobalFile(output)
    rsem_id = job.addChildJobFn(transcriptome, job_vars, disk='30 G', memory='30 G').rv()
    exon_id = job.addChildJobFn(exon_count, job_vars, disk='30 G').rv()
    return exon_id, rsem_id


def exon_count(job, job_vars):
    """
    Produces exon counts

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    uuid = input_args['uuid']
    sudo = input_args['sudo']
    # I/O
    sort_by_ref, normalize_pl, composite_bed = return_input_paths(job, work_dir, ids, 'sort_by_ref.bam',
                                                                  'normalize.pl', 'composite_exons.bed')
    # Command
    tool = 'jvivian/bedtools'
    cmd_1 = ['coverage',
             '-split',
             '-abam', docker_path(sort_by_ref),
             '-b', docker_path(composite_bed)]
    cmd_2 = ['perl',
             os.path.join(work_dir, 'normalize.pl'),
             sort_by_ref,
             composite_bed]

    popen_docker = ['docker', 'run', '-v', '{}:/data'.format(work_dir), tool]
    if sudo:
        popen_docker = ['sudo'] + popen_docker
    p = subprocess.Popen(popen_docker + cmd_1, stdout=subprocess.PIPE)
    with open(os.path.join(work_dir, 'exon_quant'), 'w') as f:
        subprocess.check_call(cmd_2, stdin=p.stdout, stdout=f)

    p1 = subprocess.Popen(['cat', os.path.join(work_dir, 'exon_quant')], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['tr', '":"', '"\t"'], stdin=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.Popen(['tr', '"-"', '"\t"'], stdin=p2.stdout, stdout=subprocess.PIPE)
    with open(os.path.join(work_dir, 'exon_quant.bed'), 'w') as f:
        subprocess.check_call(['cut', '-f1-4'], stdin=p3.stdout, stdout=f)
    # Create zip, upload to fileStore, and move to output_dir as a backup
    output_files = ['exon_quant.bed', 'exon_quant']
    tarball_files(work_dir, tar_name='exon.tar.gz', uuid=uuid, files=output_files)
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'exon.tar.gz'))


def transcriptome(job, job_vars):
    """
    Creates a bam of just the transcriptome

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    sudo = input_args['sudo']
    # I/O
    sort_by_ref, bed, hg19_fa = return_input_paths(job, work_dir, ids, 'sort_by_ref.bam',
                                                   'unc.bed', 'hg19.transcripts.fa')
    output = os.path.join(work_dir, 'transcriptome.bam')
    # Command
    parameters = ['sam-xlate',
                  '--bed', docker_path(bed),
                  '--in', docker_path(sort_by_ref),
                  '--order', docker_path(hg19_fa),
                  '--out', docker_path(output),
                  '--xgtag',
                  '--reverse']
    docker_call(tool='quay.io/ucsc_cgl/ubu:1.0--ce6807e937a7a29138f56ea1b8fc077528ee8180',
                tool_parameters=parameters, work_dir=work_dir, java_opts='-Xmx30g', sudo=sudo)
    # Write to FileStore
    ids['transcriptome.bam'] = job.fileStore.writeGlobalFile(output)
    # Run child job
    return job.addChildJobFn(filter_bam, job_vars, memory='30G', disk='30G').rv()


def filter_bam(job, job_vars):
    """
    Performs filtering on the transcriptome bam

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    cores = input_args['cpu_count']
    sudo = input_args['sudo']
    # I/O
    transcriptome_bam = return_input_paths(job, work_dir, ids, 'transcriptome.bam')
    output = os.path.join(work_dir, 'filtered.bam')
    # Command
    parameters = ['sam-filter',
                  '--strip-indels',
                  '--max-insert', '1000',
                  '--mapq', '1',
                  '--in', docker_path(transcriptome_bam),
                  '--out', docker_path(output)]
    docker_call(tool='quay.io/ucsc_cgl/ubu:1.0--ce6807e937a7a29138f56ea1b8fc077528ee8180',
                tool_parameters=parameters, work_dir=os.path.dirname(output), java_opts='-Xmx30g', sudo=sudo)
    # Write to FileStore
    ids['filtered.bam'] = job.fileStore.writeGlobalFile(output)
    # Run child job
    return job.addChildJobFn(rsem, job_vars, cores=cores, disk='30 G').rv()


def rsem(job, job_vars):
    """
    Runs RSEM to produce counts

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    cpus = input_args['cpu_count']
    sudo = input_args['sudo']
    single_end_reads = input_args['single_end_reads']
    # I/O
    filtered_bam, rsem_ref = return_input_paths(job, work_dir, ids, 'filtered.bam', 'rsem_ref.zip')
    subprocess.check_call(['unzip', '-o', os.path.join(work_dir, 'rsem_ref.zip'), '-d', work_dir])
    output_prefix = 'rsem'
    # Make tool call to Docker
    parameters = ['--quiet',
                  '--no-qualities',
                  '-p', str(cpus),
                  '--forward-prob', '0.5',
                  '--seed-length', '25',
                  '--fragment-length-mean', '-1.0',
                  '--bam', docker_path(filtered_bam)]
    if not single_end_reads:
        parameters.extend(['--paired-end'])
    parameters.extend(['/data/rsem_ref/hg19_M_rCRS_ref', output_prefix])

    docker_call(tool='quay.io/ucsc_cgl/rsem:1.1.13--be66304ff7fcd6fb0babcd1884ed289eabedc655',
                tool_parameters=parameters, work_dir=work_dir, sudo=sudo)
    os.rename(os.path.join(work_dir, output_prefix + '.genes.results'), os.path.join(work_dir, 'rsem_gene.tab'))
    os.rename(os.path.join(work_dir, output_prefix + '.isoforms.results'), os.path.join(work_dir, 'rsem_isoform.tab'))
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
    return_input_paths(job, work_dir, ids, 'rsem_gene.tab', 'rsem_isoform.tab')
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
    rseq_id, exon_id, rsem_id = flatten(output_ids)
    # Retrieve output file paths to consolidate
    # map_tar = job.fileStore.readGlobalFile(map_id, os.path.join(work_dir, 'map.tar.gz'))
    qc_tar = job.fileStore.readGlobalFile(rseq_id, os.path.join(work_dir, 'qc.tar.gz'))
    exon_tar = job.fileStore.readGlobalFile(exon_id, os.path.join(work_dir, 'exon.tar.gz'))
    rsem_tar = job.fileStore.readGlobalFile(rsem_id, os.path.join(work_dir, 'rsem.tar.gz'))
    # I/O
    out_tar = os.path.join(work_dir, uuid + '.tar.gz')
    # Consolidate separate tarballs
    with tarfile.open(os.path.join(work_dir, out_tar), 'w:gz') as f_out:
        for tar in [rsem_tar, exon_tar, qc_tar]:
            with tarfile.open(tar, 'r') as f_in:
                for tarinfo in f_in:
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        if tar == qc_tar:
                            tarinfo.name = os.path.join(uuid, 'rseq_qc', os.path.basename(tarinfo.name))
                        else:
                            tarinfo.name = os.path.join(uuid, os.path.basename(tarinfo.name))
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
        job.addChildJobFn(upload_output_to_s3, job_vars)


def upload_output_to_s3(job, job_vars):
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
    uuid_tar = return_input_paths(job, work_dir, ids, 'uuid.tar.gz')
    # Upload to S3 via boto
    conn = boto.connect_s3()
    bucket = conn.get_bucket(bucket_name)
    k = Key(bucket)
    k.key = os.path.join(bucket_dir, uuid + '.tar.gz')
    k.set_contents_from_filename(uuid_tar)


def upload_bam_to_s3(job, job_vars):
    """
    Upload bam to S3. Requires S3AM and a ~/.boto config file.
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    uuid = input_args['uuid']
    # I/O
    job.fileStore.readGlobalFile(ids['alignments.bam'], os.path.join(work_dir, 'alignments.bam'))
    bam_path = os.path.join(work_dir, 'alignments.bam')
    sample_name = uuid + '.bam'
    # Parse s3_dir to get bucket and s3 path
    s3_dir = input_args['s3_dir']
    bucket_name = s3_dir.split('/')[0]
    bucket_dir = os.path.join('/'.join(s3_dir.split('/')[1:]), 'bam_files')
    # Upload to S3 via S3AM
    s3am_command = ['s3am',
                    'upload',
                    'file://{}'.format(bam_path),
                    os.path.join('s3://', bucket_name, bucket_dir, sample_name)]
    subprocess.check_call(s3am_command)


def main():
    """
    This is a Toil pipeline for the UNC best practice RNA-Seq analysis.
    RNA-seq fastqs are combined, aligned, sorted, filtered, and quantified.

    Please read the README.md located in the same directory.
    """
    # Define Parser object and add to toil
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    # Store inputs from argparse
    inputs = {'config': args.config,
              'config_fastq': args.config_fastq,
              'input': args.input,
              'unc.bed': args.unc,
              'hg19.transcripts.fa': args.fasta,
              'composite_exons.bed': args.composite_exons,
              'normalize.pl': args.normalize,
              'output_dir': args.output_dir,
              'rsem_ref.zip': args.rsem_ref,
              'chromosomes.zip': args.chromosomes,
              'ebwt.zip': args.ebwt,
              'ssec': args.ssec,
              's3_dir': args.s3_dir,
              'sudo': args.sudo,
              'single_end_reads': args.single_end_reads,
              'upload_bam_to_s3': args.upload_bam_to_s3,
              'uuid': None,
              'sample.tar': None,
              'cpu_count': None}

    # Launch jobs
    Job.Runner.startToil(Job.wrapJobFn(download_shared_files, inputs), args)


if __name__ == "__main__":
    main()
