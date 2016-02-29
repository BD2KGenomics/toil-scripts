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
1 = Unpack/Merge fastqs
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
Docker:     wget -qO- https://get.docker.com/ | sh
Toil:       pip install toil
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
from subprocess import PIPE
import shutil
import tarfile
import multiprocessing
from urlparse import urlparse
from toil.job import Job
import logging

from toil_scripts import download_from_s3_url

logging.basicConfig(level=logging.INFO)


def build_parser():
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-c', '--config', default=None,
                       help='Path to CSV. One sample per line with the format: uuid,url (see example config). '
                            'URLs follow the format of: https://www.address.com/file.tar or file:///full/path/file.tar.'
                            ' Samples must be tarfiles that contain fastq files. \n\nAlternatively, if your data is'
                            'not tarred, you can create a CSV of one sample per line: uuid,url1,url2 .  Where url1'
                            'and url2 refer to urls to the pair of fastq files.')
    group.add_argument('-d', '--dir', default=None,
                       help='Path to directory of samples. Samples must be tarfiles that contain fastq files. '
                            'The UUID for the sample will be derived from the file basename.')
    group.add_argument('-s', '--sample_urls', nargs='+', default=None, type=str,
                       help='Sample URLs (any number). Samples must be tarfiles that contain fastq files. '
                            'URLs follow the format of: https://www.address.com/file.tar or file:///full/path/file.tar.'
                            'The UUID for the sample will be derived from the file basename.')
    group.add_argument('-g', '--genetorrent', default=None,
                       help='Path to a file with one analysis ID per line for data hosted on CGHub.')
    parser.add_argument('-k', '--genetorrent_key', default=None,
                        help='Path to a CGHub key that has access to the TCGA data being requested. An exception will'
                             'be thrown if "-g" is set but not this argument.')
    parser.add_argument('--starIndex', help='URL to download STAR Index built from HG38/gencodev23 annotation.',
                        default='https://s3-us-west-2.amazonaws.com/'
                                'cgl-pipeline-inputs/rnaseq_cgl/starIndex_hg38_no_alt.tar.gz')
    parser.add_argument('--kallistoIndex', help='URL to download Kallisto Index built from HG38 transcriptome.',
                        default='https://s3-us-west-2.amazonaws.com/'
                                'cgl-pipeline-inputs/rnaseq_cgl/kallisto_hg38.idx')
    parser.add_argument('--rsemRef', help='URL to download RSEM reference built from HG38/gencodev23.',
                        default='https://s3-us-west-2.amazonaws.com/'
                                'cgl-pipeline-inputs/rnaseq_cgl/rsem_ref_hg38_no_alt.tar.gz')
    parser.add_argument('--wiggle', default=None, action='store_true', help='Uploads a wiggle from STAR to S3.')
    parser.add_argument('--save_bam', default=None, action='store_true', help='Uploads aligned BAM to S3.')
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
    # Grab master key
    with open(key_path, 'r') as f:
        key = f.read()
    if len(key) != 32:
        raise RuntimeError('Invalid Key! Must be 32 bytes: {}'.format(key))

    key = generate_unique_key(key_path, url)
    # Create necessary headers for SSE-C encryption and download
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
    if url.startswith('s3:'):
        download_from_s3_url(file_path, url)
    else:
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
        raise RuntimeError('docker command returned a non-zero exit status: {}'.format(tool_parameters))
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
    Parses config file. Returns list of samples: [ [uuid1, url1], [uuid2, url2], ... ]
    """
    samples = []
    with open(path_to_config, 'r') as f:
        for line in f.readlines():
            if not line.isspace():
                sample = line.strip().split(',')
                assert len(sample) == 2 or len(sample) == 3, 'Error: Config file is inappropriately formatted.'
                samples.append(sample)
    return samples


def parse_genetorrent(path_to_config):
    """
    Parses genetorrent config file.  Returns list of samples: [ [id1, id1 ], [id2, id2], ... ]
    Returns duplicate of ids to follow UUID/URL standard.
    """
    samples = []
    with open(path_to_config, 'r') as f:
        for line in f.readlines():
            if not line.isspace():
                samples.append([line.strip(), line.strip()])
    return samples


def download_from_genetorrent(job, input_args, analysis_id):
    """
    Downloads a sample from CGHub via GeneTorrent.

    input_args: dict        Dictionary of input arguments
    analysis_id: str        An analysis ID for a sample in CGHub
    """
    work_dir = job.fileStore.getLocalTempDir()
    folder_path = os.path.join(work_dir, os.path.basename(analysis_id))
    sudo = input_args['sudo']
    shutil.copy(input_args['genetorrent_key'], os.path.join(work_dir, 'cghub.key'))
    parameters = ['-vv', '-c', 'cghub.key', '-d', analysis_id]
    docker_call(tool='quay.io/ucsc_cgl/genetorrent:3.8.7--9911761265b6f08bc3ef09f53af05f56848d805b',
                work_dir=work_dir, tool_parameters=parameters, sudo=sudo)
    sample = glob.glob(os.path.join(folder_path, '*tar*'))
    assert len(sample) == 1, 'More than one sample tar in CGHub download: {}'.format(analysis_id)
    return job.fileStore.writeGlobalFile(sample[0])


def partitions(l, partition_size):
    """
    >>> list(partitions([], 10))
    []
    >>> list(partitions([1,2,3,4,5], 1))
    [[1], [2], [3], [4], [5]]
    >>> list(partitions([1,2,3,4,5], 2))
    [[1, 2], [3, 4], [5]]
    >>> list(partitions([1,2,3,4,5], 5))
    [[1, 2, 3, 4, 5]]
    """
    for i in xrange(0, len(l), partition_size):
        yield l[i:i+partition_size]


# Job Functions
def parse_input_samples(job, input_args):
    """
    Parses the various input formats then launches the job batcher

    shared_ids: dict        Dictionary of fileStore IDs
    input_args: dict        Dictionary of input arguments (from main())
    """
    config = input_args['config']
    sample_dir = input_args['dir']
    sample_urls = input_args['sample_urls']
    genetorrent = input_args['genetorrent']
    samples = None
    if config:
        samples = parse_config(config)
    elif sample_dir:
        files = os.listdir(sample_dir)
        samples = [[os.path.splitext(os.path.basename(f))[0],
                    'file://' + os.path.join(sample_dir, os.path.basename(f))] for f in files]
    elif sample_urls:
        samples = [[os.path.splitext(os.path.basename(f))[0], f] for f in sample_urls]
    elif genetorrent:
        samples = parse_genetorrent(genetorrent)
    # Pass to batcher to spawn tree of jobs
    job.addChildJobFn(batcher, input_args, samples)


def batcher(job, input_args, samples):
    """
    Spawns a tree of jobs to avoid overloading the number of jobs spawned by a single parent.

    input_args: dict             Dictionary of input arguments
    samples: list/iterable       target list/iterable that spawns one job per item
    """
    num_partitions = 10
    partition_size = len(samples) / num_partitions
    if partition_size > 1:
        for partition in partitions(samples, partition_size):
            job.addChildJobFn(batcher, input_args, partition)
    else:
        for sample in samples:
            job.addChildJobFn(download_sample, input_args, sample)


def download_sample(job, input_args, sample):
    """
    Defines variables unique to a sample that are used in the rest of the pipelines

    ids: dict           Dictionary of fileStore IDS
    input_args: dict    Dictionary of input arguments
    sample: tuple       Contains uuid and sample_url
    """
    sample_input = dict(input_args)
    uuid, url = None, None
    if len(sample) == 2:
        uuid, url = sample
        sample_input['sample.tar'] = url
    if len(sample) == 3:
        uuid = sample[0]
        url = sample[1:]
        sample_input['R1.fastq'] = url[0]
        sample_input['R2.fastq'] = url[1]
    assert uuid and url, 'Issue with sample configuration retrieval: {}'.format(sample)
    # Update values unique to sample
    sample_input['uuid'] = uuid
    if sample_input['output_dir']:
        sample_input['output_dir'] = os.path.join(input_args['output_dir'], uuid)
    sample_input['cpu_count'] = multiprocessing.cpu_count()
    ids = {}
    job_vars = (sample_input, ids)
    # Download or locate local file and place in the jobStore
    if sample_input['genetorrent']:
        ids['sample.tar'] = job.addChildJobFn(download_from_genetorrent, input_args, url, disk='40G').rv()
    elif type(url) is list and len(url) == 2:
        if urlparse(url[0]) == 'file':
            ids['R1.fastq'] = job.fileStore.writeGlobalFile(urlparse(url[0]).path)
            ids['R2.fastq'] = job.fileStore.writeGlobalFile(urlparse(url[1]).path)
        else:
            if sample_input['ssec']:
                ids['R1.fastq'] = job.addChildJobFn(download_encrypted_file, sample_input, 'R1.fastq', disk='40G').rv()
                ids['R2.fastq'] = job.addChildJobFn(download_encrypted_file, sample_input, 'R1.fastq', disk='40G').rv()
            else:
                ids['R1.fastq'] = job.addChildJobFn(download_from_url, sample_input['R1.fastq'], disk='40G').rv()
                ids['R2.fastq'] = job.addChildJobFn(download_from_url, sample_input['R1.fastq'], disk='40G').rv()
    elif urlparse(url).scheme == 'file':
        ids['sample.tar'] = job.fileStore.writeGlobalFile(urlparse(url).path)
    else:
        if sample_input['ssec']:
            ids['sample.tar'] = job.addChildJobFn(download_encrypted_file, sample_input, 'sample.tar', disk='40G').rv()
        else:
            ids['sample.tar'] = job.addChildJobFn(download_from_url, sample_input['sample.tar'], disk='40G').rv()
    job.addFollowOnJobFn(static_dag_launchpoint, job_vars)


def static_dag_launchpoint(job, job_vars):
    """
    Statically define pipeline

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    input_args, ids = job_vars
    if 'sample.tar' in ids:
        a = job.wrapJobFn(process_sample_tar, job_vars, disk='70G').encapsulate()
    else:
        a = job.wrapJobFn(cutadapt, job_vars, disk='100G').encapsulate()
    b = job.wrapJobFn(consolidate_output, job_vars, a.rv(), disk='2G')
    # Take advantage of "encapsulate" to simplify pipeline wiring
    job.addChild(a)
    a.addChild(b)


def process_sample_tar(job, job_vars):
    """
    Converts sample.tar(.gz) into two fastq files.
    Due to edge conditions... BEWARE: HERE BE DRAGONS

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    # Unpack variables
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    ids['R.fastq'] = None
    # I/O
    read_from_filestore(job, work_dir, ids, 'sample.tar')
    sample_tar = os.path.join(work_dir, 'sample.tar')
    # Untar File and concat
    p = subprocess.Popen(['tar', '-xvf', sample_tar, '-C', work_dir], stderr=PIPE, stdout=PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        # Handle error if tar archive is corrupt
        if 'EOF' in stderr:
            with open(os.path.join(work_dir, 'error.txt'), 'w') as f:
                f.write(stderr)
                f.write(stdout)
            ids['error.txt'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'error.txt'))
            if input_args['s3_dir']:
                job.addChildJobFn(upload_to_s3, job_vars)
        else:
            raise subprocess.CalledProcessError
    else:
        os.remove(os.path.join(work_dir, 'sample.tar'))
        # Grab files from tarball
        fastqs = []
        for root, subdir, files in os.walk(work_dir):
            fastqs.extend([os.path.join(root, x) for x in files])
        # Check for read 1 and read 2 files
        r1 = sorted([x for x in fastqs if '_1' in x])
        r2 = sorted([x for x in fastqs if '_2' in x])
        if not r1 or not r2:
            # Check if using a different standard
            r1 = sorted([x for x in fastqs if 'R1' in x])
            r2 = sorted([x for x in fastqs if 'R2' in x])
        # Prune file name matches from each list
        if len(r1) > len(r2):
            r1 = [x for x in r1 if x not in r2]
        elif len(r2) > len(r1):
            r2 = [x for x in r2 if x not in r1]
        if not r1 or not r2:
            # Sample is assumed to be single-ended
            if fastqs[0].endswith('.gz'):
                with open(os.path.join(work_dir, 'R.fastq'), 'w') as f:
                    subprocess.check_call(['zcat'] + fastqs, stdout=f)
            elif len(fastqs) > 1:
                with open(os.path.join(work_dir, 'R.fastq'), 'w') as f:
                    subprocess.check_call(['cat'] + fastqs, stdout=f)
            else:
                shutil.move(fastqs[0], os.path.join(work_dir, 'R.fastq'))
            ids['R.fastq'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R.fastq'))
        else:
            # Sample is assumed to be paired end
            if r1[0].endswith('.gz') and r2[0].endswith('.gz'):
                with open(os.path.join(work_dir, 'R1.fastq'), 'w') as f1:
                    p1 = subprocess.Popen(['zcat'] + r1, stdout=f1)
                with open(os.path.join(work_dir, 'R2.fastq'), 'w') as f2:
                    p2 = subprocess.Popen(['zcat'] + r2, stdout=f2)
                p1.wait()
                p2.wait()
            elif len(r1) > 1 and len(r2) > 1:
                with open(os.path.join(work_dir, 'R1.fastq'), 'w') as f1:
                    p1 = subprocess.Popen(['cat'] + r1, stdout=f1)
                with open(os.path.join(work_dir, 'R2.fastq'), 'w') as f2:
                    p2 = subprocess.Popen(['cat'] + r2, stdout=f2)
                p1.wait()
                p2.wait()
            else:
                shutil.move(r1[0], os.path.join(work_dir, 'R1.fastq'))
                shutil.move(r2[0], os.path.join(work_dir, 'R2.fastq'))
            ids['R1.fastq'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1.fastq'))
            ids['R2.fastq'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R2.fastq'))
        job.fileStore.deleteGlobalFile(ids['sample.tar'])
        # Start cutadapt step
        return job.addChildJobFn(cutadapt, job_vars, disk='125G').rv()


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
    parameters = ['-a', input_args['fwd_3pr_adapter'],
                  '-m', '35']
    if ids['R.fastq']:
        input_args['single_end'] = True
        read_from_filestore(job, work_dir, ids, 'R.fastq')
        parameters.extend(['-o', '/data/R_cutadapt.fastq', '/data/R.fastq'])
    else:
        read_from_filestore(job, work_dir, ids, 'R1.fastq', 'R2.fastq')
        parameters.extend(['-A', input_args['rev_3pr_adapter'],
                           '-o', '/data/R1_cutadapt.fastq',
                           '-p', '/data/R2_cutadapt.fastq',
                           '/data/R1.fastq', '/data/R2.fastq'])
    # Call: CutAdapt
    base_docker_call = 'docker run --log-driver=none --rm -v {}:/data'.format(work_dir).split()
    if sudo:
        base_docker_call = ['sudo'] + base_docker_call
    tool = 'quay.io/ucsc_cgl/cutadapt:1.9--6bd44edd2b8f8f17e25c5a268fedaab65fa851d2'
    p = subprocess.Popen(base_docker_call + [tool] + parameters, stderr=PIPE, stdout=PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        if 'improperly paired' in stderr:
            input_args['improper_pair'] = True
            if ids['R.fastq']:
                shutil.move(os.path.join(work_dir, 'R.fastq'), os.path.join(work_dir, 'R_cutadapt.fastq'))
            else:
                shutil.move(os.path.join(work_dir, 'R1.fastq'), os.path.join(work_dir, 'R1_cutadapt.fastq'))
                shutil.move(os.path.join(work_dir, 'R2.fastq'), os.path.join(work_dir, 'R2_cutadapt.fastq'))
        else:
            logging.error('Stdout: {}\n\nStderr: {}'.format(stdout, stderr))
            raise subprocess.CalledProcessError(p.returncode, parameters, stderr)
    # Write to fileStore
    if ids['R.fastq']:
        ids['R_cutadapt.fastq'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R_cutadapt.fastq'))
    else:
        ids['R_cutadapt.fastq'] = None
        ids['R1_cutadapt.fastq'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R1_cutadapt.fastq'))
        ids['R2_cutadapt.fastq'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'R2_cutadapt.fastq'))
    # start STAR and Kallisto steps
    rsem_output = job.addChildJobFn(star, job_vars, cores=cores, disk='150G', memory='40G').rv()
    cores = min(cores, 16)
    kallisto_output = job.addChildJobFn(kallisto, job_vars, cores=cores, disk='100G').rv()
    return rsem_output, kallisto_output, input_args['improper_pair'], input_args['single_end']


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
    subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', input_args['kallisto_hg38.idx'], '-o',
                           os.path.join(work_dir, 'kallisto_hg38.idx')])
    # Retrieve files
    parameters = ['quant',
                  '-i', '/data/kallisto_hg38.idx',
                  '-t', str(cores),
                  '-o', '/data/',
                  '-b', '100']
    if ids['R_cutadapt.fastq']:
        # TODO: Look into fragment length and SD calculations
        read_from_filestore(job, work_dir, ids, 'R_cutadapt.fastq')
        parameters.extend(['--single', '-l', '200', '-s', '15', '/data/R_cutadapt.fastq'])
    else:
        read_from_filestore(job, work_dir, ids, 'R1_cutadapt.fastq', 'R2_cutadapt.fastq')
        parameters.extend(['/data/R1_cutadapt.fastq', '/data/R2_cutadapt.fastq'])
    # Call: Kallisto
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
    uuid = input_args['uuid']
    work_dir = job.fileStore.getLocalTempDir()
    # Parameters and input retrieval
    parameters = ['--runThreadN', str(cores),
                  '--genomeDir', '/data/starIndex',
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
    subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', input_args['starIndex.tar.gz'], '-o',
                           os.path.join(work_dir, 'starIndex.tar.gz')])
    if input_args['wiggle']:
        parameters.extend(['--outWigType', 'bedGraph',
                           '--outWigStrand', 'Unstranded',
                           '--outWigReferencesPrefix', 'chr'])
    if ids['R_cutadapt.fastq']:
        read_from_filestore(job, work_dir, ids, 'R_cutadapt.fastq')
        parameters.extend(['--readFilesIn', '/data/R_cutadapt.fastq'])
    else:
        read_from_filestore(job, work_dir, ids, 'R1_cutadapt.fastq', 'R2_cutadapt.fastq')
        parameters.extend(['--readFilesIn', '/data/R1_cutadapt.fastq', '/data/R2_cutadapt.fastq'])
    subprocess.check_call(['tar', '-xvf', os.path.join(work_dir, 'starIndex.tar.gz'), '-C', work_dir])
    # Call: STAR Map
    docker_call(tool='quay.io/ucsc_cgl/star:2.4.2a--bcbd5122b69ff6ac4ef61958e47bde94001cfe80',
                work_dir=work_dir, tool_parameters=parameters, sudo=sudo)
    # Write to fileStore
    ids['transcriptome.bam'] = job.fileStore.writeGlobalFile(os.path.join(work_dir,
                                                                          'rnaAligned.toTranscriptome.out.bam'))
    # Save Wiggle File
    if input_args['wiggle'] and input_args['s3_dir']:
        wiggles = [os.path.basename(x) for x in glob.glob(os.path.join(work_dir, '*.bg'))]
        # Rename extension
        for wiggle in wiggles:
            shutil.move(os.path.join(work_dir, wiggle),
                        os.path.join(work_dir, os.path.splitext(wiggle)[0] + '.bedGraph'))
        wiggles = [os.path.splitext(x)[0] + '.bedGraph' for x in wiggles]
        tarball_files(work_dir, 'wiggle.tar.gz', uuid=uuid, files=wiggles)
        ids['wiggle.tar.gz'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'wiggle.tar.gz'))
        job.addChildJobFn(upload_wiggle_to_s3, job_vars, disk='10G', cores=cores)
    if input_args['save_bam'] and input_args['s3_dir']:
        job.addChildJobFn(upload_bam_to_s3, job_vars, disk='20G', cores=cores)
    # RSEM doesn't tend to use more than 16 cores
    cores = min(cores, 16)
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
    read_from_filestore(job, work_dir, ids, 'transcriptome.bam')
    subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', input_args['rsem_ref_hg38.tar.gz'], '-o',
                           os.path.join(work_dir, 'rsem_ref_hg38.tar.gz')])
    subprocess.check_call(['tar', '-xvf', os.path.join(work_dir, 'rsem_ref_hg38.tar.gz'), '-C', work_dir])
    prefix = 'rsem'
    # Call: RSEM
    parameters = ['--quiet',
                  '--no-qualities',
                  '-p', str(cores),
                  '--forward-prob', '0.5',
                  '--seed-length', '25',
                  '--fragment-length-mean', '-1.0',
                  '--bam', '/data/transcriptome.bam',
                  '/data/rsem_ref_hg38/hg38',
                  prefix]
    if not ids['R_cutadapt.fastq']:
        parameters = ['--paired-end'] + parameters
    docker_call(tool='quay.io/ucsc_cgl/rsem:1.2.25--d4275175cc8df36967db460b06337a14f40d2f21',
                tool_parameters=parameters, work_dir=work_dir, sudo=sudo)
    os.rename(os.path.join(work_dir, prefix + '.genes.results'), os.path.join(work_dir, 'rsem_gene.tab'))
    os.rename(os.path.join(work_dir, prefix + '.isoforms.results'), os.path.join(work_dir, 'rsem_isoform.tab'))
    # Write to FileStore
    ids['rsem_gene.tab'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem_gene.tab'))
    ids['rsem_isoform.tab'] = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem_isoform.tab'))
    # Run child jobs
    return job.addChildJobFn(rsem_postprocess, job_vars).rv()


def rsem_postprocess(job, job_vars):
    """
    Parses RSEMs output to produce the separate .tab files (TPM, FPKM, counts) for both gene and isoform.
    These are two-column files: Genes and Quantifications.
    HUGO files are also provided that have been mapped from Gencode/ENSEMBLE names.

    job_vars: tuple     Tuple of dictionaries: input_args and ids
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    uuid = input_args['uuid']
    sudo = input_args['sudo']
    # I/O
    read_from_filestore(job, work_dir, ids, 'rsem_gene.tab', 'rsem_isoform.tab')
    # Convert RSEM files into individual .tab files.
    sample = input_args['uuid']
    docker_call(tool='jvivian/rsem_postprocess', tool_parameters=[sample], work_dir=work_dir, sudo=sudo)
    os.rename(os.path.join(work_dir, 'rsem_gene.tab'), os.path.join(work_dir, 'rsem_genes.results'))
    os.rename(os.path.join(work_dir, 'rsem_isoform.tab'), os.path.join(work_dir, 'rsem_isoforms.results'))
    output_files = ['rsem.genes.norm_counts.tab', 'rsem.genes.raw_counts.tab', 'rsem.genes.norm_fpkm.tab',
                    'rsem.genes.norm_tpm.tab', 'rsem.isoform.norm_counts.tab', 'rsem.isoform.raw_counts.tab',
                    'rsem.isoform.norm_fpkm.tab', 'rsem.isoform.norm_tpm.tab', 'rsem_genes.results',
                    'rsem_isoforms.results']
    # Perform HUGO gene / isoform name mapping
    genes = [x for x in output_files if 'gene' in x]
    isoforms = [x for x in output_files if 'isoform' in x]
    command = ['-g'] + genes + ['-i'] + isoforms
    docker_call(tool='jvivian/gencode_hugo_mapping', tool_parameters=command, work_dir=work_dir, sudo=sudo)
    hugo_files = [os.path.splitext(x)[0] + '.hugo' + os.path.splitext(x)[1] for x in output_files]
    # Create tarballs for outputs
    tarball_files(work_dir, tar_name='rsem.tar.gz', uuid=uuid, files=output_files)
    tarball_files(work_dir, tar_name='rsem_hugo.tar.gz', uuid=uuid, files=hugo_files)
    rsem_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem.tar.gz'))
    hugo_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'rsem_hugo.tar.gz'))
    return rsem_id, hugo_id


def consolidate_output(job, job_vars, output_ids_and_values):
    """
    Combine the contents of separate zipped outputs into one via streaming

    output_ids_and_values: tuple   Nested tuple of all the output fileStore IDs and dictionaries
    """
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve IDs
    rsem_id, hugo_id, kallisto_id, improper_pair, single_end = flatten(output_ids_and_values)
    uuid = input_args['uuid']
    cores = input_args['cpu_count']
    # Retrieve output file paths to consolidate
    rsem_tar = job.fileStore.readGlobalFile(rsem_id, os.path.join(work_dir, 'rsem.tar.gz'))
    kallisto_tar = job.fileStore.readGlobalFile(kallisto_id, os.path.join(work_dir, 'kallisto.tar.gz'))
    hugo_tar = job.fileStore.readGlobalFile(hugo_id, os.path.join(work_dir, 'rsem_hugo.tar.gz'))
    # I/O
    if improper_pair:
        uuid = 'IMPROPERLY_PAIRED.' + uuid
        input_args['uuid'] = uuid
    if single_end:
        uuid = 'SINGLE-END.' + uuid
        input_args['uuid'] = uuid
    out_tar = os.path.join(work_dir, uuid + '.tar.gz')
    # Consolidate separate tarballs into one as streams (avoids unnecessary untaring)
    with tarfile.open(os.path.join(work_dir, out_tar), 'w:gz') as f_out:
        for tar in [rsem_tar, kallisto_tar, hugo_tar]:
            with tarfile.open(tar, 'r') as f_in:
                for tarinfo in f_in:
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        if tar == rsem_tar:
                            tarinfo.name = os.path.join(uuid, 'RSEM', os.path.basename(tarinfo.name))
                        elif tar == hugo_tar:
                            tarinfo.name = os.path.join(uuid, 'RSEM', 'Hugo', os.path.basename(tarinfo.name))
                        else:
                            tarinfo.name = os.path.join(uuid, 'Kallisto', os.path.basename(tarinfo.name))
                        f_out.addfile(tarinfo, fileobj=f_in_file)
        if input_args['improper_pair']:
            with open(os.path.join(work_dir, 'WARNING.txt'), 'w') as f:
                f.write('cutadapt: error: Reads are improperly paired. Uneven number of reads in fastq pair.')
            f_out.add(os.path.join(work_dir, 'WARNING.txt'))
    # Move to output directory of selected
    if input_args['output_dir']:
        output_dir = input_args['output_dir']
        mkdir_p(output_dir)
        copy_to_output_dir(work_dir, output_dir, uuid=None, files=[uuid + '.tar.gz'])
    # Write output file to fileStore
    ids['uuid.tar.gz'] = job.fileStore.writeGlobalFile(out_tar)
    # If S3 bucket argument specified, upload to S3
    if input_args['s3_dir']:
        job.addChildJobFn(upload_to_s3, (input_args, ids), cores=cores)


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
    # Upload to S3 via boto
    conn = boto.connect_s3()
    bucket = conn.get_bucket(bucket_name)
    k = Key(bucket)
    if 'error.txt' in ids:
        read_from_filestore(job, work_dir, ids, 'error.txt')
        k.key = os.path.join(bucket_dir, uuid + '.ERROR')
        k.set_contents_from_filename(os.path.join(work_dir, 'error.txt'))
    else:
        read_from_filestore(job, work_dir, ids, 'uuid.tar.gz')
        uuid_tar = os.path.join(work_dir, 'uuid.tar.gz')
        if 'R.fastq' in ids:
            k.key = os.path.join(bucket_dir, uuid + 'single-end' + '.tar.gz')
        else:
            k.key = os.path.join(bucket_dir, uuid + '.tar.gz')
        k.set_contents_from_filename(uuid_tar)


def s3am_with_retry(cores, *args):
    retry_count = 3
    for i in xrange(retry_count):
        s3am_command = ['s3am', 'upload', '--resume', '--part-size=50M',
                        '--upload-slots={}'.format(cores),
                        '--download-slots={}'.format(cores)] + list(args)
        ret_code = subprocess.call(s3am_command)
        if ret_code == 0:
            return
        else:
            print 'S3AM failed with status code: {}'.format(ret_code)
    raise RuntimeError, 'S3AM failed to upload after {} retries.'.format(retry_count)


def upload_wiggle_to_s3(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    cores = input_args['cpu_count']
    uuid = input_args['uuid']
    # I/O
    read_from_filestore(job, work_dir, ids, 'wiggle.tar.gz')
    wiggle_tar = os.path.join(work_dir, 'wiggle.tar.gz')
    sample_name = uuid + '.wiggle.tar.gz'
    # Parse s3_dir to get bucket and s3 path
    s3_dir = input_args['s3_dir']
    bucket_name = s3_dir.split('/')[0]
    bucket_dir = os.path.join('/'.join(s3_dir.split('/')[1:]), 'wiggle_files')
    # Upload to S3 via S3AM
    s3am_with_retry(cores, 'file://{}'.format(wiggle_tar),
                    os.path.join('s3://', bucket_name, bucket_dir, sample_name))


def upload_bam_to_s3(job, job_vars):
    input_args, ids = job_vars
    work_dir = job.fileStore.getLocalTempDir()
    uuid = input_args['uuid']
    cores = input_args['cpu_count']
    key_path = input_args['ssec']
    # I/O
    read_from_filestore(job, work_dir, ids, 'transcriptome.bam')
    bam_path = os.path.join(work_dir, 'transcriptome.bam')
    sample_name = uuid + '.bam'
    # Parse s3_dir to get bucket and s3 path
    s3_dir = input_args['s3_dir']
    bucket_name = s3_dir.split('/')[0]
    bucket_dir = os.path.join('/'.join(s3_dir.split('/')[1:]), 'bam_files')
    base_url = 'https://s3-us-west-2.amazonaws.com/'
    url = os.path.join(base_url, bucket_name, bucket_dir, sample_name)
    # Generate keyfile for upload
    with open(os.path.join(work_dir, 'temp.key'), 'wb') as f_out:
        f_out.write(generate_unique_key(key_path, url))
    # Upload to S3 via S3AM
    s3am_with_retry(cores, '--sse-key-file', os.path.join(work_dir, 'temp.key'),
                    'file://{}'.format(bam_path),
                    os.path.join('s3://', bucket_name, bucket_dir, sample_name))


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
    # Sanity checks
    inputs = {'config': args.config,
              'dir': args.dir,
              'sample_urls': args.sample_urls,
              'genetorrent': args.genetorrent,
              'genetorrent_key': args.genetorrent_key,
              'starIndex.tar.gz': args.starIndex,
              'rsem_ref_hg38.tar.gz': args.rsemRef,
              'kallisto_hg38.idx': args.kallistoIndex,
              'fwd_3pr_adapter': args.fwd_3pr_adapter,
              'rev_3pr_adapter': args.rev_3pr_adapter,
              'output_dir': args.output_dir,
              'ssec': args.ssec,
              's3_dir': args.s3_dir,
              'sudo': args.sudo,
              'wiggle': args.wiggle,
              'save_bam': args.save_bam,
              'uuid': None,
              'sample.tar': None,
              'cpu_count': None,
              'improper_pair': None,
              'single_end': None}
    if args.ssec:
        assert os.path.isfile(args.ssec)
    if args.config:
        assert os.path.isfile(args.config)
    if args.dir:
        assert os.path.isdir(args.dir)
        assert os.listdir(args.dir) != []
    if args.sample_urls:
        assert args.sample_urls != []
    if args.genetorrent:
        assert os.path.isfile(args.genetorrent)
    if args.genetorrent_key:
        assert os.path.isfile(args.genetorrent_key)
    # Genetorrent key must be provided along with genetorrent option
    if args.genetorrent and not args.genetorrent_key:
        raise RuntimeError("Cannot supply -genetorrent without -genetorrent_key")
    # Program checks
    for program in ['curl', 'docker']:
        assert which(program), 'Program "{}" must be installed on every node.'.format(program)

    # Start Pipeline
    Job.Runner.startToil(Job.wrapJobFn(parse_input_samples, inputs), args)


if __name__ == '__main__':
    main()
