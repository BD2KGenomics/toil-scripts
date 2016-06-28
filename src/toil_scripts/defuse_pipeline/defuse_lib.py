import os
import tarfile
import subprocess
import time
import re
import gzip
from urlparse import urlparse
from collections import defaultdict

from toil_scripts.lib.programs import docker_call

import shutil
from pysam import Samfile

def untargz(input_targz_file, untar_to_dir):
    """
    This module accepts a tar.gz archive and untars it.

    RETURN VALUE: path to the untar-ed directory/file

    NOTE: this module expects the multiple files to be in a directory before
          being tar-ed.
    """
    assert tarfile.is_tarfile(input_targz_file), 'Not a tar file.'
    tarball = tarfile.open(input_targz_file)
    return_value = os.path.join(untar_to_dir, tarball.getmembers()[0].name)
    tarball.extractall(path=untar_to_dir)
    tarball.close()
    return return_value


def is_gzipfile(filename):
    """
    This function attempts to ascertain the gzip status of a file based on the "magic signatures" of
    the file. This was taken from the stack overflow
    http://stackoverflow.com/questions/13044562/python-mechanism-to-identify-compressed-file-type\
        -and-uncompress
    """
    assert os.path.exists(filename), 'Input {} does not '.format(filename) + \
                                     'point to a file.'
    with file(filename, 'rb') as in_f:
        start_of_file = in_f.read(3)
        if start_of_file == '\x1f\x8b\x08':
            # bam files are bgzipped and they share the magic sequence with gzip.  Pysam will error
            # if the input is gzip but not if it is a bam.
            try:
                _ = Samfile(filename)
            except ValueError:
                return True
            else:
                return False
        else:
            return False


def get_files_from_filestore(job, files, work_dir, cache=True, docker=False):
    """
    This is adapted from John Vivian's return_input_paths from the RNA-Seq pipeline.

    Returns the paths of files from the FileStore if they are not present.
    If docker=True, return the docker path for the file.
    If the file extension is tar.gz, then tar -zxvf it.

    files is a dict with:
        keys = the name of the file to be returned in toil space
        value = the input value for the file (can be toil temp file)
    work_dir is the location where the file should be stored
    cache indiciates whether caching should be used
    """
    for name in files.keys():
        outfile = job.fileStore.readGlobalFile(files[name], '/'.join([work_dir, name]), cache=cache)
        # If the file pointed to a tarball, extract it to WORK_DIR
        if tarfile.is_tarfile(outfile) and file_xext(outfile).startswith('.tar'):
            untar_name = os.path.basename(strip_xext(outfile))
            files[untar_name] = untargz(outfile, work_dir)
            files.pop(name)
            name = os.path.basename(untar_name)
        # If the file is gzipped but NOT a tarfile, gunzip it to work_dir. However, the file is
        # already named x.gz so we need to write to a temporary file x.gz_temp then do a move
        # operation to overwrite x.gz.
        elif is_gzipfile(outfile) and file_xext(outfile) == '.gz':
            ungz_name = strip_xext(outfile)
            with gzip.open(outfile, 'rb') as gz_in, open(ungz_name, 'w') as ungz_out:
                shutil.copyfileobj(gz_in, ungz_out)
            files[os.path.basename(ungz_name)] = outfile
            files.pop(name)
            name = os.path.basename(ungz_name)
        else:
            files[name] = outfile
        # If the files will be sent to docker, we will mount work_dir to the container as /data and
        # we want the /data prefixed path to the file
        if docker:
            files[name] = docker_path(files[name])
    return files


def get_file_from_s3(job, s3_url, encryption_key=None, per_file_encryption=True,
                     write_to_jobstore=True):
    """
    Downloads a supplied URL that points to an unencrypted, unprotected file on Amazon S3. The file
    is downloaded and a subsequently written to the jobstore and the return value is a the path to
    the file in the jobstore.

    :param str s3_url: URL for the file (can be s3 or https)
    :param str encryption_key: Path to the master key
    :param bool per_file_encryption: If encrypted, was the file encrypted using the per-file method?
    :param bool write_to_jobstore: Should the file be written to the job store?
    """
    work_dir = job.fileStore.getLocalTempDir()

    parsed_url = urlparse(s3_url)
    if parsed_url.scheme == 'https':
        download_url = 'S3:/' + parsed_url.path  # path contains the second /
    elif parsed_url.scheme == 's3':
        download_url = s3_url
    else:
        raise RuntimeError('Unexpected url scheme: %s' % s3_url)

    filename = '/'.join([work_dir, os.path.basename(s3_url)])
    # This is common to encrypted and unencrypted downloads
    download_call = ['s3am', 'download', '--download-exists', 'resume']
    # If an encryption key was provided, use it.
    if encryption_key:
        download_call.extend(['--sse-key-file', encryption_key])
        if per_file_encryption:
            download_call.append('--sse-key-is-master')
    # This is also common to both types of downloads
    download_call.extend([download_url, filename])
    attempt = 0
    while True:
        try:
            with open(work_dir + '/stderr', 'w') as stderr_file:
                subprocess.check_call(download_call, stderr=stderr_file)
        except subprocess.CalledProcessError:
            # The last line of the stderr will have the error
            with open(stderr_file.name) as stderr_file:
                for line in stderr_file:
                    line = line.strip()
                    if line:
                        exception = line
            if exception.startswith('boto'):
                exception = exception.split(': ')
                if exception[-1].startswith('403'):
                    raise RuntimeError('s3am failed with a "403 Forbidden" error  while obtaining '
                                       '(%s). Did you use the correct credentials?' % s3_url)
                elif exception[-1].startswith('400'):
                    raise RuntimeError('s3am failed with a "400 Bad Request" error while obtaining '
                                       '(%s). Are you trying to download an encrypted file without '
                                       'a key, or an unencrypted file with one?' % s3_url)
                else:
                    raise RuntimeError('s3am failed with (%s) while downloading (%s)' %
                                       (': '.join(exception), s3_url))
            elif exception.startswith('AttributeError'):
                exception = exception.split(': ')
                if exception[-1].startswith("'NoneType'"):
                    raise RuntimeError('Does (%s) exist on s3?' % s3_url)
                else:
                    raise RuntimeError('s3am failed with (%s) while downloading (%s)' %
                                       (': '.join(exception), s3_url))
            else:
                if attempt < 3:
                    attempt += 1
                    continue
                else:
                    raise RuntimeError('Could not diagnose the error while downloading (%s)' %
                                       s3_url)
        except OSError:
            raise RuntimeError('Failed to find "s3am". Install via "apt-get install --pre s3am"')
        else:
            break
        finally:
            os.remove(stderr_file.name)
    assert os.path.exists(filename)
    if write_to_jobstore:
        filename = job.fileStore.writeGlobalFile(filename)
    return filename


def get_file_from_cghub(job, cghub_xml, cghub_key, univ_options, write_to_jobstore=True):
    """
    This function will download the file from cghub using the xml specified by cghub_xml

    ARGUMENTS
    1. cghub_xml: Path to an xml file for cghub.
    2. cghub_key: Credentials for a cghub download operation.
    3. write_to_jobstore: Flag indicating whether the final product should be written to jobStore.

    RETURN VALUES
    1. A path to the prefix for the fastqs that is compatible with the pipeline.
    """
    work_dir = job.fileStore.getLocalTempDir()
    # Get from S3 if required
    if cghub_xml.startswith('http'):
        assert cghub_xml.startswith('https://s3'), 'Not an S3 link'
        cghub_xml = get_file_from_s3(job, cghub_xml, encryption_key=univ_options['sse_key'],
                                     write_to_jobstore=False)
    else:
        assert os.path.exists(cghub_xml), 'Could not find file: %s' % cghub_xml
    shutil.copy(cghub_xml, os.path.join(work_dir, 'cghub.xml'))
    cghub_xml = os.path.join(work_dir, 'cghub.xml')
    assert os.path.exists(cghub_key), 'Could not find file: %s' % cghub_key
    shutil.copy(cghub_key, os.path.join(work_dir, 'cghub.key'))
    cghub_key = os.path.join(work_dir, 'cghub.key')
    temp_fastqdir = os.path.join(work_dir, 'temp_fastqdir')
    os.mkdir(temp_fastqdir)
    base_parameters = ['-d',  docker_path(cghub_xml),
                       '-c', docker_path(cghub_key),
                       '-p', docker_path(temp_fastqdir)]
    attemptNumber = 0
    while True:
        # timeout increases by 10 mins per try
        parameters = base_parameters + ['-k', str((attemptNumber + 1) * 10)]
        try:
            docker_call('genetorrent', tool_parameters=parameters, work_dir=work_dir,
                        dockerhub=univ_options['dockerhub'])
        except RuntimeError as err:
            time.sleep(600)
            job.fileStore.logToMaster(err.message)
            attemptNumber += 1
            if attemptNumber == 3:
                raise
            else:
                continue
        else:
            break
    analysis_id = [x for x in os.listdir(temp_fastqdir)
                   if not (x.startswith('.') or x.endswith('.gto'))][0]
    files = [x for x in os.listdir(os.path.join(temp_fastqdir, analysis_id))
             if not x.startswith('.')]
    if len(files) == 2:
        prefixes = [os.path.splitext(x)[1] for x in files]
        if {'.bam', '.bai'} - set(prefixes):
            raise RuntimeError('This is probably not a TCGA archive for WXS or RSQ. If you are ' +
                               'sure it is, email aarao@ucsc.edu with details.')
        else:
            bamfile = os.path.join(temp_fastqdir, analysis_id,
                                   [x for x in files if x.endswith('.bam')][0])
            return bam2fastq(job, bamfile, univ_options)
    elif len(files) == 1:
        if not files[0].endswith('.tar.gz'):
            raise RuntimeError('This is probably not a TCGA archive for WXS or RSQ. If you are ' +
                               'sure it is, email aarao@ucsc.edu with details.')
        else:
            outFastqDir = os.path.join(work_dir, 'fastqs')
            os.mkdir(outFastqDir)
            fastq_file = untargz(os.path.join(temp_fastqdir, analysis_id, files[0]), outFastqDir)
            if fastq_file.endswith(('.fastq', '.fastq.gz')):
                return re.sub('_2.fastq', '_1.fastq', fastq_file)
            else:
                raise RuntimeError('This is probably not a TCGA archive for WXS or RSQ. If you ' +
                                   'are sure it is, email aarao@ucsc.edu with details.')
    else:
        raise RuntimeError('This is probably not a TCGA archive for WXS or RSQ. If you are sure ' +
                           'it is, email aarao@ucsc.edu with details.')


def file_xext(filepath):
    """
    Get the file extension wrt compression from the filename (is it tar or targz)
    :param str filepath: Path to the file
    :return str ext: Compression extension name
    """
    ext = os.path.splitext(filepath)[1]
    if ext == '.gz':
        xext = os.path.splitext(os.path.splitext(filepath)[0])[1]
        if xext == '.tar':
            ext = xext + ext
    elif ext == '.tar':
        pass # ext is already .tar
    else:
        ext = ''
    return ext


def strip_xext(filepath):
    """
    Strips the compression extension from the filename
    :param filepath: Path to compressed file.
    :return str filepath: Path to the file with the compression extension stripped off.
    """
    ext_size = len(file_xext(filepath).split('.')) - 1
    for i in xrange(0, ext_size):
        filepath = os.path.splitext(filepath)[0]
    return filepath


def docker_path(filepath):
    return os.path.join('/data', os.path.basename(filepath))


def bam2fastq(job, bamfile, univ_options):
    """
    split an input bam to paired fastqs.

    ARGUMENTS
    1. bamfile: Path to a bam file
    2. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                |- 'dockerhub': <dockerhub to use>
                +- 'java_Xmx': value for max heap passed to java
    """
    work_dir = os.path.split(bamfile)[0]
    base_name = os.path.split(os.path.splitext(bamfile)[0])[1]
    parameters = ['SamToFastq',
                  ''.join(['I=', docker_path(bamfile)]),
                  ''.join(['F=/data/', base_name, '_1.fastq']),
                  ''.join(['F2=/data/', base_name, '_2.fastq']),
                  ''.join(['FU=/data/', base_name, '_UP.fastq'])]
    docker_call(tool='picard', tool_parameters=parameters, work_dir=work_dir,
                dockerhub=univ_options['dockerhub'], java_opts=univ_options['java_Xmx'])
    first_fastq = ''.join([work_dir, '/', base_name, '_1.fastq'])
    assert os.path.exists(first_fastq)
    return first_fastq


def run_cutadapt(job, uuid, fastq1, fastq2, adapter3, adapter5, is_gzipped=False, mock=False):
    """
    This module runs cutadapt on the input RNA fastq files and then calls the RNA aligners.

    ARGUMENTS
    1. fastqs: Dict of list of input RNA-Seq fastqs
         fastqs
            +- 'tumor_rna': [<JSid for 1.fastq> , <JSid for 2.fastq>]
    2. univ_options: Dict of universal arguments used by almost all tools
         univ_options
              +- 'dockerhub': <dockerhub to use>
    3. cutadapt_options: Dict of parameters specific to cutadapt
         cutadapt_options
              |- 'a': <sequence of 3' adapter to trim from fwd read>
              +- 'A': <sequence of 3' adapter to trim from rev read>
    RETURN VALUES
    1. output_files: Dict of cutadapted fastqs
         output_files
             |- 'rna_cutadapt_1.fastq': <JSid>
             +- 'rna_cutadapt_2.fastq': <JSid>

    This module corresponds to node 2 on the tree
    """
    job.fileStore.logToMaster('Running cutadapt on %s' % uuid)
    work_dir = job.fileStore.getLocalTempDir()
    fq_extn = '.gz' if is_gzipped else ''
    inputs = {'rna_1.fastq' + fq_extn: fastq1, 'rna_2.fastq' + fq_extn: fastq2}
    inputs = get_files_from_filestore(job, inputs, work_dir, docker=True)
    parameters = ['-a', adapter3,  # Fwd read 3' adapter
                  '-A', adapter5,  # Rev read 3' adapter
                  '-m', '35',  # Minimum size of read
                  '-o', docker_path('rna_cutadapt_1.fastq'),  # Output for R1
                  '-p', docker_path('rna_cutadapt_2.fastq'),  # Output for R2
                  inputs['rna_1.fastq'],
                  inputs['rna_2.fastq']]

    docker_call(tool='aarjunrao/cutadapt:latest', parameters=parameters, work_dir=work_dir,
                inputs=inputs.keys(), outputs={'rna_cutadapt_1.fastq': None, 'rna_cutadapt_2.fastq': None},
                mock=mock)

    file_id1 = job.fileStore.writeGlobalFile('/'.join([work_dir, 'rna_cutadapt_1.fastq']))
    file_id2 = job.fileStore.writeGlobalFile('/'.join([work_dir, 'rna_cutadapt_2.fastq']))
    return file_id1, file_id2




def index_bamfile(job, bamfile, sample_options, sample_type=''):
    """
    This module indexes BAMFILE
    ARGUMENTS
    1. bamfile: <JSid for a bam file>
    2. sample_type: string of 'tumor_dna' or 'normal_dna'
    3. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    RETURN VALUES
    1. output_files: REFER output_files in run_bwa(). This module is the one is
                     the one that generates the files.
    """
    sample_name = sample_options['patient_id']
    job.fileStore.logToMaster('Running samtools-index on %s:%s' % (sample_name,
                                                                   sample_type))
    work_dir = job.fileStore.getLocalTempDir()
    in_bamfile = '_'.join([sample_name, 'fix_pg_sorted.bam'])
    input_files = {in_bamfile: bamfile}
    input_files = get_files_from_filestore(job, input_files, work_dir,
                                           docker=True)
    parameters = ['index', input_files[in_bamfile]]
    docker_call(tool='aarjunrao/samtools:latest', parameters=parameters, work_dir=work_dir)
    output_files = {in_bamfile: bamfile,
                    in_bamfile + '.bai': job.fileStore.writeGlobalFile('/'.join([work_dir,
                                                                                 in_bamfile +
                                                                                 '.bai']))}
    return output_files





def run_rsem(job, star_bam, sample_options, rsem_options):
    """
    This module will run rsem on the RNA Bam file.

    ARGUMENTS
    1. star_bams: Dict of input STAR bams
         star_bams
              +- 'rnaAligned.toTranscriptome.out.bam': <JSid>
    2. univ_options: Dict of universal arguments used by almost all tools
         univ_options
                +- 'dockerhub': <dockerhub to use>
    3. rsem_options: Dict of parameters specific to rsem
         rsem_options
              |- 'index_tar': <JSid for the rsem index tarball>
              +- 'n': <number of threads to allocate>

    RETURN VALUES
    1. output_file: <Jsid of rsem.isoforms.results>

    This module corresponds to node 9 on the tree
    """
    job.fileStore.logToMaster('Running rsem index on %s' % sample_options['patient_id'])
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {'star_transcriptome.bam': star_bam,
                   'rsem_index.tar.gz': rsem_options['index']}
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=True)
    parameters = ['--paired-end',
                  '-p', str(rsem_options['ncores']),
                  '--bam',
                  input_files['star_transcriptome.bam'],
                  '--no-bam-output',
                  '/'.join([rsem_options['index'], 'hg19']),
                  'rsem']

    result = {'rsem.isoforms.results': rsem_options['output']}
    docker_call(tool='aarjunrao/rsem:latest', parameters=parameters, work_dir=work_dir,
                mock=True, outputs=result)
    output_file = \
        job.fileStore.writeGlobalFile('/'.join([work_dir, 'rsem.isoforms.results']))
    return output_file

