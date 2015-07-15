#!/usr/bin/env python2.7
# John Vivian
# 7-15-15

"""
Tree Structure of GATK Pipeline
     0-----> 11 ---- 12
    / \
   1   2
   |   |
   3   4
   |   |
   5   6
   |   |
   7   8
   |   |
   9   10
0  = create .dict/.fai for reference genome
1,2 = samtools index
3,4 = RealignerTargetCreator
5,6 = Indel Realignment
7,8 = Base Recalibration
9,10 = Recalibrate (PrintReads)
11 = MuTect
12 = teardown / cleanup
1-10, and 12 are "Target children"
11 is a "Target follow-on", it is executed after completion of children.
=========================================================================
:Directory Structure:
local_dir = /mnt/
# For "shared" input files
shared_dir = <local_dir>/<script_name>/<UUID4>
# For files specific to a tumor/normal pair
pair_dir = <local_dir>/<script_name>/<UUID4>/<pair>/
    <pair> is defined as UUID-normal:UUID-tumor
files are uploaded to:
    s3://bd2k-<script>/<UUID4>/ if shared (.fai/.dict)
    s3://bd2k-<script>/<UUID4>/<pair> if specific to that T/N pair.
=========================================================================
:Dependencies:
curl            - apt-get install curl
samtools        - apt-get install samtools
picard-tools    - apt-get install picard-tools
boto            - pip install boto
FileChunkIO     - pip install FileChunkIO
jobTree         - https://github.com/benedictpaten/jobTree
Active Internet Connection (Boto)
"""
import argparse
import os
import shutil
import subprocess
import uuid
import errno
from jobTree.target import Target


def build_parser():
    """
    Contains arguments for the all of necessary input files
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', required=True, help="Reference Genome URL")
    parser.add_argument('-n', '--normal', required=True, help='Normal BAM URL. Format: UUID.normal.bam')
    parser.add_argument('-t', '--tumor', required=True, help='Tumor BAM URL. Format: UUID.tumor.bam')
    parser.add_argument('-p', '--phase', required=True, help='1000G_phase1.indels.hg19.sites.fixed.vcf URL')
    parser.add_argument('-m', '--mills', required=True, help='Mills_and_1000G_gold_standard.indels.hg19.sites.vcf URL')
    parser.add_argument('-d', '--dbsnp', required=True, help='dbsnp_132_b37.leftAligned.vcf URL')
    parser.add_argument('-c', '--cosmic', required=True, help='b37_cosmic_v54_120711.vcf URL')
    parser.add_argument('-g', '--gatk', required=True, help='GenomeAnalysisTK.jar')
    parser.add_argument('-u', '--mutect', required=True, help='Mutect.jar')
    parser.add_argument('-w', '--work_dir', required=True, help='Where you wanna work from? (full path please)')
    return parser


# Convenience functions used in the pipeline
def download_input(target, input_args, ids, name):
    """
    Downloads file if not present from supplied URL.
    Updates the FileStoreID to point to a file.
    :name: Key from self.input_urls.
    :returns: Path to file (work_dir path)
    :rtype: str
    """
    # Get path to file
    work_dir = input_args['work_dir']
    file_path = os.path.join(work_dir, name)

    # Create necessary directories if not present
    mkdir_p(work_dir)

    # Check if file exists, download if not present
    if not os.path.exists(file_path):
        try:
            subprocess.check_call(['curl', '-fs', input_args[name], '-o', file_path])
        except subprocess.CalledProcessError:
            raise RuntimeError('\nNecessary file could not be acquired: {}. Check input URL')
        except OSError:
            raise RuntimeError('Failed to find "curl". Install via "apt-get install curl"')

    assert os.path.exists(file_path)

    # Update FileStoreID
    target.fileStore.updateGlobalFile(ids[name], file_path)

    return file_path


def read_and_rename_global_file(target, input_args, file_store_id, new_extension='', alternate_name=None):
    """
    Given a FileStoreID, returns the filepath linked to it.
    :new_extension: Adds an extension to the file at the file_path.
    :alternate_name: A path to a filename that you want to use. (allows directory control as well as name)
    """
    work_dir = input_args['work_dir']
    name = target.fileStore.readGlobalFile(file_store_id)
    new_name = os.path.splitext(name if alternate_name is None else alternate_name)[0] + new_extension
    shutil.move(name, os.path.join(work_dir, os.path.basename(new_name)))

    # Move to work_dir so docker mount works
    # shutil.move(new_name, os.path.join(work_dir, os.path.basename(new_name)))

    return os.path.join(work_dir, os.path.basename(new_name))


def mkdir_p(path):
    """
    The equivalent of mkdir -p
    https://github.com/BD2KGenomics/bd2k-python-lib/blob/master/src/bd2k/util/files.py
    """
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def docker_path(filepath):
    return os.path.join('/data', os.path.basename(filepath))


def docker_call(input_args, tool_command, tool):
    """
    Makes subprocess call of a command to a docker container.
    Abstracts away the docker commands needed to run the tool.
    :type tool_command: str
    :type tool: str
    :param tool: the docker image for a given tool
    """
    work_dir = input_args['work_dir']
    base_docker_call = 'sudo docker run -v {}:/data'.format(work_dir)
    try:
        subprocess.check_call(base_docker_call.split() + [tool] + tool_command.split())
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status. Check error logs.')
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')


# Start of Target Functions
def start(target, input_args):
    # Create dictionary to hold empty FileStoreIDs
    symbolic_inputs = input_args.keys() + ['ref.fai', 'ref.dict', 'normal.bai', 'tumor.bai', 'normal.intervals',
                                           'tumor.intervals', 'normal.indel.bam', 'tumor.indel.bam',
                                           'normal.recal.table', 'tumor.recal.table', 'normal.bqsr.bam',
                                           'tumor.bqsr.bam', 'mutect.vcf']
    ids = {x: target.fileStore.getEmptyFileStoreID() for x in symbolic_inputs}
    tools = {'samtools': 'jvivian/samtools:1.2',
             'picard': 'jvivian/picardtools:1.113',
             'mutect': 'jvivian/mutect:1.1.7'}
    target_vars = (input_args, symbolic_inputs, ids, tools)

    target.addChildTargetFn(create_reference_index, target_vars)
    target.addChildTargetFn(create_reference_dict, target_vars)
    target.addChildTargetFn(create_normal_index, target_vars)
    target.addChildTargetFn(create_tumor_index, target_vars)

    # target.addFollowOnTargetFn(start_preprocessing, target_vars)
    # target.setFollowOnTargetFn(mutect, target_vars)


def create_reference_index(target, target_vars):
    """
    Uses Samtools to create reference index file (.fasta.fai)
    """
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve reference & store in FileStoreID
    ref_path = download_input(target, input_args, ids, 'ref.fasta')

    # Tool call
    command = 'samtools faidx {}'.format(docker_path(ref_path))
    docker_call(input_args, command, tool=tools['samtools'])

    # Update FileStoreID of output
    target.fileStore.updateGlobalFile(ids['ref.fai'], ref_path + '.fai')


def create_reference_dict(target, target_vars):
    """
    Uses Picardtools to create reference dictionary (.dict)
    """
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve reference & store in FileStoreID
    ref_path = download_input(target, input_args, ids, 'ref.fasta')

    # Tool call
    output = os.path.splitext(docker_path(ref_path))[0]
    command = 'picard-tools CreateSequenceDictionary R={} O={}.dict'.format(docker_path(ref_path), output)
    docker_call(input_args, command, tool=tools['picard'])

    # Update FileStoreID
    target.updateGlobalFile(ids['ref.dict'], os.path.splitext(ref_path)[0] + '.dict')


def create_normal_index(target, target_vars):
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve normal bam
    normal_path = download_input(target, input_args, ids, 'normal.bam')

    # Tool call
    command = 'samtools index {}'.format(docker_path(normal_path))
    docker_call(input_args, command, tool=tools['samtools'])

    # Update FileStoreID
    target.fileStore.updateGlobalFile(ids['normal.bai'], normal_path + '.bai')


def create_tumor_index(target, target_vars):
    # Unpack target variables
    input_args, symbolic_inputs, ids, tools = target_vars

    # Retrieve tumor bam
    tumor_path = download_input(target, input_args, ids, 'tumor.bam')

    # Tool call
    command = 'samtools index {}'.format(docker_path(tumor_path))
    docker_call(input_args, command, tool=tools['samtools'])

    # Update FileStoreID
    target.fileStore.updateGlobalFile(ids['tumor.bai'], tumor_path + '.bai')


if __name__ == '__main__':
    parser = build_parser()
    Target.Runner.addJobTreeOptions(parser)
    args = parser.parse_args()

    work_dir = os.path.join(str(args.work_dir),
                            'bd2k-{}'.format(os.path.basename(__file__).split('.')[0]), str(uuid.uuid4()))
    input_args = {'ref.fasta': args.reference,
                  'normal.bam': args.normal,
                  'tumor.bam': args.tumor,
                  'phase.vcf': args.phase,
                  'mills.vcf': args.mills,
                  'dbsnp.vcf': args.dbsnp,
                  'cosmic.vcf': args.cosmic,
                  'gatk.jar': args.gatk,
                  'mutect.jar': args.mutect,
                  'work_dir': work_dir}

    # Ensure user supplied URLs to files and that BAMs are in the appropriate format
    for bam in [args.normal, args.tumor]:
        if len(bam.split('/')[-1].split('.')) != 3:
            raise RuntimeError('{} BAM is not in the appropriate format: \
            UUID.normal.bam or UUID.tumor.bam'.format(str(bam).split('.')[1]))

    Target.Runner.startJobTree(Target.wrapTargetFn(start, input_args), args)
    Target.Runner.cleanup(args)
