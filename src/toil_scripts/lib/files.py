from contextlib import closing
import errno
import os
import shutil
import tarfile

from bd2k.util.files import mkdir_p

from toil_scripts.lib.urls import s3am_upload

def tarball_files(tar_name, file_paths, output_dir='.', prefix=''):
    """
    Creates a tarball from a group of files

    :param str tar_name: Name of tarball
    :param list[str] file_paths: Absolute file paths to include in the tarball
    :param str output_dir: Output destination for tarball
    :param str prefix: Optional prefix for files in tarball
    """
    with tarfile.open(os.path.join(output_dir, tar_name), 'w:gz') as f_out:
        for file_path in file_paths:
            if not file_path.startswith('/'):
                raise ValueError('Path provided is relative not absolute.')
            arcname = prefix + os.path.basename(file_path)
            f_out.add(file_path, arcname=arcname)


def copy_files(file_paths, output_dir):
    """
    Moves files from the working directory to the output directory.

    :param str output_dir: Output directory
    :param list[str] file_paths: Absolute file paths to move
    """
    for file_path in file_paths:
        if not file_path.startswith('/'):
            raise ValueError('Path provided is relative not absolute.')
        dest = os.path.join(output_dir, os.path.basename(file_path))
        shutil.copy(file_path, dest)


def copy_file_job(job, name, file_id, output_dir):
    """
    Job version of move_files for one file

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str name: Name of output file (including extension)
    :param str file_id: FileStoreID of file
    :param str output_dir: Location to place output file
    """
    work_dir = job.fileStore.getLocalTempDir()
    fpath = job.fileStore.readGlobalFile(file_id, os.path.join(work_dir, name))
    copy_files([fpath], output_dir)


def copy_to(filename, output_dir, work_dir=None):
    """
    Moves files from the working directory to the output directory.

    :param work_dir: the working directory
    :param output_dir: the output directory
    :param filenames: remaining arguments are filenames
    """
    if os.path.isabs(filename):
        origin = filename
        filename = os.path.basename(filename)
    elif work_dir is None:
        origin = os.path.abspath(filename)
    else:
        origin = os.path.join(work_dir, filename)
    dest = os.path.join(output_dir, filename)
    try:
        shutil.copytree(origin, dest)
    except OSError as e:
        if e.errno == errno.ENOTDIR:
            mkdir_p(output_dir)
            shutil.copy(origin, dest)
        else:
            raise e
    assert os.path.exists(dest)


def upload_or_move_job(job, filename, file_id, output_dir, ssec=None):
    job.fileStore.logToMaster('Writing {} to {}'.format(filename, output_dir))
    work_dir = job.fileStore.getLocalTempDir()
    filepath = job.fileStore.readGlobalFile(file_id, os.path.join(work_dir, filename))
    # are we moving this into a local dir, or up to s3?
    if output_dir.startswith('s3://'):
        s3am_upload(fpath=os.path.join(work_dir, filepath),
                    s3_dir=output_dir,
                    s3_key_path=ssec)
    elif not os.path.exists(os.path.join(output_dir, filename)):
        mkdir_p(output_dir)
        copy_to(filepath, output_dir, work_dir)
    else:
        job.fileStore.logToMaster("File already exists: {}".format(filename))


def consolidate_tarballs_job(job, fname_to_id):
    """
    Combine the contents of separate tarballs into one.
    Subdirs within the tarball will be named the keys in **fname_to_id

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param dict[str,str] fname_to_id: Dictionary of the form: file-name-prefix=FileStoreID
    :return: The file store ID of the generated tarball
    :rtype: str
    """
    work_dir = job.fileStore.getLocalTempDir()
    # Retrieve output file paths to consolidate
    tar_paths = []
    for fname, file_store_id in fname_to_id.iteritems():
        p = job.fileStore.readGlobalFile(file_store_id, os.path.join(work_dir, fname + '.tar.gz'))
        tar_paths.append((p, fname))
    # I/O
    # output_name is arbitrary as this job function returns a FileStoreId
    output_name = 'foo.tar.gz'
    out_tar = os.path.join(work_dir, output_name)
    # Consolidate separate tarballs into one
    with tarfile.open(os.path.join(work_dir, out_tar), 'w:gz') as f_out:
        for tar, fname in tar_paths:
            with tarfile.open(tar, 'r') as f_in:
                for tarinfo in f_in:
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        tarinfo.name = os.path.join(output_name, fname, os.path.basename(tarinfo.name))
                        f_out.addfile(tarinfo, fileobj=f_in_file)
    return job.fileStore.writeGlobalFile(out_tar)


def untargz(input_targz_file, untar_to_dir):
    """
    Expands a tar.gz file

    :param str input_targz_file:
    :param str untar_to_dir: path to untar-ed directory/file
    :return: path to untar-ed archive
    :rtype: str
    """
    assert tarfile.is_tarfile(input_targz_file), 'Not a tar file.'
    tarball = tarfile.open(input_targz_file)
    return_value = os.path.join(untar_to_dir, tarball.getmembers()[0].name)
    tarball.extractall(path=untar_to_dir)
    tarball.close()
    return return_value


def get_files_from_filestore(job, work_dir, input_dict):
    """
    Copies files from the file store to a work directory. Will expand tar/gzipped files.

    :param JobFunctionWrappingJob job: Toil Job instance
    :param str work_dir: current working directory
    :param dict input_dict: {filename: fileStoreID}
    :return: {filename: filepath}
    :rtype: dict
    """
    for name, fileStoreID in input_dict.iteritems():
        if not os.path.exists(os.path.join(work_dir, name)):
            file_path = job.fileStore.readGlobalFile(fileStoreID, os.path.join(work_dir, name))
        else:
            file_path = name
        if tarfile.is_tarfile(file_path):
            file_path = untargz(file_path, work_dir)
        input_dict[name] = file_path
    return input_dict
