from contextlib import closing
import os
import errno
import tarfile
import shutil


# FIXME: replace with bd2k.util.files.mkdir_p
def mkdir_p(path):
    """
    It is Easier to Ask for Forgiveness than Permission

    :param str path: path to directory to create
    """
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def tarball_files(tar_name, fpaths, work_dir='.', prefix=''):
    """
    Creates a tarball from a group of files

    :param str tar_name: Name of tarball
    :param list[str] fpaths: File paths to include in the tarball
    :param str work_dir: Current working directory
    :param str prefix: Optional prefix for files in tarball
    """
    with tarfile.open(os.path.join(work_dir, tar_name), 'w:gz') as f_out:
        for fpath in fpaths:
            arcname = prefix + os.path.basename(fpath)
            f_out.add(fpath, arcname=arcname)


def move_files(filepaths, output_dir):
    """
    Moves files from the working directory to the output directory.

    :param str output_dir: Output directory
    :param list[str] filepaths: File paths to move
    """
    for fpath in filepaths:
        dest = os.path.join(output_dir, os.path.basename(fpath))
        shutil.move(fpath, dest)


def consolidate_tarballs_job(job, **fname_to_id):
    """
    Combine the contents of separate tarballs into one.
    Subdirs within the tarball will be named the keys in **fname_to_id

    :param dict[str,str] fname_to_id: Key word arguments of the form: file_name_prefix=FileStoreID
    :return str: FileStoreID
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
