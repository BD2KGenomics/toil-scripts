import os
import subprocess
import logging
from bd2k.util.exceptions import panic

_log = logging.getLogger(__name__)


def mock_mode():
    """
    Checks whether the ADAM_GATK_MOCK_MODE environment variable is set.
    In mock mode, all docker calls other than those to spin up and submit jobs to the spark cluster
    are stubbed out and dummy files are used as inputs and outputs.
    """
    return True if int(os.environ.get('TOIL_SCRIPTS_MOCK_MODE', '0')) else False


def docker_call(tool,
                parameters=None,
                work_dir='.',
                rm=True,
                env=None,
                outfile=None,
                inputs=None,
                outputs=None,
                docker_parameters=None,
                check_output=False,
                mock=None,
                save_dir=None):
    """
    Calls Docker, passing along parameters and tool.

    :param str tool: Name of the Docker image to be used (e.g. quay.io/ucsc_cgl/samtools)
    :param list[str] parameters: Command line arguments to be passed to the tool
    :param str work_dir: Directory to mount into the container via `-v`. Destination convention is /data
    :param bool rm: Set to True to pass `--rm` flag.
    :param dict[str,str] env: Environment variables to be added (e.g. dict(JAVA_OPTS='-Xmx15G'))
    :param bool sudo: If True, prepends `sudo` to the docker call
    :param file outfile: Pipe output of Docker call to file handle
    :param list[str] inputs: A list of the input files.
    :param dict[str,str] outputs: A dictionary containing the outputs files as keys with either None
                                  or a url. The value is only used if mock=True
    :param dict[str,str] docker_parameters: Parameters to pass to docker
    :param bool check_output: When True, this function returns docker's output
    :param bool mock: Whether to run in mock mode. If this variable is unset, its value will be determined by
                      the environment variable.
    """
    from toil_scripts.lib.urls import download_url
    from toil_scripts.lib.files import copy_files_to

    if mock is None:
        mock = mock_mode()
    if parameters is None:
        parameters = []
    if inputs is None:
        inputs = []
    if outputs is None:
        outputs = {}

    for filename in inputs:
        assert(os.path.isfile(os.path.join(work_dir, filename)))

    if mock:
        for filename, url in outputs.items():
            file_path = os.path.join(work_dir, filename)
            if url is None:
                # create mock file
                if not os.path.exists(file_path):
                    f = open(file_path, 'w')
                    f.write("contents") # FIXME
                    f.close()

            else:
                file_path = os.path.join(work_dir, filename)
                if not os.path.exists(file_path):
                    outfile = download_url(url, work_dir=work_dir, name=filename)
                assert os.path.exists(file_path)
        return
    
    base_docker_call = ['docker', 'run',
                        '--log-driver=none',
                        '-v', '{}:/data'.format(os.path.abspath(work_dir))]
    if rm:
        base_docker_call.append('--rm')
    if env:
        for e, v in env.iteritems():
            base_docker_call.extend(['-e', '{}={}'.format(e, v)])
    if docker_parameters:
        base_docker_call += docker_parameters

    _log.debug("Calling docker with %s." % " ".join(base_docker_call + [tool] + parameters))

    docker_call = base_docker_call + [tool] + parameters

    try:
        if outfile:
            subprocess.check_call(docker_call, stdout=outfile)
        else:
            if check_output:
                return subprocess.check_output(docker_call)
            else:
                subprocess.check_call(docker_call)
    # Fix root ownership of output files
    except:
        # Panic avoids hiding the exception raised in the try block
        with panic():
            _fix_permissions(base_docker_call, tool, work_dir)
    else:
        _fix_permissions(base_docker_call, tool, work_dir)

    for filename in outputs.keys():
        if not os.path.isabs(filename):
            filename = os.path.join(work_dir, filename)
        assert(os.path.isfile(filename))

def _fix_permissions(base_docker_call, tool, work_dir):
    """
    Fix permission of a mounted Docker directory by reusing the tool

    :param list base_docker_call: Docker run parameters
    :param str tool: Name of tool
    :param str work_dir: Path of work directory to recursively chown
    """
    base_docker_call.append('--entrypoint=chown')
    stat = os.stat(work_dir)
    command = base_docker_call + [tool] + ['-R', '{}:{}'.format(stat.st_uid, stat.st_gid), '/data']
    subprocess.check_call(command)
