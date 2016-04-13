from toil_scripts.lib import partitions


def sample_batcher_job(job, func, samples, *args):
    """
    Spawns a tree of jobs to avoid overloading the number of jobs spawned by a single parent.

    :param function func: Function to spawn dynamically, passes one sample as first argument
    :param list samples: Array of samples to be batched
    :param list args: any arguments to be passed to the function
    """
    # num_partitions isn't exposed as an argument in order to be transparent to the user.
    # The value for num_partitions is a tested value
    num_partitions = 100
    partition_size = len(samples) / num_partitions
    if partition_size > 1:
        for partition in partitions(samples, partition_size):
            job.addChildJobFn(sample_batcher_job, func, partition, *args)
    else:
        for sample in samples:
            job.addChildJobFn(func, sample, *args)
