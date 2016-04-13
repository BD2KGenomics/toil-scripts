from toil_scripts.lib import partitions


def map_job(job, func, inputs, *args):
    """
    Spawns a tree of jobs to avoid overloading the number of jobs spawned by a single parent.
    This function is appropriate to use when batching samples greater than 1,000.

    :param function func: Function to spawn dynamically, passes one sample as first argument
    :param list inputs: Array of samples to be batched
    :param list args: any arguments to be passed to the function
    """
    # num_partitions isn't exposed as an argument in order to be transparent to the user.
    # The value for num_partitions is a tested value
    num_partitions = 100
    partition_size = len(inputs) / num_partitions
    if partition_size > 1:
        for partition in partitions(inputs, partition_size):
            job.addChildJobFn(map_job, func, partition, *args)
    else:
        for sample in inputs:
            job.addChildJobFn(func, sample, *args)
