import os


def get_server():
    """ A function for determining which server the job is currently running on.

    This is meant to be run on a calculation node.
    """

    evars = os.environ

    # Try collecting server name from LSB or SLURM environment
    server = evars.get('LSB_EXEC_CLUSTER', evars.get('SLURM_CLUSTER_NAME'))
    if server:
        local = False
        return server, local
    else:
        local = True

    # Try collecting from user defined environment variable
    server = evars.get('CLUSTER')
    if not server:
        raise ValueError(
            'Server could not be identified with user, LSB, or SLURM environment vars')

    return server, local


def get_nnodes(server=None):
    """ Get the number of nodes being used in this environment.

    This is meant to be run on a calculation node.
    """

    if not server:
        server, local = get_server()

    evars = os.environ

    if not local:
        if server == 'slac':
            nnodes = len(set(evars['LSB_HOSTS'].split()))
        elif server == 'sherlock':
            nnodes = int(evars['SLURM_NNODES'])
        elif server == 'nersc':
            nnodes = int(evars['SLURM_NNODES'])
    else:
        # If running locally, assign node count of 1
        nnodes = 1

    return nnodes
