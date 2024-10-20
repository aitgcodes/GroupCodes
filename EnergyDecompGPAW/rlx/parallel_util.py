import sys


def get_parallel(comm):
    # Parallelization settings
    for domain_par in [8, 4, 2, 1]:
        if comm.size % domain_par == 0:
            band_par = comm.size // domain_par
            break
    assert domain_par * band_par == comm.size
    return {'sl_auto': True, 'domain': domain_par, 'band': band_par,
            'augment_grids': band_par > 1}


def get_logger(outer_comm, inner_comm):
    def log(*args, **kwargs):
        if outer_comm.rank == 0:
            print('[%04d/%04d]' % (inner_comm.rank, inner_comm.size),
                  *args, **kwargs)
            sys.stdout.flush()

    return log
