cimport mpcincCcvp as Ccvp

cdef extern from "../include/mpcincfgm.h":

    cdef struct mpcinc_fgm_conf:
        int in_iter
        int warm_start

    cdef struct mpcinc_fgm:
        int *j_in
        double *u_ini
        double *u_opt
        int optvar_seqlen
        mpcinc_fgm_conf *conf

    void mpcinc_fgm_solve_problem(mpcinc_fgm *fgm)

cdef extern from "../include/mpcincfgmdynmem.h":

    mpcinc_fgm *mpcinc_fgm_allocate_solver()
    void mpcinc_fgm_setup_solver(mpcinc_fgm *fgm, Ccvp.mpcinc_cvp_prb *prb, char *fname)
