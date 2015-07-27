cdef extern from "../include/mpcinccvp.h":

    cdef enum:
        MPCINC_XR,
        MPCINC_X_K,
        MPCINC_UR,

        MPCINC_PAR_NUM

    cdef struct mpcinc_term:
        int rows
        int cols
        double *data

    cdef struct mpcinc_cvp_prb:
        mpcinc_term *g
        mpcinc_term *H
        mpcinc_term *u_lb
        mpcinc_term *u_ub


    cdef struct mpcinc_cvp:
        mpcinc_term *par[MPCINC_PAR_NUM]
        mpcinc_cvp_prb *prb

    void mpcinc_cvp_form_problem(mpcinc_cvp *cvp)

cdef extern from "../include/mpcinccvpdynmem.h":

    int mpcinc_cvp_setup_former(mpcinc_cvp *cvp, char *fname)
    mpcinc_cvp *mpcinc_cvp_allocate_former()

cdef _py2c_pardata(mpcinc_cvp *ccvp, pardata)
cdef _c2py_qp(mpcinc_cvp_prb *cprb)

