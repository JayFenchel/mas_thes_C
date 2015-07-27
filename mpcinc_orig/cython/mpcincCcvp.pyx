import numpy as np

cimport mpcincCcvp as Ccvp

cdef _py2c_pardata(Ccvp.mpcinc_cvp *ccvp, pardata):
    _py2c_term(ccvp.par[Ccvp.MPCINC_XR], pardata["xr"])
    _py2c_term(ccvp.par[Ccvp.MPCINC_X_K], pardata["x_k"])
    _py2c_term(ccvp.par[Ccvp.MPCINC_UR], pardata["ur"])


cdef _c2py_qp(Ccvp.mpcinc_cvp_prb *cprb):
    prb = dict()
    prb['g'] = _c2py_term(cprb.g)
    prb['H'] = _c2py_term(cprb.H)
    prb['u_lb'] = _c2py_term(cprb.u_lb)
    prb['u_ub'] = _c2py_term(cprb.u_ub)

    return prb

cdef _py2c_term(Ccvp.mpcinc_term *term, pydata):
    cdef int i
    cdef int j
    if pydata.shape != (term.rows, term.cols):
        print('Error: incompatible shapes')
        return

    for i in range(term.rows):
        for j in range(term.cols):
            term.data[i*term.cols + j] = pydata[i,j]
    return

cdef _c2py_term(Ccvp.mpcinc_term *term):
    a = np.zeros((term.rows, term.cols))

    cdef int i
    cdef int j
    for i in range(term.rows):
        for j in range(term.cols):
            a[i,j] = term.data[i*term.cols + j]
    return a

