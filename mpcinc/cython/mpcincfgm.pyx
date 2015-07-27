import json

from cython.view cimport array as cvarray

cimport mpcincCcvp as Ccvp
cimport mpcincCfgm as Cfgm

cdef class Solver:
    cdef Cfgm.mpcinc_fgm *fgm
    cdef Ccvp.mpcinc_cvp *cvp
    cdef cvarray u_opt
    cdef cvarray u_ini

    def __cinit__(self):
        self.fgm = Cfgm.mpcinc_fgm_allocate_solver()
        self.cvp = Ccvp.mpcinc_cvp_allocate_former()

    def __dealloc__(self):
        # TODO: free memory allocated by malloc
        pass

    property u_opt:
        def __get__(self):
          return self.u_opt

    property u_ini:
        def __get__(self):
          return self.u_ini
        def __set__(self, u_ini):
          cdef double [:]u_ini_cv = u_ini
          self.u_ini[:] = u_ini_cv

    cpdef setup_solver(self, data, fname='data.json'):
        cdef int optvar_seqlen
        # TODO: check data has consistent sizes
        with open(fname, 'w') as f:
            json.dump(data, f)

        fname = fname.encode()  # python 3, char *
        Ccvp.mpcinc_cvp_setup_former(self.cvp, fname)
        cdef Ccvp.mpcinc_cvp_prb *prb
        prb = self.cvp.prb
        Cfgm.mpcinc_fgm_setup_solver(self.fgm, prb, fname)
        optvar_seqlen = self.fgm.optvar_seqlen
        self.u_opt = <double[:optvar_seqlen]> self.fgm.u_opt
        self.u_ini = <double[:optvar_seqlen]> self.fgm.u_ini

    cpdef solve_problem(self, pardata):
        cdef double *u
        cdef Ccvp.mpcinc_cvp *ccvp = self.cvp
        Ccvp._py2c_pardata(ccvp, pardata)
        Ccvp.mpcinc_cvp_form_problem(self.cvp)
        Cfgm.mpcinc_fgm_solve_problem(self.fgm)
        return

    cpdef configure(self, int in_iter, int warm_start):
        self.fgm.conf.in_iter = in_iter
        self.fgm.conf.warm_start = warm_start

