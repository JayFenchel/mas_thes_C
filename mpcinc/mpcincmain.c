#include <stdlib.h>
#include "mpcinccvp.h"
#include "mpcinccvpdynmem.h"
#include "mpcincfgm.h"
#include "mpcincfgmdynmem.h"

int main(void) {

struct mpcinc_cvp *cvp = mpcinc_cvp_allocate_former();
struct mpcinc_fgm *fgm = mpcinc_fgm_allocate_solver();

if (mpcinc_cvp_setup_former(cvp, "data.json")) {
  return 0;
}
if (mpcinc_fgm_setup_solver(fgm, cvp->prb, "data.json")) {
  return 0;
}

cvp->par[MPCINC_X_K]->data[0] = 2.;
cvp->par[MPCINC_X_K]->data[1] = 3.;
cvp->par[MPCINC_XR]->data[0] = 5.;
cvp->par[MPCINC_XR]->data[1] = 7.;
cvp->par[MPCINC_XR]->data[2] = 5.;
cvp->par[MPCINC_XR]->data[3] = 7.;
cvp->par[MPCINC_XR]->data[4] = 11.;
cvp->par[MPCINC_XR]->data[5] = 13.;
cvp->par[MPCINC_UR]->data[0] = 17.;
cvp->par[MPCINC_UR]->data[1] = 19.;

fgm->u_ini[0] = -1000.;
fgm->u_ini[1] = 1000.;
fgm->conf->in_iter = 1;
fgm->conf->warm_start = 0;
fgm->j_in = &(fgm->conf->in_iter);
mpcinc_cvp_form_problem(cvp);
mpcinc_fgm_solve_problem(fgm);
printf("u0: %f\n", fgm->u_opt[0]);
printf("u0: %f\n", fgm->u_opt[1]);
mpcinc_fgm_solve_problem(fgm);
printf("u0: %f\n", fgm->u_opt[0]);
mpcinc_fgm_solve_problem(fgm);
printf("u0: %f\n", fgm->u_opt[0]);
fgm->conf->in_iter = 1;
fgm->conf->warm_start = 1;
fgm->j_in = &(fgm->conf->in_iter);
mpcinc_cvp_form_problem(cvp);
mpcinc_fgm_solve_problem(fgm);
printf("u0: %f\n", fgm->u_opt[0]);
mpcinc_fgm_solve_problem(fgm);
printf("u0: %f\n", fgm->u_opt[0]);
mpcinc_fgm_solve_problem(fgm);
printf("u0: %f\n", fgm->u_opt[0]);

  return 0;
}
