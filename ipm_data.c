#include "include/static_data.h"


uint32_t in_iter = 1;
uint32_t warm_start = 1;
real_t reg = 0.001;
// struct hhmpc_ipm_conf conf = {in_iter, warm_start, reg};
struct hhmpc_ipm_conf conf = {1, 1, 0.001};
struct hhmpc_ipm_P_hat P;
real_t kappa;
real_t roh;

void form_ipm(struct hhmpc_ipm *ipm, struct hhmpc_socp_prb *prb){
    ipm->conf = &conf;
    ipm->P_of_z = &P;
    ipm->kappa = &kappa;
    ipm->roh = &roh;
}