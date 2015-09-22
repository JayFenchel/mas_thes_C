#include <stdio.h>
#include <time.h>
// #include <sys/time.h>
// #include <unistd.h>

#include "include/mpcincmtxops.h"
#include "include/hhmpcsocp.h"
#include "include/mpcinccvpdynmem.h"
#include "include/hhmpcsocpdynmem.h"
#include "include/hhmpcmath.h"
#include "include/static_data.h"

#include "include/mpcincfgm.h"
#include "include/hhmpcipm.h"
#include "include/mpcincfgmdynmem.h"
#include "include/hhmpcipmdynmem.h"

#include "include/hhmpcusefull.h"


static void cmp_socp(struct hhmpc_socp *socp1, struct hhmpc_socp *socp2);

int main(void) {
    real_t time;
    clock_t begin, end;
    real_t htd[30];
    real_t wtd[30];
    real_t ut[30];
    struct hhmpc_socp socp_tmp;
    struct hhmpc_socp *socp, *socp_new;
    char *file = "test03data.json";
    /*struct mpcinc_cvp *cvp = mpcinc_cvp_allocate_former();*/
    struct hhmpc_socp *socp_old = hhmpc_socp_allocate_former();
    /*struct mpcinc_fgm *fgm = mpcinc_fgm_allocate_solver();*/
    struct hhmpc_ipm *ipm = hhmpc_ipm_allocate_solver();
    
    if (hhmpc_socp_setup_former(socp_old, file)) {
        return 0;
    }


socp_new = &socp_tmp;
     form_socp(socp_new);

socp = socp_new;
// socp = socp_old;

    
    if (hhmpc_ipm_setup_solver(ipm, socp->prb, file)) {
        return 0;
    }
/*    
    ipm->z_ini[0] = 0.9;
    ipm->z_ini[1] = 0.7;
    ipm->z_ini[2] = 0.2;
    ipm->z_ini[3] = 0.2;
    ipm->z_ini[4] = 0.5;
    ipm->z_ini[5] = 0.9;
    ipm->z_ini[6] = 0.7;
    ipm->z_ini[7] = 0.2;
    ipm->z_ini[8] = 0.2;
    ipm->z_ini[9] = 0.5;
    ipm->z_ini[10] = 0.9;
    ipm->z_ini[11] = 0.7;
    ipm->z_ini[12] = 0.2;
    ipm->z_ini[13] = 0.2;
    ipm->z_ini[14] = 0.5;
    ipm->z_ini[15] = 0.9;
    ipm->z_ini[16] = 0.7;
    ipm->z_ini[17] = 0.2;
    ipm->z_ini[18] = 0.2;
    ipm->z_ini[19] = 0.5;
    ipm->z_ini[20] = 0.9;
    ipm->z_ini[21] = 0.7;
    ipm->z_ini[22] = 0.2;
    ipm->z_ini[23] = 0.2;
    ipm->z_ini[24] = 0.5;
    ipm->z_ini[25] = 0.9;
    ipm->z_ini[26] = 0.7;
    ipm->z_ini[27] = 0.2;
    ipm->z_ini[28] = 0.2;
    ipm->z_ini[29] = 0.5;
    ipm->z_ini[30] = 0.9;
    ipm->z_ini[31] = 0.7;
    ipm->z_ini[32] = 0.2;
    ipm->z_ini[33] = 0.2;
    ipm->z_ini[34] = 0.5;
    ipm->v_ini[0] = 0.;
    ipm->v_ini[1] = 0.;
    ipm->v_ini[2] = 0.;
    ipm->v_ini[3] = 0.;
    ipm->v_ini[4] = 0.;
    ipm->v_ini[5] = 0.;
    ipm->v_ini[6] = 0.;
    ipm->v_ini[7] = 0.;
    ipm->v_ini[8] = 0.;
    ipm->v_ini[9] = 0.;
    ipm->v_ini[10] = 0.;
    ipm->v_ini[11] = 0.;
    ipm->v_ini[12] = 0.;
    ipm->v_ini[13] = 0.;
    ipm->v_ini[14] = 0.;
    ipm->v_ini[15] = 0.;
    ipm->v_ini[16] = 0.;
    ipm->v_ini[17] = 0.;
    ipm->v_ini[18] = 0.;
    ipm->v_ini[19] = 0.;
    ipm->v_ini[20] = 0.;
    ipm->v_ini[21] = 0.;
    ipm->v_ini[22] = 0.;
    ipm->v_ini[23] = 0.;
    ipm->v_ini[24] = 0.;
    ipm->v_ini[25] = 0.;
    ipm->v_ini[26] = 0.;
    ipm->v_ini[27] = 0.;
    ipm->v_ini[28] = 0.;
    ipm->v_ini[29] = 0.;
    */
    ipm->conf->in_iter = 12;
    ipm->conf->reg = .001;
    ipm->conf->warm_start = 1;
    begin = clock();
    hhmpc_socp_form_problem(socp);
//     hhmpc_socp_form_problem(socp_new);
//     cmp_socp(socp, socp_new);
    end = clock();
//     printf("begin:                         %d\n", begin);
//     printf("end:                           %d\n", end);
//     printf("clocks für form_problem: %d\n", end - begin);
    time=end - begin;
    time/=CLOCKS_PER_SEC;
    printf("Zeit für form_problem:   %f Sekunden\n", time);
//     printf("CLOCKS_PER_SEC:                %d\n", CLOCKS_PER_SEC);
    begin = clock();
    hhmpc_ipm_solve_problem(ipm);
    end = clock();
//     printf("begin:                         %d\n", begin);
//     printf("end:                           %d\n", end);
//     printf("clocks für solve_problem: %d\n", end - begin);
    time=end - begin;
    time/=CLOCKS_PER_SEC;
    printf("Zeit für solve_problem:   %f Sekunden\n", time);
//     printf("CLOCKS_PER_SEC:                %d\n", CLOCKS_PER_SEC);
//     print_mtx(ipm->z_ini, ipm->optvar_seqlen,1);
    printf("u_opt1 = %f\n", ipm->z_opt[0]);
    printf("u_opt2 = %f\n", ipm->z_opt[31]);
    printf("u_opt3 = %f\n", ipm->z_opt[62]);
    printf("u_opt4 = %f\n", ipm->z_opt[93]);
    printf("u_opt5 = %f\n", ipm->z_opt[124]);
    
//     print_mtx(ipm->z_opt, ipm->optvar_seqlen, 1);
//     
//     ipm->conf->in_iter = 8;
//     hhmpc_ipm_solve_problem(ipm);
//     printf("u_opt1 = %f\n", ipm->z_opt[0]);
//     printf("u_opt2 = %f\n", ipm->z_opt[31]);
//     printf("u_opt3 = %f\n", ipm->z_opt[62]);
//     printf("u_opt4 = %f\n", ipm->z_opt[93]);
//     printf("u_opt5 = %f\n", ipm->z_opt[124]);
// //     htd[0] = socp->prb->x_k->data[18];
// //     wtd[0] = socp->prb->x_k->data[6];
// //     ut[0] = ipm->z_opt[0];
// //     sim_next_xk(socp);
// //     htd[1] = socp->prb->x_k->data[18];
// //     wtd[1] = socp->prb->x_k->data[6];
// //     for (uint32_t i = 2; i< 30; i++){
// //   
// //     ipm->conf->in_iter = 8;
// //     hhmpc_socp_form_problem(socp);
// //     hhmpc_ipm_solve_problem(ipm);
// //     ut[i-1] = ipm->z_opt[0];
// // //     printf("u_opt1 = %f\n", ipm->z_opt[0]);
// // //     printf("u_opt2 = %f\n", ipm->z_opt[31]);
// // //     printf("u_opt3 = %f\n", ipm->z_opt[62]);
// // //     printf("u_opt4 = %f\n", ipm->z_opt[93]);
// // //     printf("u_opt5 = %f\n", ipm->z_opt[124]);
// //     sim_next_xk(socp);
// //     htd[i] = socp->prb->x_k->data[18];
// //     wtd[i] = socp->prb->x_k->data[6];
// //     }
// //     print_mtx(htd, 30, 1);
// //     print_mtx(wtd, 30, 1);
// //     printf("%f\n", htd[0]);
// //     hhmpc_socp_form_problem(socp);
    printf("%f \n", socp->constant[HHMPC_R_KL]->data[0]);
//     printf("%.15f \n", smpl_pow(E, 0.001));

//     print_mtx(socp->pmetric[HHMPC_ZR]->val->data, 155, 1);
    printf("ENDE\n");
    return 0;
}

void cmp_socp(struct hhmpc_socp *socp1, struct hhmpc_socp *socp2){
    
    printf("%d\n", mtx_cmp(socp1->prb->P->data, socp2->prb->P->data,
            socp1->prb->P->rows*socp1->prb->P->cols, 1e-10)
          );   
    printf("%d\n", mtx_cmp(socp1->prb->Psoft->data, socp2->prb->Psoft->data,
            socp1->prb->Psoft->rows*socp1->prb->Psoft->cols, 1e-10)
          );  
    
    printf("%d\n", mtx_cmp(socp1->prb->H->data, socp2->prb->H->data,
            socp1->prb->H->rows*socp1->prb->H->cols, 1e-10)
          );    
//     /*Not used.*/
//     printf("ffsoft: %d\n", mtx_cmp(socp1->prb->ffsoft->data, socp2->prb->ffsoft->data,
//             socp1->prb->ffsoft->rows, 1e-10)
//           );
    
    printf("%d\n", mtx_cmp(socp1->prb->v_ini->data, socp2->prb->v_ini->data,
            socp1->prb->v_ini->rows, 1e-10)
          );
    
    printf("%d\n", mtx_cmp(socp1->prb->z_ini->data, socp2->prb->z_ini->data,
            socp1->prb->z_ini->rows, 1e-10)
          );    
    
    printf("%d\n", mtx_cmp(socp1->prb->hsoft->data, socp2->prb->hsoft->data,
            socp2->prb->hsoft->rows, 1e-10)
          );
    
    printf("%d\n", mtx_cmp(socp1->prb->b->data, socp2->prb->b->data,
            socp1->prb->b->rows, 1e-10)
          );    
    
    printf("%d\n", mtx_cmp(socp1->prb->g->data, socp2->prb->g->data,
            socp1->prb->g->rows, 1e-10)
          );    
    
    printf("%d\n", mtx_cmp(socp1->prb->zref->data, socp2->prb->zref->data,
            socp1->prb->zref->rows, 1e-10)
          );
//     printf("%f\t%f\n", socp1->prb->zref->data[143], socp2->prb->zref->data[143]);
    
    printf("%d\n", mtx_cmp(socp1->prb->h->data, socp2->prb->h->data,
            socp1->prb->h->rows, 1e-10)
          );
//     printf("%f\t%f\n", socp1->prb->h->data[33], socp2->prb->h->data[33]);
    
}

/*// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include "include/hhmpctestfunc.h"
// #include "include/hhmpcmath.h"
// #include "include/hhmpcalg.h"
// #include "include/hhmpcusefull.h"
// 
// #define PPP     printf("help\n")
// 
// struct point makept(int x_koor, int y_koor);
// 
// struct point{
//     int x;
//     int y;
// };
// 
// struct point makept(int x, int y){
//     struct point temp;
//     temp.x = x;
//     temp.y = y;
//     return temp;
// }
// 
// int main(int argc, char *argv[]){
//     
//     real_t mtx[] = {4., 2., 0., 2., 5., 2., 0., 2., 5.};
//     real_t mtx2[9];
//     uint32_t dim = 3;
//     
//     struct point pt = makept(2, 3), *pp;
//     pp = &pt;
//     printf("Hallo Lea am Punkt %d, %d\n", pp->y, pp->x);
//     real_t test_a[] = {2, 2, 3, 3}, test_b[] = {4, 4}; 
//     printf("%f\n", mtx_out(test_a, 2, 2, test_b));
//     real_t sol[2]; 
//     real_t tests_a[] = {2, 0, 3, 1}, tests_b[] = {4, 4}; 
//     fwd_subst(sol, tests_a, 2., tests_b, 1);
//     printf("%f %f\n", sol[0], sol[1]);
//     printf("%d\n", mtx_cmp(sol, tests_b, 2, 0));
//     
//     real_t testss_a[] = {4., 1., 0, 2}, testss_b[] = {4, 4}; 
//     bwd_subst(sol, testss_a, 2., testss_b, 1);
//     printf("%f %f\n", sol[0], sol[1]);
//     printf("%f\n", smpl_sqrt(2., 2.));
//     cholesky(mtx2, mtx, dim);
//     eye(mtx, 3);
//     print_mtx(mtx, 3, 3);
//     
//     
//     
//    extern int c;
//    int i, *p1, *p2, v[LEN];
//    float *pc;
//    if (strcmp(*(argv+1), "echo") != -1){
//        printf("%d\n", strcmp(*(argv+1), "echo"));
//        printf("%d, %s\n", argc, *(argv+2));
//    }
//    
//    pc = alloc(6);
//    
//    pc[0] = 3.4;
//    *(pc+5) = 2.5;
//    pc[6] = 7.;
//    pc = alloc(2);
//    pc[0] = 4.4;
//    afree(pc-4);
//    
//    pc = alloc(2);
//    pc[0] = 4.3;
//    
//    aprint();
//    
//    i = 5;
//    p1 = &i;
//    printf("%d \n", *p1);
//    
//    for (i = 0; i < LEN; ++i){
//        printf("%d %d %d\n", i, power(2, i), power(-3, i));
//        v[i] = power(2, i);
//    }
//    
//    p2 = v;
//    for (i=0; i < LEN; ++i){
//        printf("%d ", v[i]);
//        printf("%d ", p2[i]);
//    }
//    printf("\n");
//    for (i=0; i < LEN; ++i){
//        printf("%d ", *v);
//        printf("%d ", *(v+i));
//    }
//    printf("\n");
//    for (i=0; i < LEN; ++i){
//        printf("%d ", *p2);
//        printf("%d ", *(p2+i));
//    }
//    printf("\n");
//    for (i=0; i < LEN; ++i){
//        printf("%d ", *(p2++));
//    }
//    printf("\n");
//    
//    c = SQUARE(1+1+1);
//    change(v, v);
//    for (i = 0; i < LEN; ++i)
//        printf("%d\n", c);
//    PPP;
//     return 0;
// }


// int main(/*int argc, char **argv*/
//) {
//     
//     int c, d;
//     long n[2];
//     printf(TEXT\n);
//  
// //     c = getchar() == EOF;
//     c = getchar();
//     printf("%d", c);
//     n[0] = 0;
//     n[1] = 0;
//     while ((d = getchar()) != EOF) {
//         if (d == '\n')
//             ++n[1];
//         putchar(d);
//         ++n[0];
// //         printf("%d", EOF);
// //         printf("%d", c);
// //         c = c + 1;
//             
//     }
//     printf("nc%ld nl%ld", n[0], n[1]);
//     for (n[0] = 0; getchar() != EOF; ++n[0])
//         ;
//     printf("%ld", n[0]);
//     return 0;
// }
//*/
