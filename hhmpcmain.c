#include <stdio.h>

#include "include/mpcincmtxops.h"
#include "include/hhmpcsocp.h"
#include "include/mpcinccvpdynmem.h"
#include "include/hhmpcsocpdynmem.h"

#include "include/mpcincfgm.h"
#include "include/hhmpcipm.h"
#include "include/mpcincfgmdynmem.h"
#include "include/hhmpcipmdynmem.h"

#include "include/hhmpcusefull.h"

int main(void) {
    
    /*struct mpcinc_cvp *cvp = mpcinc_cvp_allocate_former();*/
    struct hhmpc_socp *socp = hhmpc_socp_allocate_former();
    /*struct mpcinc_fgm *fgm = mpcinc_fgm_allocate_solver();*/
    struct hhmpc_ipm *ipm = hhmpc_ipm_allocate_solver();
    
    if (hhmpc_socp_setup_former(socp, "test02data.json")) {
        return 0;
    }
    if (hhmpc_ipm_setup_solver(ipm, socp->prb, "test02data.json")) {
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

    ipm->conf->in_iter = 32;
    ipm->conf->reg = .00000001;
    ipm->conf->warm_start = 1;
    
    hhmpc_socp_form_problem(socp);
    hhmpc_ipm_solve_problem(ipm);
//     print_mtx(ipm->z_ini, ipm->optvar_seqlen,1);
    printf("u_opt1 = %f\n", ipm->z_opt[0]);
    printf("u_opt2 = %f\n", ipm->z_opt[31]);
    printf("u_opt3 = %f\n", ipm->z_opt[62]);
    printf("u_opt4 = %f\n", ipm->z_opt[93]);
    printf("u_opt5 = %f\n", ipm->z_opt[124]);
    hhmpc_ipm_solve_problem(ipm);
    
    
    
    printf("u_opt1 = %f\n", ipm->z_opt[0]);
    printf("u_opt2 = %f\n", ipm->z_opt[31]);
    printf("u_opt3 = %f\n", ipm->z_opt[62]);
    printf("u_opt4 = %f\n", ipm->z_opt[93]);
    printf("u_opt5 = %f\n", ipm->z_opt[124]);
    
    printf("%f \n", socp->constant[HHMPC_R_KL]->data[0]);
    printf("ENDE\n");
    return 0;
}

void sim_next_xk()
{
    mpcinc_mtx
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
