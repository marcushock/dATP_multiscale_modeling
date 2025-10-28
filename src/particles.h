#ifndef PARTICLES_H
#define PARTICLES_H

#include "problemDefines.h"
#include <vector>

class initParticleArgs
{
public:
    std::vector< std::pair<float, float> > experimentalData;
    //---------------------------------------------
    // model reference parameters that we need to optimize
    //--------------------------------------------
    float protocol;
    float k_force_baseline;
    float k_force_drug;
    float k_plus_SR_baseline;
    float k_plus_SR_drug;
    float k_minus_SR;
    float k_xb;
    float k1_plus_ref_baseline;
    float k1_plus_ref_drug;
    float k2_plus_baseline;
    float k2_plus_drug;
    float k3_plus_baseline;
    float k3_plus_drug;
    float k4_plus_ref_baseline;
    float k4_plus_ref_drug;
    float kB_plus_ref;
    float kB_minus_ref;
    float kCa_plus_ref;
    float kCa_minus_ref;
    float percent_drug;
    float lambda;
    float gamma_B;
    float gamma_M;
    float mu_B;
    float mu_M;
    float q;
    float r;
    float x_preR;
    float x_xb;
    float conc_ADP;
    float conc_ATP;
    float conc_Pi;
    float delta_G_ATP;
    float alpha;
    float beta;
    float eta;
    float g_Cb;
    float g_Ca;
    float k_plus_SS;
    float k_minus_SS;
    float K_D; 
    float coop_N;

    // float gamma_B; // [unitless] - RU-RU cooperative coefficient
    // float gamma_M; // [unitless] - XB-RU/RU-XB coopcoefficient (Note: gamma_M = mu_B)
    // float mu_M; // [unitless] - Inter-RU XB-XB cooperative coefficient
    // float k1_plus_ref;
    // float k3_plus;
    // float k4_plus_ref;
    // float kB_plus_ref; //               - (p = plus)
    // float kB_minus_ref; //               - (m = minus)
    // float lambda; // [unitless] must be between (0,1)
    // float kCa_plus_ref;
    // float kCa_minus_ref;
    // float percent_drug;
    // float k_force;
    // float k_plus_SR_ref;
    // float k_minus_SR_ref;
    // float protocol;
    initParticleArgs(std::vector< std::pair<float, float> > experimentalData, std::vector<float> argsVector)
        :experimentalData(experimentalData),
        protocol(argsVector[0]),
        k_force_baseline(argsVector[1]),
        k_force_drug(argsVector[2]),
        k_plus_SR_baseline(argsVector[3]),
        k_plus_SR_drug(argsVector[4]),
        k_minus_SR(argsVector[5]),
        k_xb(argsVector[6]),
        k1_plus_ref_baseline(argsVector[7]),
        k1_plus_ref_drug(argsVector[8]),
        k2_plus_baseline(argsVector[9]),
        k2_plus_drug(argsVector[10]),
        k3_plus_baseline(argsVector[11]),
        k3_plus_drug(argsVector[12]),
        k4_plus_ref_baseline(argsVector[13]),
        k4_plus_ref_drug(argsVector[14]),
        kB_plus_ref(argsVector[15]),
        kB_minus_ref(argsVector[16]),
        kCa_plus_ref(argsVector[17]),
        kCa_minus_ref(argsVector[18]),
        percent_drug(argsVector[19]),
        lambda(argsVector[20]),
        gamma_B(argsVector[21]),
        gamma_M(argsVector[22]),
        mu_B(argsVector[23]),
        mu_M(argsVector[24]),
        q(argsVector[25]),
        r(argsVector[26]),
        x_preR(argsVector[27]),
        x_xb(argsVector[28]),
        conc_ADP(argsVector[29]),
        conc_ATP(argsVector[30]),
        conc_Pi(argsVector[31]),
        delta_G_ATP(argsVector[32]),
        alpha(argsVector[33]),
        beta(argsVector[34]),
        eta(argsVector[35]),
        g_Cb(argsVector[36]),
        g_Ca(argsVector[37]),
        k_plus_SS(argsVector[38]),
        k_minus_SS(argsVector[39]),
        K_D(argsVector[40]),
        coop_N(argsVector[41])
       {}


};

void init_particle(initParticleArgs & args, int replicate_number);

#endif // PARTICLES_H
