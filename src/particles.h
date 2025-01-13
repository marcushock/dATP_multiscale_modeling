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
    float k2_plus_ref_baseline;
    float k2_plus_ref_drug;
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
    float eta;
    float g_Cb;
    float g_Ca;

    // float gamma_B; // [unitless] - RU-RU cooperative coefficient
    // float gamma_M; // [unitless] - XB-RU/RU-XB coopcoefficient (Note: gamma_M = mu_B)
    // float mu_M; // [unitless] - Inter-RU XB-XB cooperative coefficient
    // float k2_plus_ref;
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
        k2_plus_ref_baseline(argsVector[7]),
        k2_plus_ref_drug(argsVector[8]),
        k3_plus_baseline(argsVector[9]),
        k3_plus_drug(argsVector[10]),
        k4_plus_ref_baseline(argsVector[11]),
        k4_plus_ref_drug(argsVector[12]),
        kB_plus_ref(argsVector[13]),
        kB_minus_ref(argsVector[14]),
        kCa_plus_ref(argsVector[15]),
        kCa_minus_ref(argsVector[16]),
        percent_drug(argsVector[17]),
        lambda(argsVector[18]),
        gamma_B(argsVector[19]),
        gamma_M(argsVector[20]),
        mu_B(argsVector[21]),
        mu_M(argsVector[22]),
        q(argsVector[23]),
        r(argsVector[24]),
        x_preR(argsVector[25]),
        x_xb(argsVector[26]),
        conc_ADP(argsVector[27]),
        conc_ATP(argsVector[28]),
        conc_Pi(argsVector[29]),
        delta_G_ATP(argsVector[30]),
        alpha(argsVector[31]),
        eta(argsVector[32]),
        g_Cb(argsVector[33]),
        g_Ca(argsVector[34])
       {}


};

void init_particle(initParticleArgs & args, int replicate_number);

#endif // PARTICLES_H
