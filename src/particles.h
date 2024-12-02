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
    float gamma_B; // [unitless] - RU-RU cooperative coefficient
    float gamma_M; // [unitless] - XB-RU/RU-XB coopcoefficient (Note: gamma_M = mu_B)
    float mu_M; // [unitless] - Inter-RU XB-XB cooperative coefficient
    float k2_plus_ref;
    float k3_plus;
    float k4_plus_ref;
    float kB_plus_ref; //               - (p = plus)
    float kB_minus_ref; //               - (m = minus)
    float lambda; // [unitless] must be between (0,1)
    float kCa_plus_ref;
    float kCa_minus_ref;
    float percent_dATP;
    float k_force;
    float k_plus_SR_ref;
    float k_minus_SR_ref;
    int protocol;
    initParticleArgs(std::vector< std::pair<float, float> > experimentalData,
                     float gamma_B,
                     float gamma_M,
                     float mu_M,
                     float k2_plus_ref,
                     float k3_plus,
                     float k4_plus_ref,
                     float kB_plus_ref,
                     float kB_minus_ref,
                     float lambda,
                     float kCa_plus_ref,
                     float kCa_minus_ref,
                     float percent_dATP,
                     float k_force,
                     float k_plus_SR_ref,
                     float k_minus_SR_ref,
                     int protocol):experimentalData(experimentalData),
        gamma_B(gamma_B),
        gamma_M(gamma_M),
        mu_M(mu_M),
        k2_plus_ref(k2_plus_ref),
        k3_plus(k3_plus),
        k4_plus_ref(k4_plus_ref),
        kB_plus_ref(kB_plus_ref),
        kB_minus_ref(kB_minus_ref),
        lambda(lambda),
        kCa_plus_ref(kCa_plus_ref),
        kCa_minus_ref(kCa_minus_ref),
        percent_dATP(percent_dATP),
        k_force(k_force),
        k_plus_SR_ref(k_plus_SR_ref),
        k_minus_SR_ref(k_minus_SR_ref),
        protocol(protocol)
       {}


};

void init_particle(initParticleArgs & args);

#endif // PARTICLES_H
