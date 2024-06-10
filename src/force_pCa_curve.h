#ifndef FORCE_PCA_CURVE_H
#define FORCE_PCA_CURVE_H
#include "particles.h"
void force_pCa_curve(initParticleArgs & args,
                     unsigned long randSeed,
                     float * ForceArrays,
                     float * Fss,
                     float * McArrays,
                     float * CArrays,
                     float * BArrays,
                     float * SRArrays,
                     int cc
                    );
#endif // FORCE_PCA_CURVE_H
