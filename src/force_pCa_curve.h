#ifndef FORCE_PCA_CURVE_H
#define FORCE_PCA_CURVE_H
#include "particles.h"
void force_pCa_curve(initParticleArgs & args,
                     unsigned long randSeed,
                     float * M3Arrays,
                     float * Fss,
                     float * M1Arrays,
                     float * CArrays,
                     float * BArrays,
                     float * SRArrays,
                     int cc
                    );
#endif // FORCE_PCA_CURVE_H
