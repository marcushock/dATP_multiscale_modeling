#ifndef FORCE_PCA_CURVE_H
#define FORCE_PCA_CURVE_H
#include "particles.h"
void force_pCa_curve(initParticleArgs & args,
                     unsigned long randSeed,
                     float * Fss,
                     float * M1Arrays,
                     float * M2Arrays,
                     float * M3Arrays,
                     float * CArrays,
                     float * BArrays,
                     float * SRArrays,
                     float * ATPaseArrays,
                     int cc
                    );
#endif // FORCE_PCA_CURVE_H
