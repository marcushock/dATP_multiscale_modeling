#ifndef PROBLEM_DEFINES_H
#define PROBLEM_DEFINES_H

#define N_RU 26 // Number of RUs // parameter defined here
#define N_S 9 // Number of states
#define MAX_REPS (400) // Max number used to repeat the simulation 32 blocks * 32 threads //150016
#define DT 20e-4f // fixed time step
#define MAX_TSTEPS (1500001) // Max number of time stepping 3 seconds
#define SEED (time(NULL))
#define RANDVAL (rand())
#endif // PROBLEM_DEFINES_H
