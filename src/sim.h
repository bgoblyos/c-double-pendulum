/* Double inclusion guard */
#ifndef SIM_H_INCLUDED
#define SIM_H_INCLUDED

/* Prceice floating poit type for the calculations.
 * The name isn't quite correct, since it's
 * most likely 80 bits, but it'll do for now.
 * It's also the same as double on Windows. */
typedef long double triple;

/* Unsigled long for storing non-negative integers (like step count) */
typedef unsigned long int ulong;

/* Structure for storing the constant parameters of the simulation. */
typedef struct {
        triple l;
        triple m;
        triple g;
} constants;

/* Stores the pendulum's state (or its derivatives) at a point in time. */
typedef struct {
        triple t1;
        triple t2;
        triple p1;
        triple p2;
} pend_state;

/* Stores the variable parameters of the simulation. */
typedef struct {
        ulong steps;
        triple dt;
        triple t;
        ulong freq;
        ulong plot_freq;
        ulong flip_length;
				constants c;
} sim_params;

/* This function runs a simulation with the given parameters and stores every
 * intermediate state in a dynamic array. It returns a pointer to this array
 * when the simulation finishes. */
pend_state *full_sim(triple theta1_0, triple theta2_0, sim_params params);

/* Runs params.flip_length^2 simulations until the lower pendulum flips over
 * and returns a dynamic matrix with the time it took for each simulation
 * (-1 if the pendulum did not flip during the simulation). */
triple **flip_matrix(sim_params params);

#endif
