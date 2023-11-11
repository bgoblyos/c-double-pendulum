#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979323846264338328

typedef long double triple;
typedef unsigned long int ulong;

typedef struct {
        triple l;
        triple m;
        triple g;
} constants;

typedef struct {
        triple t1;
        triple t2;
        triple p1;
        triple p2;
} pend_state;

typedef struct {
        ulong steps;
        triple dt;
        triple t;
        ulong freq;
        ulong plot_freq;
        ulong flip_length;
				constants c;
} sim_params;

triple d_theta_1(triple t1, triple t2, triple p1, triple p2, constants c) {
        return (6/(c.m*pow(c.l, 2)))*(2*p1-3*cosl(t1 - t2)*p2)/(16 - 9*pow(cosl(t1 - t2), 2));
}

triple d_theta_2(triple t1, triple t2, triple p1, triple p2, constants c) {
        return (6/(c.m*pow(c.l, 2)))*(8*p2-3*cosl(t1 - t2)*p1)/(16 - 9*pow(cosl(t1 - t2), 2));
}

triple d_p_1(triple t1, triple t2, triple dt1, triple dt2, constants c) {
        return -0.5 * c.m * pow(c.l, 2) * (dt1 * dt2 * sinl(t1 - t2) + 3*c.g*sinl(t1)/c.l);
}

triple d_p_2(triple t1, triple t2, triple dt1, triple dt2, constants c) {
        return -0.5 * c.m * pow(c.l, 2) * (-dt1 * dt2 * sinl(t1 - t2) + c.g*sinl(t2)/c.l);
}

triple triple_abs(triple n) {
        return (n < 0 ? -n : n);
}

/* Function for stepping the simulation. This code snippet is huge and it needs
 * to be used in both simulation types, so I decided to put it into a
 * function on its own. It's probably a good idea to enable compiler
 * optimizations, since the code uses tons of function calls that might
 * as well be inlined, but were separated for readability. */
pend_state step_sim(pend_state old, pend_state prev, constants c, triple h) {

        /* These are technically the intermediate values of the derivatives,
         * but they need the same fields as the state. */
        pend_state k[4];

        pend_state new;

        /* 1st approximation */
                k[0].t1 = d_theta_1(prev.t1, prev.t2, prev.p1, prev.p2, c);
                k[0].t2 = d_theta_2(prev.t1, prev.t2, prev.p1, prev.p2, c);
                k[0].p1 = d_p_1(prev.t1, prev.t2,
                                           (prev.t1 - old.t1)/(2*h),
                                           (prev.t2 - old.t2)/(2*h), c);
                k[0].p2 = d_p_2(prev.t1, prev.t2,
                                           (prev.t1 - old.t1)/(2*h),
                                           (prev.t2 - old.t2)/(2*h), c);

                /* 2nd approximation */
                k[1].t1 = d_theta_1(prev.t1 + h*k[0].t1, prev.t2 + h*k[0].t2,
                                                prev.p1 + h*k[0].p1, prev.p2 + h*k[0].p2, c);
                k[1].t2 = d_theta_2(prev.t1 + h*k[0].t1, prev.t2 + h*k[0].t2,
                                                prev.p1 + h*k[0].p1, prev.p2 + h*k[0].p2, c);
                k[1].p1 = d_p_1(prev.t1 + h*k[0].t1, prev.t2 + h*k[0].t2,
                                            k[0].t1, k[0].t2, c);
                k[1].p2 = d_p_2(prev.t1 + h*k[0].t1, prev.t2 + h*k[0].t2,
                                            k[0].t1, k[0].t2, c);

                /* 3rd approximation */
                k[2].t1 = d_theta_1(prev.t1 + h*k[1].t1, prev.t2 + h*k[1].t2,
                                                prev.p1 + h*k[1].p1, prev.p2 + h*k[1].p2, c);
                k[2].t2 = d_theta_2(prev.t1 + h*k[1].t1, prev.t2 + h*k[1].t2,
                                                prev.p1 + h*k[1].p1, prev.p2 + h*k[1].p2, c);
                k[2].p1 = d_p_1(prev.t1 + h*k[1].t1, prev.t2 + h*k[1].t2,
                                            k[1].t1, k[1].t2, c);
                k[2].p2 = d_p_2(prev.t1 + h*k[0].t1, prev.t2 + h*k[1].t2,
                                            k[1].t1, k[1].t2, c);

                /* 4th aproximation */
                k[3].t1 = d_theta_1(prev.t1 + 2*h*k[2].t1, prev.t2 + 2*h*k[2].t2,
                        prev.p1 + 2*h*k[2].p1, prev.p2 + 2*h*k[2].p2, c);
                k[3].t2 = d_theta_2(prev.t1 + 2*h*k[2].t1, prev.t2 + 2*h*k[2].t2,
                        prev.p1 + 2*h*k[2].p1, prev.p2 + h*k[2].p2, c);
                k[3].p1 = d_p_1(prev.t1 + 2*h*k[2].t1, prev.t2 + 2*h*k[2].t2,
                        k[2].t1, k[2].t2, c);
                k[3].p2 = d_p_2(prev.t1 + 2*h*k[2].t1, prev.t2 + 2*h*k[2].t2,
                        k[2].t1, k[2].t2, c);

                new.t1 = prev.t1 + (k[0].t1 + 2*k[1].t1 + 2*k[2].t1 + k[3].t1) * h/3;
                new.t2 = prev.t2 + (k[0].t2 + 2*k[1].t2 + 2*k[2].t2 + k[3].t2) * h/3;
                new.p1 = prev.p1 + (k[0].p1 + 2*k[1].p1 + 2*k[2].p1 + k[3].p1) * h/3;
                new.p2 = prev.p2 + (k[0].p2 + 2*k[1].p2 + 2*k[2].p2 + k[3].p2) * h/3;

                return new;
}

pend_state *full_sim(triple theta1_0, triple theta2_0, sim_params params) {

        pend_state *states =
                (pend_state*)malloc(params.steps*sizeof(pend_state));
        if (states == NULL)
                return NULL;
        /* The first two instants are the same, because
         * we take a numerical derivative later on */
        states[0].t1 = states[1].t1 = theta1_0;
        states[0].t2 = states[1].t2 = theta2_0;
        states[0].p1 = states[1].p1 = 0;
        states[0].p2 = states[1].p2 = 0;

        triple h = params.dt / 2;

        for (ulong i = 2; i < params.steps; ++i)
                states[i] = step_sim(states[i-2], states[i-1], params.c, h);

        return states;
}

triple flip_sim(triple theta1, triple theta2, sim_params params) {
        pend_state old, prev, current;
        old.t1 = prev.t1 = theta1;
        old.t2 = prev.t2 = theta2;
        prev.p1 = prev.p2 = 0;

        for (ulong i = 0; i < params.steps; ++i) {
                current = step_sim(old, prev, params.c, params.dt/2);
                if (triple_abs(current.t2) > PI)
                        return i*params.dt;
                old = prev;
                prev = current;
        }

        return -1;
}

/* This doesn't really need a long int, but on my
 * machine size_t is an unsigned long int and I
 * need to print it later, so I figure this is
 * probably a good way to not accidentally
 * break string format specifiers on other platforms. */
triple* linspace(ulong length) {
        triple step = (triple)(2*PI)/(length - 1);
        triple *result = (triple*)malloc(length*sizeof(triple));
        if (result == NULL)
                return NULL;
        for (ulong i = 0; i < length; ++i)
                result[i] = i*step - PI;
        return result;
}

triple** matrix(ulong length) {
        triple **result = (triple**)malloc(length*sizeof(triple*));
        if (result == NULL)
                return NULL;
        result[0] = (triple*)malloc(length*length*sizeof(triple));
        if (result[0] == NULL)
                return NULL;
        for (ulong i = 1; i < length; ++i)
                result[i] = result[0] + i * length;
        return result;
}

triple **flip_matrix(sim_params params) {
        triple *thetas = linspace(params.flip_length);
        triple **results = matrix(params.flip_length);
        if (thetas == NULL || results == NULL)
                return NULL;
        for (ulong i = 0; i < params.flip_length; ++i) {
                for (ulong j = 0; j < params.flip_length; ++j) {
                        results[i][j] = flip_sim(thetas[i], thetas[j], params);
                }
                printf("Row %3lu/%lu computed\n", i + 1, params.flip_length);
        }
        free(thetas);
        return results;
}
