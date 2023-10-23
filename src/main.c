#include <stdio.h>
#include <math.h>

/* This is most likely defined in math.h as well,
 * but here it is rounded to 30 decimals just to be safe */
#define PI 3.14159265358979323846264338328

/* The name isn't quite correct, since it's
 * most likely 80 bits, but it'll do for now */
typedef long double triple;

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
	size_t steps;
	triple dt;
} sim_params;


/* The following four functions are the four equations that need to be
 * numerically solved. They are separated because each will be called 4
 * times per loop.
 * Source: https://en.wikipedia.org/wiki/Double_pendulum */

triple d_theta_1(triple t1, triple t2, triple p1, triple p2, constants c) {
	return (6/(c.m*pow(c.l, 2)))*(2*p1-3*cos(t1 - t2)*p2)/(16 - 9*pow(cos(t1 - t2), 2));
}

triple d_theta_2(triple t1, triple t2, triple p1, triple p2, constants c) {
	return (6/(c.m*pow(c.l, 2)))*(8*p2-3*cos(t1 - t2)*p1)/(16 - 9*pow(cos(t1 - t2), 2));
}

triple d_p_1(triple t1, triple t2, triple dt1, triple dt2, constants c) {
	return -0.5 * c.m * pow(c.l, 2) * (dt1 * dt2 * sin(t1 - t2) + 3*c.g*sin(t1)/c.l);
}

triple d_p_2(triple t1, triple t2, triple dt1, triple dt2, constants c) {
	return -0.5 * c.m * pow(c.l, 2) * (-dt1 * dt2 * sin(t1 - t2) + c.g*sin(t2)/c.l);
}

triple triple_abs(triple n) {
	n = n < 0 ? -n : n;
	return n;
}

void simulate(triple theta1_0, triple theta2_0,
							pend_state states[], constants c,
							sim_params params, triple *flip) {

	/* The first two instants are the same, because
	 * we take a numerical derivative later on */
	states[0].t1 = states[1].t1 = theta1_0;
	states[0].t2 = states[1].t2 = theta2_0;
	states[0].p1 = states[1].p1 = 0;
	states[0].p2 = states[1].p2 = 0;

	triple h = params.dt / 2;

	/* These are technically the intermediate values of the derivatives,
	 * but they need the same fields as the state. */
	pend_state k[4]; 

	int track_flip = flip != NULL;

	for (unsigned int i = 2; i < params.steps; ++i) {
		/* It's nigh unreadable if we don't define this */
		pend_state prev = states[i-1];
		
		/* 1st approximation */
		k[0].t1 = d_theta_1(prev.t1, prev.t2, prev.p1, prev.p2, c);
		k[0].t2 = d_theta_2(prev.t1, prev.t2, prev.p1, prev.p2, c);
		k[0].p1 = d_p_1(prev.t1, prev.t2,
				           (prev.t1 - states[i-2].t1)/params.dt,
				           (prev.t2 - states[i-2].t2)/params.dt, c);
		k[0].p2 = d_p_2(prev.t1, prev.t2,
				           (prev.t1 - states[i-2].t1)/params.dt,
				           (prev.t2 - states[i-2].t2)/params.dt, c);

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

		states[i].t1 = prev.t1 + (k[0].t1 + 2*k[1].t1 + 2*k[2].t1 + k[3].t1)
				                     * params.dt/6;
		states[i].t2 = prev.t2 + (k[0].t2 + 2*k[1].t2 + 2*k[2].t2 + k[3].t2)
				                     * params.dt/6;
		states[i].p1 = prev.p1 + (k[0].p1 + 2*k[1].p1 + 2*k[2].p1 + k[3].p1)
				                     * params.dt/6;
		states[i].p2 = prev.p2 + (k[0].p2 + 2*k[1].p2 + 2*k[2].p2 + k[3].p2)
				                     * params.dt/6;
	
		if (track_flip && triple_abs(states[i].t2 > PI)) {
			*flip = i*params.dt;
			/* If we're just looking for the flipping point,
			 * there's no need to continue the simulation */
			return;
		}
	}
	if (track_flip)
		*flip = -1.0;
}

int main() {
	constants c = {1, 1, 9.81};
	sim_params params;
	params.steps = 100000;
	params.dt = 0.001;
	triple flip;
	pend_state results[100000];
	simulate(PI, PI, results, c, params, &flip);
	printf("Time to flip: %Lfs\n", flip);
	return 0;
}
