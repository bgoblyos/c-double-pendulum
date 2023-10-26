#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* This is most likely defined in math.h as well,
 * but here it is rounded to 30 decimals just to be safe */
#define PI 3.14159265358979323846264338328

/* The name isn't quite correct, since it's
 * most likely 80 bits, but it'll do for now */
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
	size_t flip_length; 
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

void full_sim(triple theta1_0, triple theta2_0, pend_state states[],
							constants c, sim_params params) {

	/* The first two instants are the same, because
	 * we take a numerical derivative later on */
	states[0].t1 = states[1].t1 = theta1_0;
	states[0].t2 = states[1].t2 = theta2_0;
	states[0].p1 = states[1].p1 = 0;
	states[0].p2 = states[1].p2 = 0;

	triple h = params.dt / 2;

	for (unsigned long int i = 2; i < params.steps; ++i)
		states[i] = step_sim(states[i-2], states[i-1], c, h);
}

void save_sim_data(pend_state *states, sim_params params) {
	unsigned long int skip = params.freq / params.plot_freq;

  FILE *upper = fopen("data/temp_upper.dat", "w");
  FILE *lower = fopen("data/temp_lower.dat", "w");
	for (unsigned long int i = 0; i < params.steps; i += skip) {
		fprintf(upper, "%Lf %Lf\n", states[i].t1, states[i].p1);
		fprintf(lower, "%Lf %Lf\n", states[i].t2, states[i].p2);
	}
	fclose(upper);
	fclose(lower);
}

void plot_phase_space(char *filename) {
	FILE *gnuplot = popen("gnuplot", "w");
	fprintf(gnuplot, "set term svg size 1920,1920 rounded background rgb 'white'\nset output \"%s\"\n", filename);
	fprintf(gnuplot, "set xlabel \"angle\"\nset ylabel \"impulse\"\n");
	fprintf(gnuplot, "plot \"data/temp_upper.dat\" u 1:2 t 'Upper' w l, \"data/temp_lower.dat\" u 1:2 t 'Lower' w l\n");
	fflush(gnuplot);
	pclose(gnuplot);
}

triple flip_sim(triple theta1, triple theta2, constants c, sim_params params) {
	pend_state old, prev, current;
	old.t1 = prev.t1 = theta1;
	old.t2 = prev.t2 = theta2;
	prev.p1 = prev.p2 = 0;

	for (unsigned long int i = 0; i < params.steps; ++i) {
		current = step_sim(old, prev, c, params.dt/2);
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
	for (ulong i = 0; i < length; ++i)
		result[i] = i*step - PI;
	return result;
}

triple** matrix(ulong length) {
	triple **result = (triple**)malloc(length*sizeof(triple*));
	result[0] = (triple*)malloc(length*length*sizeof(triple));
	for (unsigned int i = 1; i < length; ++i)
		result[i] = result[0] + i * length;
	return result;
}

/* This function takes in the flip time between
 * -1 and t, and returns a value between 0 and 1 */
triple normalize(triple in, triple max) {
	return (in + 1)/(max+1);
}

void flip_plot(triple **data, char *filename, sim_params params) {
	/* Writing the PPM file is done according to this StackOverflow answer:
	 * https://stackoverflow.com/a/4346905 */
	FILE *f = fopen(filename, "wb");
	fprintf(f, "P6\n%lu %lu 255\n", params.flip_length, params.flip_length);
	for (ulong i = 0; i < params.flip_length; i++) {
    for (ulong j = 0; j < params.flip_length; j++) {
			unsigned char col = (unsigned char)(255*normalize(data[i][j], params.t));
      /* Writing the pixel's RGB data
			 * The first parameters may be tweaked
			 * to get different color schemes      */
			fputc(col, f);         /* RED   */
      fputc(0, f);           /* GREEN */
      fputc((255-col)/5, f); /* Blue  */
    }
	}
	fclose(f);
	printf("Plot written to %s\n", filename);
}

triple **flip_matrix(sim_params params, constants c) {
	triple *thetas = linspace(params.flip_length);
	triple **results = matrix(params.flip_length);
	for (ulong i = 0; i < params.flip_length; ++i) {
		for (ulong j = 0; j < params.flip_length; ++j) {
			results[i][j] = flip_sim(thetas[i], thetas[j], c, params);
		}
		printf("Row %3lu/%lu computed\n", i + 1, params.flip_length);
	}
	return results;
}

int main() {
	constants c = {1, 1, 9.81};
	sim_params params;
	params.t = 60;
	params.flip_length = 128;
	params.freq = 1000;
	params.dt = (triple)1/params.freq;
	params.steps = (ulong)(params.t * params.freq);
	params.plot_freq = 1000;
	triple **flips = flip_matrix(params, c);
	flip_plot(flips, "data/out.ppm", params);
	return 0;
}
