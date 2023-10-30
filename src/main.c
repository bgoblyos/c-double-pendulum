#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "input.h"

/* This is most likely defined in math.h as well,
 * but here it is rounded to 30 decimals just to be safe */
#define PI 3.14159265358979323846264338328

/* The name isn't quite correct, since it's
 * most likely 80 bits, but it'll do for now.
 * It's also the same as double on Windows. */
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

/* Taken from https://stackoverflow.com/a/8465083 */
char* str_concat(const char *s1, const char *s2)
{
    char *result = malloc(strlen(s1) + strlen(s2) + 1);
    strcpy(result, s1);
    strcat(result, s2);
    return result;
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

pend_state *full_sim(triple theta1_0, triple theta2_0,
		constants c, sim_params params) {
	
	pend_state *states =
		(pend_state*)malloc(params.steps*sizeof(pend_state));
	/* The first two instants are the same, because
	 * we take a numerical derivative later on */
	states[0].t1 = states[1].t1 = theta1_0;
	states[0].t2 = states[1].t2 = theta2_0;
	states[0].p1 = states[1].p1 = 0;
	states[0].p2 = states[1].p2 = 0;

	triple h = params.dt / 2;

	for (ulong i = 2; i < params.steps; ++i)
		states[i] = step_sim(states[i-2], states[i-1], c, h);

	return states;
}

void save_sim_data(pend_state *states, sim_params params, char *fname) {
	ulong skip = params.freq / params.plot_freq;
	skip = skip < 1 ? 1 : skip;

	FILE *f = fopen(fname, "w");
	if (f == NULL) {
		printf("Could not open file for writing.\n");
		return;
	}
	for (ulong i = 0; i < params.steps; i += skip)
		fprintf(f, "%Lf, %Lf, %Lf, %Lf\n",
			states[i].t1, states[i].p1, states[i].t2, states[i].p2);
	fclose(f);
	printf("Data saved to %s\n", fname);
}

void plot_phase_space(pend_state *states, sim_params params, char *filename) {
	FILE *gnuplot;
	#ifdef _WIN32
		if (system("where gnuplot 2> nul 1> nul"))
			gnuplot = NULL;
		else
			gnuplot = _popen("gnuplot", "w");
	#else
		if (system("which gnuplot 2> /dev/null 1> /dev/null"))
			gnuplot = NULL;
		else
			gnuplot = popen("gnuplot", "w");
	#endif
	if (gnuplot == NULL) {
		printf("gnuplot could not be found, no plot will be saved.\n");
		return;
	}
	ulong skip = params.freq / params.plot_freq;
	skip = skip < 1 ? 1 : skip;
	fprintf(gnuplot, "set term svg size 1000,1000 rounded background rgb");
	fprintf(gnuplot, "'white'\nset output \"%s\"\n", filename);
	fprintf(gnuplot, "set xlabel \"angle\"\nset ylabel \"impulse\"\n");
	fprintf(gnuplot, "$dataset << EOD\n");
	for (ulong i = 0; i < params.steps; i+= skip)
		fprintf(gnuplot, "%Lf %Lf %Lf %Lf\n",
			states[i].t1, states[i].p1, states[i].t2, states[i].p2);
	fprintf(gnuplot, "EOD\n");
	fprintf(gnuplot, "plot $dataset using 1:2 t 'Upper' w l,");
	fprintf(gnuplot, "$dataset using 3:4 t 'Lower' w l\n");
	fflush(gnuplot);
	#ifdef _WIN32
		_pclose(gnuplot);
	#else
		pclose(gnuplot);
	#endif
	printf("Phase space plot saved to %s\n", filename);
}

triple flip_sim(triple theta1, triple theta2, constants c, sim_params params) {
	pend_state old, prev, current;
	old.t1 = prev.t1 = theta1;
	old.t2 = prev.t2 = theta2;
	prev.p1 = prev.p2 = 0;

	for (ulong i = 0; i < params.steps; ++i) {
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
	for (ulong i = 1; i < length; ++i)
		result[i] = result[0] + i * length;
	return result;
}

/* This function takes in the flip time between
 * -1 and t, and returns a value between 0 and 1 */
triple normalize(triple in, triple max) {
	return (in + 1)/(max+1);
}

/* Writing the PPM file is done according to this StackOverflow answer:
 * https://stackoverflow.com/a/4346905 */
void flip_plot(triple **data, char *filename, sim_params params) {
	FILE *f = fopen(filename, "wb");
	if (f == NULL) {
		printf("Failed to open %s for writing.", filename);
		return;
	}
	/* Write the magic number and the parameters of the image */
	fprintf(f, "P6\n%lu %lu 255\n", params.flip_length, params.flip_length);
	/* Iterating through the matrix */
	for (ulong i = 0; i < params.flip_length; i++) {
		for (ulong j = 0; j < params.flip_length; j++) {
			unsigned char col = (unsigned char)(255*normalize(data[i][j], params.t));
      			/* Writing the pixel's RGB data
			 * The first parameters may be tweaked
			 * to get different color schemes      */
			fputc(col, f);         /* Red   */
			fputc(0, f);           /* Green */
			fputc((255-col)/5, f); /* Blue  */
    		}
	}
	fclose(f);
	printf("Plot written to %s\n", filename);
}

void convert_plot(char *filename, char *target) {
	int has_convert = 0;
	#ifdef _WIN32
		has_convert =
			!system("where convert 2> nul 1> nul");
	#else
		has_convert =
			!system("which convert 2> /dev/null 1> /dev/null");
	#endif
	if (has_convert) {
		/* This isn't exactly elegant, I should find a better way
		 * to construct the command */
		char *command = str_concat("convert ", filename);
		command = str_concat(command, " ");
		command = str_concat(command, target);
		system(command);
		free(command);
		printf("File saved as %s\n", target);
		printf("Would you like to remove the original file? [y/N] ");
		char response = get_bool();
		if (response) {
			remove(filename);
			printf("Removed %s\n", filename);
		}
	}
	else {
		printf("convert could not be found,\
				no output will be produced.\n");
		printf("Please install ImageMagick to\
				export to other formats.\n");
	}

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
	free(thetas);
	return results;
}

char *to_dynamic(char *input) {
	char *output = (char*)malloc((strlen(input) + 1)*sizeof(char));
	strcpy(output, input);
	return output;
}

void general_setup(constants *c, sim_params *p) {
	ulong choice;
	while (1) {
		printf("\nGeneral options\n[1] m = %Lf kg\n", c->m);
		printf("[2] l = %Lf m\n[3] g = %Lf m/s^2\n", c->l, c->g);
		printf("[4] t = %Lf s\n[5] f = %lu Hz\n", p->t, p->freq);
		printf("[6] Exit\nPlease enter your choice [1-6]: ");
		choice = get_ulong(0);
		switch (choice) {
			case 1 :
				printf("Please enter new value for m [1]: ");
				c->m = get_triple(1);
				break;
			case 2 :
				printf("Please enter new value for l [1]: ");
				c->l = get_triple(1);
				break;
			case 3 :
				printf("Please enter new value for g [9.81]: ");
				c->g = get_triple(9.81);
				break;
			case 4 :
				printf("Please enter new value for t [60]: ");
				p->t = get_triple(60);
				break;
			case 5 :
				printf("Please enter new value for f [1000]: ");
				p->freq = get_ulong(1000);
				break;
			default :
				return;
		}
	}
}

void full_setup(constants *c, sim_params *p, triple *theta1, triple *theta2, char *csv_def, char *svg_def) {
	ulong choice;
	int sim_done = 0;
	char *csv_fname = to_dynamic(csv_def);
	char *svg_fname = to_dynamic(svg_def);
	
	pend_state *result;
	while (1) {
		printf("\nFull trajectory simulation options\n[1] Theta 1 = %Lf\n", *theta1);
		printf("[2] Theta 2 = %Lf\n[3] Plotting frequency: %lu Hz\n", *theta2, p->plot_freq);
		printf("[4] Run simulation\n[5] Save data to csv\n");
		printf("[6] Plot phase space\n");
		printf("[7] Exit\nPlease enter your choice [1-7]: ");
		choice = get_ulong(0);
		switch (choice) {
			case 1 :
				printf("Please enter new value for Theta 1 [0]: ");
				*theta1 = get_triple(0);
				sim_done = 0;
				break;
			case 2 :
				printf("Please enter new value for Theta 2 [0]: ");
				*theta2 = get_triple(0);
				sim_done = 0;
				break;
			case 3 :
				printf("Please enter new value for plotting frequency [1000]: ");
				p->plot_freq = get_ulong(1000);
				break;
			case 4 :
				printf("Started simulation\n");
				result = full_sim(*theta1, *theta2, *c, *p);
				sim_done = 1;
				break;
			case 5 :
				if (!sim_done) {
					printf("No up-to-date simulation found, starting it\n");
					result = full_sim(*theta1, *theta2, *c, *p);
					sim_done = 1;
				}
				printf("Enter filename for CSV file [%s]: ", csv_fname);
				csv_fname = get_fname(csv_fname);
				save_sim_data(result, *p, csv_fname);
				break;
			case 6 :
				if (!sim_done) {
					printf("No up-to-date simultion found, starting it\n");
					result = full_sim(*theta1, *theta2, *c, *p);
					sim_done = 1;
				}
				printf("Enter filename for SVG plot [%s]: ", svg_fname);
				svg_fname = get_fname(svg_fname);
				plot_phase_space(result, *p, svg_fname);
				break;
			default:
				return;
		}
	}
	free(result);
	free(csv_fname);
	free(svg_fname);
}

void flip_setup(constants *c, sim_params *p, char *ppm_def, char *img_def) {
	ulong choice;
	int sim_done = 0;
	char *ppm_fname = to_dynamic(ppm_def);
	char *img_fname = to_dynamic(img_def);
	
	triple **result;

	while (1) {
		printf("\nFlipover simulation options\n");
		printf("[1] Pixels per side: %lu\n", p->flip_length);
		printf("[2] Run simulation\n[3] Save output to PPM\n");
		printf("[4] Convert PPM to another image format\n");
		printf("[5] Exit\nPlease enter your choice [1-5]: ");
		fflush(stdin);
		choice = get_ulong(0);
		switch (choice) {
			case 1 :
				printf("Please enter new value for side length [32]: ");
				p->flip_length = get_ulong(32);
				sim_done = 0;
				break;
			case 2 :
				printf("Started simulation\n");
				result = flip_matrix(*p, *c);
				sim_done = 1;
				break;
			case 3 :
				if (!sim_done) {
					printf("No up-to-date simulation found, starting it\n");
					result = flip_matrix(*p, *c);
					sim_done = 1;
				}
				printf("Enter filename for PPM file [%s]: ", ppm_fname);
				ppm_fname = get_fname(ppm_fname);
				flip_plot(result, ppm_fname, *p);
				break;
			case 4 :
				printf("Enter filename for PPM file [%s]: ", ppm_fname);
				ppm_fname = get_fname(ppm_fname);
				printf("Enter filename for target file [%s]: ", img_fname);
				img_fname = get_fname(img_fname);
				convert_plot(ppm_fname, img_fname);
				break;
			default:
				return;
		}
	}
	free(result[0]);
	free(result);
	free(ppm_fname);
	free(img_fname);
}

int main() {
	char *csv_def = "data/sim.csv";
	char *svg_def = "data/phase_space.svg";
	char *ppm_def = "data/flip.ppm";
	char *img_def = "data/flip.png";
	/* Set default parameters */
	constants c = {1, 1, 9.81};
	triple theta1 = 0, theta2 = 0;
	sim_params params;
	params.t = 60;
	params.flip_length = 32;
	params.freq = 1000;
	params.dt = (triple)1/params.freq;
	params.steps = (ulong)(params.t * params.freq);
	params.plot_freq = 1000;
	int done = 0;
	ulong choice;
	while (!done) {
		printf("\nMain menu\n[1] General options\n");
		printf("[2] Full-trajectory simulation\n");
		printf("[3] Flipover time simulation\n");
		printf("[4] Exit\nPlease enter your choice [1-4]: ");
		choice = get_ulong(0);

		switch (choice) {
			case 1: general_setup(&c, &params); break;
			case 2:
				full_setup(&c, &params, &theta1, &theta2, csv_def, svg_def);
				break;
			case 3: 
				flip_setup(&c, &params, ppm_def, img_def);
				break;
			default:
				done = 1;
				break;
		}
	}

	return 0;
}
