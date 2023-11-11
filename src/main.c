#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "input.h"
#include "sim.h"

/* Taken from https://stackoverflow.com/a/8465083 */
char* str_concat(const char *s1, const char *s2)
{
    char *result = malloc(strlen(s1) + strlen(s2) + 1);
    strcpy(result, s1);
    strcat(result, s2);
    return result;
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
	int has_magick = 0;
	#ifdef _WIN32
		has_magick =
			!system("where magick 2> nul 1> nul");
	#else
		has_magick =
			!system("which magick 2> /dev/null 1> /dev/null");
	#endif
	if (has_magick) {
		/* This isn't exactly elegant, I should find a better way
		 * to construct the command */
		char *command = str_concat("magick ", filename);
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
		printf("magick utility could not be found.\n");
		printf("Please install ImageMagick to export to other formats.\n");
	}

}

/* This is not expected to receive huge input,
 * so I'm not handling the NULL for now. */
char *to_dynamic(char *input) {
	char *output = (char*)malloc((strlen(input) + 1)*sizeof(char));
	strcpy(output, input);
	return output;
}

void free_array(pend_state *arr) {
	if (arr != NULL)
		free(arr);
}

void free_matrix(triple **mtr) {
	if (mtr != NULL) {
		if (mtr[0] != NULL)
			free(mtr[0]);
		free(mtr);
	}
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
	
	pend_state *result = NULL;
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
				free_array(result);
				printf("Started simulation\n");
				result = full_sim(*theta1, *theta2, *c, *p);
				if (result == NULL) {
					sim_done = 0;
					printf("Failed to allocate memory for results.\n");
				}
				else
					sim_done = 1;
				break;
			case 5 :
				if (!sim_done) {
					free_array(result);
					printf("No up-to-date simulation found, starting it\n");
					result = full_sim(*theta1, *theta2, *c, *p);
					if (result == NULL) {
						printf("Failed to allocate memory for results.\n");
						break;
					}
					else
						sim_done = 1;
				}
				printf("Enter filename for CSV file [%s]: ", csv_fname);
				csv_fname = get_fname(csv_fname);
				save_sim_data(result, *p, csv_fname);
				break;
			case 6 :
				if (!sim_done) {
					free_array(result);
					printf("No up-to-date simultion found, starting it\n");
					result = full_sim(*theta1, *theta2, *c, *p);
					if (result == NULL) {
						printf("Failed to allocate memory for results.\n");
						break;
					}
					else
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
	
	triple **result = NULL;

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
				free_matrix(result);
				printf("Started simulation\n");
				result = flip_matrix(*p, *c);
				if (result == NULL) {
					printf("Failed to allocate momory for results.\n");
					sim_done = 0;
				}
				else
					sim_done = 1;
				break;
			case 3 :
				if (!sim_done) {
					free_matrix(result);
					printf("No up-to-date simulation found, starting it\n");
					result = flip_matrix(*p, *c);
					if (result == NULL) {
						printf("Failed to allocate momory for results.\n");
						break;
					}
					else
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
