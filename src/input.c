#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFF_SIZE 4096

typedef unsigned long int ulong;
typedef long double triple;

/* Strip out leading spaces from dynamic string. 
 * The input MUST be NULL-terminated, otherwise bad things happen.
 * May be extended to treat other whitespaces as well. */
size_t leading_spaces(char *string) {
	size_t offset = 0;
	/* find number of spaces */
	while (string[offset] == ' ')
		++offset;
	return offset;
}

char *get_fname(char *prev) {
	char buff[BUFF_SIZE];
	if (fgets(buff, BUFF_SIZE, stdin) == NULL)
		return prev; /* don't change anything on error */
	buff[strlen(buff) - 1] = '\0'; /* re-terminate at ending newline */
	size_t offset = leading_spaces(buff);
	size_t stripped_size = strlen(buff + offset);
	if (stripped_size == 0)
		return prev; /* do nothing on empty input */
	free(prev);
	char *result = malloc((stripped_size + 1)*sizeof(char));
	strcpy(result, buff + offset);
	return result;
}

/* Reads in integer choice from the user.
 * Heavily inspired by https://stackoverflow.com/a/58434556 */
ulong get_ulong(ulong def) {
	char buff[BUFF_SIZE];
	ulong choice = 0;

	if (fgets(buff, BUFF_SIZE, stdin) == NULL)
		return def; /* Return default on error */

	if (sscanf(buff, "%lu", &choice) == 1)
		return choice;

	return def; /* return default on incorrect input format */
}

triple get_triple(triple def) {
	char buff[BUFF_SIZE];
	triple choice;

	if (fgets(buff, BUFF_SIZE, stdin) == NULL)
		return def; /* Return default on error */

	if (sscanf(buff, "%Lf", &choice) == 1)
		return choice;

	return def; /* return default on incorrect input format */
}

int get_bool() {
	char buff[BUFF_SIZE];
	char choice;

	if (fgets(buff, BUFF_SIZE, stdin) == NULL)
		return 0; /* Return False on error */

	if (sscanf(buff, "%c", &choice) == 1)
		if (choice == 'y' || choice == 'Y')
			return 1;

	return 0; /* return False otherwise */
}
