#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFF_SIZE 4096

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
	printf("Offset: %lu\n", offset);
	size_t stripped_size = strlen(buff + offset);
	printf("Stripped size: %lu\n", stripped_size);
	if (stripped_size == 0)
		return prev; /* do nothing on empty input */
	free(prev);
	char *result = malloc((stripped_size + 1)*sizeof(char));
	strcpy(result, buff + offset);
	for (size_t i = 0; i <= strlen(result); ++i)
		printf("Character: %c\tCode: %d\n", result[i], (int)result[i]);
	return result;
}

/* Reads in integer choice from the user.
 * Heavily inspired by https://stackoverflow.com/a/58434556 */
int get_choice() {
	char buff[BUFF_SIZE];
	int choice = 0;

	if (fgets(buff, BUFF_SIZE, stdin) == NULL)
		return 0; /* Return 0 on error */

	if (sscanf(buff, "%d", &choice) == 1)
		return choice;

	return 0; /* return 0 on icorrect input format */
}
