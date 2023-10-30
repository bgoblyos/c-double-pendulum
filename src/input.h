/* function for determining the amount
 * of leading spaces in a NULL-terminated string */
size_t leading_spaces(char*);

/* Function for reading in a filename from the user,
 * stripping the input of leading spaces and writing it
 * to the supplied string if it's non-empty */
char *get_fname(char*);

/* Function for reading in a decimal choice from the user */
int get_choice(void);
