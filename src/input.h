/* function for determining the amount
 * of leading spaces in a NULL-terminated string */
size_t leading_spaces(char*);

/* Function for reading in a filename from the user,
 * stripping the input of leading spaces and writing it
 * to the supplied string if it's non-empty */
char *get_fname(char*);

/* Function for reading in an unsigned long int from the user */
unsigned long int get_ulong(unsigned int);

/* Function for reading in a long double from the user */
long double get_triple(long double);

/* Function for reading in a yes/no decision from the user */
int get_bool(void);
