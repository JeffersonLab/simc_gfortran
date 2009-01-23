/* makereg.c version 1.0.  August 1994, Allen Boozer
   Report bugs to adb2y@virginia.edu
 $Log: makereg.c,v $
 Revision 1.1  2009/01/23 13:34:01  gaskelld
 Initial revision

 Revision 1.1.22.1  2008/09/25 00:54:05  jones
 Updated for running on Fedora 8 with gfortran

 Revision 1.2  2008/09/25 00:01:29  jones
 Updated to run with gfortran compiler

 Revision 1.1.24.1  2007/09/10 21:32:47  pcarter
 Implemented changes to allow compilation on RHEL 3,4,5 and MacOSX

 Revision 1.1  1998/12/07 22:11:11  saw
 Initial setup

 Revision 1.2  1995/01/09 15:10:44  saw
 Put titles in reg calls on a continuation line

 * Revision 1.1  1994/08/26  17:47:07  saw
 * Initial revision
 *
*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>

/*
   c%%                  - Put the current line in the fortran code
   c CTPTYPE = test     - Use "test" call statements
   c CTPTYPE = parm     - Use "parm" call statements
   c CTPTYPE = event    - Use "event" call statements
   c CTPTYPE = off      - Ignore all lines until CTYPE is set to test, parm, or event
*/


#define VERSION "v1.01"
#define BUFFER_LEN 256
#define NUM_TYPES 8

#define CTPTEST  0
#define CTPPARM  1
#define CTPEVENT 2
#define CTPOFF   3

#define COMMON             -1
#define PARAMETER          -2
#define EQUIV              -3
#define NOP                -4
#define COMMON_CONTINUE    -5
#define REGISTER_CONTINUE  -6
#define MARK               -7
#define IGNORE             -8
#define SKIP               -9

/*
   Two linked lists (the register list and the common list) are used to store
   information about variables that have been declared.  The register list
   stores variables which have been registered, and the common list stores
   variables which have been defined in common blocks.  The elements of the
   linked lists are of type "node", as defined below:
*/

struct node {
  int vartype;         /* A number which represents the type of the variable */
  int action;          /* A number which tells what to do with the variable */
  int calltype;        /* Use test, parm, or event calls */
  int line_number;     /* The line number on which the variable occurs */
  char *name;          /* The name of the variable */
  char *size;          /* The size of the array, or NULL if not an array */
  char *title;         /* The title string, or NULL if no title string */
  struct node *next;   /* Ptr to the next node, or NULL if last node */
};

struct node *register_start;  /* Ptr to the first node of the register list */
struct node *common_start;    /* Ptr to the first node of the common list */

/* Variable types (as they appear when variables are declared) */
char types[NUM_TYPES][20] = {
  "logical",          "logical*4",
  "integer",          "integer*4",
  "real",             "real*4", 
  "double precision", "real*8" };

/* Variable types (as the appear in fortran call statements) */
char type_names[5][10] = {"int", "int", "real", "double", "string"};

char keywords[3][15] = {"common", "parameter", "equivalence"};
char call_names[4][10] = {"test", "parm", "event", "off"};

int variable_flags[3][5] = {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}};
int array_flags[3][5] = {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}};

FILE *input, *output, *error;
char input_filename[BUFFER_LEN];
char output_filename[BUFFER_LEN];
char error_filename[BUFFER_LEN] = {""};
char subroutine_name[BUFFER_LEN] = {""};
char command_line[BUFFER_LEN];
int current_calltype, current_line;

void eprint (char string[]);
void eprintn (char string[], int length);
void eprint_newline ();
void eprint_line (int line_num);
void memory_error ();
char *parse_array_size (char string[]);
void create (struct node *ptr, int vartype, char string[], char title[]);
int max (int a, int b);
int min (int a, int b);
void free_node (struct node *ptr);
void clear_list (struct node *ptr);
void copy (char dest[], char source[], int length);
int strncmp_i (char str1[], char str2[], int length);
struct node *find_node (struct node *start, char string[]);
void mark_node (struct node *start, char string[], int new_action);
int determine_type (char string[]);
char *skip_blanks (char string[]);
char *skip_nonblanks (char string[]);
char *find_char (char string[], char character);
int extract_text (char string[]);
void shift_left (char string[]);
void parse (char string[]);
void compare_lists ();
void write_fortran_header ();
void write_fortran_code ();
void set_call_type (char string[]);

/***************************************************************************
* Linked list functions
***************************************************************************/

/* Print a character string to "error" */
void eprint (char string[]) {
  fprintf (error, "%s", string);
  fprintf (output, "%s", string);
}

/* Print "length" characters of a character string to "error" */
void eprintn (char string[], int length) {
  char	output_buffer[BUFFER_LEN];

  copy (output_buffer, string, length);
  fprintf (error, "%s", output_buffer);
  fprintf (output, "%s", output_buffer);
}

/* Print a newline character to "error" */
void eprint_newline () {
  fprintf (error, "\n");
  fprintf (output, "\n*     ");
}

/* Print a line number to "error" */
void eprint_line (int line_num) {
  fprintf (error, "Line %d: ", line_num);
  fprintf (output, "Line %d: ", line_num);
}

/* memory_error is called if malloc returns a NULL pointer */
void memory_error () {
  printf ("Memory allocation error\n");
  fclose (input);
  fclose (output);
  fclose (error);
  exit (0);
}

/* Return a pointer to a string containing the array size */
char *parse_array_size (char string[]) {
  char output_buffer[BUFFER_LEN];
  int length;
  char *size_string, *ptr = string;

  memset (output_buffer, '\0', BUFFER_LEN);
  ptr = find_char (ptr, '(');
  while (ptr < find_char(string,')')) {
    if (strlen(output_buffer) != 0) strcat (output_buffer, "*");
    length = min (find_char(ptr,',')-ptr, find_char(ptr,')')-ptr) - 1;
    strcat (output_buffer, "(");
    if (find_char(ptr,':')-ptr < length) {
      strcat (output_buffer, "1-");
      strncat (output_buffer, ptr, find_char(ptr,':')-ptr-1);
      strcat (output_buffer, "+");
      ptr = find_char (ptr, ':');
      length = min (find_char(ptr,',')-ptr, find_char(ptr,')')-ptr) - 1;
    }
    strncat (output_buffer, ptr, length);
    strcat (output_buffer, ")");
    ptr = skip_blanks (find_char (ptr, ','));
  }
  size_string = malloc (strlen(output_buffer)+1);
  if (size_string == NULL) memory_error ();
  strcpy (size_string, output_buffer);
  return (size_string);
}

/* Add a variable to a linked list */
void create (struct node *start, int vartype, char string[], char comment[]) {
  struct node *end = start, *temp;

  if (find_node (start, string) == NULL) {
    /* Create a new node and add it to the end of the linked list */
    while (end->next != NULL) end = end->next;
    end->next = malloc (sizeof(struct node));
    if (end->next == NULL) memory_error ();
    end = end->next;
    end->vartype = vartype;
    end->action = NOP;
    end->calltype = current_calltype;
    end->line_number = current_line;
    end->name = calloc (extract_text(string)+1, sizeof(char));
    if (end->name == NULL) memory_error ();
    copy (end->name, string, extract_text(string));
    if ((find_char (string, '(') < find_char (string, ',')) &&
	(find_char (string, '(') < find_char (string, '!'))) {
      /* Variable is an array */
      if (start == common_start) {
	/* Array is defined in a common block, so print a warning */
	eprint_newline ();
	eprint ("Warning - Array size defined in common block:");
	eprint_newline ();
	eprint_line (current_line);
	eprintn (string, find_char(string, ')') - string);
	eprint_newline ();
	temp = find_node (register_start, string);
	if (temp != NULL) {
	  if (temp->vartype > 0) array_flags[temp->calltype][temp->vartype]=1;
	  free (temp->size);
	  temp->size = parse_array_size (string);
	  end->size = NULL;
	}
      }
      else end->size = parse_array_size (string);
    }
    else (end->size = NULL);
    if (comment != NULL) {
      end->title = calloc (strlen(comment), sizeof(char));
      if (end->title == NULL) memory_error ();
      copy (end->title, comment, strlen(comment)-1);
    }
    else end->title = NULL;
    end->next = NULL;
  }
}

/* Return the maximum of a and b */
int max (int a, int b) {
  return ((a < b) ? b : a);
}

/* Return the minimum of a and b */
int min (int a, int b) {
  return ((a < b) ? a : b);
}

/* Release the memory used by a node */
void free_node (struct node *ptr) {
  free (ptr->name);
  if (ptr->size != NULL) free (ptr->size);
  if (ptr->title != NULL) free (ptr->title);
  free (ptr);
}

/* Release the memory used by each node in a linked list */
void clear_list (struct node* start) {
  struct node *ptr = start->next, *old_ptr;

  while (ptr != NULL) {
    old_ptr = ptr;
    ptr = ptr->next;
    free_node (old_ptr);
  }
  start->next = NULL;
}

/* Copy length chars from "dest" to "source", terminate "dest" with a \0 */
void copy (char dest[], char source[], int length) {
  memcpy (dest, source, length);
  dest[length] = '\0';
}

/* Case insensitive string compare */
int strncmp_i (char str1[], char str2[], int length) {
  int i;

  for (i=0; i<length;i++) if (toupper(str1[i]) != toupper(str2[i])) return (1);
  return (0);
}

/* Return a pointer to the node which has "string" as it's name field */
struct node *find_node (struct node *start, char string[]) {
  int length, flag=1;
  struct node *ptr = start;

  do {
    ptr = ptr->next;
    if (ptr != NULL) {
      length = max (extract_text(string), strlen(ptr->name));
      flag = strncmp_i (ptr->name, string, length);
    }
  } while ((flag != 0) && (ptr != NULL));
  return (ptr);
}

/* Set the action of the node which has "string" as it's name field */
void mark_node (struct node *start, char string[], int new_action) {
  struct node *node_ptr;

  node_ptr = find_node (start, string);
  if (node_ptr != NULL) node_ptr->action = new_action;
}

/***************************************************************************
* Parse functions
***************************************************************************/

int determine_type (char string[]) {
  int i, length;

  for (i=0; i<NUM_TYPES; i++) {
    length = max (strlen(types[i]), skip_nonblanks(string)-string);
    if (strncmp_i(string, types[i], length) == 0) return (i/2);
  }
  if (strncmp_i(string, "character", 9) == 0) return (4);
  if (strncmp_i(string, "integer*2", 9) == 0) {
    eprint_newline ();
    eprint ("Warning - Type integer*2:");
    eprint_newline ();
    eprint_line (current_line);
    eprintn (string, strlen(string) - 1);
    eprint_newline ();
    return (NOP);
  }
  for (i=0; i<3; i++) {
    length = max (strlen(keywords[i]), extract_text(string));
    if (strncmp_i(string, keywords[i], length) == 0) return (-i-1);
  }
  return (NOP);
}

/* Return a ptr to the first nonblank character in the string */
char *skip_blanks (char string[]) {
  while (isspace (*string)) string++;
  return (string);
}

/* Return a ptr to the first blank character in the string */
char *skip_nonblanks (char string[]) {
  while ((! isspace (*string)) && (*string != '\0')) string++;
  return (string);
}

/* Return a ptr to the character after the first occurrence of "character" */
char *find_char (char string[], char character) {
  while ((*string != character) && (*string != '\0')) string++;
  if (*string == character) string++;
  return (string);
}

/* Return the number of contiguous text characters in a string */
int extract_text (char string[]) {
  int i=0;

  while ((isalnum (string[i])) || (string[i] == '_')) i++;
  return (i);
}

/* Shift a character string one character to the right */
void shift_left (char string[]) {
  char *ptr = string;

  while (*ptr != '\0') ptr++;
  *(ptr+1) = '\0';
  while (ptr != string) *ptr = *(--ptr);
}

/* Parse one line of text */
void parse (char string[]) {
  struct node *list_ptr /*, *node_ptr */;
  static int vartype = NOP, state = NOP;
  char *ptr = string, *title_ptr, *temp;

  if (((string[0] == ' ') && (strlen(string) > 6)) ||
      ((string[0] == '\t') && (strlen(string) > 2))) {
    /* The current line is not a comment */
    if (string[0] == '\t') ptr++;
    else ptr += 5;
    if ((isspace (*ptr)) || ((state != COMMON_CONTINUE) && (state != REGISTER_CONTINUE))) {
      ptr = skip_blanks (ptr);
      vartype = determine_type (ptr);
      state = vartype;
      if (strncmp_i("double precision", ptr, 16) == 0)
	ptr = skip_blanks (skip_nonblanks (ptr));
    }
    if ((vartype != NOP) && (vartype != PARAMETER) && (vartype != EQUIV)) {
      /* The current line is a list of variables */
      if (state == COMMON_CONTINUE) {
	/* The current line is the continuation of a common block */
	list_ptr = common_start;
	ptr++;
      }
      else if (state == COMMON) {
	/* The current line is a common block */
	list_ptr = common_start;
	ptr = find_char (find_char (ptr, '/'), '/');
	state = COMMON_CONTINUE;
      }
      else if (state == REGISTER_CONTINUE) {
	/* The current line is the continuation of a register statement */
	list_ptr = register_start;
	ptr++;
      }
      else {
	/* The current line contains variables to be registered */
	list_ptr = register_start;
	ptr = skip_nonblanks (ptr);
	state = REGISTER_CONTINUE;
      }
      ptr = skip_blanks (ptr);

      /* Obtain the title string */
      title_ptr = string;
      if (find_char (ptr, ',') > find_char (ptr, '!')) {
	while (*title_ptr != '!') title_ptr++;
	title_ptr = skip_blanks (title_ptr+1);
	if (*title_ptr == '\0') title_ptr = NULL;
      }
      else title_ptr = NULL;
      /* Convert ' to '' in title string */
      if (title_ptr != NULL) {
	temp = title_ptr;
	while (*temp != '\0') {
	  if (*temp == '\'') shift_left (temp++);
	  temp++;
	}
      }

      /* Add each variable to the linked list */
      while ((ptr < find_char (string, '!')) && (*ptr != '!')) {
	if ((find_char (ptr, '(') < find_char (ptr, ',')) &&
	    (find_char (ptr, '(') < find_char (ptr, '!'))) {
	  /* The variable is an array */
	  create (list_ptr, vartype, ptr, title_ptr);
	  if (vartype >= 0) array_flags[current_calltype][vartype] = 1;
	  ptr = skip_blanks (find_char (find_char (ptr, ')'), ','));
	}
	else {
	  /* The variable is not an array */
	  create (list_ptr, vartype, ptr, title_ptr);
	  if (vartype >=0) variable_flags[current_calltype][vartype] = 1;
	  ptr = skip_blanks (find_char (ptr, ','));
	}
      }
    }
    if (vartype == PARAMETER)
      /* If the line is a parameter statement, then ignore the variable */
      mark_node (register_start, skip_blanks(find_char (ptr, '(')), IGNORE);
    if (vartype == EQUIV) {
      /* If the line is an equivalence statement, then skip the variables */
      ptr = skip_blanks (find_char (ptr, '('));
      mark_node (register_start, ptr, SKIP);
      if (find_char(ptr,'(') < find_char(ptr,',')) ptr = find_char(ptr,')');
      ptr = skip_blanks (find_char (ptr, ','));
      mark_node (register_start, ptr, SKIP);
    }
  }
}

/* Compare the common list with the register list */
void compare_lists () {
  struct node *ptr, *common_node;
  
  eprint_newline ();
  eprint ("Registered, but did not occur in a common block:");
  eprint_newline ();
  ptr = register_start->next;
  while (ptr != NULL) {
    common_node = find_node (common_start, ptr->name);
    if (common_node != NULL) {
      common_node->action = MARK;
      if ((common_node->title != NULL) && (ptr->title == NULL)) {
	ptr->title = malloc (strlen(common_node->title)+1);
	if (ptr->title == NULL) memory_error ();
	strcpy (ptr->title, common_node->title);
      }
    }
    else if ((ptr->action != IGNORE) && (ptr->action != SKIP)) {
      eprint_line (ptr->line_number);
      eprint (ptr->name);
      eprint_newline ();
    }
    ptr = ptr->next;
  }

  eprint_newline ();
  eprint ("Occurred in a common block, but were not registered:");
  ptr = common_start->next;
  while (ptr != NULL) {
    if (ptr->action != MARK) {
      eprint_newline ();
      eprint_line (ptr->line_number);
      eprint (ptr->name);
    }
    ptr = ptr->next;
  }
  eprint ("\n\n");
}

/* Write a header to the output file */
void write_fortran_header () {
  time_t current_time = time (NULL);

  fprintf (output, "******************************************************");
  fprintf (output, "*************************\n");
  fprintf (output, "*     This file (%s) was generated ", output_filename);
  fprintf (output, "from %s by makereg %s\n", input_filename, VERSION);
  fprintf (output, "*     This file was created on ");
  fprintf (output, "%s", asctime (localtime (&current_time)));
  fprintf (output, "*\n");
  fprintf (output, "*     The command used to create this file was:\n");
  fprintf (output, "*     %s\n", command_line);
  fprintf (output, "*\n");
  fprintf (output, "*     Do not edit this file.\n");
  fprintf (output, "******************************************************");
  fprintf (output, "*************************\n\n");
  fprintf (output, "      subroutine %s\n\n", subroutine_name);
  fprintf (output, "      implicit none\n\n");
}

/* Write the fortran code to an output file */
void write_fortran_code () {
  struct node *ptr = register_start->next;
  int i, j;

  for (j=0; j<3; j++) {
    if (variable_flags[j][0] == 1) variable_flags[j][1] = 1;
    if (array_flags[j][0] == 1) array_flags[j][1] = 1;
  }

  fprintf (output, "      include '%s'\n\n", input_filename);
  for (j=0; j<3; j++) {
    for (i=1; i<5; i++) {
      if (variable_flags[j][i] == 1) {
	fprintf (output, "c      integer ");
	fprintf (output, "reg%s%s\n", call_names[j], type_names[i]);
	fprintf (output, "c      external ");
	fprintf (output, "reg%s%s\n", call_names[j], type_names[i]);
      }
      if (array_flags[j][i] == 1) {
	fprintf (output, "c      integer ");
	fprintf (output, "reg%s%sarray\n", call_names[j], type_names[i]);
	fprintf (output, "c      external ");
	fprintf (output, "reg%s%sarray\n", call_names[j], type_names[i]);
      }
    }
  }
  fprintf (output, "\n");
  
  /* Loop to output the reg calls */
  while (ptr != NULL) {
    if (ptr->action != IGNORE) {
      fprintf (output, "      call reg%s", call_names[ptr->calltype]);
      fprintf (output, "%s", type_names[ptr->vartype]);
      if (ptr->size == NULL)
	fprintf (output, "('%s',%s,", ptr->name, ptr->name);
      else
	fprintf (output, "array('%s',%s,%s,", ptr->name, ptr->name, ptr->size);
      if (ptr->title == NULL) fprintf (output, "0)\n");
      else {
	fprintf (output, "\n     & ");
	fprintf (output, "'%s')\n", ptr->title);
      }
    }
    ptr = ptr->next;
  }

  fprintf (output, "\n");
  fprintf (output, "      return\n");
  fprintf (output, "      end\n");
}

/* Set the call type to "test", "parm", or "event" */
void set_call_type (char string[]) {
  int i;
  char *ptr;

  if (strlen (string) != 0) {
    ptr = skip_blanks (string+1);
    if (strncmp_i(ptr, "CTPTYPE", 7) == 0) {
      while (isalpha(*ptr) && (*ptr != '\0')) ptr++;
      while ((! isalpha(*ptr)) && (*ptr != '\0')) ptr++;
      if (*ptr != '\0')
	for (i=0; i<4; i++)
	  if (strncmp_i(ptr, call_names[i], 3) == 0) current_calltype = i;
    }
  }
}

void print_usage () {
    printf ("Usage:  makereg infile [-o outfile] [-e errorfile] ");
    printf ("[-s subroutine name]\n");
    printf ("                       [-c test | parm | event]\n");
    exit (0);
}

int main (int argc, char *argv[]) {
  char buffer[BUFFER_LEN];
  struct node first_register_node, first_common_node;
  int i, j;

  first_register_node.next = NULL;
  first_common_node.next = NULL;
  register_start = &first_register_node;
  common_start = &first_common_node;

  for (i=0; i<argc; i++) {
    strcat (command_line, argv[i]);
    strcat (command_line, " ");
  }

  current_calltype = CTPTEST;
  error = stderr;
  if (argc < 2) print_usage ();
  strcpy (input_filename, argv[1]);
  strcpy (output_filename, input_filename);
  if ((strcmp (strrchr(output_filename, '.'), ".cmn")) == 0)
    strcpy (strrchr(output_filename, '.'), ".f");
  else strcat (output_filename, ".f");
  i = 2;
  while (i < argc) {
    if (strcmp(argv[i], "-o") == 0) {
      if (argc > i+1) strcpy (output_filename, argv[i+1]);
      else print_usage ();
    }
    else if (strcmp(argv[i], "-e") == 0) {
      if (argc > i+1) strcpy (error_filename, argv[i+1]);
      else print_usage ();
    }
    else if (strcmp(argv[i], "-c") == 0) {
      if (argc > i+1) {
	for (j=0; j<3; j++)
	  if (strcmp(argv[i+1], call_names[j]) == 0) current_calltype = j;
      }
      else print_usage ();
    }
    else if (strcmp(argv[i], "-s") == 0) {
      if (argc > i+1) strcpy (subroutine_name, argv[i+1]);
      else print_usage ();
    }
    i += 2;
  }
  input = fopen (input_filename, "r");
  if (input == NULL) {
    printf ("Invalid filename: %s\n", input_filename);
    print_usage ();
  }
  output = fopen (output_filename, "w");
  if (output == NULL) {
    printf ("Invalid filename: %s\n", output_filename);
    print_usage ();
  }
  if (strlen(error_filename) != 0) {
    error = fopen (error_filename, "w");
    if (error == NULL) {
      printf ("Invalid filename: %s\n", error_filename);
      print_usage ();
    }
  }
  if (strlen(subroutine_name) == 0) {
    strcpy (subroutine_name, output_filename);
    if (strrchr(subroutine_name, '.') != NULL)
      *strrchr(subroutine_name, '.') = '\0';
  }

  write_fortran_header ();
  current_line = 1;
  while (fgets (buffer, BUFFER_LEN, input) != NULL) {
    set_call_type (buffer);
    if (strncmp(buffer, "*%%", 3) == 0)
      fprintf (output, "      %s", skip_blanks(skip_nonblanks(buffer)));
    else if (current_calltype != CTPOFF) parse (buffer);
    current_line++;
  }
  compare_lists ();
  write_fortran_code ();
  clear_list (register_start);
  clear_list (common_start);
  fclose (input);
  fclose (output);
  fclose (error);
  return 0;
}
