/*-----------------------------------------------------------------------------
 * Copyright (c) 1993 Southeastern Universities Research Association,
 *                    Continuous Electron Beam Accelerator Facility
 *
 * This software was developed under a United States Government license
 * described in the NOTICE file included as part of this distribution.
 *
 * Stephen A. Wood, 12000 Jefferson Ave., Newport News, VA 23606
 * Email: saw@cebaf.gov  Tel: (804) 249-7367  Fax: (804) 249-5800
 *-----------------------------------------------------------------------------
 * 
 * Description:
 *  Utilities used by CTP
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: thUtils.c,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.3.22.1  2008/09/25 00:54:05  jones
 *   Updated for running on Fedora 8 with gfortran
 *
 *   Revision 1.4  2008/09/25 00:01:29  jones
 *   Updated to run with gfortran compiler
 *
 *   Revision 1.3.24.1  2007/09/10 21:32:47  pcarter
 *   Implemented changes to allow compilation on RHEL 3,4,5 and MacOSX
 *
 *   Revision 1.3  2003/02/21 20:55:25  saw
 *   Clean up some types and casts to reduce compiler warnings.
 *
 *   Revision 1.2  1999/11/04 20:34:07  saw
 *   Alpha compatibility.
 *   New RPC call needed for root event display.
 *   Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *   Revision 1.1  1998/12/07 22:11:14  saw
 *   Initial setup
 *
 *	  Revision 1.6  1995/04/10  15:50:27  saw
 *	  Add thSpecial to interpret directives such as #real, #int, ...
 *
 *	  Revision 1.5  1995/01/09  15:38:18  saw
 *	  Fix up potential ref to unallocated memory in thIDToken.
 *
 *	  Revision 1.4  1994/09/27  19:44:50  saw
 *	  Add fnmatch routine from BSD.  Only use when OS doesn't have it (ultrix)
 *
 *	  Revision 1.3  1994/07/21  20:40:22  saw
 *	  In thCommas, ignore commas in quotes.  Replace stderr with STDERR.
 *
 *	  Revision 1.2  1993/09/22  15:24:50  saw
 *	  Allow thIDToken to accept hex
 *
 *	  Revision 1.1  1993/05/11  18:02:05  saw
 *	  Initial revision
 *
 */

/* thUtils.c
   Routines used by thTest, thHist and thParm

*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "daVar.h"
#include "th.h"
#include "thInternal.h"
#include "thUtils.h"

daVarStatus thGetIndex(char *name, int *index, char **pptr)
/* If the name has an array index, evaluate the index and return it to
*index.  Also return a pointer to ( or [ that starts the index.
If there is no index, return a zero index, but also return the status
code S_DAVAR_NOINDEX.

I have a whole bunch of ways of dealing with these  indices.  Must be
cleaned up.
*/
{
  char *nend, *pbegin;			/* Name end */
  char *t;
  char cright;
  int istyle;

  pbegin = nend = name + strlen(name);
  if((t=strchr(name,'('))) if(t<pbegin) {
    pbegin = t;
    cright = ')';
    istyle = -1;
  }
  if((t=strchr(name,'['))) if(t<pbegin) {
    pbegin = t;  
    cright = ']';
    istyle = 0;
  }
  *pptr = pbegin;
  if(pbegin < nend) {
    char *s,*rp;
    int i;

    s = pbegin + 1;
    if((rp = strrchr(s,cright)) == 0) {
      *index = 0;
      return(S_FAILURE);
    }
    *rp = 0;
    if(thEvalImed(s,0,&i) != S_SUCCESS) {
      *index = 0;
      return(S_FAILURE);
    }
    *rp = cright;
    *index = i + istyle;
    return(S_SUCCESS);
  }
  *index = 0;
  return(S_DAVAR_NOINDEX);
}

int thComSpace(char *s,char **args)
/* Treat spaces, tabs and commas as delimiters. */
{
  char *p;
  int i;
  int comma;

  i = 0;
  p = s;
  comma = 1;
  while(*p != '\0'){
    while(isspace(*p) || (!comma && *p == ',')) {
      if(!comma) comma = (*p == ',');
      p++;
    }
    if(*p == '\0') break;
    args[i++] = p;
    while(!isspace(*p) && *p != ',' && *p != '\0') p++;
    if(*p != '\0'){
      comma = (*p == ',');
      *p++ = '\0';
    }
  }
  return(i);
}
int thCommas(char *s,char **args){
  /* commas inside of quotes are now protected */
  char *p;
  int i;
  char quotechar;
  int instring;

  args[0] = s;
  i = 1;
  p = s;
  instring = 0;
  while(*p != 0) {
    if(instring && *p == quotechar) {
      if(*(p+1) == quotechar) p++;
      else instring = 0;
    } else {
      if(*p == ','){
	*p = '\0';
	args[i++] = p+1;
      } else if (*p == QUOTECHAR1 || *p == QUOTECHAR2) {
	instring = 1;
	quotechar = *p;
      }
    }
    p++;
  }

/*  while((p=strchr(p,',')) != NULL){
    *p++ = '\0';
    args[i++] = p;
  }*/
  args[i] = 0;
  return(i);
}
char *thTokenArray(char *s,int *index)
/* Interpret the string s as an array element.  Looks for an index inside
of []'s, returning that in index.  If there is an index inside of ()'s, then
one is subtracted from the index before it is returnd.  []'s are for
c style indexing, ()'s for fortran style indexing.  A pointer to the [ or (
is returned, so that the variable name may be null terminated by the caller.
If there is an error in the balancing of the []'s or ()'s, a null is returned
to signify an error.
*/
{
  int cstyle;
  char *leftp,*rightp;
  char sindex[25];
  char cright;
  char wleft,wright;		/* Character for "other" style */
  int len;

  if((leftp=strchr(s,'('))){
    cright = ')';
    wleft = '['; wright = ']';
    cstyle = 1;
  } else if((leftp=strchr(s,'['))){
    cright = ']';
    wleft = '('; wright = ')';
    cstyle = 0;
  } else return(NULL);
  rightp = strchr(s,cright);
  if((leftp>=rightp) || strchr(s,wleft) || strchr(s,wright)) return (NULL);
  len = min(24,(rightp-leftp)-1);
  strncpy(sindex,leftp+1,len);
  sindex[len] = 0;
  *index = atol(sindex);
  if(cstyle) (*index)--;
  return(leftp);
}

thTokenType thIDToken(char *s)
/* Examines a string to determine if it is an integer, real number,
variable name, or array element.  Only works if is not an expression.
Allows hex constants as INT's */
{
  char *p;
  thTokenType typ;
  int nume;			/* Number of E characters */

  /*  printf("TYPE(%s)=",s);*/
  nume = 0;
  p = s;
  typ = TOKINT;
  while(*p != '\0'){
    if(strchr("()[]",*p)) {/*  printf("%d\n",TOKARRAY); */ return(TOKARRAY);}
    if(typ != TOKVAR){
      switch(*p)
	{
	case 'e':
	case 'E':
	  if(nume > 0 || p==s){
	    typ = TOKVAR;
	    break;
	  }
	  nume = 1;
	case '.':
	  typ = TOKFLOAT;
	  break;
	case '+':
	case '-':
	  break;
	case 'x': case 'X':
	  if(p==s || *(p-1) != '0' || nume > 0)
	    typ = TOKVAR;
	  else
	    nume = 1;
	default:
	  if(!isdigit(*p)) typ = TOKVAR;
	  break;
	}
    }
    p++;
  }
  /*printf("%d\n",typ);*/
  return(typ);
}

char *thSpaceStrip(char *s)
/* Strip leading, trailing and embedded spaces from a string.
   Modifies argument.*/
{
  char *p,*t;

  p = t = s;
  while(*s != '\0'){
    if(*s != ' ') *t++ = *s++;
    else *s++;
  }
  *t = '\0';
  return(p);

}

int thSpecial(char *line, char *default_class)
{				/* Process special commands */
  char *s,*p;
  char *command;
  char *arg;
  char *arglist[20];
  int nargs;
  char *class;
  daVarStruct *varp;
  int vartype;
  char *classlist[]={0,0};
  int i;

  s = line;
  while(isspace(*s) && *s) s++;
  if(*s != SPECIALCHAR) return(0); /* Not a special command */

  /* Split line into command and argument */
  s++;
  while(isspace(*s) && *s) s++;	/* Skip to command */
  if(!*s) return(1);		/* Empty line */
  /* s now points to the command */
  p = s+1;
  while(!isspace(*p) && *p) p++; /* Skip to end of command */
  if(!*p) return(1);		/* No argument */
  command = malloc(p-s+1);
  strncpy(command,s,p-s); command[p-s] = '\0';
  s = p;
  while(isspace(*s) && *s) s++;	/* Skip to argument */
  if(!*s) {			/* No argument */
    free(command);
    return(1);
  }
  arg = malloc(strlen(s)+1);
  strcpy(arg,s);
  nargs = thCommas(arg,arglist);
  
/* Now need to look for a . in command to see if there is a class specified */
  s = strchr(command,'.');
  if(s) {
    *s = 0;
    class = s+1;
    s = class;
    while(*s) {tolower(*s); s++;}
  } else {
    class = default_class;
  } /* Should probably use a table here */
  if(strncasecmp(command,"integer",strlen(command))==0) {
    vartype = DAVARINT;
  } else if(strncasecmp(command,"real",strlen(command))==0) {
    vartype = DAVARFLOAT;
  } else if(strncasecmp(command,"double",strlen(command))==0) {
    vartype = DAVARDOUBLE;
  } else if(strncasecmp(command,"string",strlen(command))==0) {
    vartype = DAVARSTRING;
  }
  classlist[0] = class;
  for(i=0;i<nargs;i++) {
    arglist[i] = thSpaceStrip(arglist[i]);
    thVarCreate(arglist[i],vartype,classlist,&varp);
  }
/* command and arg now contain what to do */
  free(command);
  free(arg);
  return(1);
}
thVarCreate(char *s, int vartype, char **classlist, daVarStruct **varpp)
/* Should eventually merge this with thVarResolve */
{
  int cstyle;
  char cleft,cright;
  char *leftp,*rightp;
  daVarStruct var;
  int arindex,arsize;

  *varpp = 0;
  arindex = 0;

  leftp = strchr(s,'(');
  {
    char *lb;
    lb = strchr(s,'[');
    if(leftp) {
      if(lb && lb<leftp) leftp = lb;
    } else
      leftp = lb;
  }

  if(leftp) {
    cleft = *leftp;
    *leftp = '\0';
  } 
  if(leftp){
    int cstyle;
    char *sindex;
    int indtemp;

    sindex = leftp + 1;
    
    if(cleft=='('){
      cstyle = -1;
      cright = ')';
    } else {
      cstyle = 0;
      cright = ']';
    }
    if((rightp=strrchr(sindex,cright)) == 0){
      fprintf(stderr,"Syntax error in %s\n",s);
      return(S_FAILURE);
    }
    *rightp = 0;
    if(thEvalImed(sindex,0,&indtemp)!= S_SUCCESS){
      fprintf(stderr,"Error evaluating index %s\n",sindex);
      *rightp = cright;
      return(S_FAILURE);
    }
    *rightp = cright;
    arindex = indtemp + cstyle;
    arsize = indtemp;
  } else {
    arsize = 1;
  }
  if(daVarLookupPWithClass(s,classlist,varpp) != S_SUCCESS){
    /* Doesn't exist, we can create it as we want it */
    if(strchr(s,'.')) {	/* Don't prepend a class */
      var.name = malloc(strlen(s)+1);
      strcpy(var.name,s);
    } else {
      var.name = malloc(strlen(classlist[0])
				 +strlen(s)+2);
      strcpy(var.name,classlist[0]);
      strcat(var.name,".");
      strcat(var.name,s);
    }
    if(leftp) *leftp = cleft;	/* Restore left paren */
    var.size = arsize;
    var.type = vartype;
    switch(vartype)
      {
      case DAVARINT:		/* How should we initialize the variables */
	var.varptr = malloc(var.size*sizeof(DAINT));
	break;
      case DAVARFLOAT:
	var.varptr = malloc(var.size*sizeof(DAFLOAT));
	break;
      case DAVARDOUBLE:
	var.varptr = malloc(var.size*sizeof(DADOUBLE));
	break;
      case DAVARSTRING:
	var.varptr = malloc(var.size*sizeof(char *));
	*((char *)var.varptr) = '\0';
	break;
      }
    var.opaque = 0;
    var.rhook = 0;
    var.whook = 0;
    var.flag = DAVAR_REPOINTOK | DAVAR_READONLY | DAVAR_DYNAMIC_PAR;
    var.flag = DAVAR_REPOINTOK | DAVAR_READONLY | DAVAR_DYNAMIC_PAR;
    var.title = 0;
/*    printf("Registering %s(%d) at %x\n",var.name,var.size,var.varptr);*/
    daVarRegister((int) 0,&var); /* Create the parameter */
    daVarLookupP(var.name,varpp);
    free(var.name);
  } else {
    if(leftp) *leftp = cleft;	/* Restore left paren */
    /* Already exists */
    if((*varpp)->type == vartype && (*varpp)->size == arsize) {
      /* We are OK, return */
      return(S_SUCCESS);
    } else {
      /* See if we can fix */
      if((*varpp)->flag&DAVAR_DYNAMIC_PAR) {
	if((*varpp)->size != arsize || (*varpp)->type != vartype) {
	  (*varpp)->type = vartype;
	  (*varpp)->size = arsize;
	  switch((*varpp)->type)
	    {
	    case DAVARINT:
	      (*varpp)->varptr = (void *) realloc((*varpp)->varptr,arsize * sizeof(DAINT));
	      break;
	    case DAVARFLOAT:
	      (*varpp)->varptr = (void *) realloc((*varpp)->varptr,arsize * sizeof(DAFLOAT));
	      break;
	    case DAVARDOUBLE:
	      (*varpp)->varptr = (void *) realloc((*varpp)->varptr,arsize * sizeof(DADOUBLE));
	    case DAVARSTRING:
	      (*varpp)->varptr = (void *) realloc((*varpp)->varptr,arsize * sizeof(char *));
	      *((char *)var.varptr) = '\0';
	      break;
	    }
	}
      }
    }
  }
  return(S_SUCCESS);
}
	    


#ifdef OLD
int thGetMode(char *s, enum MODES mode, enum MODES *nmode, char **blocknamep)
{
  char *args[20];
  int nargs;
  enum MODES omod e;
  static char *nulstr="";
  char *command,*last;

  omode = mode;
  *nmode = mode;
  *blocknamep = 0;

  command = s;
  while(isspace(*command)) command++;
  last = command;
  while(!isspace(*last) && *last) last++;
  if(strncasecmp(command,"begin",last-command) == 0
     || strncasecmp(command,"end",last-command) == 0){
    nargs = thComSpace(s,args);
    if(nargs < 2) {
      *nmode = M_BAD;
      return(1);
    }
    if(strcasecmp("test",args[1]) == 0)
      *nmode = M_TES;
    else if(strcasecmp("hist",args[1]) == 0
	    || strcasecmp("histogram",args[1]) == 0)
      *nmode = M_HIS;
    else if(strcasecmp("parameter",args[1]) == 0)
      *nmode = M_PAR;
    else
      *nmode = M_BAD;
    if(strcasecmp(args[0],"end") == 0){
      if(*nmode == omode)
	*nmode = M_COM;
      else
	*nmode = M_BAD;
    } else {
      if(nargs > 2)
	*blocknamep = args[2];
      else
	*blocknamep = nulstr;
    }
    return(1);
  } else
    return(0);
}

int thTokToPtr(char *token, int create, int intonly, daVarStruct *varp)
/* Find registered variable pointer for token.  If the token is an array,
   create a new entry if the array is already registered.  If the variable
   or array is not registered, register it only if create is true.  If
   intonly is true, fail if an existing variable is not integer.

   If intonly is true when the token is an element of an array that is
   not short or long, then the routine returns an error.

   If the token is a constant, then create it if it doesn't exist */
{
  thTokenType toktyp;
  daVarStruct var;
  int lookstat;

  toktyp = thIDToken(token);
  if(intonly&&create)
    if(toktyp != TOKVAR && toktyp != TOKARRAY){
      fprintf(STDERR,"Variable %s must not be a number\n",token);
      return(S_FAILURE);
    }
  
  lookstat = daVarLookup(token, &var);
  if(lookstat != S_SUCCESS)
    lookstat = thVarLookup(token, &var);
  if(lookstat != S_SUCCESS){
    if(create || toktyp != TOKVAR){
      var.name = token;
      var.title = 0;
      var.flag = DAVAR_READONLY;
      var.rhook = 0;
      var.whook = 0;
      var.opaque = 0;
      switch(toktyp)
	{
	case TOKINT:
	  {
	    register DAINT *longp;
	    longp = malloc(sizeof(DAINT));
	    *longp = atol(token);
	    var.varptr = longp;
	    var.type = DAVARINT;
	    var.size = 1;
/*	    daVarRegister(token,longp,DAVARINT,1,DAVAR_READONLY);*/
	    thVarRegister((int) 0, &var);
	  }
	  break;
	case TOKFLOAT:
	  {
	    register DAFLOAT *floatp;
	    floatp = malloc(sizeof(DAFLOAT));
	    *floatp = atof(token);
	    var.varptr = floatp;
	    var.type = DAVARFLOAT;
	    var.size = 1;
/*	    daVarRegister(token,floatp,DAVARFLOAT,1,DAVAR_READONLY);*/
	    thVarRegister((int) 0, &var);
	  }
	  break;
	case TOKVAR:
	  {
	    register DAINT *longp;
	    longp = malloc(sizeof(DAINT));
	    var.varptr = longp;
	    var.type = DAVARINT;
	    var.size = 1;
/*	    daVarRegister(token,intp,DAVARINT,1,DAVAR_READWRITE);*/
	    thVarRegister((int) 0, &var);
	  }
	  break;
	case TOKARRAY:
	  {
	    /* Only create if underlying variable exists. */
	    char *p; int index; char leftp;
	    p = thTokenArray(token,&index);
	    leftp = *p;		/* Save ( or [ character */
	    *p = 0;
	    lookstat = daVarLookup(token, &var);
	    if(lookstat != S_SUCCESS)
	      lookstat = thVarLookup(token, &var);
	    if(lookstat == S_SUCCESS && intonly && (var.type != DAVARINT)) {
	      fprintf(STDERR,
		      "Array %s must be a preregistered LONGINT array.\n",token);
	      return(S_FAILURE);
	    } else {
	      if(index >= var.size){
		fprintf(STDERR,
			"Index %d exceeds %s array length of %d\n",index,
			token,var.size);
		return(S_FAILURE);
	      }
	    }
	    *p = leftp;
	    switch(var.type)
	      {
	      case DAVARINT:
		var.varptr = (int *)(var.varptr) + index;
		var.type = DAVARINT;
		thVarRegister((int) 0, &var);
		break;
	      case DAVARFLOAT:
		var.varptr = (float *)(var.varptr) + index;
		var.type = DAVARFLOAT;
		thVarRegister((int) 0, &var);
		break;
	      default:
		*p = 0;
		fprintf(STDERR,"Illegal type %d for variable %s\n"
			,var.type,token);
		return(S_FAILURE);
	      }
	    break;
	  }
	}
      lookstat = daVarLookup(token, &var);
      if(lookstat != S_SUCCESS)
	if((lookstat = thVarLookup(token, &var)) != S_SUCCESS){
	  fprintf(STDERR,"(thToktoPtr) Bad error\n");
	  return(S_FAILURE);
	}
    } else {
      fprintf(STDERR,"Variable %s must be preregistered\n",token);
      return(S_FAILURE);
    }
  }
  daVarStructCopy(varp,&var);
  return(S_SUCCESS);
}
#endif
int thCleanLine(char *s)
{
  int blank;
  blank = 1;
  while(*s != 0){
    if(isspace(*s)) *s = ' ';	/* Remove tabs, ... */
    else if(*s == COMCHAR) {
      *s = 0;
      break;
    } else
      blank = 0;
    s++;
  }
  return(blank);
}
/*float argtoFloat(daVarStruct *x)
{
  float d;

  switch(x->type)
    {
    case DAVARINT:
      d = *(DAINT *) x->varptr;
      break;
    case DAVARFLOAT:
      d = *(DAFLOAT *) x->varptr;
      break;
    case DAVARINTP:
      d = **(DAINT **) x->varptr;
      break;
    case DAVARFLOATP:
      d = **(DAFLOAT **) x->varptr;
      break;
    }
  return(d);
}*/
int argtoInt(daVarStruct *x)
{
  DAINT l;
  DAFLOAT d;

  switch(x->type)
    {
    case DAVARINT:
      l = *(DAINT *) x->varptr;
      break;
    case DAVARFLOAT:
      d = *(DAFLOAT *) x->varptr;
      l = floatToLong(d);
      break;
    case DAVARINTP:
      l = **(DAINT **) x->varptr;
      break;
    case DAVARFLOATP:
      d = **(DAFLOAT **) x->varptr;
      l = floatToLong(d);
      break;
    }
  return(l);
}

    
/*
	strstr - public-domain implementation of standard C library function

	last edit:	02-Sep-1990	D A Gwyn

	This is an original implementation based on an idea by D M Sunday,
	essentially the "quick search" algorithm described in CACM V33 N8.
	Unlike Sunday's implementation, this one does not wander past the
	ends of the strings (which can cause malfunctions under certain
	circumstances), nor does it require the length of the searched
	text to be determined in advance.  There are numerous other subtle
	improvements too.  The code is intended to be fully portable, but in
	environments that do not conform to the C standard, you should check
	the sections below marked "configure as required".  There are also
	a few compilation options, as follows:

	#define ROBUST	to obtain sane behavior when invoked with a null
			pointer argument, at a miniscule cost in speed
	#define ZAP	to use memset() to zero the shift[] array; this may
			be faster in some implementations, but could fail on
			unusual architectures
	#define DEBUG	to enable assertions (bug detection)
	#define TEST	to enable the test program attached at the end
*/
#define ROBUST
#if !defined(__osf__) || !defined(__alpha)
#define	ZAP
#endif

#include	<stddef.h>		/* defines size_t and NULL */
#include	<limits.h>		/* defines UCHAR_MAX */

#ifdef ZAP
typedef void	*pointer;
/* Not clear why we need to do this at all */
#ifndef linux
extern pointer	memset( pointer, int, size_t );
#endif
#endif

#if defined(ultrix)
#define const	/* nothing */
#endif

#ifndef DEBUG
#define	NDEBUG
#endif
#include	<assert.h>

typedef const unsigned char	cuc;	/* char variety used in algorithm */

#define EOS	'\0'			/* C string terminator */

char *					/* returns -> leftmost occurrence,
					   or null pointer if not present */
strcasestr(const char *s1, const char *s2 )
     /*	const char	*s1;	*/	/* -> string to be searched */
     /*	const char	*s2;	*/	/* -> search-pattern string */
{
  register cuc	*t;		/* -> text character being tested */
  register cuc	*p;		/* -> pattern char being tested */
  register cuc	*tx;		/* -> possible start of match */
  register size_t	m;		/* length of pattern */
  register cuc    *top;		/* -> high water mark in text */
#if UCHAR_MAX > 255			/* too large for auto allocation */
  static				/* not malloc()ed; that can fail! */
#endif					/* else allocate shift[] on stack */
    size_t	shift[UCHAR_MAX + 1];	/* pattern shift table */
  
#ifdef ROBUST				/* not required by C standard */
  if ( s1 == NULL || s2 == NULL )
    return NULL;		/* certainly, no match is found! */
#endif
  
  /* Precompute shift intervals based on the pattern;
     the length of the pattern is determined as a side effect: */
  
#ifdef ZAP
  (void)memset( (pointer)&shift[1], 0, UCHAR_MAX * sizeof(size_t) );
#else
  {
    register unsigned char	c;
    
    c = UCHAR_MAX;
    do
      shift[c] = 0;
    while ( --c > 0 );
  }
#endif
  /* Note: shift[0] is undefined at this point (fixed later). */
  
  for ( m = 1, p = (cuc *)s2; *p != EOS; ++m, ++p )
    shift[tolower((cuc)*p)] = m;
  
  assert(s2[m - 1] == EOS);
  
  {
    register unsigned char	c;
    
    c = UCHAR_MAX;
    do
      if(!isupper(c)) shift[c] = m - shift[c];
    while ( --c > 0 );
    
    /* Note: shift[0] is still undefined at this point. */
  }
  
  shift[0] = --m; 		/* shift[EOS]; important details! */
  
  assert(s2[m] == EOS);
  
  /* Try to find the pattern in the text string: */
  
  for ( top = tx = (cuc *)s1; ; tx += shift[tolower(*(top = t))] )
    {
      for ( t = tx, p = (cuc *)s2; ; ++t, ++p )
	{
	  if ( *p == EOS )       /* entire pattern matched */
	    return (char *)tx;
	  
	  if ( tolower(*p) != tolower(*t) )
	    break;
	}
      if ( t < top)	/* idea due to ado@elsie.nci.nih.gov */
	t = top;	/* already scanned this far for EOS */
      
      do	{
	assert(m > 0);
	assert(t - tx < m);
	
	if ( *t == EOS )
	  return NULL;	/* no match */
      }
      while ( ++t - tx != m );	/* < */
    }
}

void thAddVarToList(daVarStructList **head, daVarStruct *varp)
{
  daVarStructList *nextptr;

  if(head==0) {
    fprintf(STDERR,"thAddVarToList shouldn't have been called\n");
    return;
  }
  nextptr = *head;
  while(nextptr){
    if(nextptr->varp == varp)	/* Don't add duplicates */
      return;
    head = &(nextptr->next);
    nextptr = nextptr->next;
  }
  *head = malloc(sizeof(daVarStructList));
  (*head)->next = 0;
  (*head)->varp = varp;
  return;
}
#ifdef NOFNMATCH
/*
 * Copyright (c) 1989, 1993, 1994
 *	The Regents of the University of California.  All rights reserved.
 *
 * This code is derived from software contributed to Berkeley by
 * Guido van Rossum.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *	This product includes software developed by the University of
 *	California, Berkeley and its contributors.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#if defined(LIBC_SCCS) && !defined(lint)
static char sccsid[] = "@(#)fnmatch.c	8.2 (Berkeley) 4/16/94";
#endif /* LIBC_SCCS and not lint */

/*
 * Function fnmatch() as specified in POSIX 1003.2-1992, section B.6.
 * Compares a filename or pathname to a pattern.
 */

#include "fnmatch.h"

#define	EOS	'\0'

/*static const char *rangematch __P((const char *, int, int));*/
static const char *rangematch(const char *, int, int);

int
fnmatch(pattern, string, flags)
	const char *pattern, *string;
	int flags;
{
	const char *stringstart;
	char c, test;

	for (stringstart = string;;)
		switch (c = *pattern++) {
		case EOS:
			return (*string == EOS ? 0 : FNM_NOMATCH);
		case '?':
			if (*string == EOS)
				return (FNM_NOMATCH);
			if (*string == '/' && (flags & FNM_PATHNAME))
				return (FNM_NOMATCH);
			if (*string == '.' && (flags & FNM_PERIOD) &&
			    (string == stringstart ||
			    ((flags & FNM_PATHNAME) && *(string - 1) == '/')))
				return (FNM_NOMATCH);
			++string;
			break;
		case '*':
			c = *pattern;
			/* Collapse multiple stars. */
			while (c == '*')
				c = *++pattern;

			if (*string == '.' && (flags & FNM_PERIOD) &&
			    (string == stringstart ||
			    ((flags & FNM_PATHNAME) && *(string - 1) == '/')))
				return (FNM_NOMATCH);

			/* Optimize for pattern with * at end or before /. */
			if (c == EOS)
				if (flags & FNM_PATHNAME)
					return (strchr(string, '/') == NULL ?
					    0 : FNM_NOMATCH);
				else
					return (0);
			else if (c == '/' && flags & FNM_PATHNAME) {
				if ((string = strchr(string, '/')) == NULL)
					return (FNM_NOMATCH);
				break;
			}

			/* General case, use recursion. */
			while ((test = *string) != EOS) {
				if (!fnmatch(pattern, string, flags & ~FNM_PERIOD))
					return (0);
				if (test == '/' && flags & FNM_PATHNAME)
					break;
				++string;
			}
			return (FNM_NOMATCH);
		case '[':
			if (*string == EOS)
				return (FNM_NOMATCH);
			if (*string == '/' && flags & FNM_PATHNAME)
				return (FNM_NOMATCH);
			if ((pattern =
			    rangematch(pattern, *string, flags)) == NULL)
				return (FNM_NOMATCH);
			++string;
			break;
		case '\\':
			if (!(flags & FNM_NOESCAPE)) {
				if ((c = *pattern++) == EOS) {
					c = '\\';
					--pattern;
				}
			}
			/* FALLTHROUGH */
		default:
			if (c != *string++)
				return (FNM_NOMATCH);
			break;
		}
	/* NOTREACHED */
}

static const char *
rangematch(pattern, test, flags)
	const char *pattern;
	int test, flags;
{
	int negate, ok;
	char c, c2;

	/*
	 * A bracket expression starting with an unquoted circumflex
	 * character produces unspecified results (IEEE 1003.2-1992,
	 * 3.13.2).  This implementation treats it like '!', for
	 * consistency with the regular expression syntax.
	 * J.T. Conklin (conklin@ngai.kaleida.com)
	 */
	if (negate = (*pattern == '!' || *pattern == '^'))
		++pattern;
	
	for (ok = 0; (c = *pattern++) != ']';) {
		if (c == '\\' && !(flags & FNM_NOESCAPE))
			c = *pattern++;
		if (c == EOS)
			return (NULL);
		if (*pattern == '-' 
		    && (c2 = *(pattern+1)) != EOS && c2 != ']') {
			pattern += 2;
			if (c2 == '\\' && !(flags & FNM_NOESCAPE))
				c2 = *pattern++;
			if (c2 == EOS)
				return (NULL);
			if (c <= test && test <= c2)
				ok = 1;
		} else if (c == test)
			ok = 1;
	}
	return (ok == negate ? NULL : pattern);
}
#endif
