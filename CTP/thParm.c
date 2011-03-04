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
 *  Book parameters
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: thParm.c,v $
 *   Revision 1.2  2011/03/04 20:01:28  jones
 *   Used to be %li and %ld, but that makes 8 byte result stuffed into 4 byte lval
 *
 *   Revision 1.4  2003/02/21 20:55:24  saw
 *   Clean up some types and casts to reduce compiler warnings.
 *
 *   Revision 1.3  1999/11/04 20:34:06  saw
 *   Alpha compatibility.
 *   New RPC call needed for root event display.
 *   Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *   Revision 1.16  1999/08/16 16:31:10  saw
 *   Treat numbers that start with "0x" as hex.
 *
 *   Revision 1.15  1998/09/29 18:28:47  saw
 *   We shouldn't use thIDToken to identify whether the RHS of a parameter
 *   setting line is a simple constant or an expression.  So for now all RHS's
 *   will be evaluated.  Need some thought about setting the type of new
 *   variables that get created.  This eliminates the 1998 CTP bug.
 *
 *   Revision 1.14  1995/08/03 13:54:22  saw
 *   Add thpset function to single parameter setting lines from code
 *
 *	  Revision 1.13  1995/04/10  15:41:21  saw
 *	  Handle ctp file registration (#int, #real, ...)
 *
 *	  Revision 1.12  1995/01/09  15:26:09  saw
 *	  On fprintf, indicate block type and well as name
 *
 *	  Revision 1.11  1994/08/26  13:29:37  saw
 *	  Add DAVAR_REPOINTOK to created parameter.
 *
 *	  Revision 1.10  1994/07/21  20:35:44  saw
 *	  Don't prepend parm. when creating variables that have . in them.
 *
 *	  Revision 1.9  1994/06/13  13:21:04  saw
 *	  Fix up handling of string type CTP variables.
 *
 *	  Revision 1.8  1994/06/03  18:49:54  saw
 *	  Replace stderr with STDERR
 *
 *	  Revision 1.7  1994/02/14  20:23:29  saw
 *	  Comment out debugging printf's
 *
 *	  Revision 1.6  1994/02/08  21:34:01  saw
 *	  Remove debugging statement
 *
 *	  Revision 1.5  1993/12/02  21:34:47  saw
 *	  Fully allow doubles on parm left or right hand sides
 *
 *	  Revision 1.4  1993/09/22  17:27:39  saw
 *	  Convert integer values with sscanf to allow for octal and hex.
 *
 *	  Revision 1.3  1993/09/13  20:52:34  saw
 *	  Dynamically allocated arrays allowed.  Dynamic params will automatically
 *	  float if needed.
 *
 *	  Revision 1.2  1993/05/11  17:53:56  saw
 *	  Fix header
 *
 */

/* What to do about unregistered variables?
   Register them as Int's now.  Later allow declaring of reals and arrays.

*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "daVar.h"
#include "th.h"
#include "thInternal.h"
#include "thUtils.h"
#include "cfortran.h"

#define MAXLINELENGTH 512
/*#define NULL 0*/

/* Global variables */
int thParmVarIndex;
int thParmVarType;
int thParmVarALen;
int *thParmVarLP;
float *thParmVarFP;
double *thParmVarDP;
char *thParmVarSP;
char *thParmVarName;
daVarStruct *thParmVarVarp;
int thParmVarDynamic;
/**/

char *classlist[]={PARMSTR,0};	/* Class list for parameter names */

thStatus thParmLineSet(char *line);
/*FCALLSCFUN1(INT,thLoadParameters,LOADPARM,loadparm,STRING)*/

thStatus thLoadParameters(daVarStruct *var)
/* Set the parameters as specified on the title line.
   When done, replace title line with the parameters without the values?? 
   For now, we won't modify the lines. */
{
  char *lines,*eol;
  int line_count;

  if(*((DAINT *)var->varptr) != 0) /* This block already booked */
    return(S_SUCCESS);
  *((DAINT *) var->varptr) = 1;
  lines = var->title;
  line_count = 0;
  while(*lines){
    char *lcopy;
    
    line_count++;
    eol = strchr(lines,'\n');
    if(!eol) {
      fprintf(STDERR,"L %d: Last line of parm block %s has no newline\n"
	      ,line_count,var->name);
      break;
    }
    if(*(eol+1)=='\0'){		/* This is the last line */
      if(strcasestr(lines,ENDSTR) == 0) {
	fprintf(STDERR,"L %d: Last line of parm block %s is not an END\n"
		,line_count,var->name);
      }
      break;
    }
    if(line_count == 1)
      if(strcasestr(lines,BEGINSTR) != 0){
	/*	printf("Is a begin\n");*/
	lines = eol + 1;
	continue;
      } else
	fprintf(STDERR,"First line of parm block %s is not a BEGIN\n",var->name);
    /* Ready to book the line, Add continuation lines later */
    lcopy = (char *) malloc(eol-lines+1);
    strncpy(lcopy,lines,(eol-lines));
    *(lcopy + (eol-lines)) = '\0';
/*     printf("Passing|%s|\n",lcopy);*/
    if(!thSpecial(lcopy,PARMSTR)) {
      if(thParmLineSet(lcopy)!=S_SUCCESS)
	fprintf(STDERR,"Error saving parameters on line %d\n",line_count);
    }
    free(lcopy);
    lines = eol + 1;
  }
  return(S_SUCCESS);
}
thStatus thParmLineSet(char *line)
     /* Process a line of a parameter CTP block */
{
  thTokenType toktyp;
  int vartyp;
  int vardimen;
  daVarStruct *varp;
  int i;
  char *varnam;
  int nargs;
  char *args[50];
  char *orgargs;		/* Unadulterated arguments line */

  {				/* Needs to be fixed to handle strings */
    char *s;
    int blank;
    char quotechar;
    int instring;

    s = line;
    blank = 1;
    instring = 0;
    while(*s != 0){
      if(instring && *s == quotechar) {
	if(*(s+1) == quotechar) s++;
	else instring = 0;
      } else {
	if(*s == QUOTECHAR1 || *s== QUOTECHAR2) {
	  instring = 1;
	  quotechar = *s;
	  blank = 0;
	} else if(isspace(*s)) {
	  *s = ' ';	/* Remove tabs, ... */
	} else if(*s == COMCHAR) {
	  *s = 0;
	  break;
	} else
	  blank = 0;
      }
      s++;
    }
    if(blank) return(S_SUCCESS);
    /* Now look for = and figure out what kind of variable is on left.
       If more than one number is given, left must be array. */

    s = line;
    orgargs = 0;
    if((s = strchr(s,'='))){
      *s++ = '\0';
      orgargs = (char *) malloc(strlen(s)+1);
      strcpy(orgargs,s);
      varnam = thSpaceStrip(line);
    } else {
      s = line;
      varnam = 0;
    }
    nargs = thCommas(s,args);
    for(i=0;i<nargs;i++){
      args[i] = thSpaceStrip(args[i]);/* Remove all space from the argument */
      /*      printf("%s ",args[i]);*/
    }
    if(nargs > 0)		/* If only white space after last comma, */
      if(args[nargs-1][0] == '\0') /* then don't count it as an argument */
	nargs--;
  }

  if(varnam){
    toktyp = thIDToken(varnam);
    if(toktyp != TOKVAR && toktyp != TOKARRAY){
      fprintf(STDERR,"Variable name %s can't be a number\n",varnam);
      if(orgargs) free(orgargs);
      return(S_FAILURE);
    }
    if(toktyp == TOKARRAY){
      char *p;
      p = thTokenArray(varnam,&thParmVarIndex);
      *p = 0;
    } else
      thParmVarIndex = 0;
    if(daVarLookupPWithClass(varnam,classlist,&varp) != S_SUCCESS) {
      /* Variable is not preregistered, we will automatically allocate it.
	 Later, a flag for the block will be added which will disallow
	 auto allocation.  We will allocate here an integer array of length
	 thParmVarIndex plus the number of arguments on the line.  If there
	 are subsequent lines, the array will automatically be extended below
	 since the code will see the DAVAR_DYNAMIC_PAR flag.  If floating
	 point values are found on the lines below, the array will
	 automatically be changed to floating.  (Perhaps we should just always
	 make the variables floating point.)
	 The noauto flag will  eventually be implimented to disable automatic
	 variable createion. */
      daVarStruct var;
      if(strchr(varnam,'.')) {	/* Don't prepend parm., if varname has '.'s */
	var.name = (char *) malloc(strlen(varnam)+1);
	strcpy(var.name,varnam);
      } else {
	var.name = (char *) malloc(strlen(classlist[0])
				   +strlen(varnam)+2);
	strcpy(var.name,classlist[0]);
	strcat(var.name,".");
	strcat(var.name,varnam);
      }
      var.size = thParmVarIndex + nargs;
      var.varptr = (void *) malloc(var.size*sizeof(DAINT));
      var.opaque = 0;
      var.rhook = 0;
      var.whook = 0;
      var.type = DAVARINT;
      var.flag = DAVAR_REPOINTOK | DAVAR_READONLY | DAVAR_DYNAMIC_PAR;
      var.title = 0;
      daVarRegister((int) 0,&var); /* parameter */
      daVarLookupP(var.name,&varp);
      free(var.name);
      fprintf(STDERR,"%s not registered, registering as int\n",varnam);
    }
    /*    vardimen = varp->dimension;
	  if(vardimen==0 && toktyp == TOKARRAY) {
	  fprintf(STDERR,"Variable %s not registered as an array\n",varnam);
	  return(S_FAILURE);
	  }
	  */
    vartyp = varp->type;
    thParmVarType = vartyp;
    thParmVarName = varp->name;
    thParmVarVarp = varp;
    thParmVarDynamic = (varp->flag & DAVAR_DYNAMIC_PAR);
    thParmVarALen = varp->size;
    switch(vartyp)
      {
      case DAVARINT:
	thParmVarLP = (int *) varp->varptr;
	break;
      case DAVARFLOAT:
	thParmVarFP = (float *) varp->varptr;
	break;
      case DAVARDOUBLE:
	thParmVarDP = (double *) varp->varptr;
	break;
      case DAVARSTRING:
      case DAVARFSTRING:
	thParmVarSP = (char *) varp->varptr;
	break;
      }
  }
  if(thParmVarType == DAVARINT || thParmVarType == DAVARFLOAT 
     || thParmVarType == DAVARDOUBLE) {
    for(i=0;i<nargs;i++){
      int lval;
      double dval;
      toktyp = thIDToken(args[i]);
      if(thParmVarIndex>=thParmVarALen){
	if(thParmVarDynamic) {	/* Automatically up size for dynamic pars */
/*	  printf("i=%d, thParmVarIndex=%d, thParmVarALen=%d\n",i,thParmVarIndex
		 ,thParmVarALen);
	  printf("thParmVarType=%d\n",thParmVarType);*/
	  if(thParmVarType == DAVARINT){
	    int j;
	    int *TMPP;
	    
	    thParmVarVarp->size = (thParmVarIndex + (nargs-i));
	    TMPP = (int *) malloc(thParmVarVarp->size*sizeof(DAINT));
	    for(j=0;j<thParmVarALen;j++)
	      TMPP[j] = thParmVarLP[j];
	    free(thParmVarLP);
	    thParmVarLP = TMPP;
	    thParmVarVarp->varptr = thParmVarLP;
	  } else { /*if(thParmVarType == DAVARFLOAT)*/
	    int j;	
	    float *TMPP;
	    
	    thParmVarVarp->size = (thParmVarIndex + (nargs-i));
	    TMPP = (float *) malloc(thParmVarVarp->size*sizeof(DAFLOAT));
	    for(j=0;j<thParmVarALen;j++)
	      TMPP[j] = thParmVarFP[j];
	    free(thParmVarFP);
	    thParmVarFP = TMPP;
	    thParmVarVarp->varptr = thParmVarFP;
	  }
	  thParmVarALen = thParmVarVarp->size;
	} else {
	  fprintf(STDERR,"Tried to fill past end of array %s\n",thParmVarName);
	  if(orgargs) free(orgargs);
	  return(S_FAILURE);
	}
      }
#define ALWAYSEVAL
#ifndef ALWAYSEVAL
      switch(toktyp)
	{
	case TOKINT:
	  /* Used to be %li and %ld, but that makes 8 byte result
	     stuffed into 4 byte lval */
	  if(args[i][0] == '0' && (args[i][1] == 'x' || args[i][1] == 'X')) {
	    sscanf(args[i],"%i",&lval); /* Treat as Hex */
	  } else {
	    sscanf(args[i],"%d",&lval); /* Treat as decimal */
	  }
	  dval = lval;
	  break;
	case TOKFLOAT:
	  dval = atof(args[i]);
	  lval = floatToLong(dval);
	  break;
	default:
#endif
	  if(thEvalImed(args[i],&dval,&lval) != S_SUCCESS)
	    fprintf(STDERR,"Parm: Error interpreting %s\n");
#ifndef ALWAYSEVAL
	  break;
	}
#endif
      switch(thParmVarType)
	{
	case DAVARINT:	
	  if(thParmVarDynamic) {
	    /* User must be careful, if an expression evaluated with thEvalImed ends up
	       as integer, then the type of the variable will stay as integer. */
	    if(toktyp == TOKFLOAT || (toktyp != TOKINT && dval != lval)) {
	      /* Floating point arg found */
	      int j;		/* Copy integer arry to float array */
	      thParmVarFP = (float *) malloc(thParmVarALen*sizeof(DAFLOAT));
	      for(j=0;j<thParmVarALen;j++)
		thParmVarFP[j] = thParmVarLP[j];
	      free(thParmVarLP);
	      thParmVarVarp->varptr = thParmVarFP;
	      thParmVarVarp->type = DAVARFLOAT;
	      thParmVarFP[thParmVarIndex++] = dval;
	      thParmVarType = DAVARFLOAT;
	      break;
	    }
	  }
	  thParmVarLP[thParmVarIndex++] = lval;
	  break;
	case DAVARFLOAT:
	  thParmVarFP[thParmVarIndex++] = dval;
	  break;
	case DAVARDOUBLE:
	  thParmVarDP[thParmVarIndex++] = dval;
	  break;
	}
      /*    printf("Saved args[%d] %s %d %f\n",i,args[i],lval,dval);*/
    }
  } else if(thParmVarType == DAVARSTRING || thParmVarType == DAVARFSTRING) {
    int maxlen, arglen;
    char *argptr; char *s;

    maxlen = thParmVarALen - ((thParmVarType == DAVARSTRING) ? 1 : 0);
    /* Find first non blank character after the = */
    argptr = orgargs;

    while(isspace(*argptr)) argptr++;
    if(argptr[0] == QUOTECHAR1 || argptr[0] == QUOTECHAR2){
      s = argptr+1;
      while(*s && *s != argptr[0]) s++;       /* Search for nul or matching qu
e */
      *s = 0;
      argptr++;                       /* Move to char after quote */
    }
    arglen = strlen(argptr);
    if(arglen > maxlen) arglen = maxlen;
    strncpy(thParmVarSP, argptr, arglen);
    if(thParmVarType == DAVARFSTRING) {
      while(arglen < maxlen)
	thParmVarSP[arglen++] = ' ';    /* Blank pad fortran strings */
    } else {
      thParmVarSP[arglen] = 0;
    }
  }
  if(orgargs) free(orgargs);
  return(S_SUCCESS);
}
/* Fortran routine to evaluate a line of the form parm = value */  
#ifdef NOF77extname
int thpset
#else
int thpset_
#endif
(char *A1,unsigned C1)
{
  int A0;
  char *B1;
  thStatus status;

  status = thParmLineSet((!*(int *)A1)?0:memchr(A1,'\0',C1)?A1:
		      (memcpy(B1=malloc(C1+1),A1,C1),B1[C1]='\0'
		       ,kill_trailing(B1,' ')));
  if(B1) free(B1);
  return status;
}
