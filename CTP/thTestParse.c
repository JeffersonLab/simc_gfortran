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
 *  The expression parser and stack executor for the Test Package
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: thTestParse.c,v $
 *   Revision 1.2  2011/03/04 20:01:35  jones
 *   Used to be %li and %ld, but that makes 8 byte result stuffed into 4 byte lval
 *
 *   Revision 1.4.24.1  2007/09/10 21:32:47  pcarter
 *   Implemented changes to allow compilation on RHEL 3,4,5 and MacOSX
 *
 *   Revision 1.4  2003/02/21 20:55:25  saw
 *   Clean up some types and casts to reduce compiler warnings.
 *
 *   Revision 1.3  1999/11/04 20:34:07  saw
 *   Alpha compatibility.
 *   New RPC call needed for root event display.
 *   Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *   Revision 1.19  1999/08/25 13:16:07  saw
 *   *** empty log message ***
 *
 *   Revision 1.18  1999/07/07 13:43:51  saw
 *   Don't make numbers starting with "0" be octal.  Accept 0x as hex.
 *
 *   Move thTestRHandler() into thTestParse.c
 *
 *   Revision 1.17  1999/03/01 20:00:50  saw
 *   Fix bug where a series of numbers added or subtracted in a parameter line
 *   got evaluated to be just the first number.  Improve the scientific
 *   number detection to do this.
 *
 *   Revision 1.16  1996/07/31 20:36:56  saw
 *   Support floating point for mod command.  Add trig functions.
 *
 * Revision 1.15  1995/08/03  13:56:36  saw
 * Add single argument functions
 *
 * Revision 1.14  1995/04/10  15:51:21  saw
 * thEvalImed returns INTOVF if double result is to large to convert to int.
 *
 * Revision 1.13  1995/02/14  16:53:52  saw
 * Make compatible with OSF/Alpha (64 bit pointers)
 *
 * Revision 1.12  1995/01/09  16:06:11  saw
 * Fix a short malloc for a string
 *
 * Revision 1.11  1994/11/17  18:14:21  saw
 * Strip out unary + operators when parsing expressions
 *
 * Revision 1.10  1994/11/07  14:28:34  saw
 * Add thevalchk fortran call to check for expressions.
 * Try to avoid bomb outs for bad expressions
 *
 * Revision 1.9  1994/10/03  12:41:22  saw
 * All "/" (division) has real result.  New op "//" has integerized
 * result.  thEvalImed actually gets a double from thTestExecute.
 * Added fortran interfaces to thEvalImed (itheval, ftheval, dtheval).
 *
 * Revision 1.8  1994/09/12  15:12:31  saw
 * thGetTok was missing reset of lastop on EOL
 *
 * Revision 1.7  1994/08/29  20:08:10  saw
 * Fix calculation of testscalervarname length
 *
 * Revision 1.6  1994/08/26  17:46:18  saw
 * Register test scaler results
 *
 * Revision 1.5  1994/08/26  13:36:46  saw
 * Add DAVAR_REPOINTOK to some flags
 *
 * Revision 1.4  1994/06/03  18:54:29  saw
 * Replace stderr with STDERR
 *
 * Revision 1.3  1993/12/02  21:33:36  saw
 * Fully allow use of doubles in test expressions
 *
 * Revision 1.2  1993/11/24  21:24:54  saw
 * thEvalImed now returns double instead of floating result.
 *
 *	  Revision 1.3  1993/09/22  17:51:06  saw
 *	  Allow integer constants to be octal or hex.
 *
 *	  Revision 1.2  1993/05/11  18:00:10  saw
 *	  Update header
 *
 */

/* thTestParse.c
   
   Make test result variable that are created take the type of the rhs???
   Add variable names to stack so that expressions can be recreated.
   Allow constants to be hex or octal.
   Agree on a new comment character or syntax since ! is not part
   of expressions.
   Add unary operators to executor.  Allow + to be a unary operator too.

   Need to build up a linked list of test results used in a block.  Don't
   duplicate any variables.  Print warning when a scaler test result is
   reused.

*/
/*An argument is a variable name, an array, or a number.  Numbers are not
allowed for test result.  Arrays start at 0 if []'s are used and start
at 1 if ()'s are used.  Arrays may only be used for test results if they
are already registered by the analyzer.  (May add option to declare them
in the test package.)*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#define INT_MAX 2147483647
/* limits.h is used only to get #define INT_MAX 2147483647
 * If you don't have limits.h, try #include <values.h> instead and then
 * #define INT_MAX MAXINT */
//#include <limits.h>

#include <rpc/rpc.h>
#include "daVar.h"
#include "daVarRpc.h"
#include "daVarHandlers.h"
#include "th.h"
#include "thUtils.h"
#include "thTestParse.h"
#include "thInternal.h"
#include "cfortran.h"

daVarStatus thTestRHandler(char *name, daVarStruct *varclass, any *retval);

CODE opstack[100];		/* Operator stack */
CODE typstack[100];		/* Result type stack */

typedef struct
{
  char *ops;
  int  toks[3];
} OPTABLE;

OPTABLE optable[] =
{
  {"(",{OPLP}},
  {")",{OPRP}},
  {"[",{OPLINDEXB}},
  {"]",{OPRP}},
  {"-",{OPSUB}},
  {"+",{OPADD}},
  {"<<=",{OPISLT,OPSHL,OPISLE}},
  {">>=",{OPISGT,OPSHR,OPISGE}},
  {"==",{OPEQUAL,OPISEQUAL}},
  {"!=",{OPNOT,OPISNOTEQUAL}},
  {"&&",{OPBITAND,OPLOGAND}},
  {"||",{OPBITOR,OPLOGOR}},
  {"^^",{OPBITXOR,OPLOGXOR}},
  {"*",{OPTIMES}},
  {"//",{OPDIV,OPIDIV}},
  {"%",{OPMOD}},
  {"~",{OPCOMP}},
  {",",{OPCOMMA}},
  {0,{0,0,0}}};
static char *opchars=0;

/* For Q like test package format, must be in same order as type
   types listed in typedef for thTestType. */
char *testCodes[tBAD] 
  = {"GA","PA","EQ","BI","AN","IO","EO","MA","US"};

typedef struct
{
  CODE op;
  CODE result[9];
} TYPETABLE;

TYPETABLE typetable[] =
{
  {OPLINDEX,{0,0,0,1,1,1,2,2,2}}, /* Result is same as variable */
  {OPLINDEXB,{0,0,0,1,1,1,2,2,2}}, /* being indexed */
  {OPLINDEXP,{0,0,0,1,1,1,2,2,2}}, /* Result is same as variable */
  {OPLINDEXPB,{0,0,0,1,1,1,2,2,2}}, /* being indexed */
  {OPEQUAL,{0,0,0,1,1,1,2,2,2}}, /* Set result type to LHS type */
  {OPLOGOR,{0,0,0,0,0,0,0,0,0}}, /* Result is always integer */
  {OPLOGXOR,{0,0,0,0,0,0,0,0,0}},
  {OPLOGAND,{0,0,0,0,0,0,0,0,0}},
  {OPBITOR,{0,0,0,0,0,0,0,0,0}},
  {OPBITXOR,{0,0,0,0,0,0,0,0,0}},
  {OPBITAND,{0,0,0,0,0,0,0,0,0}},
  {OPISEQUAL,{0,0,0,0,0,0,0,0,0}},
  {OPISNOTEQUAL,{0,0,0,0,0,0,0,0,0}},
  {OPISLT,{0,0,0,0,0,0,0,0,0}},
  {OPISLE,{0,0,0,0,0,0,0,0,0}},
  {OPISGT,{0,0,0,0,0,0,0,0,0}},
  {OPISGE,{0,0,0,0,0,0,0,0,0}},
  {OPSHL,{0,0,0,0,0,0,0,0,0}},
  {OPSHR,{0,0,0,0,0,0,0,0,0}},
  {OPADD,{0,2,2,2,2,2,2,2,2}},	/* Result is double unless both ops int */
  {OPSUB,{0,2,2,2,2,2,2,2,2}},
  {OPTIMES,{0,2,2,2,2,2,2,2,2}},
  {OPDIV,{2,2,2,2,2,2,2,2,2}},	/* Result always double */
  {OPIDIV,{0,0,0,0,0,0,0,0,0}},	/* Result always integer */
  {OPMOD,{0,2,2,2,2,2,2,2,2}},
  {OPNEG,{0,1,2,0,1,2,0,1,2}},	/* No lh operand, type = rh type */
  {OPNOT,{0,0,0,0,0,0,0,0,0}},	/* No lh operand, type always int */
  {OPCOMP,{0,0,0,0,0,0,0,0,0}},	/* No lh operand, type always int */
  {0,{0,0,0,0,0,0,0,0,0}}};

INTRINSIC_FUNCTIONS intrinsic_functions[] =
{
  {"abs",{0,1,2}},
  {"sqrt",{2,2,2}},
  {"exp",{2,2,2}},
  {"sin",{2,2,2}},
  {"cos",{2,2,2}},
  {"tan",{2,2,2}},
  {0,{0,0,0}}
};
char *thGetTok(char *linep, int *tokenid, char **tokstr,
	       CODE *tokval, void **tokptr, int expflag, daVarStructList **vlisthead)
/* Pass a pointer to the unscanned portion of the line.
   Returns An ID code for operators, and an operand type for operands in
   tokenid.
   Returns the string for the operand in tokstr.  (Null otherwise)
   Returns the operand value in tokval, or in tokptr if the operand is
   a pointer.
   If the operand is a function, then tokenid will be pushfunction, and
   tokval will be a the fuction id.

   Function returns pointer to remainder of the line.


   */
{
  static char string[100];
  static int lasttoktype=0;	/* Last tok was an operator */
  static CODE lastop=0;

  char *savelinep;
  char *stringp;
  char *ptr,c;
  int tindex,sindex;
  daVarStruct *varp;
  DAFLOAT f;

  /* Build up a list of characters that can start operators */
  if(opchars == 0){
    int count=0;
    int i;

    while(optable[count++].ops != 0) ;
    opchars = (char *) malloc(count);
    for(i=0;i<(count-1); i++)
      opchars[i] = optable[i].ops[0];
    opchars[count-1] = 0;
  }

  *tokstr = 0;
  *tokval = 0;
  *tokptr = 0;
  *tokenid = 0;			/* Will signify an undeclared operand */

  if(!(*linep)) {
    *tokenid = OPEOL;
    lasttoktype = 0;
    lastop = 0;
    return(0);
  }
  savelinep = linep;
  while(*linep == ' ' || *linep == '\t') linep++;
  if((ptr = strchr(opchars,*linep))) { /* Operator */
    tindex = ptr -  opchars;
    if(lasttoktype == 0 && *linep == '-') { /*  Last thing was an operator */
      *tokenid = OPNEG;		/* So the '-' must be a negative sign */
      linep++;
    } else if(lasttoktype == 0 && *linep == '+') { /* Unary plus */
      linep++;
      goto operand;
    } else if(lasttoktype == 1 && *linep == '(') {
      *tokenid = OPLINDEX;
      linep++;
    } else if(lasttoktype == 3 && *linep == '(') {
	/* How will we know when the right hand operator is the closing
	   paren of the function?  We don't need to know.  The RHP only
	   acts to determine precedence.  */
      *tokenid = OPLFARG;
      linep++;
    } else {
      linep++;
      *tokenid = optable[tindex].toks[0];
      sindex = 1;
      if(*linep) {		/* Don't search past end of line */
	while((c = optable[tindex].ops[sindex])) {
	  if(*linep == c) {
	    *tokenid = optable[tindex].toks[sindex];
	    linep++;
	    break;
	  }
	  sindex++;
	}
      }
    }
    /*  Following two lines were before last }. */
    if(*tokenid == OPRP) lasttoktype = 2; /* Minus is minus after ) or ] */
    else lasttoktype = 0;
    /* For OPLINDEX and OPLINDEXB, need to search ahead for matching ) or ]
       and check if the next operator is an = not ==).  If so, then we
       need to return OPLINDEXP or OPLINDEXPB.  */
    if(*tokenid == OPLINDEX || *tokenid == OPLINDEXB){
      char *p; char rc; int ccount=0; int bcount=0;
      if(*tokenid == OPLINDEXB) rc = ']';
      else rc = ')';
      p = linep;
      while(*p && (*p != rc || bcount || ccount)) {
	switch(*p++)
	  {
	  case '(': ccount++; break;
	  case ')': ccount--; break;
	  case '[': bcount++; break;
	  case ']': bcount--; break;
	  default: break;
	  }
      }				/* Only NULL or balanced rc terminates */
      if(*p++){			/* Search for = */
	while(*p == ' ' || *p=='\t') p++;
	if(*p == '=' && *(p+1) != '=') {
	  *tokenid += (OPLINDEXP - OPLINDEX); 
	  /* Assumes OPLINDEXBP-OPLINDEXB is the same*/
	}
      } else
	fprintf(STDERR,"thTest: Parenthesis balance problem\n");
    }
    lastop = *tokenid;
  } else {			/* Operand */
 //   int optype;
    int isnum;
    int efound;

  operand:
    lasttoktype = 1;
    stringp = string;
    /* Scan until operator or whitespace is reached */
    isnum = 1; efound = 0;
/* What a hack to check for scientific notation */
    while(*linep && *linep!=' ' && *linep!='\t' && !strchr(opchars,*linep)){
      if(*linep == 'e' || *linep == 'E') {
	if(efound) {
	  isnum = 0;
	} else {
	  if(stringp > string) {
	    efound=1;
	  } else {
	    isnum = 0;
	  }
	}
      } else if(!isdigit(*linep) && *linep != '.') isnum = 0;
      *stringp++ = *linep++;
    }
    if(isnum && efound) {	/* Exponential, scan past last digit */
      if(*linep == '-' || *linep == '+') *stringp++ = *linep++;
      while(isdigit(*linep)) {
	*stringp++ = *linep++;
      }
    }
    while(*linep == ' ' || *linep == '\t') linep++; /* Skip past whitespace */
    *stringp = 0;
    *tokstr = string;
/*    printf("token=%s\n",string);*/
    switch(thIDToken(string))
      {
      case TOKINT:
	/* Used to be %li and %ld, but that makes 8 byte result stuffed into
	   4 byte *tokval */
	if(string[0] == '0' && (string[1] == 'x' || string[1] == 'X')) {
	  sscanf(string,"%i",tokval); /* Treat as Hex */
	} else {
	  sscanf(string,"%d",tokval); /* Treat as decimal */
	}
	*tokenid = OPPUSHINT;
	break;
      case TOKFLOAT:
	f = atof(string);
	*tokval = *(DAINT *)&f;	/* Copy floating value */
	*tokenid = OPPUSHFLOAT;
	break;
      case TOKVAR:
	{
	  char **classlist;
	  thOperandType optype;

	  optype = thGetOperandType(string,linep,lastop,0);
	  classlist = thGetClassList(optype);
	  /* Probably should consistently use the same class list of
	     TEST,EVENT,PARM here */
/* If token is a result variable (and we are in non-immediate mode), and the
   variable is an integer type, then we need to add this variable to a list
   of variables for the current block.   (Probably add real  variable to
   the list anyway. ) This will allow us to acumulate scalers.  The opaque
   pointer of each variable in the list will point to the scaler array. */

	  /* First check if variable is really an intrinsic function */
	  {
	    int ifunc;
	    ifunc = 0;
	    while(intrinsic_functions[ifunc].name) {
	      if(strcasecmp(string,intrinsic_functions[ifunc].name)==0) {
		*tokenid = OPPUSHFUNCTION;
		*tokval = ifunc;
		lasttoktype = 3;
		break;
	      }
	      ifunc++;
	    }
	    if(*tokenid) break;	/* Hopefully this breaks out of case */
	  }
	  if(daVarLookupPWithClass(string,classlist,&varp) == S_SUCCESS) {
/*	    printf("Found variable %s[%s]\n",string,varp->name);*/
	    if(varp->type == DAVARFLOAT) {/* If next operator is a ( or [ */
/*	      printf("FLOAT ");*/
	      /* then push pointer instead of */
#define ISARRAYORLHS(x) (*x=='(' || *x=='[' || (*x=='=' && *(x+1)!='='))
	      *tokenid = ISARRAYORLHS(linep) ? OPPUSHPFLOAT : OPPUSHFLOATP; 
	    } else if(varp->type == DAVARDOUBLE){	/* value onto rpn stack */
/*	      printf("DOUBL ");*/
	      *tokenid = ISARRAYORLHS(linep) ? OPPUSHPDOUBLE : OPPUSHDOUBLEP;
	    } else if(varp->type == DAVARINT){	/* value onto rpn stack */
/*	      printf("INT   ");*/
	      *tokenid = ISARRAYORLHS(linep) ? OPPUSHPINT : OPPUSHINTP;
	    }
	    else {
	      fprintf(STDERR
		      ,"thTest: Variable %s[%s] must be integer, float or double\n"
		      ,string,varp->name);
	    }
/*	    *tokval = *(DAINT *)&varp->varptr;*/ /* Get the pointer */
	    *tokptr = varp->varptr;
	  } else if(*linep=='=' && (*(linep+1)!='=')){
	    /* Undefined variable is an unindexed result */
	    /* For now, create an integer variable.  Later figure out how
	       to make the variable the same type as the rhs */
	    daVarStruct var;
	    var.name = (char *) malloc(strlen(classlist[0])
				       +strlen(string)+2);
	    strcpy(var.name,classlist[0]);
	    strcat(var.name,".");
	    strcat(var.name,string);
	    var.varptr = (void *) malloc(sizeof(DAINT));
	    var.size = 1;
	    var.opaque = 0;
	    var.rhook = 0;
	    var.whook = 0;
	    var.type = DAVARINT;
	    var.flag = DAVAR_READONLY | DAVAR_REPOINTOK;
	    var.title = savelinep;
	    daVarRegister((int) 0,&var); /* Create test result */
	    daVarLookupP(var.name,&varp);
	    free(var.name);
printf("Created test result %s\n",varp->name);
	    *tokenid = OPPUSHPINT;
/*	    *tokval = *(DAINT *)&varp->varptr;*/
	    *tokptr = varp->varptr;
	  } /* else 
	       printf("%s not found\n",string);
	       }*/
	  /* If variable does not exist, caller will note that toktype and
	     tokval have not been set. */
	  if(optype == otRESULT && vlisthead){ /* Don't make scalers for */
	    thAddVarToList(vlisthead,varp); /* Variables created in */
	    if(varp->type == DAVARINT) { /* thEvalImed */
	      DAINT *sarray; int i;
	      if(varp->opaque == 0) { /* No scaler array yet */
		char *testscalervarname;
		daVarStruct *svarp; /* Pointer to scaler var struct */
		testscalervarname = /* Add the "scaler" attribute */
		  (char *) malloc(strlen(varp->name)+strlen(SCALERSTR)+2);
		strcpy(testscalervarname,varp->name);
		strcat(testscalervarname,".");
		strcat(testscalervarname,SCALERSTR);
		if(daVarLookupP(testscalervarname,&svarp) != S_SUCCESS) {
		  daVarStruct svar;
		  svar.name = testscalervarname;
		  svar.opaque = 0;
		  svar.rhook = 0;
		  svar.whook = 0;
		  svar.type = DAVARINT;
		  svar.flag = DAVAR_READONLY | DAVAR_REPOINTOK;
		  svar.varptr = (void *) malloc(varp->size*sizeof(DAINT));
		  svar.size = varp->size;
		  /* Actually not OK to repoint, but this says CTP made it */
		  svar.title = varp->name;
		  daVarRegister((int) 0, &svar);
		  daVarLookupP(svar.name,&svarp);
		}
		varp->opaque = (DAINT *) svarp->varptr;
		varp->rhook = thTestRHandler;
		free(testscalervarname);
	      }
	      sarray = varp->opaque;
	      for(i=0;i<varp->size;i++)
		sarray[i] = 0;
	    }
	  }
	}
	break;
      default:
	fprintf(STDERR,"thTest: Error understanding %s\n",string);
	break;
      }
/*    printf("token = %x\n",*tokenid);*/
    lastop = 0;
  }
  while(*linep == ' ' || *linep == '\t') linep++; /* Skip whitespace */
  return(linep);
}
  

char **thGetClassList(thOperandType optype)
{
  static char *explist[]={PARMSTR,EVENTSTR,TESTSTR,0}; /* Immediate expressions */
  static char *loglist[]={TESTSTR,PARMSTR,EVENTSTR,0}; /* Logical operand */
  static char *numlist[]={EVENTSTR,PARMSTR,TESTSTR,0}; /* Operand is a value */
  static char *resultlistp[]={TESTSTR,EVENTSTR,PARMSTR,0}; /* Operand is a result */

#define ALWAYSTESTFIRST
#ifdef ALWAYSTESTFIRST
  return(resultlistp);
#else
  switch(optype)
    {
    case otIMMED:
      return(explist);
    case otLOGIC:
      return(loglist);
    case otVALUE:
      return(numlist);
    case otRESULT:
      return(resultlistp);
    }
#endif
}

thOperandType thGetOperandType(char *soperand, char *rest, CODE lastop,
			       int expflag)
{
  if(expflag)
    return(otIMMED);
  else if(lastop == OPNOT)
    return(otLOGIC);
  else if(lastop != 0 && (lastop != OPLOGOR)
		  && (lastop != OPLOGAND) && (lastop != OPLOGXOR)
		  && (lastop != OPEQUAL) && (lastop != OPCOMMA)
		  && (lastop != OPLP))
    return(otVALUE);
  else {
    /* This is really ugly code to determine if the operand
       is a result, logical operand, or numerical operand from the
       surrounding operators.  The last operator is known, but it must
       search ahead for the next operator.  This code should be  burried
       in a subroutine. */
    
    char *p;
    p = rest;
    if(*p == '(' || *p == '[') {
      int ccount=0; int bcount=0;
      if(*p++ == '(') ccount++; else bcount++;
      while(*p && (bcount || ccount)){
	/*	      printf("%c(%d,%d)\n",*p,ccount,bcount);*/
	switch(*p++) {
	case '(': ccount++; break; case ')': ccount--; break;
	case '[': bcount++; break; case ']': bcount--; break;
	default: break;
	}
      }
      /*	    printf("pos=%c, %d %d ",*p,ccount,bcount);*/
    }
    while(*p == ' ' || *p =='\t') p++;
#define ISLOG(x,y) (*x==y && *(x+1)==y)
    /*	  printf(", Nextchar=%c:  ",*p);*/
    if(*p=='=' && *(p+1)!='=') {
      return(otRESULT);
    } else if((ISLOG(p,'|') || ISLOG(p,'&') || ISLOG(p,'^')
	       || *p=='\0' || *p==',' || *p == ')'))
      return(otLOGIC);
  }
  return(otVALUE);
}

CODE thGetResultType(CODE operator, CODE leftoptype, CODE rightoptype)
{
  /* For a given operator, determine the data type of the result a given
     combination of the types of the lh and rh operands.
     Assumes that only types 0, 1, or 2 are allowed. */
  
  int lrindex;
  int i;

  if(leftoptype < 0 || leftoptype > 2 || rightoptype < 0 || rightoptype > 2) {
    fprintf(STDERR,"thTest: Illegal operand type %x %x\n",leftoptype,rightoptype);
    return(0);
  }
  lrindex = (leftoptype * 3) + rightoptype;
  for(i=0; typetable[i].op; i++) { /* Do Linear search for the operator */
    if(operator == typetable[i].op) {
      return(typetable[i].result[lrindex]);
    }
  }
  fprintf(STDERR,"Operator %x not found in result type table\n",operator);
  return(0);
}

thStatus thEvalImed(char *line, DADOUBLE *d, DAINT *i)
/* ImmedOBiately evaluate the expression in line.  Will internally evaluate to
   a float, and then pass back both the float and interized values. */
{
  CODEPTR codehead, codenext, codelimit, codelastop;
  int codesize;
#define RDOUBLE
#ifdef RDOUBLE
  double result;
#else
  float result;			/* Should    change to double */
#endif
  thStatus retcode;

/*  printf("%s=",line);*/
  codesize = 10+2*strlen(line);
  codehead = codenext = (CODEPTR) malloc(sizeof(CODE)*codesize);
  codelimit = codehead + codesize;
#ifdef RDOUBLE
  *codenext++ = OPPUSHPDOUBLE;
#ifdef USEMEMCPY
  {
    void *resultp;
    resultp = &result;
  memcpy(((void **)codenext)++, (void *) &resultp, sizeof(void *));
  }
#else
  *((void **) codenext) = (void *) &result; /*phil*/
  codenext = (CODEPTR) (void **) ((void **)codenext +1);
#endif
/*  printf("%x\n",codenext);*/
#else
  *codenext++ = OPPUSHPFLOAT;	/* Should change to double */
  *((void **) codenext)++ = (void *) &result;
#endif
  retcode = S_SUCCESS;
  if(thBookaTest(line,&codehead,&codenext,&codelimit,&codelastop,0)!=S_SUCCESS) {
    fprintf(STDERR,"Failure interpreting expression |%s|\n",line);
    result = 0.0;
    retcode = S_FAILURE;
  } else {
    int exptype;
    CODE lastop;
#if 0
    printf("%x-%x=%d\n",codenext,codehead,codenext-codehead);
    {
      CODEPTR code;
      for(code=codehead;code < codenext; code++)
	if(code==codelastop) printf("* %x\n",*code);
	else printf("  %x\n",*code);
    }
#endif
    codenext = codelastop;
    exptype = *codenext++ & OPRESTYPEMASK;
    lastop = *codelastop & OPCODEMASK;
    if(lastop == OPPUSHPINT || lastop == OPPUSHINTP) {
      codenext = (CODEPTR) (DAINT **) ((DAINT **)codenext + 1);/*phil*/
    } else if(lastop == OPPUSHINT) {
      if(exptype == OPRDOUBLE) {
      codenext = (CODEPTR) (DADOUBLE **) ((DADOUBLE **)codenext + 1);/*phil*/
      } else {			/* Assume ints, floats have size */
        codenext = (CODEPTR) (DAINT *) ((DAINT *)codenext + 1);/*phil*/
      }
    }
#ifdef RDOUBLE
    *codenext++ = OPEQUAL | 0x202 | (exptype<<4);
#else
    *codenext++ = OPEQUAL | 0x101 | (exptype<<4);
#endif
#ifdef RDOUBLE
    *codenext++ = OPEOL | (OPRDOUBLE<<4);
#else
    *codenext++ = OPEOL;
#endif
    if(thExecuteCode("IMMED",codehead,codenext)!=S_SUCCESS){
      fprintf(STDERR,"Failure evaluating expression |%s|\n",line);
      result = 0.0;
      retcode = S_FAILURE;
    }
  }
/*  printf("%f\n",result);*/
  free(codehead);
  if(d) *d = result;
  if(i) {
    if(result>=INT_MAX || result <=-INT_MAX) {
      if(retcode==S_SUCCESS)
	retcode=S_INTOVF;
    } else {
      *i = floatToLong(result);
    }
  }
  return(retcode);
}
  
thStatus thBookaTest(char *line, CODEPTR *codeheadp, CODEPTR *codenextp,
		     CODEPTR *codelimitp, CODEPTR *codelastop, daVarStructList **vlisthead)
/* if expflag != 0, still treat as an expression even if there is no
   equal sign in the line.
   Return codes:
     S_SUCCESS = Line OK
     S_FAILURE = Line not executable
*/
{
  /*  int type;*/
  char *args[20];
  int nargs;
  thTokenType toktyp;
  daVarStruct var, *varp;
  thTestType test_type;
  int forcefloat;
  int iarg;
  char *token;
  CODEPTR codenext;
  int index;			/* Used for index into arrays */
  thStatus status;
  int expflag;

  if(codelastop) expflag = 1; else expflag = 0;
  status = S_SUCCESS;
  if(*codenextp + 2*strlen(line) > *codelimitp) {
    CODEPTR src,dst,newhead;
    int newsize;
/*    printf("Increasing the size of the code stack from %d ",
	   *codelimitp-*codeheadp);*/
    src = *codeheadp;
    newsize = max((*codelimitp-*codeheadp)+CODEGROWSIZE
		  ,(*codenextp-*codeheadp)+2*strlen(line));
    newhead = dst = (CODEPTR) malloc(sizeof(CODE)*newsize);
    while(src < *codenextp) *dst++ = *src++;
    if(*codeheadp) free(*codeheadp);
    *codelimitp = newhead + newsize;
    *codeheadp = newhead;
    *codenextp = *codenextp + (dst - src);
    
    /*printf("to %d, using %d\n",*codelimitp-*codeheadp,*codenextp - *codeheadp);*/
  }
  codenext = *codenextp;

/*  printf("Booking \"%s\"\n",line);*/
  if(strchr(line,'=')||expflag) {
    char *linep;
    int TOKEN,TOKCOMP;
    char *tokstr; CODE tokval;
    void *tokptr;
    CODE *osp, *tsp, opcode;
    CODE rightoptype,leftoptype,resultoptype;

    osp = opstack;		/* Stack of pending operators */
    *osp = '\0';

    tsp = typstack;		/* Stack of Current result type */
				/* Like the stack in the executor but only */
				/* contains the data types */
    linep = line;
    do {
      /* Get tokens until there are no more (last token will be OPEOL) */
      linep = thGetTok(linep,&TOKEN, &tokstr, &tokval, &tokptr, expflag, vlisthead);
      if(tokstr) {		/* Operand */
/*	printf("Operand %s |",tokstr);*/
	if(codelastop) *codelastop = codenext; /* HACK for thImmed: Save ptr to last operator */
	if(TOKEN) {
	  if(tokptr == 0) {	/* Value operand - 4 bytes */
	    *codenext++ = TOKEN;	/* String not put on stack at moment */
	    *codenext++ = tokval;
	  } else {		/* Pointer operand - maybe 8 bytes */
	    *codenext++ = TOKEN;
#ifdef USEMEMCPY
	    memcpy(((void **)codenext)++,&tokptr,sizeof(void *));
#else
	    *(void **)codenext = tokptr;/*phil*/
            codenext = (CODEPTR) (void **) ((void **)codenext +1);
#endif
	  }
	  /* If TOKEN is push function, then tokval is an index into a list of
	     functions.  We put this index on tsp instead of the result type. */
	  if(TOKEN==OPPUSHFUNCTION) {
	    *tsp++ = tokval;
	  } else {
	    *tsp++ = TOKEN & OPRESTYPEMASK;
	  }
	} else {
	  fprintf(STDERR,"thTest: Unregistered variable %s\n",tokstr);
          status = S_TH_UNREG;
	  *codenext++ = OPPUSHINT;
	  *codenext++ = 0;
	  *tsp++ = OPPUSHINT & OPRESTYPEMASK;
	}
      } else {			/* Operator */
	switch(TOKEN)
	  {
	  case 0:
	    fprintf(STDERR,"thTest: Bad token\n");
	    return(S_FAILURE);
	    break;
	  case OPLP:
	    *++osp = TOKEN;
	    break;
	  default:
/*	    printf("OSP:");
	    {CODE *sp; for(sp=opstack;sp<=osp;sp++)
	       printf("%x:",*sp);}
	    printf("\n");
*/
	    /* Generate code for all operators of equal or higher precedence
	       that are pending on the operator stack. */
	    if((TOKEN & OPGROUPMASK) == OPLINDEXGROUP)
	      TOKCOMP = 0xFFFFFFF; /* Nothing higher in precedence */
	    else
	      TOKCOMP = TOKEN & OPPRECMASK;
	    while((*osp & OPPRECMASK) >= TOKCOMP){
/*	      if((*osp & OPPRECMASK) == OPLINDEX){*/
	      if((*osp & OPGROUPMASK) == OPLINDEXGROUP){
		if(TOKEN == OPRP) {
		  if(*osp == OPLFARG) TOKEN = OPRFARG;
		  else TOKEN = OPRINDEX; /* Break from case */
		}
		TOKCOMP = 0xFFFFFFF; /* Terminate osp rundown */
	      }
	      rightoptype = *--tsp;
	      leftoptype = ((*osp & OPPRECMASK) == OPUNARY) ? 0 : (*--tsp);
	      /* If the Operator is "evaluate function", we need to find out
		 what the function is so that we can get the correct
		 result type.  leftoptype should be an index into
		 "intrinsic_functions".  We can use that and rightoptype
		 to look up the resulttype. */
	      if(*osp==OPLFARG) {
		resultoptype = 
		  intrinsic_functions[leftoptype].result[rightoptype];
	      } else {
		resultoptype = thGetResultType(*osp,leftoptype,rightoptype);
	      }
	      opcode = *osp--;
	      opcode |= (leftoptype << 8) | (rightoptype << 4)
		| resultoptype;
	      if(codelastop) if((opcode&&OPCODEMASK) !=OPEOL) *codelastop = codenext; /* HACK for thImmed: Save ptr to last operator */
	      *codenext++ = opcode;
	      *tsp++ = resultoptype; /* Keep a rpn stack of the data type */
	    }
	    if(TOKEN == OPRINDEX || TOKEN == OPRFARG) break; /* No clean up needed */

	    if(TOKEN == OPRP) {
	      if(*osp == OPLP) osp--; /* ) removes matching ( */
	      else {
		fprintf(STDERR,"Right paren not matched by left\n");
		return(S_FAILURE);
	      }
	    } else if(TOKEN == OPEOL || TOKEN == OPCOMMA) {
	      if(codelastop) if(TOKEN==OPCOMMA) *codelastop = codenext; /* HACK for thImmed: Save ptr to last operator */
	      *codenext++ = TOKEN | (*--tsp) << 4; /* Leave type in Right type field */
	    } else {
	      *++osp = TOKEN;
	    }
	    break;
	  }
      }
      /* Token processed */
    } while (linep);
/* Check that stacks are OK.  Need to add some clean up of allocated memory. */
    if(tsp != typstack) {
      fprintf(STDERR,"%d items left on type stack\n",tsp-typstack);
      return(S_FAILURE);
    }
    if(osp != opstack) {
      fprintf(STDERR,"%d items left on operand stack\n",osp-opstack);
      return(S_FAILURE);
    }
  } else {			/* Old style test lines */
    int i;
    nargs = thCommas(line,args);
    for(i=0;i<nargs;i++){
      args[i] = thSpaceStrip(args[i]); /* Remove all space from the argument */
    }
    
    if(nargs <= 1) return(S_FAILURE);
    
    {				/* Interpret the test type. */
      
      for(test_type=0;test_type<tBAD;test_type++){
	if(testCodes[test_type][0] == toupper(args[1][0]) &&
	   testCodes[test_type][1] == toupper(args[1][1])) break;
      }
      if(test_type == tBAD) return(S_FAILURE);
      /*    printf("%s\n",testCodes[test_type]);*/
    }
    if(test_type == tGATE || test_type == tEQ) {
      forcefloat = 1;
    } else forcefloat = 0;
    for(iarg=2;iarg<nargs;iarg++){
      DAFLOAT f;		/* Should do double  here */
      token = args[iarg];
      toktyp = thIDToken(token);
      switch((toktyp = thIDToken(token)))
	{
	case TOKINT:
	  *codenext++ = PUSHI;
	  if(forcefloat) {
	    f = atof(token);
	    *codenext++ = *(DAINT *)&f;
	  } else {
	    DAINT i;
	    /* Used to be %li and %ld, but that makes 8 byte result
	       stuffed into 4 byte i */
	    if(token[0] == '0' && (token[1] == 'x' || token[1] == 'X')) {
	      sscanf(token,"%i",&i); /* Treat as Hex */
	    } else {
	      sscanf(token,"%d",&i); /* Treat as decimal */
	    }
	    *codenext++ = i;
	  }
	  break;
	case TOKFLOAT:		/* Should Do all floats as doubles */
	  *codenext++ = PUSHI;
	  if(forcefloat) {
	    f = atof(token);
	    *codenext++ = *(DAINT *)&f;
	  } else {
	    *codenext++ = (DAINT) floatToLong(atof(token));
	  }
	  break;
	case TOKARRAY:
	case TOKVAR:
	  {
	    char *p; int index; char leftp;
	    if(toktyp == TOKARRAY) {
	      p = thTokenArray(token,&index);
	      leftp = *p; *p = 0;	/* Save ( or [ then null terminate */
	    } else
	      index = 0;
	    if(daVarLookup(token,&var)!=S_SUCCESS) {
	      fprintf(STDERR,"(thTest) %s not registered\n",token);
	      *codenext++ = PUSHI;
	      if(forcefloat) {
		f = 0.0;
		*codenext++ = *(DAINT *)&f;
	      } else
		*codenext++ = 0;
	    } else {
	      if(forcefloat)
		if(var.type == DAVARINT)
		  *codenext++ = PUSHITOFS; /* Push converting to float and skip */
		else if(var.type == DAVARFLOAT)
		  *codenext++ = PUSHS;
		else
		  *codenext++ = PUSHI; /* Push immediate */
	      else
		if(var.type == DAVARINT)
		  *codenext++ = PUSHS; /* Push and skip */
		else if(var.type == DAVARFLOAT)
		  *codenext++ = PUSHFTOIS;
		else
		  *codenext++ = PUSHI; /* Push immediate */
	      if(toktyp == TOKARRAY)
		*p = leftp;
	      if(var.type == DAVARINT || var.type == DAVARFLOAT) {
		*(void **)codenext = ((DAINT *) var.varptr+index);/*phil*/
                codenext = (CODEPTR) (void **) ((void **)codenext + 1);
		*((void **)codenext) = (void *) malloc(sizeof(token)+1);
		strcpy((char *) *(void **)codenext,token);
                codenext = (CODEPTR) (void **) ((void **)codenext + 1);
	      } else {
		if(forcefloat) {
		  f = 0.0;
		  *codenext++ = *(DAINT *)&f;
		} else
		  *codenext++ = 0;
	      }
	    }
	  }
	  break;
	}
    }
    *codenext++ = test_type;	/* Operation to do on pushed args */
    *codenext++ = nargs-2;	/* # of args pushed on stack for this op */
    
    /* Now push test result on stack */
    *codenext++ = POPS;
    
    token = args[0];
    toktyp = thIDToken(token);
    index = 0;
    switch((toktyp = thIDToken(token)))
      {
      case TOKINT:
      case TOKFLOAT:
	fprintf(STDERR,"(thTest) Test result must be a variable name\n");
	return(S_FAILURE);	/* No test is added to program */
      case TOKARRAY:
	/* First check if variable with index has been already registered
	   perhaps from a previous booking of histograms */
	if(daVarLookup(token,&var) != S_SUCCESS){
	  char *p; char leftp;
	  p = thTokenArray(token,&index);
	  leftp = *p; *p = 0;	/* Save ( or [ then null terminate */
	  if(daVarLookup(token,&var) != S_SUCCESS){
	    fprintf(STDERR,
		    "(thTest) Arrays %s must be registered\n",token);
	    return(S_FAILURE);
	  }
	  *p = leftp;		/* Restore the left ( or [ */
	  if(index >= var.size) {
	    fprintf(STDERR,
		    "(thTest) Array size for %s exceeded\n",token);
	    return(S_FAILURE);
	  }
	  if(var.type != DAVARINT) {
	    fprintf(STDERR,
		    "(thTest) Array %s must be of integer*4\n",token);
	    return(S_FAILURE);
	  }
	  var.varptr = (DAINT *) var.varptr + index;
	  var.name = token;
	  var.opaque = 0;
	}
	var.title = token;	/* Eventually be the input line */
	break;
      case TOKVAR:
	if(daVarLookup(token,&var)!=S_SUCCESS) {
	  var.name = token;
	  var.varptr = (void *) malloc(sizeof(DAINT));
	  var.opaque = 0;
	  var.rhook = 0;
	  var.whook = 0;
	  var.type = DAVARINT;
	  var.flag = DAVAR_READONLY | DAVAR_REPOINTOK;
/* Do I need to set the size to 1 here??? */
	}
	var.title = token;
	break;
      }
    daVarRegister((int) 0, &var);	/* Create or replace variable */
    *(void **)codenext = ((DAINT *) var.varptr);/*phil*/
    codenext = (CODEPTR) (void **) ((void **)codenext + 1);
    /* Save the token string for future reference */
    *((void **)codenext) = ((void *) malloc(strlen(token)+1));
    strcpy((char *) *(void **)codenext,token);
    codenext = (CODEPTR) (void **) ((void **)codenext + 1);
  }
  *codenextp = codenext;
  return(status);
  
}
int thevalchk_(char *A1,unsigned C1)
/* Check if an expression is valid.  Return's zero if valid */
{
  int A0;
  char *B1;
  thStatus status;

  status = thEvalImed((!*(int *)A1)?0:memchr(A1,'\0',C1)?A1:
		      (memcpy(B1=malloc(C1+1),A1,C1),B1[C1]='\0'
		       ,kill_trailing(B1,' ')),0,0);
  if(B1) free(B1);
  return(status);
}

int itheval_(char *A1,unsigned C1)
{
  int A0;
  char *B1;
  DAINT i;
  double d;
  thStatus status;

  status = thEvalImed((!*(int *)A1)?0:memchr(A1,'\0',C1)?A1:
		      (memcpy(B1=malloc(C1+1),A1,C1),B1[C1]='\0'
		       ,kill_trailing(B1,' ')),0,&i);
  if(B1) free(B1);
  return i;
}
double dtheval_(char *A1,unsigned C1)
{
  char *B1;
  double d;
  thStatus status;

  status = thEvalImed((!*(int *)A1)?0:memchr(A1,'\0',C1)?A1:
		      (memcpy(B1=malloc(C1+1),A1,C1),B1[C1]='\0'
		       ,kill_trailing(B1,' ')),&d,0);
  if(B1) free(B1);
  return d;
}
float ftheval_(char *A1,unsigned C1)
{
  char *B1;
  DAINT i;
  double d;
  float f;
  thStatus status;

  status = thEvalImed((!*(int *)A1)?0:memchr(A1,'\0',C1)?A1:
		      (memcpy(B1=malloc(C1+1),A1,C1),B1[C1]='\0'
		       ,kill_trailing(B1,' ')),&d,0);
  if(B1) free(B1);
  f = d;
  return f;
}

daVarStatus thTestRHandler(char *name, daVarStruct *varclass, any *retval)
/* The default Read handler */
{
  daVarStruct *varp;
  char *attribute;
  daVarStatus status;
  int index;

  status = daVarAttributeFind(name, varclass, &varp, &attribute, &index);
  status = daVarRegRatr(varp, attribute, index, retval);
  if(status == S_SUCCESS) {
    if(strcasecmp(attribute,DAVAR_RATR) == 0){
      retval->any_u.s = realloc(retval->any_u.s,strlen(retval->any_u.s)
				+strlen(TH_SCALER) + 2);
      strcat(retval->any_u.s,TH_SCALER);
      strcat(retval->any_u.s,"\n");
    }
  } else {
    if(strcasecmp(attribute,TH_SCALER) == 0){
      int i;
      if(varp->opaque){
	retval->valtype = DAVARINT_RPC;
	if(index == DAVAR_NOINDEX) {
	  retval->any_u.i.i_len = varp->size;
	  retval->any_u.i.i_val = (int *) malloc(varp->size*sizeof(int));
	  for(i=0;i<varp->size;i++) {
	    retval->any_u.i.i_val[i] = ((DAINT *)varp->opaque)[i];
	  }
	} else {
	  retval->any_u.i.i_len = 1;
	  retval->any_u.i.i_val = (int *) malloc(sizeof(int));
	  retval->any_u.i.i_val[0] = ((DAINT *)varp->opaque)[index];
	}
      } else {
	retval->valtype = DAVARERROR_RPC;
	retval->any_u.error = S_SUCCESS;
      }
    }
  }
  /* A special handler would check more attributes if status != SUCCESS */
  return(status);
}
