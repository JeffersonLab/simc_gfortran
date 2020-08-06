/*-----------------------------------------------------------------------------
 * Copyright (c) 1995,1996 Southeastern Universities Research Association,
 *                         Continuous Electron Beam Accelerator Facility
 *
 * This software was developed under a United States Government license
 * described in the NOTICE file included as part of this distribution.
 *
 * Stephen A. Wood, 12000 Jefferson Ave., Newport News, VA 23606
 * Email: saw@cebaf.gov  Tel: (804) 249-7367  Fax: (804) 249-5800
 *-----------------------------------------------------------------------------
 * 
 * Description:
 *  Code for dealing with blocks in groups.
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: thGroup.c,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.5.6.1  2008/09/25 00:54:05  jones
 *   Updated for running on Fedora 8 with gfortran
 *
 *   Revision 1.6  2008/09/25 00:01:29  jones
 *   Updated to run with gfortran compiler
 *
 *   Revision 1.5.8.1  2007/09/10 21:32:47  pcarter
 *   Implemented changes to allow compilation on RHEL 3,4,5 and MacOSX
 *
 *   Revision 1.5  2004/07/08 20:05:37  saw
 *   Use dummy fortran ctp tree routines when ROOTSYS not defined.
 *
 *   Revision 1.4  2004/07/07 18:15:27  saw
 *   Consistenly use thtreeexeg
 *
 *   Revision 1.3  2004/07/02 20:11:07  saw
 *   Make fortran tree group function names sane
 *
 *   Revision 1.2  1999/11/04 20:34:05  saw
 *   Alpha compatibility.
 *   New RPC call needed for root event display.
 *   Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *   Revision 1.1  1998/12/07 22:11:12  saw
 *   Initial setup
 *
 *   Revision 1.2  1996/07/31 20:31:52  saw
 *   Book blocks in the order they appear in the input files.
 *
 *   Revision 1.1  1996/01/30 15:35:16  saw
 *   Initial revision
 *
 */
#include <stdio.h>
#include "daVar.h"
#include "th.h"
#include "thInternal.h"
#include "thUtils.h"
#include "cfortran.h"
#include "thGroup.h"

thStatus thBookGroup(char *name);
/* Fortran Interface */
FCALLSCFUN0(INT,thBook,THBOOK,thbook)
FCALLSCFUN1(INT,thBookGroup,THGBOOK,thgbook,STRING)
FCALLSCFUN1(INT,thExecuteGroup,THGEXE,thgexe,STRING)
FCALLSCFUN1(INT,thClearGroup,THGCLR,thgclr,STRING)
FCALLSCFUN1(INT,thClearScalersGroup,THGCLS,thgcls,STRING)
FCALLSCFUN1(INT,thIncrementScalersGroup,THGINS,thgins,STRING)

struct thHook {
  char *type;
  thStatus (*book)();		/* Hooks to the appropriate routines */
  thStatus (*execute)();
  thStatus (*clear)();
  thStatus (*clearScalers)();
  thStatus (*incrementScalers)();
  thStatus (*ctpwrite)();
  thStatus (*ctpclose)();
};
typedef struct thHook thHook;
/* V indicates these calls take a variable pointer, not a name */
/* Eventually calls without V may go away, and these will all be renamed */
/* How many of these c routines are "advertised */
thHook thHooks[] = {
  {PARMSTR   ,thLoadParameters,0,0,0,0,0,0},
  /*{GETHITSTR ,thBookGethits,thExecuteaGethitBlock,0,0,0,0,0},*/
  {TESTSTR   ,thBookTests,thExecuteTestsV,thClearTestFlagsV,thClearTestScalersV
     ,thIncTestScalersV,0,0},
  /*{HISTSTR,thBookHists,thExecuteHistsV,thClearHistsV,0,0,0,0},*/
  /*{UHISTSTR,thBookHists,0,0,0,0,0,0},*/
  {TREESTR,thBookTree,thFillTreeV,thClearTreeV,0,0,thWriteTreeV,thCloseTreeV},
  {REPORTSTR,thBookReports,0,0,0,0,0,0},
  {0,0,0,0,0,0,0,0}};

int thGroupClassesSet=0;

thStatus thSetGroupClasses();

void thInitGroupOpaque(char *name, thGroupOpaque *opqptr)
{
  daVarStatus status;
  daVarStruct *varclass;
  thGroupOpaque *classopqptr;

  if(!thGroupClassesSet) {
    status = thSetGroupClasses();
    thGroupClassesSet = 1;
  }
  status = daVarClassFind(name,&varclass);
  classopqptr = varclass->opaque;
  opqptr->blocklist = 0;
  opqptr->type = classopqptr->type;
  opqptr->book = classopqptr->book;
  opqptr->execute = classopqptr->execute;
  opqptr->clear = classopqptr->clear;
  opqptr->clearScalers = classopqptr->clearScalers;
  opqptr->incrementScalers = classopqptr->incrementScalers;
  opqptr->ctpwrite = classopqptr->ctpwrite;
  opqptr->ctpclose = classopqptr->ctpclose;
  /* Find out what group type we are copy the information out of it's
     opaque block */
}
thStatus thSetGroupClasses()
{
/*  thBookList *booklist;
  thBookList **booknext;*/
  thGroupOpaque *opqptr;
  char *classname;
  daVarStruct var;
  int i;

  for(i=0;thHooks[i].type;i++){
    /* Need to decide later to do automatic booking by groups 
     What do I mean by that ^^^ ?*/
    classname = (char *) malloc(strlen(GROUPSTR)+strlen(thHooks[i].type)+2);
    strcpy(classname,GROUPSTR);
    strcat(classname,".");
    strcat(classname,thHooks[i].type);
    
    var.name = classname;
    var.title = 0;/*""*/
    var.type = DAVARINT;
    var.varptr = (DAINT *) malloc(sizeof(DAINT));
    *((DAINT *) var.varptr) = i; /* Booking order */
    var.size = 1;
    var.whook = 0;		/* Need handlers? */
    var.rhook = 0;
    var.flag = DAVAR_READWRITE | DAVAR_REPOINTOK;
//    var.opaque = (void *)opqptr = (thGroupOpaque *) malloc(sizeof(thGroupOpaque));/*phil*/
    var.opaque = (void *) (opqptr = (thGroupOpaque *) (void *)
            ((thGroupOpaque *) malloc(sizeof(thGroupOpaque))));
    opqptr->blocklist = 0;
    opqptr->type = (char *) malloc(strlen(thHooks[i].type) + 1);
    strcpy(opqptr->type,thHooks[i].type);
    opqptr->book = thHooks[i].book;
    opqptr->execute = thHooks[i].execute;
    opqptr->clear = thHooks[i].clear;
    opqptr->clearScalers = thHooks[i].clearScalers;
    opqptr->incrementScalers = thHooks[i].incrementScalers;
    opqptr->ctpwrite = thHooks[i].ctpwrite;
    opqptr->ctpclose = thHooks[i].ctpclose;
    if(daVarRegister((int) 0, &var) == S_FAILURE){
      fprintf(STDERR,"Failed to register %s\n",var.name);
      return(S_FAILURE);
    }
/*    printf("Registered %s\n",var.name);*/
  }
  return(S_SUCCESS);
}
thStatus thBookGroup(char *group)
{
  daVarStruct *varp, *bvar;
  daVarStructList *blocklist;
  thStatus stat;
  thStatus (*hook)();

  if(daVarLookupP(group, &varp) != S_SUCCESS){
/*    fprintf(STDERR,"Failed to find %s\n",group);*/
    return(S_FAILURE);
  }
  hook = ((thGroupOpaque *)varp->opaque)->book;
  blocklist = ((thGroupOpaque *)varp->opaque)->blocklist;
  while(blocklist) {
    bvar = blocklist->varp;
    if(!bvar->varptr) {	/* Only book those not already booked */
      bvar->varptr = (DAINT *) malloc(sizeof(DAINT));
      *((DAINT *) bvar->varptr) = 0; /* Initializec execution counter */
      if(hook) {
/*	printf("Booking %s\n",bvar->name);*/
	stat = (hook)(bvar);
      } else {
/*	printf("No booking required for %s\n",bvar->name);*/
      }
    }
    blocklist = blocklist->next;
  }
  return(S_SUCCESS);
}

#define MAKEGSUB(SUBNAME,QSUBNAME,ELEMENT) \
thStatus SUBNAME(char *group)\
{\
  daVarStruct *varp, *bvar;\
  daVarStructList *blocklist;\
  thStatus stat;\
  thStatus (*hook)();\
\
  if(daVarLookupP(group, &varp) != S_SUCCESS){\
/*    fprintf(STDERR,"(%s) Failed to find %s\n",QSUBNAME,group);*/\
    return(S_SUCCESS);\
  }\
  hook = ((thGroupOpaque *)varp->opaque)->ELEMENT;\
  blocklist = ((thGroupOpaque *)varp->opaque)->blocklist;\
  if(hook) {\
    while(blocklist) {\
      bvar = blocklist->varp;\
      stat = (hook)(bvar);\
      blocklist = blocklist->next;\
    }\
  }\
  return(S_SUCCESS);\
}
MAKEGSUB(thExecuteGroup,"thExecuteGroup",execute)
MAKEGSUB(thClearGroup,"thClearGroup",clear)
MAKEGSUB(thClearScalersGroup,"thClearScalersGroup",clearScalers)
MAKEGSUB(thIncrementScalersGroup,"thIncrementScalersGroup",incrementScalers)
MAKEGSUB(thWriteGroup,"thWriteGroup",ctpwrite)
MAKEGSUB(thCloseGroup,"thCloseGroup",ctpclose)

thStatus thBook()
/* Book all the parameter, tests and histograms and anything else that is
   in the booking order list.  Books groups in alphabetical order.  Within
   groups booking is done in the order that the blocks appear in the CTP
   files.

1/30/96 New behaviour.  Just book the "all" group for each class.  That way
booking order will be precisely as appears in the CTP files.
*/
{
  int i;
  thStatus status;

  status = S_FAILURE;	/* Return failure if everything fails to book */

  for(i=0;thHooks[i].type;i++){
    char *prefix; char **glist0; char **glist; int count;

#if 0
    prefix = (char *) malloc(strlen(GROUPSTR) + strlen(thHooks[i].type) + 3);
#else
    prefix = (char *) malloc(strlen(GROUPSTR) + strlen(thHooks[i].type) 
			     + strlen(ALLGRPSTR) + 3);
#endif
    strcpy(prefix,GROUPSTR);
    strcat(prefix,".");
    strcat(prefix,thHooks[i].type);
    strcat(prefix,".");
    
#if 0
/* Make a list of all groups in this class */
    
    daVarList(prefix,&glist0,&count);
    glist = glist0;
    while(count-- > 0){
/*      printf("Booking group %s\n",*glist);*/
      if(thBookGroup(*glist) == S_SUCCESS){
	status = S_SUCCESS;
      } else {
	fprintf(STDERR,"Failed to book %s\n",*glist);
      }
      glist++;
    }
#else
    strcat(prefix,ALLGRPSTR);
/*    printf("Booking %s\n",prefix);*/
    if(thBookGroup(prefix) == S_SUCCESS){
	status = S_SUCCESS;
    } else {
/*      fprintf(STDERR,"Failed to book %s\n",prefix);*/
    }
#endif
    free(prefix);
  }
  return(status);
}

#define MAKEFSUB(SUBNAME,TYPESTR,CSUBNAME) \
int SUBNAME(char *A1,unsigned C1)\
{\
  int A0;\
  char *B1=0;\
  int newsize;\
  static char *full_name=0;\
  static int full_name_size;\
\
  newsize = strlen(GROUPSTR) + strlen(TYPESTR) + C1 + 3;\
  if(!full_name) {\
    full_name_size = newsize;\
    full_name = (char *) malloc(full_name_size);\
  } else {\
    if(newsize > full_name_size) {\
      full_name = realloc(full_name,newsize);\
    }\
  }\
  strcpy(full_name,GROUPSTR);\
  strcat(full_name,".");\
  strcat(full_name,TYPESTR);\
  strcat(full_name,".");\
  strcat(full_name,(!*(int *)A1)?0:memchr(A1,'\0',C1)?A1:\
		      (memcpy(B1=malloc(C1+1),A1,C1),B1[C1]='\0'\
		       ,kill_trailing(B1,' ')));\
  A0 = CSUBNAME(full_name);\
  if(B1) free(B1);\
  return(A0);\
}

#ifdef AbsoftUNIXFortran
MAKEFSUB(thtstexeg,TESTSTR,thExecuteGroup)
MAKEFSUB(thtstclrg,TESTSTR,thClearGroup)
MAKEFSUB(thtstclsg,TESTSTR,thClearScalersGroup)
MAKEFSUB(thtstinsg,TESTSTR,thIncrementScalersGroup)
MAKEFSUB(thhstexeg,HISTSTR,thExecuteGroup)
MAKEFSUB(thgethitg,GETHITSTR,thExecuteGroup)
MAKEFSUB(thtreeexeg,TREESTR,thExecuteGroup)
MAKEFSUB(thtreecloseg,TREESTR,thCloseGroup)
MAKEFSUB(thtreewriteg,TREESTR,thWriteGroup)
#else
MAKEFSUB(thtstexeg_,TESTSTR,thExecuteGroup)
MAKEFSUB(thtstclrg_,TESTSTR,thClearGroup)
MAKEFSUB(thtstclsg_,TESTSTR,thClearScalersGroup)
MAKEFSUB(thtstinsg_,TESTSTR,thIncrementScalersGroup)
MAKEFSUB(thhstexeg_,HISTSTR,thExecuteGroup)
MAKEFSUB(thgethitg_,GETHITSTR,thExecuteGroup)
MAKEFSUB(thtreeexeg_,TREESTR,thExecuteGroup)
MAKEFSUB(thtreecloseg_,TREESTR,thCloseGroup)
MAKEFSUB(thtreewriteg_,TREESTR,thWriteGroup)
#endif

/*
#define MAKEFCALL(SUBNAME,CLASS) \

MAKEFCALL(thhstexeg,hist)
*/
