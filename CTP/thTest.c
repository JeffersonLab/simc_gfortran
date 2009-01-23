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
 *  Book tests, execute or operate on tests.
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: thTest.c,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.2  1999/11/04 20:34:07  saw
 *   Alpha compatibility.
 *   New RPC call needed for root event display.
 *   Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *   Revision 1.1  1998/12/07 22:11:13  saw
 *   Initial setup
 *
 *   Revision 1.7  1996/07/31 20:30:17  saw
 *   (SAW) List line number of booking errors.  Add code in support of groups.
 *
 *	  Revision 1.6  1995/04/10  15:42:35  saw
 *	  Handle ctp file registration (#int, #real, ...)
 *	  No defined test blocks is not an error in thTestExecute
 *
 *	  Revision 1.5  1995/01/09  15:55:25  saw
 *	  On fprintf, indicate block type and well as name.
 *
 *	  Revision 1.4  1994/08/31  13:05:02  saw
 *	  Add code for WALK_CLEAR_FLAGS to thWalkTree
 *
 *	  Revision 1.3  1994/06/03  18:51:22  saw
 *	  Replace stderr with STDERR
 *
 *	  Revision 1.2  1993/05/11  17:56:37  saw
 *	  Fix header
 *
 */

/* 
   Need to add     means of accessing scalers.  Two ways.  1 Register the
   scaler array, 2 make a call with the name and a pointer to an array.
*/

/*An argument is a variable name, an array, or a number.  Numbers are not
allowed for test result.  Arrays start at 0 if []'s are used and start
at 1 if ()'s are used.  Arrays may only be used for test results if they
are already registered by the analyzer.  (May add option to declare them
in the test package.)*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "daVar.h"
#include "th.h"
#include "thInternal.h"
#include "thUtils.h"
#include "thTestParse.h"
#include "cfortran.h"
#define max(a,b) 		(a<b ? b : a)
#define min(a,b) 		(a>b ? b : a)

#define MAXLINELENGTH 512

struct thTBOpaque {
  CODEPTR code;
  CODEPTR end;
  daVarStructList *vlisthead;
};
typedef struct thTBOpaque thTBOpaque;

struct thTBlockList {
  struct thTBlockList *next;
  char *blockname;		/* Block name with out the "block.test"  */
  daVarStruct *var;		/* Pointer to variable that describes block */
/* Block execution counter, pointer to list of test names and results */
};
typedef struct thTBlockList thTBlockList;

thTBlockList *thTBlockListP = NULL;

/* End of newstuff */

thStatus thBookTests(daVarStruct *var)
{
  char *lines,*eol;
  int line_count;
  thTBOpaque *opqptr;
  CODEPTR codehead, codenext, codelimit;

  /* Need to zero out the count.  Make sure variable for count is
     created */

/*  printf("Booking tests %s\n",var->name);*/

  if(var->opaque == 0) {
    opqptr = (thTBOpaque *) malloc(sizeof(thTBOpaque));
    var->opaque = (void *) opqptr;
    codehead = codenext = (CODEPTR) malloc(sizeof(CODE)*CODESTARTSIZE);
    codelimit =  codehead + CODESTARTSIZE;
    opqptr->vlisthead = 0;
  } else {
    daVarStructList *next,*this;

    opqptr = (thTBOpaque *) var->opaque;
    codehead = codenext = opqptr->code;
    codelimit = opqptr->end;
    this = opqptr->vlisthead;
    while(this) {
      next = this->next;
      free(this);
      this = next;
    }
    opqptr->vlisthead = 0;
  }

  lines = var->title;
  line_count = 0;
  while(*lines){
    char *lcopy;

    line_count++;
    eol = strchr(lines,'\n');
    if(!eol) {
      fprintf(STDERR,"L %d: Last line of test block %s has no newline\n"
	      ,line_count,var->name);
      break;
    }
/*    {
      *eol = 0;
      printf("L %d:%s\n",line_count,lines);
      *eol = '\n';
    }*/
      
/*    printf("Next char = %d\n",*(eol+1));*/
    if(*(eol+1)=='\0'){		/* This is the last line */
      if(strcasestr(lines,ENDSTR) == 0) {
	fprintf(STDERR,"L %d: Last line of test block %s is not an END\n"
		,line_count,var->name);
/*	fprintf(STDERR,"%s",var->title);*/
      }
      break;
    }
    if(line_count == 1)
      if(strcasestr(lines,BEGINSTR) != 0){
/*	printf("Is a begin\n");*/
	lines = eol + 1;
	continue;
      } else
	fprintf(STDERR,"First line of test block %s is not a BEGIN\n",var->name);
    /* Ready to book the line, Add continuation lines later */
    lcopy = (char *) malloc(eol-lines+1);
    strncpy(lcopy,lines,(eol-lines));
    *(lcopy + (eol-lines)) = '\0';
/*    printf("Passing|%s|\n",lcopy);*/
    if(!thSpecial(lcopy,TESTSTR)) {
      if(!thCleanLine(lcopy)){
	if(thBookaTest(lcopy,&codehead,&codenext,&codelimit,0,
		       &(opqptr->vlisthead))==S_SUCCESS) {
	  /*	{
		daVarStructList *walk;
		walk = opqptr->vlisthead;
		while(walk){
		walk = walk->next;
		}*/
	  /*      {
		  CODEPTR code;
		  for(code=codehead;code < codenext; code++)
		  printf("  %x\n",*code);
		  }*/
	} else {
	  fprintf(STDERR,"(%s): Test booking error in line %d\n",var->name,line_count);
	}
      }
    }
    free(lcopy);
    lines = eol+1;
  }

  {				/*  */
    CODEPTR src,dst;
    dst = opqptr->code = (CODEPTR) malloc(sizeof(CODE)*(codenext-codehead));
    src = codehead;
    while(src < codenext) *dst++ = *src++;
    opqptr->end = dst;
    free(codehead);
/*    printf("%x %x\n",opqptr->code,opqptr->end);*/
  }
  /* Update internal table of test blocks.  */
  {
    thTBlockList *thisblock,*nextblock,**lastblockp;
    nextblock = thTBlockListP;
    lastblockp = &thTBlockListP;
    thisblock = thTBlockListP;
    while(thisblock){
      if((strcmp(thisblock->var->name,var->name)) == 0){
	/* Replacing a block with a new definition */
	fprintf(STDERR,"Replacing %s with new definition\n",var->name);
	if(thisblock->var != var){
	  fprintf(STDERR,"ERR: Same name, different var pointer\n");
	}
	break;
      }
      lastblockp = &thisblock->next;
      thisblock = thisblock->next;
    }
    if(!thisblock){		/* Create entry for New block */
      char *s;
      int i;
      *lastblockp = thisblock = (thTBlockList *) malloc(sizeof(thTBlockList));
      thisblock->var = var;
      thisblock->next = (thTBlockList *) NULL;
      /* Get the name without the block.test on it */
      s = var->name;		/* If name doesn't fit pattern, use whole */
      if(strcasestr(var->name,BLOCKSTR)==var->name){
	i = strlen(BLOCKSTR) + 1;
	if(strcasestr((var->name + i),TESTSTR)==(var->name + i)){
	  i += strlen(TESTSTR);
	  if(*(var->name + i) == '.'){
	    s += i + 1;
	  }
	}
      }
      thisblock->blockname = (char *) malloc(strlen(s) + 1);
      strcpy(thisblock->blockname,s);
    }
  }
  return(S_SUCCESS);
}

thStatus thClearTestFlagsV(daVarStruct *var)
{
  thTBOpaque *opqptr;
  thStatus status;
  daVarStructList *vlist;

  opqptr = (thTBOpaque *) var->opaque;
  vlist = opqptr->vlisthead;
  while(vlist){
    daVarStruct *varp;
    varp = vlist->varp;
    if(varp->type == DAVARINT) { /* Only clear when flag is integer */
      int i;
      for(i=0;i<varp->size;i++)
	((DAINT *) varp->varptr)[i] = 0;
    }
    vlist = vlist->next;
  }
}
thStatus thClearTestScalersV(daVarStruct *var)
{
  thTBOpaque *opqptr;
  thStatus status;
  daVarStructList *vlist;

  opqptr = (thTBOpaque *) var->opaque;
  vlist = opqptr->vlisthead;
  while(vlist){
    daVarStruct *varp;
    varp = vlist->varp;
    if(varp->type == DAVARINT) { /* Only clear when flag is integer */
      int i;
      for(i=0;i<varp->size;i++)
	((DAINT *) varp->opaque)[i] = 0;
    }
    vlist = vlist->next;
  }
}
thStatus thIncTestScalersV(daVarStruct *var)
{
  thTBOpaque *opqptr;
  thStatus status;
  daVarStructList *vlist;

  opqptr = (thTBOpaque *) var->opaque;
  vlist = opqptr->vlisthead;
  while(vlist){
    daVarStruct *varp;
    varp = vlist->varp;
    if(varp->type == DAVARINT) { /* Only clear when flag is integer */
      int i;
      for(i=0;i<varp->size;i++)
	((DAINT *) varp->opaque)[i] += (((DAINT *) varp->varptr)[i] != 0);

    }
    vlist = vlist->next;
  }
}
thStatus thExecuteTestsV(daVarStruct *var)
{
  thTBOpaque *opqptr;
  thStatus status;

  CODEPTR codehead, codenext, codelimit;

  opqptr = (thTBOpaque *) var->opaque;
  status = thExecuteCode(var->name,opqptr->code,opqptr->end);
  if(status == S_SUCCESS)
    (*((DAINT *)var->varptr))++; /* Increment block counter */
  return(status);
}

thStatus thWalkTree(char *block_name, WALKOP walkop)
{
  thTBlockList *thisblock,*nextblock,**lastblockp;
  int *result;
  thTBOpaque *opqptr;

  if(block_name)
    if(*block_name=='\0')
      block_name = 0;
    else
      lastblockp = &thTBlockListP; /* If no block specified, do default group */

  if(thTBlockListP == 0){
    return(S_SUCCESS);		/* No tests defined */
  } else {
    thisblock = thTBlockListP;
    while(thisblock){
      if(block_name)
	if(strcasecmp(block_name,thisblock->blockname)!=0) {
	  lastblockp = &thisblock->next;
	  thisblock = thisblock->next;
	  continue;
	}
      if(walkop == WALK_DISPLAY)
	fprintf(STDERR,"           /%s\n",thisblock->blockname);
      else {
	opqptr = thisblock->var->opaque;
	if(walkop == WALK_EXECUTE) {
	  thExecuteCode(thisblock->var->name,opqptr->code,opqptr->end);
	  (*((DAINT *)thisblock->var->varptr))++; /* Increment block counter */
	} else if(walkop == WALK_CLEAR_SCALERS 
		  || walkop == WALK_INCREMENT_SCALERS
		  || walkop == WALK_CLEAR_FLAGS) {
	  daVarStructList *vlist;
	  vlist = opqptr->vlisthead;
	  while(vlist) {
	    daVarStruct *varp;
	    varp = vlist->varp;
	    if(varp->type == DAVARINT) {
	      int i;
	      if(walkop == WALK_CLEAR_SCALERS) {
		if(varp->opaque)
		  for(i=0;i<varp->size;i++)
		    ((DAINT *) varp->opaque)[i] = 0;
	      } else if(walkop == WALK_INCREMENT_SCALERS) {
		if(varp->opaque)
		  for(i=0;i<varp->size;i++)
		    ((DAINT *) varp->opaque)[i]
		      += (((DAINT *) varp->varptr)[i] != 0);
	      } else {		/* WALK_CLEAR_FLAGS */
		if(varp->type == DAVARINT) {
		  for(i=0;i<varp->size;i++)
		    ((DAINT *) varp->varptr)[i] = 0;
		}
	      }
	    }
	    vlist = vlist->next;
	  }
	} else
	  fprintf(STDERR,"Unimplimented WALK code\n");
      }
      nextblock = thisblock->next;
      if(block_name) {
	if(strcasecmp(block_name,thisblock->blockname)==0) {
	  if(walkop == WALK_REMOVE){
	    *lastblockp = nextblock;
	    free(thisblock);
	  }
	  return(S_SUCCESS);
	  }
      } else if (walkop == WALK_REMOVE){
	free(thisblock);
      }
      thisblock = nextblock;
    }
    if(block_name) {
      fprintf(STDERR,"Test block %s not found\n",block_name);
      return(S_FAILURE);
    }
    if(walkop == WALK_REMOVE)
      thTBlockListP = (thTBlockList *) NULL;
  }
  return(S_SUCCESS);
}
    
/* Fortran callable versions of the various test walk routines */

#ifdef NOF77extname
int thtstexe()
#else
int thtstexe_()
#endif
{
  int A0;
  A0 = thWalkTree(0,WALK_EXECUTE);
  return A0;
}
#ifdef NOF77extname
int thtstexeb
#else
int thtstexeb_
#endif
(char *A1,unsigned C1)
{
  int A0;
  char *B1;
  A0 = thWalkTree((!*(int *)A1)?0:memchr(A1,'\0',C1)?A1:
                  (memcpy(B1=malloc(C1+1),A1,C1),B1[C1]='\0'
                   ,kill_trailing(B1,' ')),WALK_EXECUTE);
  if(B1) free(B1);
  return A0;
}
#ifdef NOF77extname
int thtstdis()
#else
int thtstdis_()
#endif
{
  int A0;
  A0 = thWalkTree(0,WALK_DISPLAY);
  return A0;
}
#ifdef NOF77extname
int thtstdisb
#else
int thtstdisb_
#endif
(char *A1,unsigned C1)
{
  int A0;
  char *B1;
  A0 = thWalkTree((!*(int *)A1)?0:memchr(A1,'\0',C1)?A1:
                  (memcpy(B1=malloc(C1+1),A1,C1),B1[C1]='\0'
                   ,kill_trailing(B1,' ')),WALK_DISPLAY);
  if(B1) free(B1);
  return A0;
}
#ifdef NOF77extname
int thtstclr()
#else
int thtstclr_()
#endif
{
  int A0;
  A0 = thWalkTree(0,WALK_CLEAR_FLAGS);
  return A0;
}
#ifdef NOF77extname
int thtstclrb
#else
int thtstclrb_
#endif
(char *A1,unsigned C1)
{
  int A0;
  char *B1;
  A0 = thWalkTree((!*(int *)A1)?0:memchr(A1,'\0',C1)?A1:
                  (memcpy(B1=malloc(C1+1),A1,C1),B1[C1]='\0'
                   ,kill_trailing(B1,' ')),WALK_CLEAR_FLAGS);
  if(B1) free(B1);
  return A0;
}
#ifdef NOF77extname
int thtstins()
#else
int thtstins_()
#endif
{
  int A0;
  A0 = thWalkTree(0,WALK_INCREMENT_SCALERS);
  return A0;
}
#ifdef NOF77extname
int thtstinsb
#else
int thtstinsb_
#endif
(char *A1,unsigned C1)
{
  int A0;
  char *B1;
  A0 = thWalkTree((!*(int *)A1)?0:memchr(A1,'\0',C1)?A1:
                  (memcpy(B1=malloc(C1+1),A1,C1),B1[C1]='\0'
                   ,kill_trailing(B1,' ')),WALK_INCREMENT_SCALERS);
  if(B1) free(B1);
  return A0;
}
#ifdef NOF77extname
int thtstcls()
#else
int thtstcls_()
#endif
{
  int A0;
  A0 = thWalkTree(0,WALK_CLEAR_SCALERS);
  return A0;
}
#ifdef NOF77extname
int thtstclsb
#else
int thtstclsb_
#endif
(char *A1,unsigned C1)
{
  int A0;
  char *B1;
  A0 = thWalkTree((!*(int *)A1)?0:memchr(A1,'\0',C1)?A1:
                  (memcpy(B1=malloc(C1+1),A1,C1),B1[C1]='\0'
                   ,kill_trailing(B1,' ')),WALK_CLEAR_SCALERS);
  if(B1) free(B1);
  return A0;
}
thStatus thExecuteTests(char *block_name)
{
  return(thWalkTree(block_name,WALK_EXECUTE));
}
thStatus thClearTestFlags(char *block_name)
{
  return(thWalkTree(block_name,WALK_CLEAR_FLAGS));
}
thStatus thClearTestScalers(char *block_name)
{
  return(thWalkTree(block_name,WALK_CLEAR_SCALERS));
}
thStatus thIncTestScalers(char *block_name)
{
  return(thWalkTree(block_name,WALK_INCREMENT_SCALERS));
}
