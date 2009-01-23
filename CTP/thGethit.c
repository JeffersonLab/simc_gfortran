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
 *  Retrieve data from Hall C Engine hit lists (gen_datastructures.cmn).
 *    
 * Author: Stephen A. Wood, CEBAF, Hall C
 *
 * Revision History:
 *   $Log: thGethit.c,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.2  1999/11/04 20:34:05  saw
 *   Alpha compatibility.
 *   New RPC call needed for root event display.
 *   Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *   Revision 1.6  1996/08/01 01:31:14  saw
 *   Change argument of thExecuteaGethitBlock from block pointer to a var
 *
 * Revision 1.5  1995/04/10  15:52:52  saw
 * No defined gethit blocks is not an error in thExecuteGethits
 *
 * Revision 1.4  1995/01/13  16:26:01  saw
 * Add missing return(S_SUCCESS) calls to various routines
 *
 * Revision 1.3  1994/09/27  19:26:21  saw
 * Remove linux dependencies
 *
 * Revision 1.2  1994/08/26  12:30:23  saw
 * Register result and test variables as needed
 *
 * Revision 1.1  1994/07/21  18:48:12  saw
 * Initial revision
 *
 */
/*
Any gethit block can only operate on one "data structure".

begin gethit wc
missing_value=-1 ! Value to give when hit not found
hitcount=numhits ! Pointer to word holding hit count
matchlist=plane  ! First array to match with
matchlist=counter! Second array to match with
valuelist=data   ! Default list containing data values

xsin1,txsin1:,1,4 ! Use default array   !Makes hit.xscin1, and test.txsin1
xsin2,txsin2:data2,1,5
end gethist wc
*/
#include <stdio.h>
#include <string.h>
#include "daVar.h"
#include "th.h"
#include "thInternal.h"
#include "thUtils.h"
#include "cfortran.h"

struct thGethitList {
  struct thGethitList *next;
  daVarStruct *src;		/* Array to get data from when hit is found */
  DAINT *coord;			/* Coordinate values to match*/
  daVarStruct *dest,*test;	/* Variable to put value of hit detector
				   and a flag for hit found/not found */
  int destindex, testindex;
};
typedef struct thGethitList thGethitList;
  
struct thGethitOpaque {
  DAINT *vnhits;		/* CTP variable with number of hits */
  int ncoord;			/* Number of coordinates that must match */
  DAINT **coord_array;		/* Point to list of arrays containing hit
				   coordinates */
  thGethitList *hlisthead;	/* List of gethit definitions */
  daVarStruct *srcdefault;	/* Default data array */
  DADOUBLE default_real;	/* Default value to give when hit not found */
  DAINT default_int;		/* Int of default_real */
};

typedef struct thGethitOpaque thGethitOpaque;

struct thGBlockList {
  struct thGBlockList *next;
  char *blockname;		/* Block name without the block.gethit */
  daVarStruct *var;		/* Pointer to variable that describes block */
};
typedef struct thGBlockList thGBlockList;

thGBlockList *thGBlockListP = NULL;
char *hitsourceclasslist[]={EVENTSTR,0};

/* Function prototypes */
thStatus thGethitPreamble(char *line, thGethitOpaque *opqptr,int  *icoord);
char *thGetline(char *lines,char **lcopy);
thStatus thBookaGethit(char *line, thGethitOpaque *opqptr
		       ,thGethitList **thGethitNext);
thStatus thExecuteaGethitBlock(daVarStruct *var);
thStatus thExecuteaGethit(thGBlockList *block);

thStatus thBookGethits(daVarStruct *var)
{
  char *lines,*eol,*lcopy;
  thGethitOpaque *opqptr;
  thGethitList **thGethitNext;
  int ncoord,icoord;
  int line_count,mode;
  char *blockname;

  blockname = var->name;	/* If name doesn't fit pattern, use whole */
/*  printf("Booking gethit block %s\n",blockname);*/
  if(strcasestr(var->name,BLOCKSTR)==var->name){
    int i;
    i = strlen(BLOCKSTR) + 1;
    if(strcasestr((var->name + i),GETHITSTR)==(var->name + i)){
      i += strlen(GETHITSTR);
      if(*(var->name + i) == '.'){
	blockname += i + 1;
      }
    }
  }
      
  if(var->opaque == 0) {
    opqptr = (thGethitOpaque *) malloc(sizeof(thGethitOpaque));
    var->opaque = (void *) opqptr;
  } else {
    daVarStructList *next,*this;

    opqptr = (thGethitOpaque *) var->opaque;
    free(opqptr->coord_array);
    /* Walk down optptr->hlisthead freeing the individual gethit structures */
  }
  /* Initialize/Clear the Gethit structure */
  opqptr->vnhits = (DAINT *) 0;
  opqptr->ncoord = 0;
  opqptr->coord_array = (DAINT **) 0;
  opqptr->srcdefault = (daVarStruct *) 0;
  opqptr->hlisthead = (thGethitList *) 0;
  opqptr->default_real = opqptr->default_int = 0;

  lines = var->title;
  line_count = 0;
  ncoord=0;

  while(*lines) {		/* First Pass */
				/* Count the number of arrays to match */

    line_count++;
    lines = thGetline(lines,&lcopy);
    if(!thCleanLine(lcopy)){
      if(strcasestr(lcopy,"matchlist") && strchr(lcopy,'=')) {
	ncoord++;
      }
      if(strchr(lcopy,':')) {	/* First gethit definition */
	break;
      }
    }
  }
  opqptr->ncoord = ncoord;
  opqptr->coord_array = (DAINT **) malloc(ncoord*sizeof(DAINT *));

  mode = 0;			/* setup mode */
  thGethitNext = (thGethitList **) &opqptr->hlisthead;
  line_count = 0;
  icoord = 0;

  lines = var->title;
  while(*lines) {		/* Second Pass */
    line_count++;
    lines = thGetline(lines,&lcopy);
    if(!thCleanLine(lcopy)){
      if(strchr(lcopy,':'))
	mode = 1;
      if(mode) {		/* Defining the individual gethits */
	if(thBookaGethit(lcopy,opqptr,thGethitNext)==S_SUCCESS){
	  thGethitNext = &((*thGethitNext)->next);
	} else {
	  fprintf(STDERR,"Gethit booking error in line %d\n",line_count);
	}
      } else {			/* Preamble */
	if(thGethitPreamble(lcopy,opqptr,&icoord)!=S_SUCCESS){
	  fprintf(STDERR,"Gethit booking error in line %d\n",line_count);
	}
      }
    }
  }
  /* Update internal table of gethit blocks */
  {
    thGBlockList *thisblock,*nextblock,**lastblockp;
    nextblock = thGBlockListP;
    lastblockp = &thGBlockListP;
    thisblock = thGBlockListP;
    while(thisblock){
      if((strcasecmp(thisblock->var->name,var->name)) == 0){
	/* Replacing a block with a new definition */
	fprintf(stderr,"Replacing %s with new definition\n",var->name);
	if(thisblock->var != var){
	  fprintf(stderr,"ERR: Same name, different var pointer\n");
	}
	break;
      }
      lastblockp = &thisblock->next;
      thisblock = thisblock->next;
    }
    if(!thisblock){		/* Create entry for New block */
      *lastblockp = thisblock = (thGBlockList *) malloc(sizeof(thGBlockList));
      thisblock->var = var;
      thisblock->next = (thGBlockList *) NULL;
      thisblock->blockname = (char *) malloc(strlen(blockname) + 1);
      strcpy(thisblock->blockname,blockname);
    }
  }
/*  printf("Returning from booking Gethit's\n");*/
  return(S_SUCCESS);
}
char *thGetline(char *lines,char **lcopy)
/* Pull out the characters from lines up to a newline and copy them into
   a new array.  Return the pointer to that array.  If a line is missing
   a newline, return an error. */
{
  static char *line_copy=0;
  static line_copy_size=0;
  char *next;
  int len;
  char *eol;

  eol = strchr(lines,'\n');
  if(!eol) {
    len = strlen(lines);
    next = lines + len;
  } else {
    len = (eol-lines);
    next = eol+1;
  }
  if(!line_copy) {
    line_copy = (char *) malloc(len+1);
    line_copy_size = len+1;
  } else {
    if(len >= line_copy_size) {
      line_copy = (char *) realloc(line_copy,len+1);
      line_copy_size = len+1;
    }
  }
  strncpy(line_copy,lines,len);
  line_copy[len] = '\0';
  *lcopy = line_copy;
  return(next);
}
thStatus thGethitPreamble(char *line, thGethitOpaque *opqptr,int  *icoord){
  char *equal;
  char *arg;

/*  printf("Processing preamble line %x %s\n",line,line);*/
  if(!(equal=strchr(line,'=')))
    return(S_SUCCESS);
  arg = equal + 1;
  *equal = 0;
  arg = thSpaceStrip(arg);
/*  printf("line=%s,arg=%s\n",line,arg); */
  if(strcasestr(line,"missing_value")) {
    opqptr->default_real = atof(arg);
    opqptr->default_int = floatToLong(opqptr->default_real);
  } else {
    daVarStruct *varp;
/*    printf("line=%s,arg=%x %s\n",line,arg,arg);*/
    if(daVarLookupPWithClass(arg,hitsourceclasslist,&varp)!=S_SUCCESS){
      fprintf(STDERR,"(thGethitPreamble )Variable %s not found\n",arg);
      return(S_FAILURE);	/* Variable not found */
    }
/*    printf("line=%s,arg=%x %s\n",line,arg,arg);*/
/*    printf("Searching for keyword in %s\n",line);*/
    if(strcasestr(line,"valuelist")) {
      opqptr->srcdefault = varp;
    } else {
      if(varp->type != DAVARINT) {
	return(S_FAILURE);	/* Not an integer array */
      }
      if(strcasestr(line,"hitcount")) {
	opqptr->vnhits = (DAINT *) varp->varptr;
/*	printf("Setting vnhits to %s\n",varp->name);*/
      } else if(strcasestr(line,"matchlist")) {
/*	printf("Coordinate %d=%s\n",*icoord,varp->name);*/
	opqptr->coord_array[(*icoord)++] = (DAINT *) varp->varptr;
      } else {
	return(S_FAILURE);
      }
    }
  }
  return(S_SUCCESS);
}
thStatus thBookaGethit(char *line, thGethitOpaque *opqptr, thGethitList **thGethitNext)
{
  char *colon;
  int ndestarg;			/* Number of destination args */
  char *destargs[20];
  int nsrcarg;			/* Number of destination args */
  char *srcargs[20];
  daVarStruct *destp, *testp;
  int destind,testind;		/* Indexes into destination and test arrays */
  daVarStruct *source;
  DAINT *coord;
  int icoord;
  thGethitList *Gethit;
  thStatus status;

  if(!(colon=strchr(line,':'))) {
    return(S_SUCCESS);
  }
  
  *colon = '\0';
  colon++;

  ndestarg = thCommas(line,destargs);
  nsrcarg = thCommas(colon,srcargs);

  /* Parse the source side first so that we know the data type of the
     source in case we need to register the destination */
  
  /*  printf("%s found as %s at %x\n",destargs[0],destp->name,destp->varptr);*/
  if(nsrcarg < opqptr->ncoord+1) {
    fprintf(STDERR,"Insufficient arguments after :\n");
    return(S_FAILURE);
  }
  srcargs[0] = thSpaceStrip(srcargs[0]);
  if(strlen(srcargs[0]) == 0) {
    source = opqptr->srcdefault;
  } else {
    if(daVarLookupPWithClass(srcargs[0],hitsourceclasslist,*source)
       != S_SUCCESS) {
      fprintf(STDERR,"%s not registered\n",srcargs[0]);
      return(S_FAILURE);		/* Destination not registered */
    }
  }

  destargs[0] = thSpaceStrip(destargs[0]);
  if((status=thVarResolve(destargs[0],&destp,&destind,2,source->type)) != S_SUCCESS){
    return(S_FAILURE);
    /* ASAP we must change this to register variables as they are needed */
      /* If the variable exists, then we also must check to make sure that
	 the requested index does not exceed the size of the array.
	 a new thVarResolve should also increase the size of the array if
	 it was created by CTP */
  }
  if(ndestarg > 1){
    destargs[1] = thSpaceStrip(destargs[1]);
    if(thVarResolve(destargs[1],&testp,&testind,3,DAVARINT) != S_SUCCESS){
      return(S_FAILURE);		/* Test flag not registered */
      /* ASAP we must change this to register variables as they are needed */
      /* If the variable exists, then we also must check to make sure that
	 the requested index does not exceed the size of the array.
	 a new thVarResolve should also increase the size of the array if
	 it was created by CTP */
      
    }
  } else {
    testp = 0;
    testind = 0;
  }

  coord = (DAINT *) malloc(opqptr->ncoord*sizeof(DAINT));
  for(icoord=0;icoord<opqptr->ncoord;icoord++) { 
    srcargs[icoord+1] = thSpaceStrip(srcargs[icoord+1]);
    if(thEvalImed(srcargs[icoord+1],0,&coord[icoord]) != S_SUCCESS) {
      fprintf(STDERR,"Error evaluating %s\n",srcargs[icoord+1]);
      free(coord);
      return(S_FAILURE);
    }
  }
/* Everything obtained from line now */
  Gethit = *thGethitNext = (thGethitList *) malloc(sizeof(thGethitList));
  Gethit->next = (thGethitList *) NULL;
  Gethit->src = source;
  Gethit->coord = coord;
  Gethit->dest = destp;
  Gethit->destindex = destind;
  Gethit->test = testp;
  Gethit->testindex = testind;

  return(S_SUCCESS);
}
thStatus thExecuteGethits(char *block_name){
  thGBlockList *thisblock;
  

  if(block_name) if(*block_name=='\0') block_name = 0;

  if(thGBlockListP == 0){
    return(S_SUCCESS);		/* No gethits defined */
  } else {
    thisblock = thGBlockListP;
    while(thisblock){
      if(block_name)
	if(strcasecmp(block_name,thisblock->blockname)!=0){
	  thisblock = thisblock->next;
	  continue;
	}
      thExecuteaGethitBlock(thisblock->var);
      (*((DAINT *)thisblock->var->varptr))++; /* Increment block counter */
      if(block_name)
	if(strcasecmp(block_name,thisblock->blockname)==0) return(S_SUCCESS);
      thisblock = thisblock->next;
    }
  }
  return(S_SUCCESS);
}
  
thStatus thExecuteaGethitBlock(daVarStruct *var)
     /* Execute a gethit block */
{
  thGethitOpaque *opqptr;
  thGethitList *thisgethit;
  int nhits;
  int ihit;
  int icoord,ncoord;
  DAINT **coord_array;
  double dval;
  int ival;
  
/*  opqptr = block->var->opaque;*/	/* Structure that describes the gethits */
  opqptr = var->opaque;
/*  printf("opqptr=%x\n",opqptr);*/
  nhits = *(opqptr->vnhits);
  ncoord = opqptr->ncoord;
  coord_array = opqptr->coord_array;

  thisgethit = opqptr->hlisthead;
/*  printf("%d hits this event\n",nhits);*/
  while(thisgethit){		/* Inefficient algorithm */
/*    printf("Getting %s\n",thisgethit->dest->name);*/
    if(thisgethit->test) {
      /* Assume that test is integer */
      *((DAINT *) thisgethit->test->varptr + thisgethit->testindex) = FALSE;
    }
    for(ihit=0;ihit<nhits;ihit++){
      for(icoord=0;icoord<ncoord;icoord++){
/*	printf("Hit=%d, Coord=%d %d %d\n",ihit,icoord
	       ,coord_array[icoord][ihit],thisgethit->coord[icoord]);*/
	if((coord_array[icoord])[ihit] != thisgethit->coord[icoord]) break;
      }
      if(icoord >= ncoord){	/* All coordinates matched */
	int srctype;
/*	printf("Matched at hit %d, %d %d\n",ihit,icoord,ncoord);*/
	if(thisgethit->test) 
	  *((DAINT *) thisgethit->test->varptr + thisgethit->testindex) = TRUE;
	/* Need to grab the data value, stuff it into variable */
	srctype = thisgethit->src->type;
	switch(srctype)
	  {
	  case DAVARINT:
	    ival = ((DAINT *)thisgethit->src->varptr)[ihit];
	    switch(thisgethit->dest->type)
	      {
	      case DAVARINT:
		((DAINT *)thisgethit->dest->varptr)[thisgethit->destindex] = 
		  ival;
		break;
	      case DAVARFLOAT:
		((DAFLOAT *)thisgethit->dest->varptr)[thisgethit->destindex] = 
		  ival;
		break;
	      case DAVARDOUBLE:
		((DADOUBLE *)thisgethit->dest->varptr)[thisgethit->destindex]
		  = ival;
		break;
	      }
	    break;
	  case DAVARFLOAT:
	  case DAVARDOUBLE:
	    if(srctype == DAVARFLOAT)
	      dval = ((DAFLOAT *)thisgethit->src->varptr)[ihit];
	    else 
	      dval = ((DADOUBLE *)thisgethit->src->varptr)[ihit];
	    switch(thisgethit->dest->type)
	      {
	      case DAVARINT:
		((DAINT *)thisgethit->dest->varptr)[thisgethit->destindex] = 
		  floatToLong(dval);
		break;
	      case DAVARFLOAT:
		((DAFLOAT *)thisgethit->dest->varptr)[thisgethit->destindex] = 
		  dval;
		break;
	      case DAVARDOUBLE:
		((DADOUBLE *)thisgethit->dest->varptr)[thisgethit->destindex]
		  = dval;
		break;
	      }
	    break;
	  }
	break;			/* Only match one hit for now */
      }
    }
    thisgethit = thisgethit->next;
  }
  return(S_SUCCESS);
}
int thgethit_()
{
  int A0;
  A0 = thExecuteGethits(0);
  return A0;
}
int thgethitb_(char *A1,unsigned C1)
{
  int A0;
  char *B1;
  A0 = thExecuteGethits((!*(int *)A1)?0:memchr(A1,'\0',C1)?A1:
		      (memcpy(B1=malloc(C1+1),A1,C1),B1[C1]='\0'
		       ,kill_trailing(B1,' ')));
  if(B1) free(B1);
  return A0;
}
