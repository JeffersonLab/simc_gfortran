/*-----------------------------------------------------------------------------
 * Copyright (c) 1993,1994 Southeastern Universities Research Association,
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
 *  Books histograms.  Constants or array indices in configuration line may
 *  be expressions.  The expressions get evaluated at run time.
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: thHist.c,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.3  2003/02/21 20:55:24  saw
 *   Clean up some types and casts to reduce compiler warnings.
 *
 *   Revision 1.2  1999/11/04 20:34:06  saw
 *   Alpha compatibility.
 *   New RPC call needed for root event display.
 *   Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *   Revision 1.1  1998/12/07 22:11:12  saw
 *   Initial setup
 *
 *   Revision 1.16  1996/07/31 20:33:53  saw
 *   Show line number for booking errors.
 *
 *   Revision 1.15  1996/01/30 15:37:37  saw
 *   Add thExecuteHistsV to execute the block contained by a var.  Add thClearHists.
 *
 *	  Revision 1.14  1995/04/10  15:54:18  saw
 *	  No defined hist blocks is not an error in thExecuteHists
 *
 *	  Revision 1.13  1995/01/13  16:25:49  saw
 *	  Add missing return(S_SUCCESS) calls to various routines
 *
 *	  Revision 1.12  1995/01/09  15:53:45  saw
 *	  On fprintf, indicate block type and well as name.  Fix some memory leaks
 *	  and shorted mallocs.
 *
 *	  Revision 1.11  1994/09/27  20:14:17  saw
 *	  Remove some linux dependencies.  Define dummy hbook routines.
 *
 *	  Revision 1.10  1994/08/26  13:17:44  saw
 *	  Add option to thVarResolve to register unknown variables when requested
 *
 *	  Revision 1.9  1994/07/11  18:40:07  saw
 *	  (SAW) Move thHistOpaque structure here from thInternals.h
 *
 *	  Revision 1.8  1994/06/28  19:40:04  saw
 *	  Fix infinite loops when uhist blocks were encountered in thhstexe calls
 *
 *	  Revision 1.7  1994/06/13  13:05:01  saw
 *	  Update ifdef's so this can compile under linux (no f77)
 *
 *	  Revision 1.6  1994/04/05  19:08:18  saw
 *	  Use HEXIST to make sure we don't book existing histograms.
 *
 *	  Revision 1.5  1993/12/02  21:26:12  saw
 *	  Allow histogram filling from doubles (REAL*8)
 *
 *	  Revision 1.4  1993/11/24  21:32:31  saw
 *	  Floating argument of thEvalImed is now double.
 *
 *	  Revision 1.3  1993/07/09  18:45:45  saw
 *	  Fix detection of 2d user histograms
 *
 *	  Revision 1.2  1993/05/10  20:46:54  saw
 *	  Fix header
 *
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "daVar.h"
#include "th.h"
#include "thInternal.h"
#include "thUtils.h"
#include "hbook.h"
#ifdef NOHBOOK
int hexist_(){return(0);};
float hi_(){};
float hij_(){};
float hx_(){};
float hxy_(){};
float hie_(){};
float hxe_(){};
float hif_(){};
float hmax_(){};
float hmin_(){};
float hsum_(){};
float hstati_(){};
float hspfun_(){};
float hrndm1_(){};
float hrndm2_(){};
void hbook1_(int *id, char *chtitle, int *nx, float *xmi, float *xma,
	     float *vmx,...){};
void hbook2_(int *id, char *chtitle, int *nx, float *xmi, float *xma,
	     int *ny, float *ymi, float *yma, float *vmx,...){};
void hf1_(int *id, float *x, float *weight,...){};
void hf2_(int *id, float *x, float *y, float *weight,...){};
void hdelet_(int *id,...){};
void hlimit_(){};
void hrfile_(){};
void hrout_(){};
void hrend_(){};
/*
void HF1(){};
void HF2(){};
void HDELET(){};
void HBOOK1(){};
void HBOOK2(){};
int HEXIST(){};
void HRFILE(){};
void HREND(){};
*/
#endif

extern daVarStatus thHistRHandler();

struct thHistSpecList {		/* Opaque structure for blocks */
  daVarStruct *varname;		/* name is "hist.xxx", varptr points to ID
				   title is optional title
				   opaque points nd, and x,y,test vars */
  struct thHistSpecList *next;
};
typedef struct thHistSpecList thHistSpecList;
				   
typedef enum {HISTCLASS, UHISTCLASS} hclasstype;

struct thHistOpaque {		/* Opaque structure for histogram definition */
  int nd;			/* Number of dimensions */
  daVarStruct *x; int xindex;
  daVarStruct *y; int yindex;
  daVarStruct *test; int testindex;
  daVarStruct *weight; int weightindex;
  /* Include here the limits? */
  /* Include the pointer to the variable of how much to fill by. */
};
typedef struct thHistOpaque thHistOpaque;

struct thHBlockList {
  struct thHBlockList *next;
  char *blockname;		/* Block name without the "block.hist" */
  hclasstype classtype;
  daVarStruct *var;		/* Varptr points to # times called
				   Title is code from file.
				   opaque is pointer to hist speclist */
};
typedef struct thHBlockList thHBlockList;

thHBlockList *thHBlockListP;	/* Pointer to list of hist blocks */

thStatus thBookaHist(char *line, thHistSpecList **thHistNext);
thStatus thExecuteaHist(thHistSpecList *Hist);
thStatus thRemoveHists(char *block_name);
thStatus thVarResolve(char *s, daVarStruct **var,int *index,int datatest,int newvartype);

FCALLSCFUN1(INT,thGetHistID,THGETID,thgetid,STRING)
FCALLSCFUN1(INT,thHistAliasWrite,THWHALIAS,thwhalias,STRING)
     
/*
enum attribute_value {HNAME, XSOURCE, YSOURCE, NBINS, NBINX, NBINY, XLOW,XHIGH
, YLOW, YHIGH, TEST, TITLE, WEIGHT};
char *attribute_string[] = {"name","xsource","ysource","nbins","nbinx","nbiny"
,"xlow","xhigh","ylow","yhigh","title","test","weight"};
*/
int thMaxID=0;			/* Maximum Hbook ID used so far */
char *datasourcelist[]={EVENTSTR,TESTSTR,0}; /* Immediate expressions */
char *testflaglist[]={TESTSTR,EVENTSTR,PARMSTR,0}; /* Logical operand */
char *histclasslist[]={UHISTSTR,HISTSTR,0}; /* Hist name classes */
     
thStatus thBookHists(daVarStruct *var)
{
  char *lines,*eol;
  int line_count;
  thHistSpecList **thHistNext;
  char *blockname;
  hclasstype classtype;

#ifndef NOHBOOK
  thHistZeroLastId();
#endif
  
  {
    int i;
/*    printf("In bookhists\n");*/
    /* Get the name without the block.hist on it */
    blockname = var->name;	/* If name doesn't fit pattern, use whole */
    if(strcasestr(var->name,BLOCKSTR)==var->name){
      i = strlen(BLOCKSTR) + 1;
      if(strcasestr((var->name + i),HISTSTR)==(var->name + i)){
	classtype = HISTCLASS;
	i += strlen(HISTSTR);
	if(*(var->name + i) == '.'){
	  blockname += i + 1;
	}
      } else if(strcasestr((var->name + i),UHISTSTR)==(var->name + i)){
	classtype = UHISTCLASS;
	i += strlen(UHISTSTR);
	if(*(var->name + i) == '.'){
	  blockname += i + 1;
	}
      }
    }
  }
/*  printf("Booking histogram block %s\n",blockname);*/

  if(var->opaque) thRemoveHists(blockname);

  thHistNext = (thHistSpecList **) &var->opaque;
  lines = var->title;
  line_count = 0;
  while(*lines){
    char *lcopy;

    line_count++;
    eol = strchr(lines,'\n');
    if(!eol) {
      fprintf(stderr,"L %d: Last line of hist block %s has no newline\n"
	      ,line_count,var->name);
      break;
    }
    if(*(eol+1)=='\0'){		/* This is the last line */
      if(strcasestr(lines,ENDSTR) == 0)
	fprintf(stderr,"L %d: Last line of hist block %s is not an END\n"
		,line_count,var->name);
      break;
    }
    if(line_count == 1)
      if(strcasestr(lines,BEGINSTR) != 0){
/*	printf("Is a begin\n");*/
	lines = eol + 1;
	continue;
      } else
	fprintf(stderr,"First line of hist block %s is not a BEGIN\n",var->name);
    /* Ready to book the line, Add continuation lines later */
    lcopy = (char *) malloc(eol-lines+1);
    strncpy(lcopy,lines,(eol-lines));
    *(lcopy + (eol-lines)) = '\0';
/*    printf("Passing|%s|\n",lcopy);*/
    if(!thCleanLine(lcopy)){
      if(thBookaHist(lcopy,thHistNext)==S_SUCCESS){
	thHistNext = &((*thHistNext)->next);
      } else {
	fprintf(stderr,"(%s): Hist booking error in line %d\n",var->name,line_count);
      }
    }
    free(lcopy);
    lines = eol+1;
  }
  /* Update internal table of hist blocks. */
  {
    thHBlockList *thisblock,*nextblock,**lastblockp;
    nextblock = thHBlockListP;
    lastblockp = &thHBlockListP;
    thisblock = thHBlockListP;
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
      *lastblockp = thisblock = (thHBlockList *) malloc(sizeof(thHBlockList));
      thisblock->var = var;
      thisblock->next = (thHBlockList *) NULL;
      thisblock->blockname = (char *) malloc(strlen(blockname) + 1);
      thisblock->classtype = classtype;
      strcpy(thisblock->blockname,blockname);
    }
  }
/*  printf("Returning from booking hists\n");*/
  return(S_SUCCESS);
}

thStatus thBookaHist(char *line, thHistSpecList **thHistNext)
{
  int nargs, nd, n, id;
  char *long_title;
  daVarStruct *(varp[2]),*testp,*weightp;
  int vind[2],tind,wind;
  thHistOpaque *histpars;
  char *args[20];
  /*  int type;
      thTokenType toktyp;
      int nd,n;
      */
  int userflag;			/* 0 for normal hist, 1 for user hist */
  
  {
    char *s;
    s = line;
    long_title = 0;
    while(*s != 0) {
      if(*s == HTITLECHAR){
	*s++ = 0;
	long_title = (char *) malloc(strlen(s)+1);
	strcpy(long_title,s);	/* Make sure to free this  */
	break;
      } else s++;
    }
  }
  if(!long_title) {
    long_title = (char *) malloc(strlen(line)+1);
    strcpy(long_title,line);
  }
  /* All between # and title char, comment char, or EOL is the weight */
  {
    char *s;
    if((s=strchr(line,WEIGHTCHAR))) {
      *s++=0;
      s = thSpaceStrip(s);
      if(thVarResolve(s,&weightp,&wind,0,0) != S_SUCCESS) {
	free(long_title);
	return(S_FAILURE);
      }
    } else {
      weightp = 0;
      wind = 0;
    }
  }
  
  nargs = thCommas(line,args);
  args[0] = thSpaceStrip(args[0]);
  
  userflag = 0;
  if(nargs >= 4 && nargs <=6) {
    nd = 1;
    if(nargs == 4) userflag = 1;
  } else if(nargs == 7) {
    nd = 2;
    userflag = 1;
  } else if(nargs >= 9 && nargs <= 10){
    nd = 2;
  } else {
    fprintf(stderr,"Incorrect number of arguments\n");
    free(long_title);
    return(S_FAILURE);
  }
  
  varp[1] = 0; vind[1] = 0;
  testp = 0; tind = 0;
  if(userflag) {
    varp[0] = 0; vind[0] = 0;
  } else {
    for(n=0;n<nd;n++){		/* Interpret the data sources */
      args[n+1] = thSpaceStrip(args[n+1]);
      if(thVarResolve(args[n+1],&varp[n],&vind[n],0,0) != S_SUCCESS) {
	free(long_title);
	return(S_FAILURE);
      }
    }				/* Data sources now defined */
    if((nd==1&&nargs==6)||(nd==2&&nargs==10)){
      args[nargs-1] = thSpaceStrip(args[nargs-1]);
      if(thVarResolve(args[nargs-1],&testp,&tind,1,0) != S_SUCCESS) {
	free(long_title);
	return(S_FAILURE);
      }
      if(testp->type != DAVARINT){
	fprintf(stderr,"Test flag %s must be an integer\n",args[nargs-1]);
	free(long_title);
	return(S_FAILURE);
      }
    } else {
/*      printf("Test is NULL\n");*/
      testp = NULL;
    }
  }
  /* Find/Create the variable to hold hist def stuff */
  {
    char *name;
    daVarStruct var;
    thHistSpecList *Hist;

    name = (char *) malloc(strlen(HISTSTR)+strlen(args[0])+2);
    strcpy(name,HISTSTR);
    strcat(name,".");
    strcat(name,args[0]);
    if(daVarLookup(name,&var)!=S_SUCCESS){
      var.name = name;
      var.size = 1;
      var.varptr = (void *) malloc(sizeof(DAINT));
      *((DAINT *) var.varptr) = ++thMaxID;
      var.type = DAVARINT;
      var.flag = DAVAR_READONLY | DAVAR_REPOINTOK;
      if(userflag)
	var.opaque = 0;
      else
	var.opaque = (void *) malloc(sizeof(thHistOpaque));
      var.whook = 0;
#ifndef NOHBOOK
      var.rhook = thHistRHandler;
#endif
      var.title = long_title;
/*      printf("Registering %s\n",var.name);*/
      daVarRegister((int) 0,&var); /* Create variable for histogram */
    }
    /* Make sure that we don't use an existing histogram id */
    while(HEXIST((id= *((DAINT *) var.varptr))) != 0) {
      *((DAINT *) var.varptr) = ++thMaxID;
      }
    id = *((DAINT *) var.varptr);
    if(!userflag) {
      histpars = var.opaque;
    }
    /* Perhaps we don't want to put the user hists in this list?? */
    Hist = *thHistNext = (thHistSpecList *) malloc(sizeof(thHistSpecList));
    Hist->next = (thHistSpecList *) NULL;
    if(daVarLookupP(name,&Hist->varname)!=S_SUCCESS){
      fprintf(stderr,"This can't happen\n");
      free(long_title);
      return(S_FAILURE);
    }
    free(name);
  }
  if(!userflag) {
    histpars->nd = nd;
    histpars->x = varp[0];
    histpars->xindex = vind[0];
    histpars->y = varp[1];
    histpars->yindex = vind[1];
    histpars->test = testp;
    histpars->testindex = tind;
    histpars->weight = weightp;
    histpars->weightindex = wind;
  }

  /* Data sources and test result now interpreted */
  {
    int nbinx,nbiny;
    double xlow,xhigh,ylow,yhigh;
    int ixargoffset, iyargoffset;

    if(userflag) {
      ixargoffset = 1;
    } else {
      ixargoffset = 1 + nd;
    }
    iyargoffset = ixargoffset + 3;
    if((thEvalImed(args[ixargoffset],0,&nbinx) != S_SUCCESS) ||
       (thEvalImed(args[ixargoffset+1],&xlow,0) != S_SUCCESS) ||
       (thEvalImed(args[ixargoffset+2],&xhigh,0) != S_SUCCESS))
      fprintf(stderr,"Error intrepreting histogram arguments\n");
    if(nd == 2){
      thEvalImed(args[iyargoffset],0,&nbiny);
      thEvalImed(args[iyargoffset+1],&ylow,0);
      thEvalImed(args[iyargoffset+2],&yhigh,0);
    }
    
    if(nd==1){
      HBOOK1(id,long_title,nbinx,(float) xlow,(float) xhigh,0.0);;
    } else {
      HBOOK2(id,long_title,nbinx,(float) xlow,(float) xhigh,nbiny,
	     (float) ylow,(float) yhigh,0.0);;
    }
  }
  /* Need to add to alias file */
  free(long_title);
  return(S_SUCCESS);
}
thStatus thVarResolve(char *s, daVarStruct **varpp, int *index,
		      int datatest, int newvartype)
/* Interpret the string as a variable name with a possible index.
   Pass the index off the thEvalImed to evaluate it.

Interpret the string s as an array element.  Looks for an index inside
of []'s, returning that in index.  If there is an index inside of ()'s, then
one is subtracted from the index before it is returnd.  []'s are for
c style indexing, ()'s for fortran style indexing.  A pointer to the [ or (
is returned, so that the variable name may be null terminated by the caller.
If there is an error in the balancing of the []'s or ()'s, a null is returned
to signify an error.
datatest values
      0   -   data source
      1   -   test flag
      2   -   data destination (create float if doesn't exist)
      3   -   test flag (create if it doens't exist)
*/
{
  int cstyle;
  char cleft,cright;
  char *leftp,*rightp;
  char **classlist;
  daVarStruct var;

  *varpp = 0;
  *index = 0;

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

  if(datatest & 1)
    classlist = testflaglist;
  else
    classlist = datasourcelist;
  if(daVarLookupPWithClass(s,classlist,varpp) != S_SUCCESS){
    if(datatest & 2) { /* Create the variable (as an int) */
      if(strchr(s,'.')) {	/* Don't prepend a class */
	var.name = (char *) malloc(strlen(s)+1);
	strcpy(var.name,s);
      } else {
	var.name = (char *) malloc(strlen(classlist[0])
				   +strlen(s)+2);
	strcpy(var.name,classlist[0]);
	strcat(var.name,".");
	strcat(var.name,s);
      }
      var.size = 1;
      var.type = newvartype;
      switch(newvartype)
	{
	case DAVARINT:
	  var.varptr = (void *) malloc(sizeof(DAINT));
	  *((DAINT *)var.varptr) = 0;
	  break;
	case DAVARFLOAT:
	  var.varptr = (void *) malloc(sizeof(DAFLOAT));
	  *((DAFLOAT *)var.varptr) = 0.0;
	  break;
	case DAVARDOUBLE:
	  var.varptr = (void *) malloc(sizeof(DADOUBLE));
	  *((DADOUBLE *)var.varptr) = 0.0;
	  break;
	}
      var.opaque = 0;
      var.rhook = 0;
      var.whook = 0;
      var.flag = DAVAR_REPOINTOK | DAVAR_READONLY | DAVAR_DYNAMIC_PAR;
      var.title = 0;
      printf("Registering %s at %d\n",var.name,var.varptr);
      daVarRegister((int) 0,&var); /* Create the parameter */
      daVarLookupP(var.name,varpp);
      free(var.name);
    } else {
      fprintf(stderr,"Variable %s must be registered\n",s);
      if(leftp) *leftp = cleft;
      return(S_FAILURE);
    }
  }
    
  if(leftp){
    int cstyle;
    char *sindex;
    int indtemp;

    *leftp = cleft;
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
    *index = indtemp + cstyle;
  }
  if(datatest & 2) {		/* See if array needs to be larger */
    if(*index >= (*varpp)->size) {
      (*varpp)->size = *index+1;
      switch((*varpp)->type)
	{
	case DAVARINT:
	  (*varpp)->varptr = (void *) realloc((*varpp)->varptr,(*index+1) * sizeof(DAINT));
	  break;
	case DAVARFLOAT:
	  (*varpp)->varptr = (void *) realloc((*varpp)->varptr,(*index+1) * sizeof(DAFLOAT));
	  break;
	case DAVARDOUBLE:
	  (*varpp)->varptr = (void *) realloc((*varpp)->varptr,(*index+1) * sizeof(DADOUBLE));
	  break;
	}
    }
  }
  return(S_SUCCESS);
}
  
thStatus thExecuteHistsV(daVarStruct *var){
  thHistSpecList *thishist;
  
  thishist = var->opaque;  
  while(thishist){
    thExecuteaHist(thishist);
    thishist = thishist->next;
  }
  (*((DAINT *)var->varptr))++; /* Increment block counter */
  return(S_SUCCESS);
}
thStatus thClearHistsV(daVarStruct *var){
  thHistSpecList *thishist;
  
  thishist = var->opaque;  
  while(thishist){
    HRESET(*(DAINT *) thishist->varname->varptr, thishist->varname->title);
    thishist = thishist->next;
  }
  (*((DAINT *)var->varptr)) = 0; /* Increment block counter */
  return(S_SUCCESS);
}

thStatus thExecuteHists(char *block_name){
  thHistSpecList *thishist;
  thHBlockList *thisblock;

  if(block_name) if(*block_name=='\0') block_name = 0;

  if(thHBlockListP == 0){
    return(S_SUCCESS);		/* No hists defined */
  } else {
    thisblock = thHBlockListP;
    while(thisblock){
      if(block_name)
	if(strcasecmp(block_name,thisblock->blockname)!=0){
	  thisblock = thisblock->next;
	  continue;
	}
      if(thisblock->classtype == UHISTCLASS){
	thisblock = thisblock->next;
	continue;
      }
      thishist = thisblock->var->opaque;
      while(thishist){
	thExecuteaHist(thishist);
	thishist = thishist->next;
      }
      (*((DAINT *)thisblock->var->varptr))++; /* Increment block counter */
      if(block_name)
	if(strcasecmp(block_name,thisblock->blockname)==0) return(S_SUCCESS);
      thisblock = thisblock->next;
    }
    if(block_name) {
      fprintf(STDERR,"Hist block %s not found\n",block_name);
      return(S_FAILURE);
    }
  }
  return(S_SUCCESS);
}
thStatus thExecuteaHist(thHistSpecList *Hist)
{
  float vals[2];

  int id,nd;
  daVarStruct *(varp[2]);
  int indxy[2];
  thHistOpaque *opqptr;
  double weight;

/*  nd = Hist->nd;*/
  opqptr = Hist->varname->opaque;
  if(!opqptr) return(S_SUCCESS); /* Assume it is a user filled histogram */
  nd = opqptr->nd;

  varp[0] = opqptr->x; indxy[0] = opqptr->xindex;
  varp[1] = opqptr->y; indxy[1] = opqptr->yindex;
  
  for(id=0;id<nd;id++){
    switch(varp[id]->type)
      {
      case DAVARINT:
	vals[id] = *((DAINT *) varp[id]->varptr + indxy[id]);
	break;
      case DAVARFLOAT:
	vals[id] = *((DAFLOAT *) varp[id]->varptr + indxy[id]);
	break;
      case DAVARDOUBLE:
	vals[id] = *((DADOUBLE *) varp[id]->varptr + indxy[id]);
	break;
      }
  }
  if(opqptr->weight) {
    switch(opqptr->weight->type)
      {
      case DAVARINT:
	weight = *((DAINT *) opqptr->weight->varptr + opqptr->weightindex);
	break;
      case DAVARFLOAT:
	weight = *((DAFLOAT *) opqptr->weight->varptr + opqptr->weightindex);
	break;
      case DAVARDOUBLE:
	weight = *((DADOUBLE *) opqptr->weight->varptr + opqptr->weightindex);
	break;
      }
  } else {
    weight = 1.0;
  }
  if((opqptr->test ? *((DAINT *) opqptr->test->varptr + opqptr->testindex) : 1)) {
    if(nd==1){
      HF1(*(DAINT *) Hist->varname->varptr,vals[0],weight);
/*      printf("Filling %s at %f\n",Hist->varname->name,vals[0]);*/
    }
    else
      HF2(*(int *) Hist->varname->varptr,vals[0],vals[1],weight);
  }
  return(S_SUCCESS);
}
thStatus thClearHists(char *block_name){
  thHistSpecList *thishist;
  thHBlockList *thisblock;

  if(block_name) if(*block_name=='\0') block_name = 0;

  if(thHBlockListP == 0){
    return(S_SUCCESS);		/* No tests defined */
  } else {
    thisblock = thHBlockListP;
    while(thisblock){
      if(block_name)
	if(strcasecmp(block_name,thisblock->blockname)!=0){
	  thisblock = thisblock->next;
	  continue;
	}
/*      if(thisblock->classtype == UHISTCLASS){
	      thisblock = thisblock->next;
	      continue;
	    }*/
      thishist = thisblock->var->opaque;
      while(thishist){
	HRESET(*(DAINT *) thishist->varname->varptr, thishist->varname->title);
	thishist = thishist->next;
      }
      (*((DAINT *)thisblock->var->varptr))=0; /* Increment block counter */
      if(block_name)
	if(strcasecmp(block_name,thisblock->blockname)==0) return(S_SUCCESS);
      thisblock = thisblock->next;
    }
  }
  return(S_SUCCESS);
}
thStatus thRemoveHists(char *block_name){
  thHistSpecList *thishist,*nexthist;
  thHBlockList *thisblock,*nextblock,**lastblockp;

  if(block_name)
    if(*block_name=='\0')
      block_name = 0;
    else
      lastblockp = &thHBlockListP;

  if(thHBlockListP == 0){
    return(S_FAILURE);		/* No histograms defined */
  } else {
    thisblock = thHBlockListP;
    if(!block_name)
      HDELET(0);		/* Free all space */
    while(thisblock) {
      if(block_name)
	if(strcasecmp(block_name,thisblock->blockname)!=0){
	  lastblockp = &thisblock->next;
	  thisblock = thisblock->next;
	  continue;
	}
      thishist = thisblock->var->opaque;
      while(thishist){
	nexthist = thishist->next;
	if(block_name) {
	  HDELET(*(DAINT *) thishist->varname->varptr);
	}
	free(thishist);
	thishist = nexthist;
      }
      nextblock = thisblock->next;
      if(block_name){
	if(strcasecmp(block_name,thisblock->blockname)==0) {
	  *lastblockp = nextblock;
	  free(thisblock);
	}
	return(S_SUCCESS);
      } else {
	free(thisblock);
      }
      thisblock = nextblock;
    }
    if(block_name) {
      fprintf(stderr,"Hist block %s not found\n",block_name);
      return(S_FAILURE);
    }
    thHBlockListP = (thHBlockList *) NULL;
  }
  return(S_SUCCESS);
}
int thGetHistID(char *name)
{
  daVarStruct *varp;

  if(daVarLookupPWithClass(name,histclasslist,&varp) == S_SUCCESS){
    if(varp->type = DAVARINT){
      return((int) *((DAINT *) varp->varptr));
    }
  }
  return(0);			/* No hist ID found */
}
thStatus thHistAliasWrite(char *fname)
{
  FILE *OFILE;
  thHBlockList *thisblock;
  thHistSpecList *thishist;

  if(thHBlockListP == 0) {
    return(S_FAILURE);		/* No tests defined */
  }

  if((OFILE = fopen(fname,"w")) == NULL) {
    fprintf(stderr,"(thHistAliasWrite) Failed to open %s for write\n",fname);
    return(S_FAILURE);
  }
  thisblock = thHBlockListP;
  while(thisblock){
    thishist = thisblock->var->opaque;
    while(thishist){
      char *name;
      name = strchr(thishist->varname->name,'.');
      if(name) {
	name++;
	fprintf(OFILE,"alias/create %s %d\n",name
		,*((DAINT *) thishist->varname->varptr));
      } else {
	fprintf(stderr,"No . in histogram name %s\n",thishist->varname->name);
      }
      thishist = thishist->next;
    }
    thisblock = thisblock->next;
  }
  fclose(OFILE);
  return(S_SUCCESS);
}
      
int thhstexe_()
{
  int A0;
  A0 = thExecuteHists(0);
  return A0;
}
/* This routine needs to append block.hist in front.  Ultimately it
should probably only put this in front if there isn't a block. in the
name already??  */
int thhstexeb_(char *A1,unsigned C1)
{
  int A0;
  char *B1;
  A0 = thExecuteHists((!*(int *)A1)?0:memchr(A1,'\0',C1)?A1:
                      (memcpy(B1=malloc(C1+1),A1,C1),B1[C1]='\0'
                       ,kill_trailing(B1,' ')));
  if(B1) free(B1);
  return A0;
}
