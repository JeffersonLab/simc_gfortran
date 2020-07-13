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
 *  Read CTP configuration files saving them as the titles of registered 
 *  block variables.  thBook books all unbooked blocks.
 *
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: thLoad.c,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.5  2004/07/13 15:04:53  saw
 *   Get count of groups on begin line correct when there are other attributes
 *
 *   Revision 1.4  2004/07/08 20:06:45  saw
 *   CTP Root Trees can exist in input files even if root trees not compiled in.
 *
 *   Revision 1.3  2003/02/21 20:55:24  saw
 *   Clean up some types and casts to reduce compiler warnings.
 *
 *   Revision 1.2  1999/11/04 20:34:06  saw
 *   Alpha compatibility.
 *   New RPC call needed for root event display.
 *   Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *   Revision 1.17  1999/08/25 13:16:06  saw
 *   *** empty log message ***
 *
 *   Revision 1.16  1999/03/01 19:56:12  saw
 *   Work on INPUT_FILE open/close problem
 *
 *   Revision 1.15  1996/07/31 20:35:19  saw
 *   Book all groups instead of default groups.  (Forces booking order to
 *   be same as input files)
 *
 *   Revision 1.14  1996/01/30 15:38:32  saw
 *   Add groups
 *
 *	  Revision 1.13  1995/08/03  13:52:42  saw
 *	  Add a missing return statement
 *
 *	  Revision 1.12  1995/02/04  19:08:29  saw
 *	  Fix problem in finding quotes around include filenames.
 *
 *	  Revision 1.11  1995/01/09  15:22:51  saw
 *	  Add include files, fix assorted paren placement problems that were only
 *	  apparent on ultrix and linux
 *
 *	  Revision 1.10  1994/09/27  19:41:03  saw
 *	  Remove linux dependencies
 *
 *	  Revision 1.9  1994/08/26  13:20:24  saw
 *	  Add DAVAR_REPOINTOK to some flags
 *
 *	  Revision 1.8  1994/07/21  18:43:58  saw
 *	  Add gethit and report to list of valid block classes
 *
 *	  Revision 1.7  1994/06/13  13:08:27  saw
 *	  getblock: fix bug in reading of non standard block types
 *
 *	  Revision 1.6  1994/06/03  18:48:41  saw
 *	  Replace stderr with STDERR
 *
 *	  Revision 1.5  1993/10/18  18:28:28  saw
 *	  Fix blank lines causing failure on HP.
 *
 *	  Revision 1.4  1993/09/10  17:37:07  saw
 *	  *** empty log message ***
 *
 *	  Revision 1.3  1993/09/10  17:32:59  saw
 *	  Fix up improper use of blockname in getblock after freeing blockname.
 *
 *	  Revision 1.2  1993/05/10  21:18:56  saw
 *	  Fix header
 *
 */
/*
   Read in lines of text and save them in a long string.  Test, histogram
   booking and parameter lines will be saved in a long carriage return
   separated string.
*/

/* Scheme.  Call a routine with an open file.  It reads the file, ignoreing
   input until a "begin" is found.  It then reads until an end line is found.
   It then returns a string that contains the whole contents of the file
   between the begin and the end.

This is loading phase.  When loading is done, user must call a book all
routine.

*/

#include <stdio.h>
#include "daVar.h"
#include "th.h"
#include "thInternal.h"
#include "thUtils.h"
#include "cfortran.h"
#include "thGroup.h"

#define MAXLINELENGTH 512
#define MAXTOKS 20
#define max(a,b) 		(a<b ? b : a)
#define INCLUDESTR "#include"

/*extern thStatus thBookGethits(daVarStruct *var);
extern thStatus thBookTests(daVarStruct *var);
extern thStatus thBookHists(daVarStruct *var);
extern thStatus thBookReports(daVarStruct *var);
extern thStatus thBookGroup(daVarStruct *var);
extern daVarStatus thWHandler();
extern daVarStatus thRHandler();*/

char *types[]={PARMSTR, GETHITSTR, TESTSTR, HISTSTR, UHISTSTR,
               TREESTR,
               REPORTSTR, 0};
thStatus (*hooks[])()={thLoadParameters, /*thBookGethits,*/ thBookTests,
		       /*thBookHists, thBookHists,*/
                       thBookTree,
                       thBookReports, 0};

char *qualifiers[]={"obey", "read", "write", 0};
int qualflags[]={DAVAR_OBEYMF, DAVAR_READONLY, DAVAR_READWRITE, DAVAR_OBEYMF};

struct thBookList {
  char *classname;		/* Variable name containing the whook */
  struct thBookList *next;
};
typedef struct thBookList thBookList;

thBookList *thBookListP=NULL;

thStatus thSetBlockClasses();

/* Fortran Interface */
FCALLSCFUN0(INT,thOBook,THOBOOK,thobook)
FCALLSCFUN1(INT,thLoad,THLOAD,thload,STRING)

struct thFdList {
  FILE *fd;
  struct thFdList *next;
};
typedef struct thFdList thFdList;
thFdList *thFdListP=NULL;

FILE *myopen(char *fname);
FILE *myclose(FILE *fd);
void mycloseall();
FILE *do_include(char *fname);

thStatus thOBook()
/* Book all the parameter, tests, and histograms and anything else that is
   in the booking order list.
   This is the pre-group booking routine.  Within each class it booked blocks
   in alphabetical order.
*/
{
  thBookList *booklist;

  /* Need to find a way to make sure these are called by anything else
     that might do booking.  Perhaps need a thInit routine. */
/*  printf("In thBook\n");*/
  if (!thBookListP) {		/* Booking order not defined */
    if(thSetBlockClasses() != S_SUCCESS){
      fprintf(STDERR,"Failed to set the th Block Class handlers\n");
      return(S_FAILURE);
    }
  }

  booklist = thBookListP;
  while(booklist){
    char *prefix; char **blist0; char **blist; int count;
    daVarStruct var, *bvar;

    if(daVarLookup(booklist->classname,&var) != S_SUCCESS){
      fprintf(STDERR,"Failed to find class variable %s\n",booklist->classname);
      return(S_FAILURE);
    }
    
    prefix = (char *) malloc(strlen(booklist->classname) + 2);
    strcpy(prefix,booklist->classname);
    strcat(prefix,".");
    
    daVarList(prefix,&blist0,&count);
    blist = blist0;
    while(count-- > 0){
      thStatus stat;

      /* printf("Booking %s\n",*blist); */
      if(daVarLookupP(*blist,&bvar) != S_SUCCESS){
	fprintf(STDERR,"Failed to find %s\n",*blist);
	return(S_FAILURE);
      }
      if(!bvar->varptr) {	/* Only book those not already booked */
	bvar->varptr = (DAINT *) malloc(sizeof(DAINT));
	*((DAINT *) bvar->varptr) = 0; /* Initializec execution counter */
	stat = ((thStatus (*)()) var.opaque)(bvar);
      }
      blist++;
    }
    
    daVarFreeList(blist0);

    booklist = booklist->next;
  }
  return(S_SUCCESS);
}
thStatus thSetBlockClasses()
{
  thBookList *booklist;
  thBookList **booknext;
  daVarStruct var;
  int i;

  booknext = &thBookListP;
  for(i=0;types[i];i++){
    booklist = *booknext = (thBookList *) malloc(sizeof(thBookList));
    booklist->classname = (char *) 
      malloc(strlen(BLOCKSTR) + strlen(types[i]) + 2);
    strcpy(booklist->classname,BLOCKSTR);
    strcat(booklist->classname,".");
    strcat(booklist->classname,types[i]);
    booklist->next = (thBookList *) NULL;
    booknext = &booklist->next;
    
    var.name = booklist->classname;
    var.title = 0;/*""*/
    var.type = DAVARINT;
    var.varptr = (DAINT *) malloc(sizeof(DAINT));
    *((DAINT *) var.varptr) = i; /* Booking order */
    var.size = 1;
    var.whook = thWHandler;
    var.rhook = thRHandler;
    var.flag = DAVAR_READWRITE | DAVAR_REPOINTOK;
    var.opaque = (void *) hooks[i];
    if(daVarRegister((int) 0, &var) == S_FAILURE){
      fprintf(STDERR,"Failed to register %s\n",var.name);
      return(S_FAILURE);
    }
  }
  return(S_SUCCESS);
}

thStatus thLoad(char *fname)
/* Open the file and read in all the blocks of test, histogram, or parameter
   code.  Put the contents of each block into the title field of variables
   with the names block.TYPE.NAME.
   May want to change this later to take a file descriptor, not a file
   name (so that pipes can be used.)
   */
{
  FILE *INPUT_FILE;		/* Descriptor file */
  char *varname, *vartitle, **grouplist;
  int qualid,ig;

  if((INPUT_FILE = myopen(fname))==NULL) {
    fprintf(STDERR,"(thLoad) Failed to open %s\n",fname);
    return(S_FAILURE);
  }
  while(getblock(&INPUT_FILE,&varname,&vartitle,&qualid
		 ,&grouplist)==S_SUCCESS){
    daVarStruct var,*varp;

    if(daVarLookup(varname,&var) == S_SUCCESS){	/* Variable exists */
      if(var.type != DAVARINT){
	fprintf(STDERR,"Type for variable %s must be integer\n",varname);
	return(S_FAILURE);
      }
      /* Free and zero the varptr to indicate that this block has not been
	 booked. */
      if(var.varptr) {
	free(var.varptr);
	var.varptr = 0;
      }
      /*      printf("XX: Found %s\n",varname);*/
    } else {
      var.name = varname;
      var.type = DAVARINT;
      var.varptr = 0;		/* Null pointer means not yet booked */
      var.size = 1;
      var.opaque = 0;		/* Booking routine may add opaque data */
      var.rhook = 0;
      var.whook = 0;
      var.flag = DAVAR_READONLY | DAVAR_REPOINTOK;
    }
    var.title = vartitle;	/* The lines to book test, hists, pars etc. */
    var.flag = qualflags[qualid] | DAVAR_REPOINTOK;
/*    printf("XX: Registering %s\n",varname);*/
    if(daVarRegister((int) 0, &var) == S_FAILURE){ /* Create block desciptor */
      fprintf(STDERR,"Failure to register %s\n",varname);
      fclose(INPUT_FILE);
      return(S_FAILURE);
    }
    daVarLookupP(var.name,&varp);
    /* Go attach this block to each group, create var for group
       if it doesn't exist, check that block is not in group before
       adding it at the end.  Do we do the group stuff before or after
       the block is fully read in? */
    for(ig=0;grouplist[ig];ig++) {
      daVarStruct varg;
      daVarStructList *blist; /* Block list */
      thGroupOpaque *opqptr;
      
      if(daVarLookup(grouplist[ig],&varg) != S_SUCCESS){ /* Variable exists */
	varg.name = grouplist[ig];
	varg.type = DAVARINT;
	varg.varptr = (void *) malloc(sizeof(DAINT));
	*(DAINT *)varg.varptr = 0;
	varg.size = 1;
	varg.opaque = (void *) malloc(sizeof(thGroupOpaque));
	opqptr = (thGroupOpaque *) varg.opaque;
	thInitGroupOpaque(varg.name,opqptr); /* Fill the elelments of the structure */
	varg.rhook = 0;
	varg.whook = 0; /* Can put something neat here */
	varg.flag = DAVAR_READONLY | DAVAR_REPOINTOK;
	varg.title = 0;	/* Would be nice to put group descriptions here */
      } else
	opqptr = (thGroupOpaque *) varg.opaque;
      blist = opqptr->blocklist;
      thAddVarToList(&blist, varp);
      if(!opqptr->blocklist) {	/* Create group list var if it doesn't exist */
	opqptr->blocklist = blist;
	if(daVarRegister((int) 0, &varg) == S_FAILURE){
	  fprintf(STDERR,"Failure to register %s\n",varg.name);
	  fclose(INPUT_FILE);
	  return(S_FAILURE);
	}
      }
/*      printf("X %s\n",grouplist[ig]);*/
    }
  }
  /* Are we covering up something.  Should we be able to get to here with
     a null value for INPUT_FILE?  Does it mean we have left a file open? */
  if(INPUT_FILE) fclose(INPUT_FILE); 

#if 0
  {
    daVarStruct varg;
    daVarStructList *blist;

    if(daVarLookup("group.test.all",&varg)==S_SUCCESS) {
      printf("All Test group blocks:\n");
      blist = ((thGroupOpaque *)varg.opaque)->blocklist;
      while(blist) {
	printf("  %s\n",blist->varp->name);
	blist = blist->next;
      }
    }
    if(daVarLookup("group.hist.all",&varg)==S_SUCCESS) {
      printf("All Histogram group blocks:\n");
      blist = ((thGroupOpaque *)varg.opaque)->blocklist;
      while(blist) {
	printf("  %s\n",blist->varp->name);
	blist = blist->next;
      }
    }
    if(daVarLookup("group.parm.all",&varg)==S_SUCCESS) {
      printf("All Parameter group blocks:\n");
      blist = ((thGroupOpaque *)varg.opaque)->blocklist;
      while(blist) {
	printf("  %s\n",blist->varp->name);
	blist = blist->next;
      }
    }
  }
#endif
  return(S_SUCCESS);
}
int getblock(FILE **TEST_FILE, char **varname, char **vartitle, int *qualid
	     ,char ***grouplist)
/*
  Scan the file until a begin line is found.  When "begin xxx" is found, read
  all the lines up until end xxx is found, stuffing them into a line pointed
  to by vartitle.

  This code is really ugly.  Can't I make it better?
*/
{

  static char *start=0;
  static char *last=0;
  static char *next=0;
  char oline[MAXLINELENGTH];
  char line[MAXLINELENGTH];
  static char *blockname=0;
  static char **groups=0;
  char *toks[MAXTOKS];
  int toklens[3], ntoks, itok, qualifier;
  char *stat;
  char *s;
  static char *thistype=0;
  int endflag, len;
  static int maxgroups=0;

  
  int i;
  
  static int EOF_PENDING=0;	/* For cases where EOF happened in a file */

  /* Scan for a begin line */

  if(EOF_PENDING){
    EOF_PENDING = 0;
    if(start) free(start);
    start = last = 0;
    return(S_FAILURE);		/* Should be another code for EOF */
  }
  toks[0] = line;		/* Protect against blank lines */
  /* Scan for "begin" statement or #include */
 keep_looking_for_begin:
  while((stat=fgets(oline,MAXLINELENGTH,*TEST_FILE))){
    strcpy(line,oline);
    s = line;
    while(*s && *s != '\n' && *s != COMCHAR) s++;
    *s = 0;
    {				/* Break into tokens */
      s = line;
      while(*s && isspace(*s)) s++; /* Skip white space to first character */
      ntoks = 0;
      while(*s){
/*	printf("|%s|\n",s);*/
	toks[ntoks++] = s;
	while(*s && !isspace(*s)) s++; /* Skip over token */
	if(*s){
	  *s++ = 0;
	  while(*s && isspace(*s)) s++;	/* Skip white space to next */
	}
      }
    }
    itok = 0;
    if(strcasecmp(toks[itok++],BEGINSTR) == 0) break;
    if(strcasecmp(toks[0],INCLUDESTR) == 0) *TEST_FILE = do_include(toks[1]);
  }

  if(stat == NULL) {
    *TEST_FILE = myclose(*TEST_FILE);
    if(*TEST_FILE) {		/* Not the true end of file */
      goto keep_looking_for_begin;
    } else {
      if(start) free(start);
      start = last = 0;
      return(S_FAILURE);		/* Should be another code for EOF */
    }
  }

  /* Begin found at this point */


  for(i=0;qualifiers[i];i++){
    if(strncasecmp(toks[itok],
		      qualifiers[i],
		      strlen(qualifiers[i])) == 0){
      itok++;
      break;
    }
  }
  *qualid = i;

  if(thistype)
    thistype = (char *) realloc(thistype,strlen(toks[itok])+1);
  else
    thistype = (char *) malloc(strlen(toks[itok])+1);


  for(i=0;types[i];i++){
    if(strncasecmp(toks[itok],types[i],strlen(types[i])) == 0){
      break;
    }
  }
  if(types[i] == 0) {
    fprintf(STDERR,"Unknown block type %s\n",toks[itok]);
    strcpy(thistype,toks[itok]);
    /* Need provision for user defined block types? */
  } else {
    strcpy(thistype,types[i]);
  }
  itok++;
  {
    int len;
    len = strlen(BLOCKSTR) + strlen(thistype) + strlen(toks[itok]) + 3;
    if(blockname)
      blockname = (char *) realloc(blockname,len);
    else
      blockname = (char *) malloc(len);
  }
  strcpy(blockname,BLOCKSTR);
  strcat(blockname,".");
  strcat(blockname,thistype);
  strcat(blockname,".");
  strcat(blockname,toks[itok]);

/*  printf("|%s|\n",blockname);*/
  /* Look for the group id */
  itok++;
  if(!groups) {			/* Initialize groups array */
    groups = (char **) malloc(3*(sizeof(char *)));
    groups[0] = 0;
    groups[1] = 0;
    groups[2] = 0;
    maxgroups = 2;
  }
  /* Put it in the all group */
  if(groups[0]) free(groups[0]);
  groups[0] = (char *) malloc(strlen(GROUPSTR)+strlen(thistype)
			      +strlen(ALLGRPSTR)+3);
  strcpy(groups[0],GROUPSTR);
  strcat(groups[0],".");
  strcat(groups[0],thistype);
  strcat(groups[0],".");
  strcat(groups[0],ALLGRPSTR);

  if(ntoks-itok+1 > maxgroups) {
    int ig;
    groups = (char **) realloc(groups,(ntoks-itok+2)*(sizeof(char *)));
    for(ig=maxgroups+1;ig<ntoks-itok+2;ig++) {
      groups[ig] = 0;		/* Zero new elements, */
    }	/* Makeing sure realloc will work */ 
    maxgroups = ntoks-itok+1;
  }
  if(ntoks > itok) {
    int ig=1;			/* Skip over all group */
    for(;itok < ntoks;itok++) {		/* Scan for the group keyword */
      if(strncasecmp(toks[itok],GROUPSTR,strlen(GROUPSTR))==0) {
	char *grpnam;
	if(grpnam = strchr(toks[itok],'=')) {
	  grpnam++;
	  if(groups[ig]) free(groups[ig]);
#if 0
	  /* For now, groups can only contain a single block type */
	  if(strchr(grpnam,'.')) { /* Class assigned already */
	    groups[ig] =  (char *) malloc(strlen(GROUPSTR)+
					  strlen(toks[itok])+2);
	    strcpy(groups[ig],GROUPSTR);
	    strcat(groups[ig],".");
	    strcat(groups[ig],grpnam);
	  } else { /* Prepend block type */
#endif
	    
	    groups[ig] =  (char *) malloc(strlen(GROUPSTR)+
					  strlen(grpnam)+strlen(thistype)+3);
	    strcpy(groups[ig],GROUPSTR);
	    strcat(groups[ig],".");
	    strcat(groups[ig],thistype);
	    strcat(groups[ig],".");
	    strcat(groups[ig],grpnam);
#if 0
	  }
#endif
	}	/* else { ignore if = is missing } */
	ig++;
      }
    }
    if(groups[ig]) {		/* Null terminate list */
      free(groups[ig]);
      groups[ig] = 0;
    }
  } else {			/* No groups specified, use default group */
    if(groups[1]) free(groups[1]);
    groups[1] = (char *) malloc(strlen(GROUPSTR)+strlen(thistype)
				+strlen(DEFAULTGRPSTR)+3);
    strcpy(groups[1],GROUPSTR);
    strcat(groups[1],".");
    strcat(groups[1],thistype);
    strcat(groups[1],".");
    strcat(groups[1],DEFAULTGRPSTR);
    if(groups[2]) {		/* Null terminate list */
      free(groups[2]);
      groups[2] = 0;
    }
  }

  strcpy(line,oline);		/* Save the line */
  endflag = 0;
    
  next = start;			/* Initialize */
  while(1){
    len = strlen(line);
    if(next+len >= last){
      int newsize;
/*      char *src,*dst,*newstart;*/
      
/*      src = start;*/
      newsize = max((last-start) + MAXLINELENGTH,(next-start) + len + 1);
/*      newstart = dst = (char *) malloc(newsize);
      while(src < next) *dst++ = *src++;
      if(start) free(start);
      last = newstart + newsize;
      start = newstart;*/
      if(start) {
	char *newstart;
	newstart = (char *) realloc(start,newsize);
	next = next + (newstart-start);
	start = newstart;
      }  else {
	start = (char *) malloc(newsize);
	next = start;
      }
/*      next = next + (dst-src);*/
    }
    strcpy(next,line);
    next += len;

    if(endflag) break;
  reread:
    if((stat=fgets(line,MAXLINELENGTH,*TEST_FILE))){
      s = line;
      while(*s && isspace(*s)) s++;
      if(strncasecmp(s,INCLUDESTR,strlen(INCLUDESTR)) == 0) {
	*TEST_FILE = do_include(s + strlen(INCLUDESTR));
	goto reread;		/* Get another line */
      }
      if(strncasecmp(s,ENDSTR,strlen(ENDSTR))!=0) continue;
      s += strlen(ENDSTR)+1;
      while(*s && isspace(*s)) s++;
      if(strncasecmp(s,thistype,strlen(thistype))!=0) continue;
/* Don't require block name to be restated
      s += strlen(thistype)+1;
      while(*s && isspace(*s)) s++;
      if(strncasecmp(s,BLOCKSTR,strlen(BLOCKSTR))!=0) continue;
*/
      endflag = 1;
    } else {			/* EOF */
      *TEST_FILE = myclose(*TEST_FILE);
      if(!*TEST_FILE) {		/* Test for true end of file */
	strcpy(line,"end ");
	strcat(line,s);
	strcat(line," block\n");
	endflag = 1;			/* Terminate after writing end block */
      } else {
	goto reread;		/* Get another line */
      }
    }
  }

  *varname = blockname;
  *vartitle = start;
  *grouplist = groups;

  return(S_SUCCESS);
}

FILE *myopen(char *fname)
{
  /* Opens  a new file, saving the last file descriptor on a stack. */
  FILE *fd;
  thFdList *next;

  fd = fopen(fname,"r");
  if(fd) {
    next = thFdListP;
    thFdListP = (thFdList *) malloc(sizeof(thFdList));
    thFdListP->fd = fd;
    thFdListP->next = next;
  }
  return(fd);
}
  
FILE *myclose(FILE *fd)
{
  /* Close a file.  If that fd is at top of stack, remove it from stack */
  thFdList *next;

  if(thFdListP->fd == fd) {
    next = thFdListP->next;
    free(thFdListP);
    thFdListP = next;
    fclose(fd);
    if(thFdListP) return(thFdListP->fd);
    else return(0);
  } else {
    fclose(fd);
    return(0);
    /* Should probably return some kind of error because fd was not what
       was expected. */
  }
}
void mycloseall()
{
  /* Close all files on stack */
  thFdList *next;

  while(thFdListP) {
    fclose(thFdListP->fd);
    next = thFdListP->next;
    free(thFdListP);
    thFdListP = next;
  }
}

FILE *do_include(char *include_string)
{
  /* Skip over whitespace:
     Is first character a quote or double quote, if so, interpret string a
     filename.  Else interpret string as a variable and evaluate it.  If it
     is not of string type, assume it is a filename. */
  char *fstring,*fname;
  char *s;
  char quotechar;
  int len;
  FILE *fd;

  s = include_string;
  while(*s && isspace(*s)) s++;	/* Skip whitespace */
  quotechar = 0;
  if(*s == QUOTECHAR1 || *s == QUOTECHAR2) {
    quotechar = *s++;
  }
  fstring = s;
  while(*s && *s != quotechar && !isspace(*s) && (quotechar || *s != COMCHAR)){
    s++;
  }
  len = s-fstring;
  fname = (char *) malloc(len+1);
  strncpy(fname,fstring,len);
  fname[len] = '\0';
  if(quotechar) {
    fd = myopen(fname);
  } else {			/* TODO: Interpret as CTP variable */
    fd = myopen(fname);
  }
  if(!fd) {
    fprintf(STDERR,"(thLoad) Failed to open '%s'\n",fname);
    if(thFdListP) fd = thFdListP->fd;
    fd = 0;
  }
  free(fname);
  return(fd);
}
