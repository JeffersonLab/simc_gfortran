/*-----------------------------------------------------------------------------
 * Copyright (c) 1994 Southeastern Universities Research Association,
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
 *  Report generator
 *
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 * $Log: thReport.c,v $
 * Revision 1.1  2009/01/23 13:34:01  gaskelld
 * Initial revision
 *
 * Revision 1.3  2002/07/31 20:07:48  saw
 * Add files for ROOT Trees
 *
 * Revision 1.2  1999/11/04 20:34:06  saw
 * Alpha compatibility.
 * New RPC call needed for root event display.
 * Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 * Revision 1.11  1999/03/01 19:56:52  saw
 * Add Absoft Fortran stuff
 *
 * Revision 1.10  1997/05/30 14:06:17  saw
 * Fix some memory leaks
 *
 * Revision 1.9  1995/05/09 17:03:13  saw
 * Have threp and threpa return S_FAILURE when files can't be opened.
 *
 * Revision 1.8  1995/04/10  15:40:22  saw
 * Force to floating point format on integer overflow from thEvalImed.
 * Add a missing error return to thReportFromVar.
 *
 * Revision 1.7  1994/11/07  14:27:09  saw
 * Bug fixes
 *
 * Revision 1.2  1994/06/14  21:13:28  saw
 * Add fortran calls.  Strip trailing spaces off of fortran strings.
 *
 * Revision 1.1  1994/06/13  13:27:25  saw
 * Initial revision
 *
*/

#define REPORTSTR "report"
#define VLISTSTR "vlist."
#define SPECIALTYPE 9999
#include <stdio.h>
#include "daVar.h"
#include "th.h"
#include "thInternal.h"
#include "cfortran.h"

char *reportclasslist[]={TESTSTR,PARMSTR,EVENTSTR,0};

extern thStatus thReportFd(char *varname, FILE *output);
thStatus thReportFromVar(daVarStruct *bvar,FILE *output);

#ifdef AbsoftUNIXFortran
int threp
#else
int threp_
#endif
(char *repname, char *filename
	    ,unsigned l_repname, unsigned l_filename){
  int A0;
  char *BR=0;
  char *BF=0;
  FILE *fd;

  fd = fopen(((!*(int *)filename)?0:memchr(filename,'\0',l_filename)?filename:
		 (memcpy(BF=(char *) malloc(l_filename+1),filename,l_filename)
		  ,BF[l_filename]='\0',kill_trailing(BF,' '))),"w");
/* Need to add an error message if fopen fails */
  if(fd)
    A0 = thReportFd(((!*(int *)repname)?0:memchr(repname,'\0',l_repname)
		     ?repname:(memcpy(BR=(char *) malloc(l_repname+1)
				      ,repname,l_repname)
			       ,BR[l_repname]='\0',kill_trailing(BR,' '))),fd);
  else {
    A0 = S_FAILURE;
    fprintf(stderr,"Failed to open file %s\n",BF);
  }
  if(fd) fclose(fd);
  if(BF) free(BF);
  if(BR) free(BR);
  return(A0);
}
#ifdef AbsoftUNIXFortran
int threpa
#else
int threpa_
#endif
(char *repname, char *filename
	    ,unsigned l_repname, unsigned l_filename){
  int A0;
  char *BR=0;
  char *BF=0;
  FILE *fd;

  fd = fopen(((!*(int *)filename)?0:memchr(filename,'\0',l_filename)?filename:
		 (memcpy(BF=(char *) malloc(l_filename+1),filename,l_filename)
		  ,BF[l_filename]='\0',kill_trailing(BF,' '))),"a");
/* Need to add an error message if fopen fails */
  if(fd)
    A0 = thReportFd(((!*(int *)repname)?0:memchr(repname,'\0',l_repname)
		     ?repname:(memcpy(BR=(char *) malloc(l_repname+1)
				      ,repname,l_repname)
			       ,BR[l_repname]='\0',kill_trailing(BR,' '))),fd);
  else
    A0 = S_FAILURE;
  if(fd) fclose(fd);
  if(BF) free(BF);
  if(BR) free(BR);
  return(A0);
}
  
thStatus thBookReports(daVarStruct *var)
{
  /* Dummy routine.  No booking action needed for reports */
  return(S_SUCCESS);
}
thStatus thReportFd(char *block_name, FILE *output)
{
  static daVarStruct *bvar;
  static char *fullvarname;

  fullvarname = malloc(strlen(BLOCKSTR)+strlen(REPORTSTR)
		       +strlen(block_name)+3);
  strcpy(fullvarname,BLOCKSTR);
  strcat(fullvarname,".");
  strcat(fullvarname,REPORTSTR);
  strcat(fullvarname,".");
  strcat(fullvarname,block_name);
  
  if(daVarLookupP(fullvarname,&bvar) != S_SUCCESS){
    fprintf(stderr,"Failed to find %s\n",fullvarname);
    free(fullvarname);
    return(S_FAILURE);
  }
  free(fullvarname);
  return(thReportFromVar(bvar,output));
}

thStatus thReportFromVar(daVarStruct *bvar,FILE *output)
{
  char *lines,*eol,*s;
  int line_count;
  double dval;	char *dfmt="%f";
  int lval; char *lfmt="%d";
  char *sval; char *sfmt="%s";
  daVarStruct *varp;
  int dtype;
  int expression;
  
  lines = bvar->title;
/*  printf("%s",lines);*/
  (*((DAINT *) bvar->varptr))++; /* Number of times report has been printed */

  line_count = 0;
  while(*lines){
    line_count++;
    eol = strchr(lines,'\n');
    if(!eol) {
      fprintf(STDERR,"L %d: Last line of %s has no newline\n"
	      ,line_count,bvar->name);
    }
    if(*(eol+1)=='\0'){		/* This is the last line */
      if(strcasestr(lines,ENDSTR) == 0)
	fprintf(STDERR,"L %d: Last line of %s is not an END\n"
		,line_count,bvar->name);
      else
	break;
    }
    if(line_count == 1)
      if(strcasestr(lines,BEGINSTR) != 0){
	lines = eol + 1;
	continue;
      } else
	fprintf(STDERR,"First line of %s is not a BEGIN\n",bvar->name);
    /* Ready to process the line, Add continuation lines later */
    s = lines;
    while(s <= eol){
      if(*s != '{') {
	if(*s == '\\') {
	  if(*++s != '\n') fputc(*s++,output);/* Escape chars and newlines */
	  else s++;
	} else fputc(*s++,output);
      } else {
	char *eptr;		/* Expression pointer */
	char *fptr;		/* Format pointer */
	char *eend;		/* Pointer to end of expression */
	char *fend;		/* Pointer to end of format pointer */
	char ceend, cfend;

	eptr = ++s;
	fptr = 0; eend = 0; fend = 0;
	while(s < eol){
	  if(*s == ':') {
	    eend = s; ceend = *s;
	    *s++ = 0;		/* Replace : with null to terminate expr */
	    fptr = s;
	    while(s < eol) {
	      if(*s == '}') {
		fend = s; cfend = *s;
		*s++ = 0;		/* s points to next character to print */
		break;
	      }
	      s++;
	    }
	    if(!fend) {fend = eol; cfend = *eol;}
	    break;
	  } else if (*s == '}') {
	    eend = s; ceend = *s;
	    *s++ = 0;		/* s points to next character to print */
	    break;
	  } else s++;
	}
	if(!eend) {eend = eol; ceend = *eol;}

	/* Look first for atomic values or variable names (could be a string) */
	dtype = 0;		/* Data type of variable or expression */
	expression = 0;
	if(strncasecmp(eptr,VLISTSTR,strlen(VLISTSTR)) == 0) {
	  dtype = SPECIALTYPE;
	} else {
	  char *p;
	  p = eptr;
	  if(strchr(".0123456789",*p)) expression = 1;
	  while(*p && !expression){
	    if(strchr("*/-+%=!&~,<>[]()",*p++)) expression=1;
	  }
	  if(!expression) if(daVarLookupPWithClass(eptr,reportclasslist,&varp)==S_SUCCESS) {
	    switch(varp->type)
	      {
	      case DAVARSTRING:
		sval = varp->varptr;
		dtype = DAVARSTRING;
		break;
	      case DAVARFSTRING:	/* Fortran string */	
		sval = (char *) malloc(varp->size+1);
		strncpy(sval,varp->varptr,varp->size); /* Copy and  */
		sval[varp->size] = '\0'; /* Null terminate */
		{
		  char *strip;

		  strip = sval+varp->size-1;
		  while(*strip==' ' && strip>sval) strip--;
		  *(strip+1) = '\0';
		}
		printf("%s:",sval);
		dtype = DAVARFSTRING;
		break;
	      default:
		break;		/* Let thEvalImed handle other types */
	      }
	  }
	}
	if(dtype == 0){
	  int stat;
	  stat = thEvalImed(eptr,&dval, &lval);
	  if(stat==S_SUCCESS || stat==S_INTOVF){
	    /* printf("%f %d\n",dval,lval);*/
	    if(stat==S_INTOVF) {
	      dtype = DAVARDOUBLE;
	    } else if(dval == (double)  lval) {
	      dtype = DAVARINT;
	    } else {
	      dtype = DAVARDOUBLE;
	    }
	  } /* else dtype stays zero and everything between {}'s is printed */
	}
	if(fptr&&(dtype==DAVARINT || dtype==DAVARDOUBLE)) { 
	  char *sf;
	  sf = fptr;
	  while(*sf) {
	    if(*sf == '\\'){
	      if(*(sf+1)) sf++;
	    } else if(*sf=='%') {
	      if(*++sf) {
		if(*sf != '%') break;
	      }
	    }
	    sf++;
	  }
	  while(*sf){
	    if(strchr("eEfgG",*sf)) {
	      dtype = DAVARDOUBLE;
	      break;
	    } else if(strchr("diouxXDOUcp",*sf)) { 
	      dtype = DAVARINT;
	      break;
	    } else if(*sf=='s') {
	      fprintf(STDERR,"(thReport): %s not allowed with numbers\n");
	      break;
	    }
	    sf++;
	  }
	}
	switch(dtype) {
	case DAVARINT:
	  if(!fptr) fptr = lfmt;
	  fprintf(output,fptr,lval);
	  break;
	case DAVARDOUBLE:
	  if(!fptr) fptr = dfmt;
	  fprintf(output,fptr,dval);
	  break;
	case DAVARSTRING:
	case DAVARFSTRING:
	  if(!fptr) fptr = sfmt;
	  fprintf(output,fptr,sval);
	  if(dtype==DAVARFSTRING) free(sval);
	  break;
	case SPECIALTYPE:
	  thReportSpecial(output,eptr);
	  break;
	default:
	  fprintf(output,"{%s",eptr);
	  if(fptr) {
	    fprintf(output,":%s}",fptr);
	  } else {
	    fprintf(output,"}",fptr);
	  }
	}
	if(eend) {*eend = ceend;} /*s = eend+1;*/
	if(fend) {*fend = cfend;} /*s = fend+1;}*/
      }
    }
    lines = eol + 1;
  }
  return(S_SUCCESS);
}
thReportSpecial(FILE *output, char *eptr)
{
  char *vptr;			/* Pointer to pattern to match */
  char **vlist;
  daVarStruct *varp;
  int count;
  int i;
  char ftemp[20];
  int len;

  vptr = eptr + strlen(VLISTSTR);
  daVarList(vptr,&vlist,&count);
  for(i=0;i<count;i++) {
    if(daVarLookupP(vlist[i],&varp)==S_SUCCESS) {
      if(varp->flag & DAVAR_REPOINTOK)
	fputc('*',output);	/* CTP created */
      else
	fputc(' ',output);	/* Explicitely registered */
      if(strchr(vptr,'*') || strchr(vptr,'?')) {
	fprintf(output,"%-30.30s",vlist[i]);
      } else {
	fprintf(output,"%-30.30s",vlist[i]+strlen(vptr));
      }
      if(varp->size == 1){
	switch(varp->type)
	  {
	  case DAVARINT:
	    fprintf(output,"%-12i",*((DAINT *)(varp->varptr)));
	    break;
	  case DAVARFLOAT:
	    fprintf(output,"%-12f",*((DAFLOAT *)(varp->varptr)));
	    break;
	  case DAVARDOUBLE:
	    fprintf(output,"%-12lf",*((DADOUBLE *)(varp->varptr)));
	    break;
	  case DAVARSTRING:
	  case DAVARFSTRING:
	    fprintf(output,"%1.1s           ",((char *)(varp->varptr)));
	    break;
	  }
	fprintf(output," %-36.36s",varp->title);
      } else {
	switch(varp->type)
	  {
	  case DAVARINT:
	  case DAVARFLOAT:
	  case DAVARDOUBLE:
 	    if(varp->type==DAVARINT) fputs("I*4",output);
	    else if(varp->type==DAVARFLOAT) fputs("R*4",output);
	    else if(varp->type==DAVARDOUBLE) fputs("R*8",output);
	    sprintf(ftemp,"(%d)",varp->size);
	    fputs(ftemp,output);
	    len = strlen(ftemp)+3;
	    while(len++<12) fputc(' ',output);
	    fprintf(output," %-36.36s",varp->title);
	    break;
	  case DAVARSTRING:
	    fprintf(output,"%-48.48s",((char *)(varp->varptr)));
	    break;
	  case DAVARFSTRING:
	    len = varp->size;
	    if(len>49) len=48;
	    sprintf(ftemp,"%%-%d.%ds",len,len);
	    fprintf(output,ftemp,((char *)(varp->varptr)));
	    break;
	  }
      }
    }
    fprintf(output,"\n");
  }
  daVarFreeList(vlist);
}
