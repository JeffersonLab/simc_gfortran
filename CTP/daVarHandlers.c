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
 *  Handlers for RPC services.
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: daVarHandlers.c,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.1  1998/12/07 22:11:09  saw
 *   Initial setup
 *
 *	  Revision 1.7  1995/01/09  15:59:53  saw
 *	  Move location of kill_trailing function in the file
 *
 *	  Revision 1.6  1994/11/07  14:13:49  saw
 *	  Finish name change by include daVarHandlers.h instead of daVarServ.h
 *
 *	  Revision 1.5  1994/10/16  21:37:57  saw
 *	  Add RPC support for Fortran strings
 *	  Change name from daVarServ to daVarHandlers
 *
 *	  Revision 1.4  1994/06/03  21:04:45  saw
 *	  Add RPC support for doubles
 *
 *	  Revision 1.3  1993/08/12  15:03:49  saw
 *	  Add #include <rpc/rpc.h>
 *
 *	  Revision 1.2  1993/05/10  21:05:05  saw
 *	  Fix description
 *
 *	  Revision 1.1  1993/05/10  21:02:55  saw
 *	  Initial revision
 *
*/

/*
If an index is put on a variable for a read, just that array element is read.
If an index is put on a variable for a write, that array element is used as
the starting point for however long a list of values is contained in the
write call.  If the write  would go past the end of the array, nothing is
written and an error is returned.

A write of a string to a value (of the type string) will be allowed if the
string to be written has a length <= to the allocated size specified in
the ->size.

If a string is written to a title, the existing string will be free'd and
the new string written.  The new string may be any size.

*/

#include <string.h>
#include <rpc/rpc.h>

#include "daVar.h"
#include "daVarRpc.h"
#include "daVarHandlers.h"

char *daVarRAtrList[]={"value","title","size","flag","type","watr","ratr",0};
char *daVarWAtrList[]={"value","title",0};


daVarStatus daVarClassFind(char *name, daVarStruct **varp);
daVarStatus daVarAttributeFind(char *name, daVarStruct *varclass,
			       daVarStruct **varp, char **attribute,
			       int *index);
void daVarRhandler(char *name, daVarStruct *varclass, any *retval);
daVarStatus daVarWhandler(char *name,daVarStruct *varclass,any *setval);
char *daVarMakeRAtrList();
char *daVarMakeWAtrList();
void daVarCopyAlist(char *nllist,char **list);

static char *kill_trailing(char *s, char t)
{
  char *e; 
  e = s + strlen(s);
  if (e>s) {                           /* Need this to handle NULL string.*/
    while (e>s && *--e==t);            /* Don't follow t's past beginning. */
    e[*e==t?0:1] = '\0';               /* Handle s[0]=t correctly.       */
  }
  return s;
}

daVarStatus daVarWriteVar(char *name, any *setval)
{
  daVarStruct *varclass;

  if(daVarClassFind(name,&varclass)!=S_SUCCESS){
    return(S_FAILURE);

  }
  if(varclass->whook){
    return((varclass->whook)(name,varclass,setval));
  } else {
    return(daVarWhandler(name,varclass,setval));
  }
}


void daVarReadVar(char *name, any *retval)
{
  daVarStruct *varclass;

  if(daVarClassFind(name,&varclass)!=S_SUCCESS){
    retval->valtype = DAVARERROR_RPC;
    retval->any_u.error = S_FAILURE;
    return;
  }
  if(varclass->rhook){
    (varclass->rhook)(name,varclass,retval);
  } else {
    daVarRhandler(name,varclass,retval);
  }
}

daVarStatus daVarWhandler(char *name,daVarStruct *varclass,any *setval)
/* The default Write handler */
{
  daVarStruct *varp;
  char *attribute;
  daVarStatus status;
  int index;

  status = daVarAttributeFind(name, varclass, &varp, &attribute, &index);
  if(status == S_SUCCESS)
    status = daVarRegWatr(varp, attribute, index, setval);
  /* A special handler would check more attributes if status != SUCCESS */
  return(status);
}

void daVarRhandler(char *name, daVarStruct *varclass, any *retval)
/* The default Read handler */
{
  daVarStruct *varp;
  char *attribute;
  daVarStatus status;
  int index;

  status = daVarAttributeFind(name, varclass, &varp, &attribute, &index);
  status = daVarRegRatr(varp, attribute, index, retval);
  /* A special handler would check more attributes if status != SUCCESS */
  return;
}

daVarStatus daVarRegRatr(daVarStruct *varp, char *attribute
			 ,int index, any *retval)
/* Regular Read Attribute routine.  Returns failure if attribute is not
   of the standard set. */
{
  int i;

  if(*attribute == '\0' || strcasecmp(attribute,DAVAR_VALUE) == 0){
    if(varp->type == DAVARINT){
      retval->valtype = DAVARINT_RPC;
      if(index ==  DAVAR_NOINDEX) {
	retval->any_u.i.i_len = varp->size;
	retval->any_u.i.i_val = (int *) malloc(varp->size*sizeof(int));
	for(i=0;i<varp->size;i++) {
	  retval->any_u.i.i_val[i] = ((DAINT *)varp->varptr)[i];
	  /*	printf("%d %d\n",i,retval->any_u.i.i_val[i]);*/
	}
      } else {
	retval->any_u.i.i_len = 1;
	retval->any_u.i.i_val = (int *) malloc(sizeof(int));
	retval->any_u.i.i_val[0] = ((DAINT *)varp->varptr)[index];
      }
    } else if(varp->type == DAVARFLOAT){
      retval->valtype = DAVARFLOAT_RPC;
      if(index ==  DAVAR_NOINDEX) {
	retval->any_u.r.r_len = varp->size;
	retval->any_u.r.r_val = (float *) malloc(varp->size*sizeof(float));
	for(i=0;i<varp->size;i++) {
	  retval->any_u.r.r_val[i] = ((DAFLOAT *)varp->varptr)[i];
	  /*	printf("%d %f\n",i,retval->any_u.r.r_val[i]);*/
	}
      } else {
	retval->any_u.r.r_len = 1;
	retval->any_u.r.r_val = (float *) malloc(sizeof(float));
	retval->any_u.r.r_val[0] = ((DAFLOAT *)varp->varptr)[index];
      }
    } else if(varp->type == DAVARDOUBLE){
      /* Return a float type for now since doubles don't work with our RPC */
      retval->valtype = DAVARFLOAT_RPC;
      if(index ==  DAVAR_NOINDEX) {
	retval->any_u.r.r_len = varp->size;
	retval->any_u.r.r_val = (float *) malloc(varp->size*sizeof(float));
	for(i=0;i<varp->size;i++) {
	  retval->any_u.r.r_val[i] = ((DADOUBLE *)varp->varptr)[i];
/*	  printf("%d %f\n",i,retval->any_u.r.r_val[i]);*/
	}
      } else {
	retval->any_u.r.r_len = 1;
	retval->any_u.r.r_val = (float *) malloc(sizeof(float));
	retval->any_u.r.r_val[0] = ((DADOUBLE *)varp->varptr)[index];
      }
    } else if(varp->type == DAVARSTRING && index == DAVAR_NOINDEX){
      /* indices to strings not supported */
      retval->valtype = DAVARSTRING_RPC;
      retval->any_u.s = (char *) malloc(strlen((char *)varp->varptr) + 1);
      strcpy(retval->any_u.s,(char *)varp->varptr);
      /*      printf("%s\n",retval->any_u.s);*/
    } else if(varp->type == DAVARFSTRING && index == DAVAR_NOINDEX){
      retval->valtype = DAVARSTRING_RPC;
      retval->any_u.s = (char *) malloc(varp->size+1);
      strncpy(retval->any_u.s,(char *)varp->varptr,varp->size);
      retval->any_u.s[varp->size] = '\0';
      kill_trailing(retval->any_u.s,' ');
    } else {
      retval->valtype = DAVARERROR_RPC;
      retval->any_u.error = S_SUCCESS;
    }
    return(S_SUCCESS);
  } else if(strcasecmp(attribute,DAVAR_TITLE) == 0){
    retval->valtype = DAVARSTRING_RPC;
    retval->any_u.s = (char *) malloc(strlen(varp->title) + 1);
    strcpy(retval->any_u.s,varp->title);
    /*    printf("%s\n",retval->any_u.s);*/
    return(S_SUCCESS);
  } else if(strcasecmp(attribute,DAVAR_SIZE) == 0){
    retval->valtype = DAVARINT_RPC;
    retval->any_u.i.i_len = 1;
    retval->any_u.i.i_val = (int *) malloc(sizeof(int));
    retval->any_u.i.i_val[0] = varp->size;
    return(S_SUCCESS);
  } else if(strcasecmp(attribute,DAVAR_TYPE) == 0){
    retval->valtype = DAVARINT_RPC;
    retval->any_u.i.i_len = 1;
    retval->any_u.i.i_val = (int *) malloc(sizeof(int));
    retval->any_u.i.i_val[0] = varp->type;
    return(S_SUCCESS);
  } else if(strcasecmp(attribute,DAVAR_FLAG) == 0){
    retval->valtype = DAVARINT_RPC;
    retval->any_u.i.i_len = 1;
    retval->any_u.i.i_val = (int *) malloc(sizeof(int));
    retval->any_u.i.i_val[0] = varp->flag;
    return(S_SUCCESS);
  } else if(strcasecmp(attribute,DAVAR_WATR) == 0){
    retval->valtype = DAVARSTRING_RPC;
    retval->any_u.s = daVarMakeWAtrList();
    return(S_SUCCESS);
  } else if(strcasecmp(attribute,DAVAR_RATR) == 0){
    retval->valtype = DAVARSTRING_RPC;
    retval->any_u.s = daVarMakeRAtrList();
    return(S_SUCCESS);
  } else {
    retval->valtype = DAVARERROR_RPC;
    retval->any_u.error = S_DAVAR_UNKATTR; /* Why success??? */
    return(S_DAVAR_UNKATTR);
  }
}  

char *daVarMakeRAtrList()
     /* Caller must free up the allocated memory */
{
  static int listsize=0;
  char *nllist;

  if(listsize==0) listsize=daVarGetAListSize(daVarRAtrList);
  nllist = (char *) malloc(listsize+1);
  daVarCopyAlist(nllist,daVarRAtrList);	/* Copy to \n delineated list */
  return(nllist);
}
char *daVarMakeWAtrList()
     /* Caller must free up the allocated memory */
{
  static int listsize=0;
  char *nllist;

  if(listsize==0) listsize=daVarGetAListSize(daVarWAtrList);
  nllist = (char *) malloc(listsize+1);
  daVarCopyAlist(nllist,daVarWAtrList);
  return(nllist);
}
int daVarGetAListSize(char **list)
{
  int len;
  len = 0;
  while(*list) len += strlen(*list++)+1;
  return(len);
}
void daVarCopyAlist(char *nllist,char **list)
{
  char *nlp;
  char *lp;

  nlp = nllist;
  while(lp=*list++){
    while(*nlp = *lp++) nlp++;
    *nlp++ = '\n';
  }
  *nlp = '\0';
}

daVarStatus daVarRegWatr(daVarStruct *varp, char *attribute
			 , int index, any *setval)
/* Regular Write Attribute routine.  Returns failure if attribute is not
   of the standard set. */
{
  int i;


  if(*attribute == '\0' || strcasecmp(attribute,DAVAR_VALUE) == 0){
    if(index == DAVAR_NOINDEX) index = 0;
    if(setval->valtype == DAVARINT_RPC) {
      if(varp->size >= setval->any_u.i.i_len+index){
	if(varp->type == DAVARINT){
	  for(i=0;i<setval->any_u.i.i_len;i++){
	    ((DAINT *)varp->varptr)[i+index] = setval->any_u.i.i_val[i];
	  }
	} else if(varp->type == DAVARFLOAT){
	  for(i=0;i<setval->any_u.i.i_len;i++){
	    ((DAFLOAT *)varp->varptr)[i+index] = setval->any_u.i.i_val[i];
	  }
	} else if(varp->type == DAVARDOUBLE){
	  for(i=0;i<setval->any_u.i.i_len;i++){
	    ((DADOUBLE *)varp->varptr)[i+index] = setval->any_u.i.i_val[i];
	  }
	} else
	  return(S_DAVAR_ILLCONV);
      } else
	return(S_DAVAR_TOOMANY);
      } else if(setval->valtype == DAVARFLOAT_RPC) {
	if(varp->size >= setval->any_u.r.r_len+index){
	  if(varp->type == DAVARINT){
	    for(i=0;i<setval->any_u.r.r_len;i++){
	      ((DAINT *)varp->varptr)[i+index] = floatToLong(setval->any_u.r.r_val[i]);
	    }
	  } else if(varp->type == DAVARFLOAT){
	    for(i=0;i<setval->any_u.r.r_len;i++){
	      ((DAFLOAT *)varp->varptr)[i+index] = setval->any_u.r.r_val[i];
	    }
	  } else if(varp->type == DAVARDOUBLE){
	    for(i=0;i<setval->any_u.r.r_len;i++){
	      ((DADOUBLE *)varp->varptr)[i+index] = setval->any_u.r.r_val[i];
	    }
	  } else
	    return(S_DAVAR_TOOMANY);
	  
	} else
	return(S_DAVAR_ILLCONV);
      }
#ifdef DOUBLERPC
      else if(setval->valtype == DAVARDOUBLE_RPC) {
	if(varp->size >= setval->any_u.d.d_len+index){
	  if(varp->type == DAVARINT){
	    for(i=0;i<setval->any_u.d.d_len;i++){
	      ((DAINT *)varp->varptr)[i+index] = floatToLong(setval->any_u.d.d_val[i]);
	    }
	  } else if(varp->type == DAVARFLOAT){
	    for(i=0;i<setval->any_u.d.d_len;i++){
	      ((DAFLOAT *)varp->varptr)[i+index] = setval->any_u.d.d_val[i];
	    }
	  } else if(varp->type == DAVARDOUBLE){
	    for(i=0;i<setval->any_u.d.d_len;i++){
	      ((DADOUBLE *)varp->varptr)[i+index] = setval->any_u.d.d_val[i];
	    }
	  } else
	    return(S_DAVAR_TOOMANY);
	} else
	  return(S_DAVAR_ILLCONV);
      }
#endif	  
      else if(setval->valtype == DAVARSTRING_RPC) {
	if(varp->type == DAVARSTRING) {
	  if(varp->size < strlen(setval->any_u.s)+index) {
	    return(S_DAVAR_TOOMANY); /* Maybe we should truncate strings??? */
	  }
	  strcpy(((char *) varp->varptr)+index, setval->any_u.s);
	} else if(varp->type == DAVARFSTRING) {
	  int in_len;
	  int blanks;
	  char *p;
	  in_len = strlen(setval->any_u.s);
/*	  printf("|%s|\n",setval->any_u.s);*/
	  if(varp->size < in_len+index) {
	    return(S_DAVAR_TOOMANY);
	  }
	  strncpy(((char *) varp->varptr)+index,
		  setval->any_u.s,in_len);
	  blanks = varp->size - index - in_len;
	  p = ((char *) varp->varptr) + index + in_len;
	  while(blanks-- > 0) *p++ = ' ';	/* Blank pad */
	} else
	  return(S_DAVAR_ILLCONV);
      } else
	return(S_FAILURE);
  } else if(strcasecmp(attribute,DAVAR_TITLE) == 0){
    /* Index ignored for title attribute */
/*    printf("setval->valtype = %d\n",setval->valtype);*/
    if(setval->valtype == DAVARSTRING_RPC) {
      varp->title = (char *) realloc(varp->title,strlen(setval->any_u.s) + 1);
/*      printf("Writing string %s\n",setval->any_u.s);*/
      strcpy(varp->title,setval->any_u.s);
    } else
      return(S_DAVAR_ILLCONV);
  } else {
    return(S_DAVAR_UNKATTR);
  }
  return(S_SUCCESS);
}  

daVarStatus daVarClassFind(char *name, daVarStruct **varp)
/* Pass ptr to name of full variable.
   Name is of form a.b.c.d.e..., searches for variable "a", then "a.b"
   and so on until a registered variable is found. 
   varp is the   pointer to the structure of the registered variable.

   Must have write access to the string name (even though it will be left
   unchanged.)

*/
{
  char c, *s, *t, *end;
  
  s = name;

  end = s + strlen(s);
  if((t=strchr(s,'('))) if(t<end) end = t; /* ) (to keep c-mode happy) */
  if((t=strchr(s,'['))) if(t<end) end = t; /* ] (to keep c-mode happy) */
  c = *end;
  *end = '\0';

  while((t=strchr(s,'.'))){
    *t = 0;
    if(daVarLookupP(name, varp) == S_SUCCESS) {
      *t = '.';
      *end = c;
      return(S_SUCCESS);
    }
    *t = '.';
    s = t + 1;
  }
  if(daVarLookupP(name, varp) == S_SUCCESS) {
    *end = c;
    return(S_SUCCESS);
  }
  *end = c;
  return(S_FAILURE);
}
daVarStatus daVarAttributeFind(char *name, daVarStruct *varclass, 
			daVarStruct **varp, char **attribute, int *index)
     /* Return DAVAR_NOINDEX for index if none specified.   */
{
  char *atptr;			/* Pointer to attribute string */
  daVarStatus status;
  char *end, *pptr, parenchar;

  end = name + strlen(name);
  if((status=thGetIndex(name, index, &pptr)) == S_SUCCESS) {
    parenchar = *pptr;
  } else if(status == S_DAVAR_NOINDEX) {
    parenchar = *pptr;
    *index = DAVAR_NOINDEX;
  } else			/* Failed */
    return(S_FAILURE);
  
  *pptr = '\0';
  atptr = pptr;
  

  if(strcasecmp(name,varclass->name) == 0){ /* Lookup alread done */
    *varp = varclass;
    *attribute = end;		/* Point to null */
  } else if(daVarLookupP(name, varp) == S_SUCCESS){ /* Fullly qualified name?*/
    *attribute = end;		/* Point to null */
  } else {			/* Start striping attributes off pptr */
    int clnamlen;
    clnamlen = strlen(varclass->name);
    atptr--;
    *varp = varclass;
    while((atptr - name) > clnamlen) { /* Stop before "class" name */
      if(*atptr == '.'){
	*atptr = '\0';
	if(daVarLookupP(name, varp) == S_SUCCESS){
	  *atptr = '.';
	  break;
	  *attribute = atptr + 1;
	} else
	  *atptr-- = '.';
      } else
	atptr--;
    }
    *attribute = atptr + 1;
  }
/*  *pptr = parenchar;  Leave name null terminated at the left paren */
  return(S_SUCCESS);
}  
