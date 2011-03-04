#if defined(__osf__) || defined(__LP64__)
#define BIT64
#endif
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
 *  Routines used by a client to retrieve CTP variables
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *  $Log: thClient.c,v $
 *  Revision 1.2  2011/03/04 20:01:51  jones
 *  Add check for 64bit by looking for LP64
 *
 *  Revision 1.5  2003/02/21 20:55:24  saw
 *  Clean up some types and casts to reduce compiler warnings.
 *
 *  Revision 1.4  1999/11/04 20:34:05  saw
 *  Alpha compatibility.
 *  New RPC call needed for root event display.
 *  Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *  Revision 1.6  1999/08/25 13:16:06  saw
 *  *** empty log message ***
 *
 *  Revision 1.5  1999/03/01 19:53:05  saw
 *  Add OSF stuff
 *
 * Revision 1.4  1995/08/03  13:50:52  saw
 * Add SGI compatibility
 *
 * Revision 1.3  1994/11/07  14:09:42  saw
 * Bug fixes in thGetList_test and callback server code.
 *
 * Revision 1.2  1994/10/17  17:07:28  saw
 * Add thGetList_test and the callback service davar_readmultiple_test_cb_1
 *
 * Revision 1.1  1994/09/27  19:19:09  saw
 * Initial revision
 *
 */
#include <stdio.h>
#include <string.h>
#include <rpc/rpc.h>
#include "daVar.h"
#include "daVarRpc.h"
#include "cfortran.h"

int thCreateList();        /* Move to some  include file */
int thAddToList(int handle, char *pattern);
int thRemoveFromList(int handle, char *pattern);
#ifdef BIT64
int thGetList(int handle, int client);
int thGetList_test(int handle, int client, char *test_condition,
		    int max_time_wait, int max_event_wait);
int thImportVars(char *pattern, int client);
#else
int thGetList(int handle, CLIENT *clnt);
int thGetList_test(int handle, CLIENT *clnt, char *test_condition,
		    int max_time_wait, int max_event_wait);
int thImportVars(char *pattern, CLIENT *clnt);
#endif
int thPrintList(int handle);
#ifdef BIT64
int myClntCreate(char *host, int prog, int vers, char *proto);
#endif

FCALLSCFUN0(INT,thCreateList,THCRLIST,thcrlist);
FCALLSCFUN2(INT,thAddToList,THADDLIST,thaddlist,INT,STRING);
FCALLSCFUN2(INT,thRemoveFromList,THREMLIST,thremlist,INT,STRING);
FCALLSCFUN2(INT,thGetList,THGETLIST,thgetlist,INT,INT);
FCALLSCFUN5(INT,thGetList_test,THCGETLIST,thcgetlist,INT,INT,STRING,INT,INT);
FCALLSCFUN1(INT,thPrintList,THPRTLIST,thprtlist,INT);
/* Don't really understand the following.  What about ultrix?
   This is probably because of the _ in clnt_create */
#ifndef __osf__
FCALLSCFUN4(INT,clnt_create,CLNT_CREATE,clnt_create,STRING,INT,INT,STRING);
#else
#ifdef BIT64
FCALLSCFUN4(INT,myClntCreate,CLNT_CREATE,clnt_create,STRING,INT,INT,STRING);
#else
FCALLSCFUN4(INT,clnt_create,CLNT_CREATE_,clnt_create_,STRING,INT,INT,STRING);
#endif
#endif

struct thNameNode {
  char *name;
  struct thNameNode *next;
};
typedef struct thNameNode thNameNode;

struct thNameList {
  thNameNode *namehead;
  int nnames;			/* Number of names in list */
  int rpc_made;
  NAMELIST rpc;		/* List in form for RPC argument */
};
typedef struct thNameList thNameList;

TESTNAMELIST *pending_arg=0;
int pending_flag=0;
int callback_result;

#ifdef BIT64
#define MAXHANDLES 10
static thNameList *handle_list[MAXHANDLES]={0,0,0,0,0,0,0,0,0,0};
static CLIENT *clnt_list[MAXHANDLES]={0,0,0,0,0,0,0,0,0,0};
#endif

#ifdef BIT64
/* Keep client pointers in an array so that we can return a 32 bit "handle"
   instead of 64 bit.  We probably won't need more than one, but we have an
   array to hold 10 client pointers just for the heck of it.
*/
int myClntCreate(char *host, int prog, int vers, char *proto)
{
  CLIENT *clnt;
  int client;

  clnt = clnt_create(host, prog, vers, proto);
  for(client=0;client<MAXHANDLES;client++){
    if(clnt_list[client]==0) {
      clnt_list[client] = clnt;
      return(client+1);
    }
  }
  return(0);
}
#endif
int thCreateList(){
  /* Create a handle for a list of variables */
  thNameList *list;
#ifdef BIT64
  int ihandle;
#endif

  list = (thNameList *) malloc(sizeof(thNameList));;
  list->namehead = 0;
  list->nnames = 0;
  list->rpc_made = 0;
  list->rpc.NAMELIST_len = 0;
  list->rpc.NAMELIST_val = 0;
#ifdef BIT64
  for(ihandle=0;ihandle<MAXHANDLES;ihandle++){
    if(handle_list[ihandle]==0) {
      handle_list[ihandle] = list;
      printf("cr: handle_list[%d]=%x\n",ihandle,handle_list[ihandle]);
      return(ihandle+1);
    }
  }
  free(list);
  return(0);
#else
  return((int) list);
#endif
}
int thAddToList(int handle, char *pattern)
/* Add registered variables to a list of variables to get from a server
   Return the number of variables added
   Should we check for duplicates?  Not now.  No harm in duplicates.
   thRemoveFromList will remove all duplicates though. */
{
  thNameList *list;
  thNameNode *next,*this;
  char **vlist;			/* List of characters matching pattern */
  int count;			/* Number of variables being added */
  int i;

#ifdef BIT64
  list = (thNameList *) handle_list[handle-1];
#else
  list = (thNameList *) handle;
#endif
  daVarList(pattern,&vlist,&count);
  for(i=0;i<count;i++) {
    next = list->namehead;		/* The current list */
    this = list->namehead = (thNameNode *) malloc(sizeof(thNameNode)); /* New name */
    this->name = (char *) malloc(strlen(vlist[i])+1);
    strcpy(this->name,vlist[i]);
    this->next = next;		/* Attach rest of list */
  }
  list->nnames += count;	/* Perhaps I should just count when needed ? */
  list->rpc_made = 0;		/* RPC format list now out of date. */
  return(list->nnames);
}
int thRemoveFromList(int handle, char *pattern)
{
  thNameList *list;
  thNameNode *this,**thisp;
  char **vlist;			/* List of characters matching pattern */
  int count;			/* Number of variables being removed */
  int nremove;
  int i;

#ifdef BIT64
  list = (thNameList *) handle_list[handle-1];
#else
  list = (thNameList *) handle;
#endif
  daVarList(pattern,&vlist,&count);
  nremove = 0;
  for(i=0;i<count;i++) {
    this = list->namehead;	/* Start of list */
    thisp = &list->namehead;	/* Pointer to next field of previous name */
    while(this) {
      if(strcasecmp(this->name,vlist[i])==0) { /* Remove matching string */
	*thisp = this->next;
	free(this->name);
	free(this);
	this = *thisp;
	nremove++;
      } else {
        thisp = &(this->next);
	this = this->next;
      }
    }
  }
  list->nnames -= nremove;	/* Perhaps I should just count when needed ? */
  return(list->nnames);
}
int thPrintList(int handle)
     /* For debugging */
{
  thNameList *list;
  thNameNode *next,*this;
  int count;

#ifdef BIT64
  list = (thNameList *) handle_list[handle-1];
#else
  list = (thNameList *) handle;
#endif
  this = list->namehead;
  printf("Variables attached to handle %d\n",handle);
  count++;
  while(this) {
    printf("%s\n",this->name);
    count++;
    this = this->next;
  }
  return(count);
}
#ifdef BIT64
int thGetList(int handle, int client)
#else
int thGetList(int handle, CLIENT *clnt)
#endif
     /* Returns 0 for total success, -1 for total failure, a positive number
	for the number of variables that didn't work */
{
  thNameList *list;
  thNameNode *next,*this;
  int i;
  RVALLIST *vals;
  int nerrors;
#ifdef BIT64
  CLIENT *clnt;

  clnt = clnt_list[client-1];
#endif

/*
  printf("sizeof(handle)=%d\n",sizeof(handle));
  printf("thGetLIst: handle=%d\n",handle);
  printf("thGhandle_list[0]=%x\n",handle_list[0]);
  printf("handle_list[%d-1]=%x\n",handle,handle_list[handle-1]);
*/
#ifdef BIT64
  list = (thNameList *) handle_list[handle-1];
#else
  list = (thNameList *) handle;
#endif
/*  printf("list->nnames=%d\n",list->nnames);*/
  if(!list->rpc_made) {
    if(list->rpc.NAMELIST_len == 0) {
      list->rpc.NAMELIST_len = list->nnames;
      list->rpc.NAMELIST_val = (char **) malloc(list->rpc.NAMELIST_len
						*sizeof(char *));
    } else if (list->rpc.NAMELIST_len != list->nnames) {
      list->rpc.NAMELIST_len = list->nnames;
      list->rpc.NAMELIST_val = (char **)
	realloc(list->rpc.NAMELIST_val,list->rpc.NAMELIST_len*sizeof(char *));
    }
    this = list->namehead;
    for(i=0;(i<list->nnames && this);i++){
      list->rpc.NAMELIST_val[i] = this->name;
      this = this->next;
    }
  }
  nerrors = 0;
  if(vals = davar_readmultiple_1(&(list->rpc),clnt)) {
    this = list->namehead;
/*    printf("list->rpc.NAMELIST_len=%d\n",list->rpc.NAMELIST_len);*/
    for(i=0;(i<list->rpc.NAMELIST_len && this);i++){
/*      printf("%s\n",this->name);*/
      if(vals->RVALLIST_val[i].valtype != DAVARERROR_RPC){
	if(daVarWriteVar(this->name,&(vals->RVALLIST_val[i])) != S_SUCCESS)
	  nerrors++;
      } else {
	nerrors++;
      }
      this = this->next;
    }
  } else {
    nerrors = -1;
  }
  return(nerrors);
}
#if 0
#ifdef BIT64
int thPutList(int handle, int client)
#else
int thPutList(int handle, CLIENT *clnt)
#endif
     /* Returns 0 for total success, -1 for total failure, a positive number
	for the number of variables that didn't work */
{
  thNameList *list;
  thNameNode *next,*this;
  int i;
  WVALLIST vals;
  int nerrors;
#ifdef BIT64
  CLIENT *clnt;

  clnt = clnt_list[client-1];
#endif

  /* Create the write structure */
#ifdef BIT64
  list = (thNameList *) handle_list[handle-1];
#else
  list = (thNameList *) handle;
#endif
/*  printf("list->nnames=%d\n",list->nnames);*/
  if(!list->rpc_made) {
    if(list->rpc.NAMELIST_len == 0) {
      list->rpc.NAMELIST_len = list->nnames;
      list->rpc.NAMELIST_val = (char **) malloc(list->rpc.NAMELIST_len
						*sizeof(char *));
    } else if (list->rpc.NAMELIST_len != list->nnames) {
      list->rpc.NAMELIST_len = list->nnames;
      list->rpc.NAMELIST_val = (char **)
	realloc(list->rpc.NAMELIST_val,list->rpc.NAMELIST_len*sizeof(char *));
    }
    this = list->namehead;
    for(i=0;(i<list->nnames && this);i++){
      list->rpc.NAMELIST_val[i] = this->name;
      this = this->next;
    }
  }
  nerrors = 0;
  if(vals = davar_readmultiple_1(&(list->rpc),clnt)) {
    this = list->namehead;
/*    printf("list->rpc.NAMELIST_len=%d\n",list->rpc.NAMELIST_len);*/
    for(i=0;(i<list->rpc.NAMELIST_len && this);i++){
/*      printf("%s\n",this->name);*/
      if(vals->RVALLIST_val[i].valtype != DAVARERROR_RPC){
	if(daVarWriteVar(this->name,&(vals->RVALLIST_val[i])) != S_SUCCESS)
	  nerrors++;
      } else {
	nerrors++;
      }
      this = this->next;
    }
  } else {
    nerrors = -1;
  }
  return(nerrors);
}
#endif
#if 0
int tsize=0;
struct timeval timeout;
int servone(int wait)
/* Need to move something that does this into CTP proper */
{
  fd_set readfdset;
  extern int errno;
 
  timeout.tv_sec = wait;
  timeout.tv_usec = 1;

#ifdef hpux
  if(!tsize) tsize = NFDBITS;
#else
  if(!tsize) tsize = getdtablesize(); /* how many descriptors can we have */
#endif

  readfdset = svc_fdset;
  switch(select(tsize, &readfdset, (fd_set *) NULL, (fd_set *) NULL,
		&timeout)) {
  case -1:
    if (errno == EBADF) break;
    perror("select failed");
    break;
  case 0:
    /* perform other functions here if select() timed-out */
    break;
  default:
    svc_getreqset(&readfdset);
  }
}
#endif
#ifdef BIT64
int thGetList_test(int handle, int client, char *test_condition,
		    int max_time_wait, int max_event_wait)
#else
int thGetList_test(int handle, CLIENT *clnt, char *test_condition,
		    int max_time_wait, int max_event_wait)
#endif
     /* Returns 0 for total success, -1 for total failure, a positive number
	for the number of variables that didn't work */
{
  thNameList *list;
  thNameNode *next,*this;
  int i;
  int *status;
  TESTNAMELIST *arg;
  int servret;
#ifdef BIT64
  CLIENT *clnt;

  clnt = clnt_list[client-1];
#endif

  /* Can return some kind of error if pending_arg is not zero */
#ifdef BIT64
  list = (thNameList *) handle_list[handle-1];
#else
  list = (thNameList *) handle;
#endif
/*  printf("list->nnames=%d\n",list->nnames);*/
  if(!list->rpc_made) {
    if(list->rpc.NAMELIST_len == 0) {
      list->rpc.NAMELIST_len = list->nnames;
      list->rpc.NAMELIST_val = (char **) malloc(list->rpc.NAMELIST_len
						*sizeof(char *));
    } else if (list->rpc.NAMELIST_len != list->nnames) {
      list->rpc.NAMELIST_len = list->nnames;
      list->rpc.NAMELIST_val = (char **)
	realloc(list->rpc.NAMELIST_val,list->rpc.NAMELIST_len*sizeof(char *));
    }
    this = list->namehead;
    for(i=0;(i<list->nnames && this);i++){
      list->rpc.NAMELIST_val[i] = this->name;
      this = this->next;
    }
  }
  arg = (TESTNAMELIST *) malloc(sizeof(TESTNAMELIST));
  arg->test_condition = (char *) malloc(strlen(test_condition)+1);
  strcpy(arg->test_condition,test_condition);
  arg->max_time_wait = max_time_wait;
  arg->max_event_wait = max_event_wait;
  arg->prog = DAVARSVR;
  arg->vers = DAVARVERS+1;
  arg->NAMELISTP = &list->rpc;
  pending_arg = arg;
  pending_flag = 1;

  if(!(status = davar_readmultiple_test_1(arg,clnt)))
    return(-1);

  /* Now wait for the incoming network call */

  servret = 1;
  while(pending_flag && servret > 0) /* Wait for timeout, completion or failur*/
    servret = daVarServOnce(arg->max_time_wait+10); /* Will wait double?? */
  if(servret == 0) callback_result = -2;	/* Timeout */
  else if(servret == -1) callback_result = -3;

  free(arg->test_condition);
  free(arg);
  pending_arg = 0;

  return(callback_result);
}
#ifdef BIT64
int thImportVars(char *pattern, int client)
#else
int thImportVars(char *pattern, CLIENT *clnt)
#endif
     /* Returns 0 for total success, -1 for total failure, a positive number
	for the number of variables that didn't work */
{
  WVALLIST *vals;
  int count;

  thNameList *list;
  thNameNode *next,*this;
  int i;
  int nerrors;
#ifdef BIT64
  CLIENT *clnt;

  clnt = clnt_list[client-1];
#endif

  /* need to initialize the hash tables */

  if(!(vals = davar_readpatternmatch_1(&pattern,clnt))) {
    return(-1);			/* Failed */
  }

  count = vals->WVALLIST_len;
  /*printf("daVarImportVars got %d variables matching %s\n",count,pattern);*/
  nerrors = 0;
  for(i=0;i<count;i++) {
    char *name;
    int valtype;
    daVarStruct var;

    name = vals->WVALLIST_val[i].name;
    /*printf("%d: %s\n",i,name);*/
    /* Don't do anything if it already exists */
    if(daVarWriteVar(name,vals->WVALLIST_val[i].val) == S_SUCCESS) continue;
    var.type = vals->WVALLIST_val[i].val->valtype;
    if(var.type == DAVARERROR_RPC){
      printf("Error getting %s\n",name);
      nerrors++;
      continue;
    }
    var.name = name;
    /*printf("Vartype = %d\n",var.type);*/
    switch(var.type)
      {
      case DAVARINT_RPC:
	var.size = vals->WVALLIST_val[i].val->any_u.i.i_len;
	/*printf("size=%d\n",var.size);*/
	var.varptr = (void *) malloc(var.size*sizeof(DAINT));
	break;
      case DAVARFLOAT_RPC:
	var.size = vals->WVALLIST_val[i].val->any_u.r.r_len;
	/*printf("size=%d\n",var.size);*/
	var.varptr = (void *) malloc(var.size*sizeof(DAFLOAT));
	break;
      case DAVARDOUBLE_RPC:
	var.size = vals->WVALLIST_val[i].val->any_u.d.d_len;
	/*printf("size=%d\n",var.size);*/
	var.varptr = (void *) malloc(var.size*sizeof(DADOUBLE));
	break;
      case DAVARSTRING_RPC:
	var.size = strlen(vals->WVALLIST_val[i].val->any_u.s) + 1;
	/*printf("size=%d\n",var.size);*/
	var.varptr = malloc(var.size);
	break;
      }
    var.opaque = 0;
    var.rhook = 0;
    var.whook = 0;
    var.flag = DAVAR_REPOINTOK | DAVAR_DYNAMIC_PAR;
    var.flag = DAVAR_REPOINTOK | DAVAR_DYNAMIC_PAR;
    var.title = 0;
    daVarRegister((int) 0, &var);
    /*    free(var.name);*/
    if(daVarWriteVar(name,vals->WVALLIST_val[i].val) != S_SUCCESS) {
      printf("daVarWriteVar of %s should have worked\n",name);
      nerrors++;
    }
  }
  return(nerrors);
}
int *davar_readmultiple_test_cb_1(RVALLIST *vals, CLIENT *clnt)
{
  static int result;

  TESTNAMELIST *argp;
  thNameNode *next,*this;
  int i;

  if(pending_arg) argp = pending_arg;
  else {
    pending_flag = 0;
    return(&result);		/* What error code ?? */
  }

  callback_result = 0;
  if(argp->NAMELISTP->NAMELIST_len == vals->RVALLIST_len) {
    for(i=0;(i<argp->NAMELISTP->NAMELIST_len);i++){
/*      printf("%s\n",this->name);*/
      if(vals->RVALLIST_val[i].valtype != DAVARERROR_RPC){
	if(daVarWriteVar(argp->NAMELISTP->NAMELIST_val[i]
			 ,&(vals->RVALLIST_val[i])) != S_SUCCESS)
	  callback_result++;
      } else {
	callback_result++;
      }
    }
  } else if (vals->RVALLIST_len>0) {
    printf("Lengths: %d %d",argp->NAMELISTP->NAMELIST_len,vals->RVALLIST_len);
    callback_result = -1;
  } else {
    callback_result = -2;       /* Server send timeout signal */
  }
  pending_flag = 0;
  return(&result);
}
