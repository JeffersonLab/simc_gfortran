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
 *  Special test and histogram handlers for RPC services.
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: thHandlers.c,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.3.22.1  2008/09/25 00:54:05  jones
 *   Updated for running on Fedora 8 with gfortran
 *
 *   Revision 1.4  2008/09/25 00:01:29  jones
 *   Updated to run with gfortran compiler
 *
 *   Revision 1.3.24.1  2007/09/10 21:32:47  pcarter
 *   Implemented changes to allow compilation on RHEL 3,4,5 and MacOSX
 *
 *   Revision 1.3  2003/02/21 20:55:24  saw
 *   Clean up some types and casts to reduce compiler warnings.
 *
 *   Revision 1.2  1999/11/04 20:34:06  saw
 *   Alpha compatibility.
 *   New RPC call needed for root event display.
 *   Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *   Revision 1.7  1999/08/25 13:16:06  saw
 *   *** empty log message ***
 *
 *   Revision 1.6  1999/07/07 13:43:58  saw
 *   Move thTestRHandler() into thTestParse.c
 *
 *   Revision 1.5  1996/08/01 01:31:56  saw
 *   Have thRHandler return a status
 *
 *	  Revision 1.4  1995/01/09  15:41:11  saw
 *	  Change "linux" ifdef's to NOHBOOK.
 *
 *	  Revision 1.3  1994/10/16  21:42:21  saw
 *	  Change an include file name from daVarServ.h daVarHandlers.h
 *
 *	  Revision 1.2  1993/08/12  14:57:39  saw
 *	  Add #include <rpc/rpc.h>
 *
 *	  Revision 1.1  1993/05/10  21:06:46  saw
 *	  Initial revision
 *
 */
#include <string.h>
#include <rpc/rpc.h>

#include "daVar.h"
#include "daVarRpc.h"
#include "daVarHandlers.h"
#include "th.h"
#include "thInternal.h"

#ifndef NOHBOOK
#include "hbook.h"
#endif

int thLastIdRhandled;

daVarStatus thWHandler(char *name,daVarStruct *varclass,any *setval)
/* The write handler used by block.test, block.hist and block.parm */
{
  daVarStruct *varp;
  char *attribute;
  daVarStatus status;
  int index;

  status = daVarAttributeFind(name, varclass, &varp, &attribute, &index);
  if(status == S_SUCCESS) {
    status = daVarRegWatr(varp, attribute, index, setval);
    if(strcasecmp(attribute,DAVAR_TITLE) == 0 && status == S_SUCCESS){
      status = ((daVarStatus (*)()) varclass->opaque)(varp);
    }
  }
  return(status);
}
daVarStatus thRHandler(char *name, daVarStruct *varclass, any *retval)
/* The default Read handler */
{
  daVarStruct *varp;
  char *attribute;
  daVarStatus status;
  int index;

  status = daVarAttributeFind(name, varclass, &varp, &attribute, &index);
  status = daVarRegRatr(varp, attribute, index, retval);
  /* scaler attribute a synonym for the value which holds the block counter */
  if(status == S_SUCCESS) {
    if(strcasecmp(attribute,DAVAR_RATR) == 0){
      retval->any_u.s = realloc(retval->any_u.s,strlen(retval->any_u.s)
				+strlen(TH_SCALER) + 2);
      strcat(retval->any_u.s,TH_SCALER);
      strcat(retval->any_u.s,"\n");
    }
  } else {
    if(strcasecmp(attribute,TH_SCALER) == 0){
      retval->valtype = DAVARINT_RPC;
      retval->any_u.i.i_len = 1;
      retval->any_u.i.i_val = (int *) malloc(sizeof(int));
      retval->any_u.i.i_val[0] = ((DAINT *)varp->varptr)[0];
    }
  }
  return(status);
}

#ifndef NOHBOOK
void thHistZeroLastId()
{
  thLastIdRhandled = 9999999;
  return;
}
daVarStatus thHistRHandler(char *name, daVarStruct *varclass, any *retval)
     /* Read Handler for Histograms */
{
  daVarStruct *varp;
  char *attribute;
  daVarStatus status;
  int index;
  static int NX,NY,NWT,LOC ; static float XMI,XMA,YMI,YMA;
/*  thHistOpaque *hopq;*/

  status = daVarAttributeFind(name, varclass, &varp, &attribute, &index);
  status = daVarRegRatr(varp, attribute, index, retval);
  if(status == S_SUCCESS) {
    if(strcasecmp(attribute,DAVAR_RATR) == 0){
      retval->any_u.s = realloc(retval->any_u.s,strlen(retval->any_u.s)
				+ strlen(TH_ND) + strlen(TH_NX) 
				+ strlen(TH_NY) + strlen(TH_XMI)
				+ strlen(TH_XMA) + strlen(TH_YMI)
				+ strlen(TH_YMA) + strlen(TH_CONTEN) + 9);
      strcat(retval->any_u.s,TH_ND); strcat(retval->any_u.s,"\n");
      strcat(retval->any_u.s,TH_NX); strcat(retval->any_u.s,"\n");
      strcat(retval->any_u.s,TH_NY); strcat(retval->any_u.s,"\n");
      strcat(retval->any_u.s,TH_XMI); strcat(retval->any_u.s,"\n");
      strcat(retval->any_u.s,TH_XMA); strcat(retval->any_u.s,"\n");
      strcat(retval->any_u.s,TH_YMI); strcat(retval->any_u.s,"\n");
      strcat(retval->any_u.s,TH_YMA); strcat(retval->any_u.s,"\n");
      strcat(retval->any_u.s,TH_CONTEN); strcat(retval->any_u.s,"\n");
    }
  } else {
    char chtitle[80];
    
    retval->valtype = DAVARERROR_RPC;
    retval->any_u.error = S_DAVAR_UNKATTR;
    if(thLastIdRhandled != *((DAINT *) varp->varptr)) {
      thLastIdRhandled = *((DAINT *) varp->varptr);
      /*HGIVE(*((DAINT *) varp->varptr),chtitle,NX,XMI,XMA,NY,YMI,YMA
	,NWT,LOC);*/
    }
    if(strcasecmp(attribute,TH_ND) == 0){
      retval->valtype = DAVARINT_RPC;
      retval->any_u.i.i_len = 1;
      retval->any_u.i.i_val = (int *) malloc(varp->size*sizeof(int));
      retval->any_u.i.i_val[0] = (NY == 0 ? 1 : 2);
    } else if(strcasecmp(attribute,TH_NX) == 0){
      retval->valtype = DAVARINT_RPC;
      retval->any_u.i.i_len = 1;
      retval->any_u.i.i_val = (int *) malloc(varp->size*sizeof(int));
      retval->any_u.i.i_val[0] = NX;
    } else if(strcasecmp(attribute,TH_NY) == 0){
      retval->valtype = DAVARINT_RPC;
      retval->any_u.i.i_len = 1;
      retval->any_u.i.i_val = (int *) malloc(varp->size*sizeof(int));
      retval->any_u.i.i_val[0] = NY;
    } else if(strcasecmp(attribute,TH_XMI) == 0){
      retval->valtype = DAVARFLOAT_RPC;
      retval->any_u.r.r_len = 1;
      retval->any_u.r.r_val = (float *) malloc(varp->size*sizeof(float));
      retval->any_u.r.r_val[0] = XMI;
    } else if(strcasecmp(attribute,TH_XMA) == 0){
      retval->valtype = DAVARFLOAT_RPC;
      retval->any_u.r.r_len = 1;
      retval->any_u.r.r_val = (float *) malloc(varp->size*sizeof(float));
      retval->any_u.r.r_val[0] = XMA;
    } else if(strcasecmp(attribute,TH_YMI) == 0){
      retval->valtype = DAVARFLOAT_RPC;
      retval->any_u.r.r_len = 1;
      retval->any_u.r.r_val = (float *) malloc(varp->size*sizeof(float));
      retval->any_u.r.r_val[0] = YMI;
    } else if(strcasecmp(attribute,TH_YMA) == 0){
      retval->valtype = DAVARFLOAT_RPC;
      retval->any_u.r.r_len = 1;
      retval->any_u.r.r_val = (float *) malloc(varp->size*sizeof(float));
      retval->any_u.r.r_val[0] = YMA;
    } else if(strcasecmp(attribute,TH_CONTEN) == 0){
      int size;
      retval->valtype = DAVARFLOAT_RPC;
      size = NX;
      if(NY != 0) size *= NY;
      retval->any_u.r.r_len = size;
      retval->any_u.r.r_val = (float *)malloc(size*sizeof(float));
      /* Next line gives warning "assignment of read-only location */
      char tmpstring[] = "HIST";
      /*HUNPAK(thLastIdRhandled,retval->any_u.r.r_val,tmpstring,(int) 1);*/
    }
  }
  return(status);
}
#endif
#ifdef ROOTTREE
daVarStatus thTreeRHandler(char *name, daVarStruct *varclass, any *retval)
/* The default Read handler */
{
  return(thRHandler(name, varclass, retval));
}
#endif
