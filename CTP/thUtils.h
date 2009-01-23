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
 *  prototypes and defs for when thUtils.c routines are used.
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: thUtils.h,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.1.22.1  2008/09/25 00:54:05  jones
 *   Updated for running on Fedora 8 with gfortran
 *
 *   Revision 1.2  2008/09/25 00:01:28  jones
 *   Updated to run with gfortran compiler
 *
 *   Revision 1.1.24.1  2007/09/10 21:32:47  pcarter
 *   Implemented changes to allow compilation on RHEL 3,4,5 and MacOSX
 *
 *   Revision 1.1  1998/12/07 22:11:14  saw
 *   Initial setup
 *
 *	  Revision 1.2  1995/08/03  15:02:59  saw
 *	  Declare thSpecial for #real, #int, ... directives
 *
 *	  Revision 1.1  1993/05/11  17:40:51  saw
 *	  Initial revision
 *
 */

#ifndef _TH_UTILS_H
#define _TH_UTILS_H

#ifndef _DAVAR_H
#include "daVar.h"
#endif

#define min(a,b) (a>b ? b : a)

typedef enum {TOKINT, TOKFLOAT, TOKVAR, TOKARRAY} thTokenType;

daVarStatus thGetIndex(char *name, int *index, char **pptr);
int thComSpace(char *s,char **args);
int thCommas(char *s,char **args);
char *thTokenArray(char *s,int *index);
thTokenType thIDToken(char *s);
char *thSpaceStrip(char *s);
int thSpecial(char *line, char *default_class);
#ifdef OLD
*int thGetMode(char *s, enum MODES mode, enum MODES *nmode, char **blocknamep);
int thTokToPtr(char *token, int create, int intonly, daVarStruct *varp);
#endif
int thCleanLine(char *s);
/*float argtoFloat(daVarStruct *x);*/
int argtoInt(daVarStruct *x);
char *strcasestr(const char *s1, const char *s2);

struct daVarStructList {
  struct daVarStructList *next;
  daVarStruct *varp;
};

typedef struct daVarStructList daVarStructList;

void thAddVarToList(daVarStructList **head,daVarStruct *varp);
#endif
