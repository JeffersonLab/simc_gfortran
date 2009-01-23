/*-----------------------------------------------------------------------------
 * Copyright (c) 1993-1996 Southeastern Universities Research Association,
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
 *  Structures, defs, and prototypes for internal use by CTP.
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: thInternal.h,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.3  2004/07/08 20:06:09  saw
 *   Always have CTP Tree routines (perhaps dummy) available
 *
 *   Revision 1.2  1999/11/04 20:34:06  saw
 *   Alpha compatibility.
 *   New RPC call needed for root event display.
 *   Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *   Revision 1.8  1999/08/25 13:16:06  saw
 *   *** empty log message ***
 *
 *   Revision 1.7  1999/03/01 19:54:48  saw
 *   Add weighted histograms
 *
 *	  Revision 1.6  1996/01/30  15:44:49  saw
 *	  Add prototypes of calls needed by groups, add DEFAULTGRPSTR and ALLGRPSTR
 *
 *	  Revision 1.5  1995/08/03  13:53:43  saw
 *	  Add parameters for #, real, integer, double
 *
 *	  Revision 1.4  1994/09/27  20:07:54  saw
 *	  Add SCALERSTR definition
 *
 *	  Revision 1.3  1994/07/21  20:48:27  saw
 *	  Add "gethit" and "report" string definitions
 *
 *	  Revision 1.2  1993/11/22  21:21:07  saw
 *	  Add QUOTECHARS definitions.
 *
 *	  Revision 1.1  1993/05/11  17:37:30  saw
 *	  Initial revision
 *
 */

#ifndef _TH_INTERNAL_H
#define _TH_INTERNAL_H

#ifndef _DAVAR_H
#include "daVar.h"
#endif

thStatus thLoadParameters(daVarStruct *var);
thStatus thBookGethits(daVarStruct *var);
thStatus thBookTests(daVarStruct *var);
thStatus thBookHists(daVarStruct *var);
thStatus thBookReports(daVarStruct *var);
thStatus thExecuteTestsV(daVarStruct *var);
thStatus thClearTestFlagsV(daVarStruct *var);
thStatus thClearTestScalersV(daVarStruct *var);
thStatus thIncTestScalersV(daVarStruct *var);
thStatus thExecuteHistsV(daVarStruct *var);
thStatus thClearHistsV(daVarStruct *var);
thStatus thExecuteaGethitBlock(daVarStruct *var);
thStatus thBookTree(daVarStruct *var);
thStatus thFillTreeV(daVarStruct *var);
thStatus thClearTreeV(daVarStruct *var);
thStatus thWriteTreeV(daVarStruct *var);
thStatus thCloseTreeV(daVarStruct *var);

extern daVarStatus thWHandler();
extern daVarStatus thRHandler();

#define BEGINSTR "begin"
#define ENDSTR "end"
#define BLOCKSTR "block"
#define GROUPSTR "group"
#define DEFAULTGRPSTR "default"
#define ALLGRPSTR "all"
#define COMCHAR ';'
#define SPECIALCHAR '#'

/* Double and single quote chars */
#define QUOTECHAR1 0x22
#define QUOTECHAR2 0x27
#define PARMSTR "parm"
#define TESTSTR "test"
#define HISTSTR "hist"
#define UHISTSTR "uhist"
#define NTUPLESTR "ntup"
#define EVENTSTR "event"
#define GETHITSTR "gethit"
#define REPORTSTR "report"
#define TREESTR "tree"
#define SCALERSTR "scaler"
#define REAL "real"
#define INTEGER "integer"
#define DOUBLE "double"

#define HTITLECHAR '$'
#define WEIGHTCHAR '#'

#define TH_SCALER "scaler"
#define TH_ND "nd"
#define TH_NX "nx"
#define TH_NY "ny"
#define TH_XMI "xmi"
#define TH_XMA "xma"
#define TH_YMI "ymi"
#define TH_YMA "yma"
#define TH_CONTEN "conten"
#define TH_X "x"
#define TH_Y "y"
#define TH_TEST "test"

#endif

