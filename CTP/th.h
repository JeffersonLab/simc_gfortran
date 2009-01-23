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
 *  Include file with prototypes for CTP routines.
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: th.h,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.2  1999/11/04 20:34:05  saw
 *   Alpha compatibility.
 *   New RPC call needed for root event display.
 *   Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *   Revision 1.6  1999/08/25 13:16:05  saw
 *   *** empty log message ***
 *
 *   Revision 1.5  1996/01/30 15:39:45  saw
 *   Add new prototypes: group calls, class calls used by groups
 *
 *	  Revision 1.4  1995/08/03  13:47:48  saw
 *	  Add Error code for integer overflow
 *
 *	  Revision 1.3  1994/11/07  14:40:22  saw
 *	  Add error code S_TH_UNREG
 *
 *	  Revision 1.2  1994/09/27  20:34:46  saw
 *	  Define STDERR to be stdout
 *
 *	  Revision 1.1  1993/05/11  15:10:08  saw
 *	  Initial revision
 *
 */

#ifndef _TH_H
#define _TH_H
/* Variable registration */

typedef int thStatus;		/* Return status type */

/* General booking routines */
thStatus thLoad(char *fname);
thStatus thOBook();
thStatus thBook();


/* Test package routines */
typedef enum {
  WALK_DISPLAY, WALK_CLEAR_FLAGS, WALK_CLEAR_SCALERS, WALK_EXECUTE,
  WALK_REMOVE, WALK_INCREMENT_SCALERS} WALKOP;

thStatus thWalkTree(char *block_name, WALKOP walkop);
thStatus thExecuteTests(char *block_name);
thStatus thClearTestFlags(char *block_name);
thStatus thClearTestScalers(char *block_name);
thStatus thIncTestScalers(char *block_name);
/*#define thExecuteTests(block_name) thWalkTree(block_name,WALK_EXECUTE)*/
/*#define thIncTestScalers(block_name) \
  thWalkTree(block_name,WALK_INCREMENT_SCALERS)*/
#define thDisplayTests(block_name) thWalkTree(block_name,WALK_DISPLAY)
#define thRemoveTests(block_name) thWalkTree(block_name,WALK_REMOVE)
/*#define thClearTestFlags(block_name) \
  thWalkTree(block_name,WALK_CLEAR_FLAGS)*/
/*#define thClearTestScalers(block_name) \
  thWalkTree(block_name,WALK_CLEAR_SCALERS)*/
thStatus thExecuteGroup(char *group_name);
thStatus thClearGroup(char *group_name);
thStatus thClearScalersGroup(char *group_name);
thStatus thIncrementScalersGroup(char *group_name);
thStatus thWriteGroup(char *group_name);
thStatus thCloseGroup(char *group_name);

/* Histogram package routines */
thStatus thExecuteHists(char *block_name);
thStatus thClearHists(char *block_name);
int thGetHistID(char *name);
thStatus thHistAliasWrite(char *fname);

#ifndef S_SUCCESS
#define S_SUCCESS 0
#define S_FAILURE -1
#endif
#define S_TH_UNREG 1            /* Unregistered variable in test expression */
#define S_INTOVF 2

#endif
#define STDERR stdout
