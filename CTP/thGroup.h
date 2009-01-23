/*-----------------------------------------------------------------------------
 * Copyright (c) 1996 Southeastern Universities Research Association,
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
 *  Include file with group opaque structure and internal group calls
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: thGroup.h,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.2  1999/11/04 20:34:05  saw
 *   Alpha compatibility.
 *   New RPC call needed for root event display.
 *   Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *   Revision 1.2  1999/08/25 13:16:06  saw
 *   *** empty log message ***
 *
 *   Revision 1.1  1996/01/30 15:42:43  saw
 *   Initial revision
 *
 */
struct thGroupOpaque {		/* Opaque structure for group definitions */
  daVarStructList *blocklist;
  char *type;			/* i.e. hist, test, gethit, ... */
  thStatus (*book)();		/* Hooks to the appropriate routines */
  thStatus (*execute)();
  thStatus (*clear)();
  thStatus (*clearScalers)();
  thStatus (*incrementScalers)();
  thStatus (*ctpwrite)();
  thStatus (*ctpclose)();
};
typedef struct thGroupOpaque thGroupOpaque;

void thInitGroupOpaque(char *name, thGroupOpaque *opqptr);
/*void thInitGroupOpaque(char *name, (thGroupOpaque *) opqptr);*/


