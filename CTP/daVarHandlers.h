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
 *  Header for RPC handlers
 *	
 * Author:  Stephen A. Wood, CEBAF Hall C
 *
 * Revision History:
 *  $Log: daVarHandlers.h,v $
 *  Revision 1.1  2009/01/23 13:34:01  gaskelld
 *  Initial revision
 *
 *  Revision 1.2  1999/11/04 20:34:04  saw
 *  Alpha compatibility.
 *  New RPC call needed for root event display.
 *  Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *  Revision 1.4  1999/08/25 13:16:05  saw
 *  *** empty log message ***
 *
 *  Revision 1.3  1994/11/07 14:31:34  saw
 *  Add pending callback requests structure definition
 *
 *	  Revision 1.2  1993/05/11  17:34:58  saw
 *	  Fix $Log: daVarHandlers.h,v $
 *	  Fix Revision 1.1  2009/01/23 13:34:01  gaskelld
 *	  Fix Initial revision
 *	  Fix
 *	  Fix Revision 1.2  1999/11/04 20:34:04  saw
 *	  Fix Alpha compatibility.
 *	  Fix New RPC call needed for root event display.
 *	  Fix Start of code to write ROOT trees (ntuples) from new "tree" block
 *	  Fix
 *	  Fix Revision 1.4  1999/08/25 13:16:05  saw
 *	  Fix *** empty log message ***
 *	  Fix
 *	  Fix Revision 1.3  1994/11/07 14:31:34  saw
 *	  Fix Add pending callback requests structure definition
 *	  Fix
 *
 */

void daVarReadVar(char *name, any *retval);
daVarStatus daVarWriteVar(char *name, any *retval);
daVarStatus daVarRegRatr(daVarStruct *varp, char *attribute
			 ,int index, any *retval);
daVarStatus daVarRegWatr(daVarStruct *varp, char *attribute
			 ,int index, any *setval);

#define DAVAR_VALUE "value"
#define DAVAR_TITLE "title"
#define DAVAR_SIZE "size"
#define DAVAR_FLAG "flag"
#define DAVAR_TYPE "type"
#define DAVAR_WATR "watr"
#define DAVAR_RATR "ratr"

#define DAVAR_NOINDEX -12345678

/* Structures for holding information about call back requests */

struct daVarCallBackList {
  struct sockaddr_in *sock_in;	/* Caller socket info */
  struct daVarCallBackList *next;
  struct TESTNAMELIST *list;
  time_t start_time;
};
typedef struct daVarCallBackList daVarCallBackList;
