/*-----------------------------------------------------------------------------
 * Copyright (c) 1991,1992 Southeastern Universities Research Association,
 *                         Continuous Electron Beam Accelerator Facility
 *
 * This software was developed under a United States Government license
 * described in the NOTICE file included as part of this distribution.
 *
 * CEBAF Data Acquisition Group, 12000 Jefferson Ave., Newport News, VA 23606
 * Email: coda@cebaf.gov  Tel: (804) 249-7101  Fax: (804) 249-7363
 *-----------------------------------------------------------------------------
 * 
 * Description:
 *     header file for CODA readout language symbol hash table
 *	
 * Author:  Jie Chen, CEBAF Data Acquisition Group
 *
 * Revision History:
 *   $Log: daVarHash.h,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.1  1998/12/07 22:11:09  saw
 *   Initial setup
 *
*	  Revision 1.1  94/03/15  12:53:09  12:53:09  heyes (Graham Heyes)
*	  Initial revision
*	  
 *	  
 */
#ifndef _crl_hash_h
#define _crl_hash_h

#define CTPHASH
#ifndef CTPHASH
typedef struct _symbol{
  char *var_name;
  int  var_type;      /*0: integer, 1: unsigned long */
}CrlSymbol;
#else
typedef void *CrlSymbol;
#endif
typedef struct _SLOT_ENTRY
{
 CrlSymbol         crlSymbol;
 struct _SLOT_ENTRY *next;
}symbolEntry;

extern void crlHashCreate(symbolEntry **hash_table_head);
extern int  crlHashAdd(CrlSymbol symbol,symbolEntry **hash_table_head);
extern int  crlHashDelete(CrlSymbol symbol,symbolEntry **hash_table_head);
CrlSymbol *crlHashFind(CrlSymbol symbol,symbolEntry **hash_table_head);
extern void crlHashWalk(symbolEntry **hash_table_head,void (*action)());
/*extern int  crlHashDestroy();*/
/*extern void crlAddSymbols();*/
/*extern void isSymbolFound();*/

#define TABLE_SIZE 2053


#endif
