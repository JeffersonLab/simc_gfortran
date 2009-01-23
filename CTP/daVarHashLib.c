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
 *     CODA readout language symbol hashtable related routines
 *	
 * Author:  Jie Chen, CEBAF Data Acquisition Group
 *
 * Revision History:
 * $Log: daVarHashLib.c,v $
 * Revision 1.1  2009/01/23 13:34:01  gaskelld
 * Initial revision
 *
 * Revision 1.2  1999/11/04 20:34:04  saw
 * Alpha compatibility.
 * New RPC call needed for root event display.
 * Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 * Revision 1.1  1996/07/31 20:33:15  saw
 * Initial revision
 *
 *   $Log: daVarHashLib.c,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.2  1999/11/04 20:34:04  saw
 *   Alpha compatibility.
 *   New RPC call needed for root event display.
 *   Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *   Revision 1.1  1996/07/31 20:33:15  saw
 *   Initial revision
 *
*	  Revision 1.1  94/03/15  12:53:11  12:53:11  heyes (Graham Heyes)
*	  Initial revision
*	  
 *	  
 */
#define CASE_INSENSITIVE

#include <stdio.h>
#include <stdlib.h>
#include "daVarHash.h"
#include "daVar.h"

/* hash table size, change to a different prime number if necessary */
/* Some primes: 97 257 1031 2053 4099 */

/*static symbolEntry *hash_table_head[TABLE_SIZE];*/
static void crlSymbol_copy();
static void crlSymbol_delete();
static int  hashFunction();

#ifdef CASE_INSENSITIVE
#define STRCMP daVarComp
#define TOLOWER tolower
#else
#define STRCMP strcmp
#define TOLOWER
#endif

#ifdef _NO_STRDUP
char *strdup (s)
char *s;
{
  char *p = (char *)malloc ((strlen(s) + 1)*sizeof(char));
  if (p == 0) {
    fprintf (stderr, "Cannot allocate memory for strdup\n");
    exit (1);
  }
  strcpy (p, s);
  return p;
}
#endif

/****************************************************************
 *           void crlHashCreate()                               *
 * Hash table creation routine. Hash table is organized in      *
 * such a way to eliminate the collision                        *
 ***************************************************************/
void crlHashCreate(symbolEntry **hash_table_head)
{
  register int i;
   
  for (i = 0; i < TABLE_SIZE; i++){
    hash_table_head[i]=(symbolEntry *)malloc(sizeof(symbolEntry));
    if(hash_table_head[i] == NULL){
      fprintf(stderr,"Cannot allocate memory for crl hash table\n");
      exit(1);
    }
    else{
      hash_table_head[i]->crlSymbol = 0;
      hash_table_head[i]->next=(symbolEntry *)0;
    }
  }
}

/*************************************************************
 *              int crlHashAdd()                             *
 * Add an item to the hashTable, return 1 on success         *
 * return 0 on failure                                       *
 ************************************************************/
int crlHashAdd (CrlSymbol symbol,symbolEntry **hash_table_head)
{
  int hashvalue,i;
  symbolEntry *ptr,*txt;

  hashvalue = hashFunction(symbol, TABLE_SIZE);
/*  printf ("Hash value for is %d\n",hashvalue);*/
  for(ptr = hash_table_head[hashvalue]->next; ptr != NULL; ptr = ptr->next){
    if((STRCMP(ptr->crlSymbol,symbol)) == 0)
      break;
    else
      ;
  }
  if(ptr !=NULL)
    return 0;
  else{
    txt=(symbolEntry *)malloc(sizeof(symbolEntry));
    if(txt == NULL){
      fprintf(stderr,"Cannot allocate memory for this cns entry.\n");
      exit(1);
    }
/*    crlSymbol_copy(&txt->crlSymbol,symbol);*/
    txt->crlSymbol = symbol;
    txt->next = hash_table_head[hashvalue]->next;
    hash_table_head[hashvalue]->next=txt;
    return 1;
  }
}

/************************************************************
 *            int crlHashDelete()                           *
 * delete a single item from the Hash Table                 *
 * return 0 on success, return 1 on failure                 *
 ***********************************************************/
int crlHashDelete(CrlSymbol symbol,symbolEntry **hash_table_head)
{
  int hashvalue;
  symbolEntry *ptr,*qtr;

  hashvalue=hashFunction (symbol, TABLE_SIZE);
  qtr = hash_table_head [hashvalue];
  for(ptr = hash_table_head[hashvalue]->next; ptr!=NULL; ptr=ptr->next){
    if((STRCMP(ptr->crlSymbol, symbol)) == 0)
      break;
    else
      qtr=qtr->next;
  }

  if(ptr !=NULL){
    qtr->next=ptr->next;
    ptr->next=NULL;
    /* free memory */
/*    crlSymbol_delete (&(ptr->crlSymbol));*/
    free (ptr);
    return 0;
  }
  else
    return 1;
}

/**************************************************************
 *            int crlHashFind()                               *
 * Find out whether a particular item which has key 'key'     *
 * existed in the hash Tabel, return item                     *
 * return 0 for failure                                       *
 *************************************************************/
CrlSymbol *crlHashFind(CrlSymbol symbol,symbolEntry **hash_table_head)
{
  int hashvalue;
  symbolEntry *ptr;

  hashvalue = hashFunction (symbol, TABLE_SIZE);
/*  printf ("Hash value inside find for is %d\n", hashvalue);*/
  for(ptr = hash_table_head[hashvalue]->next; ptr != NULL; ptr = ptr->next){
    if((STRCMP(ptr->crlSymbol, symbol)) == 0)
      break;
    else
      ;
  }
  if(ptr != NULL){
    return(&(ptr->crlSymbol));
  }
  else
    return 0;
}

/**********************************************************
 *            int crlHashDestroy()                        *
 * Destroy hashTable and free memory                      *
 * return 0 on success, return 1 on failure               *
 *********************************************************/
#if 0
int crlHashDestroy()
{
  int i;
  symbolEntry *ptr,*qtr;

  for(i=0;i < TABLE_SIZE; i++){
    ptr = hash_table_head[i];
    while(ptr != NULL){
      qtr = ptr->next;
      /* free all memory */
      crlSymbol_delete (&(ptr->crlSymbol));
      free(ptr);
      /* update pointer information */
      ptr = qtr;
    }
  }
  return 0;
}
#endif


/**********************************************************
 *           int hashFunction()                           *
 * return hash value according to the char key            *
 * case sensitive                                         *
 *********************************************************/
static int hashFunction(s,slot_num)
daVarStruct *s;
int slot_num;
{
  char *p;
  unsigned h=0,g;
   
/*  printf("H(%s)=",s->name);*/
  for(p = s->name; *p !='\0'; p++) {
    h= (h<<4) + (TOLOWER(*p));
    if(g = h & 0xf0000000) {
      h=h^(g>>24);
      h=h^g;
    }
  }
/*  printf("%d\n",h%slot_num);*/
  return h%slot_num;
}
extern void crlHashWalk(symbolEntry **hash_table_head,void (*action)())
{
  int i;
  symbolEntry *ptr;

  for(i=0;i < TABLE_SIZE; i++){
    ptr = hash_table_head[i]->next;
    while(ptr != NULL){
      (*action) (ptr->crlSymbol);
      ptr = ptr->next;
    }
  }
}    

/**********************************************************
 *        void crlSymbol_copy()                           *
 * locally used structure copy routine                    *
 * copy expid struct from id1 to id2                      *
 * id2 memory must be allocated before use this routine   *
 **********************************************************/
/*static void crlSymbol_copy(id2,id1) 
CrlSymbol *id2,*id1;
{
  id2->var_name = strdup(id1->var_name);
  id2->var_type = id1->var_type;
}*/

/**********************************************************
 *        void crlSymbol_delete()                         *
 * locally used structure delete routine                  *
 * free memory pointed by id                              *
 **********************************************************/
#if 0
static void crlSymbol_delete(id)
CrlSymbol *id;
{
  free (id->var_name);
}
#endif
/*********************************************************
 *        void crlAddSymbols()                           *
 * Description:                                          *
 *    add list of var names to symbol hashtable          *
 ********************************************************/
#if 0
void crlAddSymbols(var_list,type)
char        *var_list;
int         type;
{
  char *p, *q,temp[80];
  int  status;
  CrlSymbol crlSymbol;

  if(strchr(var_list,',') == NULL){  /* single var name */
    crlSymbol.var_name = strdup(var_list);
    crlSymbol.var_type = type;
    status = crlHashAdd(crlSymbol.var_name, &crlSymbol);
  }
  else{
    p = var_list;
    q = temp;
    while(*p != '\n' && *p != '\0'){
      if(*p == ','){
	*q = '\0';
	p++;
	crlSymbol.var_name = strdup(temp);
	crlSymbol.var_type = type;
	status = crlHashAdd(crlSymbol.var_name,&crlSymbol);
	q = temp;
      }
      else{
	*q = *p;
	q++; p++;
      }
    }
  }
}
#endif
/*****************************************************************
 *        void isSymbolFound()                                   *
 * Description:                                                  *
 *    Check to see whether a var name existed in the hashTable   *
 ****************************************************************/
#if 0
void isSymbolFound(symbol)
char        *symbol;
{
  int status = 0;
  
  status = crlHashFind(symbol);
  if(status != 0){
    fprintf(stderr,"Error: Undefined symbol \"%s\"\n", symbol);
    fprintf(stderr,"Cannot continue, Quit.\n");
    exit(1);
  }
}
#endif
