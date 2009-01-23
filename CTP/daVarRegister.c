/*-----------------------------------------------------------------------------
 * Copyright (c) 1992 Southeastern Universities Research Association,
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
 *	C and Fortran routines for registering variables to be used by
 *      the test, histogram and parameter packages.
 *	
 * Author:  Stephen Wood, CEBAF HALL C
 *
 * Revision History:
 *   $Log: daVarRegister.c,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.3  2003/02/21 20:55:24  saw
 *   Clean up some types and casts to reduce compiler warnings.
 *
 *   Revision 1.2  1999/11/04 20:34:04  saw
 *   Alpha compatibility.
 *   New RPC call needed for root event display.
 *   Start of code to write ROOT trees (ntuples) from new "tree" block
 *
 *   Revision 1.13  1999/08/25 13:16:05  saw
 *   *** empty log message ***
 *
 *   Revision 1.12  1999/03/01 19:51:32  saw
 *   Add Absoft Fortran stuff
 *
 *   Revision 1.11  1997/05/29 18:56:25  saw
 *   Lock changes before adding Absoft(Linux) compatibility
 *
 *	  Revision 1.10  1996/07/31  20:37:53  saw
 *	  Use hash table for name storage.
 *
 *	  Revision 1.9  1994/09/27  20:20:53  saw
 *	  Remove linux dependencies, allow  wild cards in daVarList
 *
 *	  Revision 1.8  1994/08/24  14:27:00  saw
 *	  Have daVarLookupPWithClass return S_DAVAR_UNKNOWN if var not found
 *
 *	  Revision 1.7  1994/06/03  20:59:26  saw
 *	  Replace stderr with STDERR
 *
 *	  Revision 1.6  1994/02/10  21:58:33  saw
 *	  Change node variable name to nd to not conflict with node type.
 *
 *	  Revision 1.5  1994/02/10  18:34:05  saw
 *	  Small fixes for SGI compatibility
 *
 *	  Revision 1.4  1993/11/24  21:37:39  saw
 *	  Add fortran calls for registering double (REAL *8) variable type.
 *
 *	  Revision 1.3  1993/11/22  20:09:42  saw
 *	  Add REGPARMSTRING fortran call for new Fortran string type DAVARFSTRING
 *
 *	  Revision 1.2  1993/08/12  19:58:10  saw
 *	  On HPUX don't use native tsearch.
 *
 *	  Revision 1.1  1993/05/10  20:05:09  saw
 *	  Initial revision
 *
 *
 *  18-dec-92 saw Original version
 *
 *
 * Routines available to general users:
 * --------
 *              daVarRegister(int flag, daVarStruct *args)
 *              daVarLookup(char *name, daVarStruct *results)
 *
 * Routines available to "friendly" packages (e.g.) RPC service routines
 * --------
 *              daVarLookupP(char *name, daVarStruct **results)
 *              daVarList(char ***listp)
 *              daVarFreeList(char **list)
 *
 *
 */
#include "cfortran.h"
#include "daVar.h"

#define USEHASH
#ifdef USEHASH
#include "daVarHash.h"
#else
#if !defined(ultrix)
#define NO_TSEARCH
#endif

#include <stdio.h>
#ifdef NOFNMATCH
#include "fnmatch.h"		/* For non POSIX systems */
#else
#include <fnmatch.h>
#endif

/* Stuff for Tsearch routines*/
typedef struct node_t
{
    daVarStruct *key;
    struct node_t *left, *right;
} node;

#ifdef NO_TSEARCH
typedef enum { preorder, postorder, endorder, leaf } VISIT;
node *mytsearch(void *key, node **rootp, int (* compar)());
node *mytfind(void *key, node **rootp, int (* compar)());
void mytwalk();
#else
#include <search.h>
#define mytsearch tsearch
#define mytfind tfind
#define mytwalk twalk
#endif
#endif

#ifdef USEHASH
symbolEntry *hash_table[TABLE_SIZE];
int hashNotInited=1;
#else
node *daVarRoot=0;
#endif
int daVarCount=0;		/* Used by daVarList */
char **daVarListGlob;
char *daVarListPattern;
int (*daVarListCompFunction)();
int daVarListPattern_length;

/* Local prototypes */
int daVarComp(daVarStruct *item1, daVarStruct *item2);

/* Code */
int daVarRegister(int flag, daVarStruct *args)
/* Should accept a title arg of zero and create a null string in that
   case.
*/
{
  daVarStruct search, *new, **searchresult;
  int fullnamelen;

  if(flag != 0) {
    fprintf(STDERR,
	    "(daVarRegister) Only zero allowed for flag argument now.\n");
    return(S_FAILURE);
  }

  search.name = args->name;
/*  printf("Searching for %s\n",args->name);*/
#ifdef USEHASH
  if(hashNotInited) {crlHashCreate(hash_table); hashNotInited = 0;}
  if(searchresult = (daVarStruct **) crlHashFind((CrlSymbol) &search,hash_table)) {
#else
  if(searchresult = (daVarStruct **) mytfind(&search,&daVarRoot,daVarComp)){
#endif
    fprintf(STDERR,
	    "(daVar) Replacing definition of variable \"%s\" in table\n",
	    args->name);
    free((*searchresult)->title);
    if(args->title) {
      if(((*searchresult)->title = (char *) malloc(strlen(args->title)+1))
	 == NULL)
	return(S_FAILURE);
      strcpy((*searchresult)->title,args->title);
    } else {
      if(((*searchresult)->title = (char *) malloc(1))
	 == NULL)
	return(S_FAILURE);
      (*searchresult)->title[0] = '\0';
    }
    (*searchresult)->varptr = args->varptr;
    (*searchresult)->size = (args->size<=0) ? 1 : args->size;
    (*searchresult)->type = args->type;
    (*searchresult)->flag = args->flag;
    (*searchresult)->rhook = args->rhook;
    (*searchresult)->whook = args->whook;
    (*searchresult)->opaque = args->opaque;
    return(S_DAVAR_REPLACED);
  } else {
    if((new = (daVarStruct *) malloc(sizeof(daVarStruct))) == NULL)
      return(S_FAILURE);
    if((new->name =  (char *) malloc(strlen(args->name)+1)) == NULL)
       return(S_FAILURE);
    strcpy(new->name,args->name);


    if(args->title) {
      if((new->title =  (char *) malloc(strlen(args->title)+1)) == NULL)
	return(S_FAILURE);
      strcpy(new->title,args->title);
    } else {
      if((new->title =  (char *) malloc(1)) == NULL)
	return(S_FAILURE);
      new->title[0] = '\0';
    }
    new->type = args->type;
    new->varptr = args->varptr;
    new->size = (args->size<=0) ? 1 : args->size;
    new->flag = args->flag;
    new->rhook = args->rhook;
    new->whook = args->whook;
    new->opaque = args->opaque;

#ifdef USEHASH
    if(crlHashAdd((CrlSymbol) new, hash_table))
#else    
    if(mytsearch((void *) new,&daVarRoot,daVarComp))
#endif
      return(S_SUCCESS);
    else
      return(S_FAILURE);
  }
}


int daVarLookup(char *name, daVarStruct *result)
{
  daVarStruct search, **searchresult;
  static char *namel=0;		/* Pointers to  static space for copies of */
  static int namelsize=0;
  static char *titlel=0;	/* the name and title pointers */
  static int titlelsize=0;
  int len;

  search.name = name;
#ifdef USEHASH
  if(searchresult = (daVarStruct **) crlHashFind((CrlSymbol) &search,hash_table)) {
#else
  if(searchresult = (daVarStruct **) mytfind(&search,&daVarRoot,daVarComp)){
#endif

    len=strlen((*searchresult)->name);
    if(len >= namelsize) {
      if(namel) free(namel);
      namel = (char *) malloc(len+1);
      namelsize = len+1;
    }
    strcpy(namel,(*searchresult)->name);
    result->name = namel;

    len=strlen((*searchresult)->title);
    if(len >= titlelsize) {
      if(titlel) free(titlel);
      titlel = (char *) malloc(len + 1);
      titlelsize = len+1;
    }
    strcpy(titlel,(*searchresult)->title);
    result->title = titlel;

    result->type = (*searchresult)->type;
    result->varptr = (*searchresult)->varptr;
    result->size = (*searchresult)->size;
    result->opaque = (*searchresult)->opaque;
    result->rhook = (*searchresult)->rhook;
    result->whook = (*searchresult)->whook;
    return(S_SUCCESS);
  } else
    return(S_DAVAR_UNKNOWN);
}
int daVarStrcmp(register char *s1, register char *s2)
{
  while(toupper(*s1) == toupper(*s2++))
    if(*s1++ == '\0')
      return(0);
  return(toupper(*s1) - toupper(*--s2));
}
int daVarFnmatch(register char *pattern, register char *s, register int n)
{
  return(fnmatch(pattern,s,0));
}
int daVarStrncmp(register char *s1, register char *s2, register int n)
{
  while(toupper(*s1) == toupper(*s2++))
    if(*s1++ == '\0' || (--n) <= 0)
      return(0);
  return(toupper(*s1) - toupper(*--s2));
}

int daVarComp(daVarStruct *item1, daVarStruct *item2)
/* Do case insensitive comparisons of keys */
{
  return(daVarStrcmp(item1->name,item2->name));
}

int daVarLookupP(char *name, daVarStruct **varstructptr)
{
  daVarStruct search, **searchresult;

  search.name = name;
#ifdef USEHASH
  if(searchresult = (daVarStruct **) crlHashFind((CrlSymbol) &search,hash_table)) {
#else
  if(searchresult = (daVarStruct **) mytfind(&search,&daVarRoot,daVarComp)){
#endif
    *varstructptr = *searchresult;
    return(S_SUCCESS);
  } else
    return(S_DAVAR_UNKNOWN);
}

daVarLookupPWithClass(char *name, char **prefixlist, daVarStruct **varp)
{ 
  int namlen,namtrylen;
  char *namtry;

  namlen = strlen(name);
  if(daVarLookupP(name,varp)==S_SUCCESS) return(S_SUCCESS);
  namtrylen = namlen + 10;
  namtry = (char *) malloc(namtrylen+1);
  while(*prefixlist){
    int thislen;
    thislen = strlen(*prefixlist) + namlen + 1;
    if(thislen > namtrylen) {
      namtrylen = thislen;
      namtry = (char *) realloc(namtry,namtrylen);
    }
    strcpy(namtry,*prefixlist);
    strcat(namtry,".");
    strcat(namtry,name);
    if(daVarLookupP(namtry,varp)==S_SUCCESS) {
      free(namtry);
      return(S_SUCCESS);
    }
    prefixlist++;
  }
  free(namtry);
  return(S_DAVAR_UNKNOWN);	/* Variable not registered */
}

void daVarCount_node
#ifdef USEHASH
(void *entry)
{
#else
(node *nd,VISIT order, int level)
{
  if(order==postorder || order == leaf)
#endif
  daVarCount++;
}

void daVarList_node
#ifdef USEHASH
(void *entry)
{
#else
(node *nd,VISIT order, int level)
{
  if(order==postorder || order == leaf)
#endif
    {
      char *name;

      name = ((daVarStruct *) entry)->name;
      if(daVarListPattern_length == 0 ||
	 (daVarListCompFunction)(daVarListPattern,name,daVarListPattern_length) == 0)
	daVarListGlob[daVarCount++] = name;
    }
}


int daVarList(char *pattern, char ***listp, int *count)
/* User is not allowed to muck with the strings pointed to in the list
   because they are the actual strings in the tables. */
{

  if(strchr(pattern,'*') || strchr(pattern,'?')) {
    daVarListCompFunction = daVarFnmatch;
  } else {
    daVarListCompFunction = daVarStrncmp;
  }
  if(pattern) {
    daVarListPattern = pattern;
    daVarListPattern_length = strlen(daVarListPattern);
  } else
    daVarListPattern_length = 0;
  daVarCount = 0;
#ifdef USEHASH
  crlHashWalk(hash_table,daVarCount_node);
#else
  mytwalk(daVarRoot,daVarCount_node);/* Should make list only big enough
					for what matches */
#endif

  if((*listp = daVarListGlob = 
      (char **) malloc((daVarCount)*sizeof(char *))) == NULL)
    return(S_FAILURE);
  daVarCount = 0;
#ifdef USEHASH
  crlHashWalk(hash_table,daVarList_node);
#else
  mytwalk(daVarRoot,daVarList_node);
#endif
  *count = daVarCount;
  return(S_SUCCESS);
}
#ifndef USEHASH
daVarPrint_node(node *nd,VISIT order, int level)
{
  char *name,*title;

  if(order==postorder || order == leaf) {
    name = ((daVarStruct *) nd->key)->name;
    title = ((daVarStruct *) nd->key)->title;
    printf("XX: %s %s %x %x\n",name,title,nd,nd->key);
  }
}

int daVarPrint()
{
  mytwalk(daVarRoot,daVarPrint_node);
  return(S_SUCCESS);
}
#endif
int daVarFreeList(char **list)
/* Free's up the list of variables in listp */
{
  int i;

  free(list);
  return(S_SUCCESS);
}

/* Fortran entry points */

#define LENDEFARRAY int *size,
#define LENDEFSCALER
#define LENARGARRAY *size
#define LENARGSCALER 1

#define MAKEFSUB(SUBNAME,CLASS,TYPENAME,DATYPE,ARRAY) \
int SUBNAME(char *name, TYPENAME *vptr, LENDEF##ARRAY char *title\
	      ,unsigned l_name, unsigned l_title)\
{\
  int A0;\
  daVarStruct args;\
  char *BN=0;\
  char *BT=0;\
  char *BF = 0;\
\
  BF = malloc(strlen(CLASS)+l_name+1);\
  strcpy(BF,CLASS);\
  args.name = strcat(BF,((!*(int *)name)?0:memchr(name,'\0',l_name)?name:\
			 (memcpy(BN=(char *) malloc(l_name+1),name,l_name)\
			  ,BN[l_name]='\0',kill_trailing(BN,' '))));\
  args.title = ((!*(int *)title)?0:memchr(title,'\0',l_title)?title:\
		 (memcpy(BT=(char *) malloc(l_title+1),title,l_title)\
		  ,BT[l_title]='\0',kill_trailing(BT,' ')));\
  args.size = LENARG##ARRAY;\
  args.varptr = (void *) vptr;\
  args.flag = DAVAR_READWRITE;\
  args.type = DATYPE;\
  args.opaque = 0;\
  args.rhook = 0;\
  args.whook = 0;\
  A0 = daVarRegister((int) 0, &args);\
  if(BF) free(BF);\
  if(BN) free(BN);\
  if(BT) free(BT);\
  return(A0);\
}

/* Can't figure out a more clever way */
#ifdef AbsoftUNIXFortran
MAKEFSUB(regreal,"",float,DAVARFLOAT,SCALER)
MAKEFSUB(regdouble,"",double,DAVARDOUBLE,SCALER)
MAKEFSUB(regint,"",int,DAVARINT,SCALER)
MAKEFSUB(regrealarray,"",float,DAVARFLOAT,ARRAY)
MAKEFSUB(regdoublearray,"",double,DAVARDOUBLE,ARRAY)
MAKEFSUB(regintarray,"",int,DAVARINT,ARRAY)

MAKEFSUB(regparmreal,"parm.",float,DAVARFLOAT,SCALER)
MAKEFSUB(regparmdouble,"parm.",double,DAVARDOUBLE,SCALER)
MAKEFSUB(regeventreal,"event.",float,DAVARFLOAT,SCALER)
MAKEFSUB(regeventdouble,"event.",double,DAVARDOUBLE,SCALER)
MAKEFSUB(regparmint,"parm.",int,DAVARINT,SCALER)
MAKEFSUB(regeventint,"event.",int,DAVARINT,SCALER)
MAKEFSUB(regparmrealarray,"parm.",float,DAVARFLOAT,ARRAY)
MAKEFSUB(regparmdoublearray,"parm.",double,DAVARDOUBLE,ARRAY)
MAKEFSUB(regeventrealarray,"event.",float,DAVARFLOAT,ARRAY)
MAKEFSUB(regeventdoublearray,"event.",double,DAVARDOUBLE,ARRAY)
MAKEFSUB(regparmintarray,"parm.",int,DAVARINT,ARRAY)
MAKEFSUB(regeventintarray,"event.",int,DAVARINT,ARRAY)

MAKEFSUB(regtestint,"test.",int,DAVARINT,SCALER)
MAKEFSUB(regtestintarray,"test.",int,DAVARINT,ARRAY)
#else
MAKEFSUB(regreal_,"",float,DAVARFLOAT,SCALER)
MAKEFSUB(regdouble_,"",double,DAVARDOUBLE,SCALER)
MAKEFSUB(regint_,"",int,DAVARINT,SCALER)
MAKEFSUB(regrealarray_,"",float,DAVARFLOAT,ARRAY)
MAKEFSUB(regdoublearray_,"",double,DAVARDOUBLE,ARRAY)
MAKEFSUB(regintarray_,"",int,DAVARINT,ARRAY)

MAKEFSUB(regparmreal_,"parm.",float,DAVARFLOAT,SCALER)
MAKEFSUB(regparmdouble_,"parm.",double,DAVARDOUBLE,SCALER)
MAKEFSUB(regeventreal_,"event.",float,DAVARFLOAT,SCALER)
MAKEFSUB(regeventdouble_,"event.",double,DAVARDOUBLE,SCALER)
MAKEFSUB(regparmint_,"parm.",int,DAVARINT,SCALER)
MAKEFSUB(regeventint_,"event.",int,DAVARINT,SCALER)
MAKEFSUB(regparmrealarray_,"parm.",float,DAVARFLOAT,ARRAY)
MAKEFSUB(regparmdoublearray_,"parm.",double,DAVARDOUBLE,ARRAY)
MAKEFSUB(regeventrealarray_,"event.",float,DAVARFLOAT,ARRAY)
MAKEFSUB(regeventdoublearray_,"event.",double,DAVARDOUBLE,ARRAY)
MAKEFSUB(regparmintarray_,"parm.",int,DAVARINT,ARRAY)
MAKEFSUB(regeventintarray_,"event.",int,DAVARINT,ARRAY)

MAKEFSUB(regtestint_,"test.",int,DAVARINT,SCALER)
MAKEFSUB(regtestintarray_,"test.",int,DAVARINT,ARRAY)
#endif

/* Entry points for String registration.  Do entry points for anything other
than parmameters make sense? */
#ifdef AbsoftUNIXFortran
int regparmstring
#else
int regparmstring_
#endif
(char *name, char *sptr, char *title
		    ,unsigned l_name, unsigned l_sptr, unsigned l_title)
{
  int A0;
  daVarStruct args;
  char *BN=0;

  char *BT=0;
  char *BF = 0;

  BF = malloc(5+l_name+1);
  strcpy(BF,"parm.");
  args.name = strcat(BF,((!*(int *)name)?0:memchr(name,'\0',l_name)?name:
			 (memcpy(BN=(char *) malloc(l_name+1),name,l_name)
			  ,BN[l_name]='\0',kill_trailing(BN,' '))));
  args.title = ((!*(int *)title)?0:memchr(title,'\0',l_title)?title:
		 (memcpy(BT=(char *) malloc(l_title+1),title,l_title)
		  ,BT[l_title]='\0',kill_trailing(BT,' ')));
  args.size = l_sptr;
  args.varptr = (void *) sptr;
  args.flag = DAVAR_READWRITE;
  args.type = DAVARFSTRING;
  args.opaque = 0;
  args.rhook = 0;
  args.whook = 0;
  A0 = daVarRegister((int) 0, &args);
  if(BF) free(BF);
  if(BN) free(BN);
  return(A0);
}

#ifdef NO_TSEARCH
/*
 * Tree search generalized from Knuth (6.2.2) Algorithm T just like
 * the AT&T man page says.
 *
 * The node_t structure is for internal use only, lint doesn't grok it.
 *
 * Written by reading the System V Interface Definition, not the code.
 *
 * Totally public domain.
 */
/*LINTLIBRARY*/

/*
#include <search.h>

typedef struct node_t
{
    char	  *key;
    struct node_t *left, *right;
}
node;
*/

node *mytsearch(key, rootp, compar)
/* find or insert datum into search tree */
void 	*key;			/* key to be located */
register node	**rootp;	/* address of tree root */
int	(*compar)();		/* ordering function */
{
    register node *q;

    if (rootp == (struct node_t **)0)
	return ((struct node_t *)0);
    while (*rootp != (struct node_t *)0)	/* Knuth's T1: */
    {
	int r;

	if ((r = (*compar)(key, (*rootp)->key)) == 0)	/* T2: */
	    return (*rootp);		/* we found it! */
	rootp = (r < 0) ?
	    &(*rootp)->left :		/* T3: follow left branch */
	    &(*rootp)->right;		/* T4: follow right branch */
    }
    q = (node *) malloc(sizeof(node));	/* T5: key not found */
    if (q != (struct node_t *)0)	/* make new node */
    {
	*rootp = q;			/* link new node to old */
	q->key = key;			/* initialize new node */
	q->left = q->right = (struct node_t *)0;
    }
    return (q);
}

node *mytdelete(key, rootp, compar)
/* delete node with given key */
char	*key;			/* key to be deleted */
register node	**rootp;	/* address of the root of tree */
int	(*compar)();		/* comparison function */
{
    node *p;
    register node *q;
    register node *r;
    int cmp;

    if (rootp == (struct node_t **)0 || (p = *rootp) == (struct node_t *)0)
	return ((struct node_t *)0);
    while ((cmp = (*compar)(key, (*rootp)->key)) != 0)
    {
	p = *rootp;
	rootp = (cmp < 0) ?
	    &(*rootp)->left :		/* follow left branch */
	    &(*rootp)->right;		/* follow right branch */
	if (*rootp == (struct node_t *)0)
	    return ((struct node_t *)0);	/* key not found */
    }
    r = (*rootp)->right;			/* D1: */
    if ((q = (*rootp)->left) == (struct node_t *)0)	/* Left (struct node_t *)0? */
	q = r;
    else if (r != (struct node_t *)0)		/* Right link is null? */
    {
	if (r->left == (struct node_t *)0)	/* D2: Find successor */
	{
	    r->left = q;
	    q = r;
	}
	else
	{			/* D3: Find (struct node_t *)0 link */
	    for (q = r->left; q->left != (struct node_t *)0; q = r->left)
		r = q;
	    r->left = q->right;
	    q->left = (*rootp)->left;
	    q->right = (*rootp)->right;
	}
    }
    free((struct node_t *) *rootp);	/* D4: Free node */
    *rootp = q;				/* link parent to new node */
    return(p);
}

static void trecurse(root, action, level)
/* Walk the nodes of a tree */
register node	*root;		/* Root of the tree to be walked */
register void	(*action)();	/* Function to be called at each node */
register int	level;
{
    if (root->left == (struct node_t *)0 && root->right == (struct node_t *)0)
	(*action)(root, leaf, level);
    else
    {
	(*action)(root, preorder, level);
	if (root->left != (struct node_t *)0)
	    trecurse(root->left, action, level + 1);
	(*action)(root, postorder, level);
	if (root->right != (struct node_t *)0)
	    trecurse(root->right, action, level + 1);
	(*action)(root, endorder, level);
    }
}

void mytwalk(root, action)		/* Walk the nodes of a tree */
node	*root;			/* Root of the tree to be walked */
void	(*action)();		/* Function to be called at each node */
{
    if (root != (node *)0 && action != (void(*)())0)
	trecurse(root, action, 0);
}

/* mytsearch.c ends here */
/*
 * Tree search generalized from Knuth (6.2.2) Algorithm T just like
 * the AT&T man page says.
 *
 * The node_t structure is for internal use only, lint doesn't grok it.
 *
 * Written by reading the System V Interface Definition, not the code.
 *
 * Totally public domain.
 */
/*LINTLIBRARY*/
/*
#include <search.h>

typedef struct node_t
{
    char	  *key;
    struct node_t *left, *right;
} node;
*/

node *mytfind(key, rootp, compar)
/* find a node, or return 0 */
void		*key;		/* key to be found */
register node	**rootp;	/* address of the tree root */
int		(*compar)();	/* ordering function */
{
    if (rootp == (struct node_t **)0)
	return ((struct node_t *)0);
    while (*rootp != (struct node_t *)0)	/* T1: */
    {
	int r;
	if ((r = (*compar)(key, (*rootp)->key)) == 0)	/* T2: */
	    return (*rootp);		/* key found */
	rootp = (r < 0) ?
	    &(*rootp)->left :		/* T3: follow left branch */
	    &(*rootp)->right;		/* T4: follow right branch */
    }
    return (node *)0;
}
#endif
