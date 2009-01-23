/*-----------------------------------------------------------------------------
 * Copyright (c) 1999      Thomas Jefferson National Accelerator Facility
 *
 * This software was developed under a United States Government license
 * described in the NOTICE file included as part of this distribution.
 *
 * Stephen A. Wood, 12000 Jefferson Ave., Newport News, VA 23606
 * Email: saw@jlab.org  Tel: (758) 269-7367  Fax: (757) 269-5235
 *-----------------------------------------------------------------------------
 * 
 * Description:
 *  Book ROOT trees.
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: thTree.c,v $
 *   Revision 1.1  2009/01/23 13:34:01  gaskelld
 *   Initial revision
 *
 *   Revision 1.6.4.1  2008/09/25 00:54:05  jones
 *   Updated for running on Fedora 8 with gfortran
 *
 *   Revision 1.7  2008/09/25 00:01:29  jones
 *   Updated to run with gfortran compiler
 *
 *   Revision 1.6.6.1  2007/09/10 21:32:47  pcarter
 *   Implemented changes to allow compilation on RHEL 3,4,5 and MacOSX
 *
 *   Revision 1.6  2005/02/22 16:54:46  saw
 *   Clean up some diagnostic printfs
 *
 *   Revision 1.5  2004/07/09 20:44:11  saw
 *   Can now put a test on a tree block
 *
 *
 *   Revision 1.1.16.2  2004/07/09 20:41:50  saw
 *   Can now put a test on a tree block
 *
 *   Revision 1.1.16.1  2004/07/09 14:12:11  saw
 *   Add ability for CTP to make ROOT Trees
 *
 *   Revision 1.4  2004/07/08 20:07:00  saw
 *   Supply dummy routines when ROOTSYS not defined
 *
 *   Revision 1.3  2004/07/07 18:16:30  saw
 *   Use properly exported names from thRootStuff.cpp
 *
 *   Revision 1.2  2004/07/02 18:46:29  saw
 *   Update ugly cpp routine for gcc 3.2.3.  Need to find a better way to
 *   reference C++ routines.
 *
 *   Revision 1.1  2002/07/31 20:07:48  saw
 *   Add files for ROOT Trees
 *
 *   Revision 1.1  1999/08/25 13:16:07  saw
 *   *** empty log message ***
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "daVar.h"
#include "th.h"
#include "thInternal.h"
#include "thUtils.h"

#ifdef ROOTTREE
extern daVarStatus thTreeRHandler();

struct thLeafList {		/* Variable and index list */
  struct thLeafList *next;
  char *name;
  int leafsize;			/* Number of bytes in leaf */
  char leaftype;		/* Single letter with leaf type */
  daVarStruct *varp; int index;
};
typedef struct thLeafList thLeafList;

struct thTreeBranchList {
  struct thTreeBranchList *next;
  char *branchname;		/* Block name without the "block.hist" */
  void *evstruct;		/* Event structure to fill with data */
  struct thLeafList *leaflistp;
};
typedef struct thTreeBranchList thTreeBranchList;

struct thTreeOpaque {		/* Opaque structure for histogram definition */
  void *file;
  /*  void *file Need to point to structure that has file info.  File structures will also
  be on a linked list so that we a new file is opened, you can see if it exists.  Wait, we
should just use davars.  Make a new class for files? ? 
 */
  void *treeptr;		/* Pointer to tree object */
  daVarStruct *test; int testindex;
  thTreeBranchList *branchlistP;
};
typedef struct thTreeOpaque thTreeOpaque;

struct thRBlockList {
  struct thRBlockList *next;
  char *blockname;		/* Block name without the "block.tree" */
  daVarStruct *var;		/* Varptr points to # times called
				   Title is code from file.
				   opaque is pointer to hist speclist */
};
typedef struct thRBlockList thRBlockList;

thRBlockList *thRBlockListP;	/* Pointer to list of tree blocks */

thStatus thBookaBranch(thTreeOpaque *tree, char *line, thTreeBranchList **thBranchNext);
     /*thStatus thExecuteaHist(thHistSpecList *Hist);*/
thStatus thRemoveTree(char *block_name);

thStatus thBookTree(daVarStruct *var)
{
  char *lines,*eol;
  int line_count;
  thTreeBranchList **thBranchNext;
  char *blockname;
  thTreeOpaque *treedef;
  int i;
  int isopen;
  char *lbuf=0;

  /*  thHistZeroLastId(); */
  
  /* printf("In booktrees\n");*/
  /* Get the name without the block.test on it */
  blockname = var->name;	/* If name doesn't fit pattern, use whole */
  if(strcasestr(var->name,BLOCKSTR)==var->name){
    i = strlen(BLOCKSTR) + 1;
    if(strcasestr((var->name + i),TREESTR)==(var->name + i)){
      i += strlen(TREESTR);
      if(*(var->name + i) == '.'){
	blockname += i + 1;
      }
    }
  }

  /*printf("Booking tree %s\n",blockname);*/

  if(var->opaque) thRemoveTree(blockname);

  /* We assume for now that thRemoveTree completely removed the opaque definition)
   */
  treedef = var->opaque = (thTreeOpaque *) malloc(sizeof(thTreeOpaque));
  thBranchNext = (thTreeBranchList **) &treedef->branchlistP;
  
  lines = var->title;
  line_count = 0;
  isopen = 0;
  while(*lines){
    char *lcopy;

    line_count++;
    eol = strchr(lines,'\n');
    if(!eol) {
      fprintf(stderr,"L %d: Last line of hist block %s has no newline\n"
	      ,line_count,var->name);
      break;
    }
    if(*(eol+1)=='\0'){		/* This is the last line */
      if(strcasestr(lines,ENDSTR) == 0)
	fprintf(stderr,"L %d: Last line of tree block %s is not an END\n"
		,line_count,var->name);
      break;
    }
    if(line_count == 1) {
      char *fname=0;
      if(strcasestr(lines,BEGINSTR) !=0) {
	char *p;
	if(p = strcasestr(lines,"file=")) {
	  p += 5;
	  /* If " or ', grab to next matching char */
	  /* other interpret as variable.  But interpret as file if variable not found */
	  if(*p == QUOTECHAR1 || *p == QUOTECHAR2) {
	    char *s; int len;
	    s = p+1;
	    while(*s && *s != *p && *s !='\n') s++;
	    len = (s - p) - 1;
	    fname = malloc(len+1);
	    strncpy(fname,p+1,len);
	    fname[len] = '\0';
	  } else {		/* Probably a variable */
	    char *varname=0; char *s; int len; int index;
	    daVarStruct *varp;
	    s = p;
	    while(*s && !isspace(*s) && *s !='\n') s++;
	    len = (s-p);
	    varname = malloc(len+1);
	    strncpy(varname,p,len);
	    varname[len] = '\0';
	    if(thVarResolve(varname,&varp,&index,1,0)==S_SUCCESS) {
	      /*printf("%s,type=%d, size=%d\n",varp->name,varp->type,varp->size);*/
	      if(varp->type == DAVARSTRING) {
		fname = malloc(strlen((char *)varp->varptr)+1);
		strcpy(fname,(char *)varp->varptr);
	      } else if(varp->type == DAVARFSTRING) {
		fname = malloc(varp->size+1);
		strncpy(fname,(char *)varp->varptr,varp->size);
		fname[varp->size] = '\0';
		p = fname;
		while(*p && !isspace(*p)) p++;
		*p = '\0';	/* Null terminate at first blank */
	      }
	      /*printf("|%s|\n", fname);*/
	    }
	    if(!fname) {
	      fname = malloc(len+1);
	      strncpy(fname,p,len);
	      fname[len] = '\0';
	    }
	  }
	}
	if(p = strcasestr(lines,"test=")) {
	  /* RHS must be a variable */
	  char *varname=0; char *s; int len, testindex;
	  daVarStruct *testp;
	  p += 5;
	  s = p;
	  while(*s && !isspace(*s) && *s !='\n') s++;
	  len = (s-p);
	  varname = (char *) malloc(len+1);
	  strncpy(varname,p,len);
	  varname[len] = '\0';
	  if(thVarResolve(varname,&testp,&testindex,1,0) != S_SUCCESS) {
	    return(S_FAILURE); /* Test flag not registered */
      /* ASAP we must change this to register variables as they are needed */
      /* If the variable exists, then we also must check to make sure that
	 the requested index does not exceed the size of the array.
	 a new thVarResolve should also increase the size of the array if
	 it was created by CTP */
	  }
	  treedef->test = testp;
	  treedef->testindex = testindex;
	} else {
	  treedef->test = 0; /* No test, always true */
	}
      }
      
      if(fname) {
	printf("Opening Root file %s\n",fname);
	treedef->file = (void *) thRoot_TFile(fname);
	free(fname);
      } else {
	/*printf("Opening Root file %s\n","ctp.tree");*/
	treedef->file = (void *) thRoot_TFile("ctp.root");
      }
    
      /*printf("Call to TTree(\"%s\",\"title\") goes here\n",blockname);*/
      treedef->treeptr = (void *) thRoot_TTree(blockname);

      if(strcasestr(lines,BEGINSTR) != 0){
	/*	printf("Is a begin\n");*/
	lines = eol + 1;
	continue;
      } else
	fprintf(stderr,"First line of tree block %s is not a BEGIN\n",var->name);
    }
    /* Ready to book the line, Add continuation lines later */
    lcopy = (char *) malloc(eol-lines+1);
    strncpy(lcopy,lines,(eol-lines));
    *(lcopy + (eol-lines)) = '\0';
    if(!thCleanLine(lcopy)){
      if(strchr(lcopy,'=')) {	/* Start of a new branch */
	if(lbuf) {		/* Do we have a pending branch */
	  /*printf("Passing 1 |%s|\n",lbuf);*/
	  if(thBookaBranch(treedef,lbuf,thBranchNext)==S_SUCCESS){
	    thBranchNext = &((*thBranchNext)->next);
	  } else {
	    fprintf(stderr,"(%s): Tree booking error in line %d\n",var->name,line_count);
	  }
	  free(lbuf);
	  lbuf = 0;
	}
      }
      if(lbuf) {		/* Append */
	char *lastcomma;
	char *firstcomma;
	int addcomma=0; 
	lastcomma = lbuf + strlen(lbuf) - 1;
	while(*lastcomma == ' ') lastcomma--;
	if(*lastcomma != ',' && *lastcomma != '=') lastcomma = 0;
	firstcomma = lcopy;
	while(*firstcomma == ' ') firstcomma++;
	if(*firstcomma != ',') firstcomma = 0;
	if(firstcomma && lastcomma) {
	  *firstcomma = ' ';
	} else if (!firstcomma && !lastcomma) {
	  addcomma = 1;
	}
	lbuf = realloc(lbuf,strlen(lbuf) + strlen(lcopy) + 2);
	if(addcomma) strcat(lbuf,",");
	strcat(lbuf,lcopy);
      } else {			/* Make new lbuf */
	lbuf = malloc(strlen(lcopy)+1);
	strcpy(lbuf,lcopy);
      }
    }
    free(lcopy);
    lines = eol+1;
  }
  if(lbuf) {			/* Do the last branch we were building */
    /*printf("Passing 2 |%s|\n",lbuf);*/
    if(thBookaBranch(treedef,lbuf,thBranchNext)==S_SUCCESS){
      thBranchNext = &((*thBranchNext)->next);
    } else {
      fprintf(stderr,"(%s): Tree booking error in line %d\n",var->name,line_count);
    }
    free(lbuf);
  }
  /* Update internal table of trees. */
  {
    thRBlockList *thisblock,*nextblock,**lastblockp;
    nextblock = thRBlockListP;
    lastblockp = &thRBlockListP;
    thisblock = thRBlockListP;
    while(thisblock){
      if((strcasecmp(thisblock->var->name,var->name)) == 0){
	/* Replacing a block with a new definition */
	fprintf(stderr,"Replacing %s with new definition\n",var->name);
	if(thisblock->var != var){
	  fprintf(stderr,"ERR: Same name, different var pointer\n");
	}
	break;
      }
      lastblockp = &thisblock->next;
      thisblock = thisblock->next;
    }
    if(!thisblock){		/* Create entry for New block */
      *lastblockp = thisblock = (thRBlockList *) malloc(sizeof(thRBlockList));
      thisblock->var = var;
      thisblock->next = (thRBlockList *) NULL;
      thisblock->blockname = (char *) malloc(strlen(blockname) + 1);
      strcpy(thisblock->blockname,blockname);
    }
  }
  /*printf("Returning from booking a tree\n");*/
  return(S_SUCCESS);
}

thStatus thBookaBranch(thTreeOpaque *treedef, char *line, thTreeBranchList **thBranchNext)
     /* Interpret a branch def of the form branch=leaf1,leaf2,leaf3,... */
     /* For now require the "branch=" part */
{

  /*  char *long_title;*/
  int n,nleafs;
  int lenbrancharg;
  char *brancharg;
  char *sleafs,*branchname;
  thTreeBranchList *Branch;
  thLeafList **LeafNext;
  thLeafList *thisleaf;
  daVarStruct *varp;
  int vind;
  char *args[100];

  /*printf("In thBookaBranch\n");*/
  if(!(sleafs = strchr(line,'='))) {
    return(S_FAILURE);
  }    
  *sleafs=0;
  sleafs++;			/* Pointer to list of leaves */
  nleafs = thCommas(sleafs,args);
  if(nleafs <=0) {
    return(S_FAILURE);
  }
  Branch = *thBranchNext = (thTreeBranchList *) malloc(sizeof(thTreeBranchList));
  Branch->next = (thTreeBranchList *) NULL;
  branchname = thSpaceStrip(line);
  Branch->branchname = malloc(strlen(branchname)+1);
  LeafNext = (thLeafList **) &Branch->leaflistp;
  strcpy(Branch->branchname,branchname);
  lenbrancharg = 0;
  for(n=0;n<nleafs;n++) {
    char *nameptr;
    /* Need to look for $name here.  name will be the name given to root */
    args[n] = thSpaceStrip(args[n]);
    /*printf("Leaf %s\n",args[n]);*/
    if(nameptr=strchr(args[n],'$')) *nameptr++=0;
    if(thVarResolve(args[n],&varp,&vind,0,0) == S_SUCCESS) {
      char *p, snum[25];
      /*printf("Index=%d\n",vind);*/
      thisleaf = *LeafNext = (thLeafList *) malloc(sizeof(thLeafList));
      /*printf("thisleaf = %x\n",thisleaf);*/
      thisleaf->next = (thLeafList *) NULL;
      thisleaf->varp = varp;
      thisleaf->index = vind;
      /* Pick a good name */
      if(nameptr) {
	thisleaf->name = (char *) malloc(strlen(nameptr)+1);
	strcpy(thisleaf->name,nameptr);
      } else {
	if(p=strpbrk(args[n],"()[]")) {
	  sprintf(snum,"%d",vind+1);
	  thisleaf->name = (char *) malloc(strlen(args[n])+strlen(snum)+2);
	  strncpy(thisleaf->name,args[n],p-args[n]);
	  thisleaf->name[p-args[n]] = '\0';
	  strcat(thisleaf->name,"_");
	  strcat(thisleaf->name,snum);
	} else {
	  thisleaf->name = (char *) malloc(strlen(args[n])+1);
	  strcpy(thisleaf->name,args[n]);
	}
      }
      LeafNext = &((*LeafNext)->next);
      lenbrancharg += strlen(args[n]) + 3;
    }	else {
      fprintf(stderr,"Bad variable %s\n",args[n]);
    }
  }
  /* Walk down the leaf list and build the Branch call argument */
  /* What do I do about leaf names with subscripts? */
  thisleaf = Branch->leaflistp;
  brancharg = malloc(lenbrancharg+10);
  brancharg[0] = '\0';
  while(thisleaf) {
    /*printf("thisleaf = %x\n",thisleaf);
      printf("Adding %s to branchlist\n",thisleaf->name);*/
    strcat(brancharg,thisleaf->name);
    if(thisleaf->varp->type == DAVARINT) {
      strcat(brancharg,"/I");
      thisleaf->leafsize=4;
    } else if(thisleaf->varp->type == DAVARFLOAT) {
      strcat(brancharg,"/F");
      thisleaf->leafsize=4;
    } else if(thisleaf->varp->type == DAVARDOUBLE) {
      strcat(brancharg,"/D");
      thisleaf->leafsize=8;
    } else {
      fprintf(stderr,"Variable %s has unknown type\n");
    }

    if(thisleaf->next) {
      strcat(brancharg,":");
    }
    thisleaf = thisleaf->next;
  }

  /* Reserve enough space as if they were all double */
  Branch->evstruct = (void *) malloc(lenbrancharg*sizeof(double));

  /*
       * leaflist is the concatenation of all the variable names and types
         separated by a colon character :
         The variable name and the variable type are separated by a slash (/).
         The variable type may be 0,1 or 2 characters. If no type is given,
         the type of the variable is assumed to be the same as the previous
         variable. If the first variable does not have a type, it is assumed
         of type F by default. The list of currently supported types is given below:
            - C : a character string terminated by the 0 character
            - B : an 8 bit signed integer (Char_t)
            - b : an 8 bit unsigned integer (UChar_t)
            - S : a 16 bit signed integer (Short_t)
            - s : a 16 bit unsigned integer (UShort_t)
            - I : a 32 bit signed integer (Int_t)
            - i : a 32 bit unsigned integer (UInt_t)
            - F : a 32 bit floating point (Float_t)
            - D : a 64 bit floating point (Double_t)
  */

  printf("Branch=%s  Leafs=%s\n",Branch->branchname,brancharg);
  thRoot_Branch(treedef->treeptr,Branch->branchname,(Branch->evstruct),brancharg);

  free(brancharg);
  printf("Exiting book a branch\n");
  return(S_SUCCESS);
}
  
thStatus thFillTreeV(daVarStruct *var){
  thTreeOpaque *treedef;
  thTreeBranchList *thisbranch;
  thLeafList *thisleaf;
  void *structp;
  /* printf("Executing Tree %s\n",var->name);*/
  treedef = ((thTreeOpaque *)(var->opaque));
  thisbranch = treedef->branchlistP;
  if(! (treedef->test ? *((DAINT *) treedef->test->varptr
			  + treedef->testindex) : 1)) {
    return(S_SUCCESS);  /* Test was false */
  }
  while(thisbranch) {
    structp = thisbranch->evstruct;
    /*    printf("Filling branch %s at %x\n",thisbranch->branchname,structp);*/
    thisleaf = thisbranch->leaflistp;
    while(thisleaf) {
      if(thisleaf->varp->type == DAVARINT) {
	*(DAINT *)(structp) = *((DAINT *)thisleaf->varp->varptr + thisleaf->index);/*phil*/
        structp = (void *) (DAINT *) ((DAINT *)structp + 1);
	/*	printf("   %s=%d\n",thisleaf->name,*((DAINT *)thisleaf->varp->varptr
		+ thisleaf->index));*/
      } else if(thisleaf->varp->type == DAVARFLOAT) {
	*(DAFLOAT *)(structp) = *((DAFLOAT *)thisleaf->varp->varptr + thisleaf->index);/*phil*/
        structp = (void *) (DAFLOAT *) ((DAFLOAT *)structp + 1);
	/*	printf("   %s=%f\n",thisleaf->name,*((DAFLOAT *)thisleaf->varp->varptr
		+ thisleaf->index));*/
      } else if(thisleaf->varp->type == DAVARDOUBLE) {
	*(DADOUBLE *)(structp) = *((DADOUBLE *)thisleaf->varp->varptr + thisleaf->index);/*phil*/
        structp = (void *) (DADOUBLE *) ((DADOUBLE *)structp + 1);
	/*	printf("   %s=%lf\n",thisleaf->name,*((DADOUBLE *)thisleaf->varp->varptr
		+ thisleaf->index));*/
      }
      thisleaf = thisleaf->next;
    }
    thisbranch = thisbranch->next;
  }
  thRoot_Fill(treedef->treeptr);

  (*((DAINT *)var->varptr))++; /* Increment block counter */
  return(S_SUCCESS);
}
thStatus thClearTreeV(daVarStruct *var){
  
  /* printf("Clearing Tree %s\n",var->name); */

  (*((DAINT *)var->varptr)) = 0; /* Increment block counter */
  return(S_SUCCESS);
}
thStatus thWriteTreeV(daVarStruct *var){
  
  /* printf("Writing Tree %s\n",var->name); */
  thRoot_Write(((thTreeOpaque *)(var->opaque))->file);

  (*((DAINT *)var->varptr)) = 0; /* Increment block counter */
  return(S_SUCCESS);
}
thStatus thCloseTreeV(daVarStruct *var){
  
  /*printf("Closing Tree %s\n",var->name);*/
  thRoot_Close(((thTreeOpaque *)(var->opaque))->file);

  (*((DAINT *)var->varptr)) = 0; /* Increment block counter */
  return(S_SUCCESS);
}

thStatus thRemoveTree(char *treename) {
  printf("Dummy routine to remove tree %s\n",treename);
  return(S_SUCCESS);
}

/*
int thtreewrite_()
{
  thRoot_Write();
}
*/

#else

thStatus thBookTree(daVarStruct *var) {
  return(S_SUCCESS);
}

thStatus thFillTreeV(daVarStruct *var) {
  return(S_SUCCESS);
}
thStatus thClearTreeV(daVarStruct *var) {
  return(S_SUCCESS);
}
thStatus thWriteTreeV(daVarStruct *var) {
  return(S_SUCCESS);
}
thStatus thCloseTreeV(daVarStruct *var) {
  return(S_SUCCESS);
}

#endif
