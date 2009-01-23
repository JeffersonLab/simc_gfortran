//
// thRootStuff.cpp
//    C++ wrapper routines to interface between CTP and Root libraries
//
// $Log: thRootStuff.cpp,v $
// Revision 1.1  2009/01/23 13:34:01  gaskelld
// Initial revision
//
// Revision 1.3  2005/02/22 16:25:51  saw
// Make sure next pointer is zeroed in rootfilelist
//
// Revision 1.2  2004/07/07 18:16:55  saw
// use extern "C" to export names needed in thTree.c
//
// Revision 1.1  2002/07/31 20:07:48  saw
// Add files for ROOT Trees
//
// Revision 1.1  1999/08/25 13:16:07  saw
// *** empty log message ***
//
#include <stdio.h>
#include <stdlib.h>

#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TRandom.h>
#include <TTree.h>

TROOT CTP("CTP","CTP Histograms and trees");

TTree *tree;
TFile *hfile;

struct thRootFileList {
  char *filename;
  TFile *tfile;
  int count;
  struct thRootFileList *next;
};
typedef struct thRootFileList thRootFileList;

thRootFileList *rootfilelistp=0;
  
extern "C" void *thRoot_TFile(char *filename);

extern "C" void *thRoot_TTree(char *treename);

extern "C" void thRoot_Branch(TTree *tree, char *branchname, void *structp, char *brancharg);

extern "C" void thRoot_Fill(TTree *tree);

extern "C" void thRoot_Write(thRootFileList *file);

extern "C" void thRoot_Close(thRootFileList *file);

void *thRoot_TFile(char *filename)
{
  thRootFileList *thisfile,**lastp;
  thisfile = rootfilelistp;
  lastp = &rootfilelistp;
  while(thisfile) {
    if(strcmp(thisfile->filename,filename)==0) {
      thisfile->count++;
      return((void *) thisfile);
    }
    lastp = &(thisfile->next);
    thisfile = thisfile->next;
  }
  /* Need to check if file has been opened */
  printf("Tfile(\"%s\",\"RECREATE\",\"CTP ROOT file with trees\")\n",filename);
  *lastp = (thRootFileList *) malloc(sizeof(thRootFileList));
  thisfile = *lastp;
  thisfile->tfile = new TFile(filename,"RECREATE","CTP ROOT file with trees");
  thisfile->count = 1;
  thisfile->filename = (char *)malloc(strlen(filename)+1);
  thisfile->next = (thRootFileList *) 0;
  strcpy(thisfile->filename,filename);
  return((void *) thisfile);
}

void *thRoot_TTree(char *treename)
{
  TTree *tree;

  /* Perhaps Check if a tree exists by this name?? */
  printf("new TTree(\"%s\",\"CTP ROOT tree\")\n",treename);
  tree = new TTree(treename,"CTP ROOT tree");
  return((void *)tree);
}

void thRoot_Branch(TTree *tree, char *branchname, void *structp, char *brancharg)
{
  tree->Branch(branchname,structp,brancharg);
}

void thRoot_Fill(TTree *tree)
{
  tree->Fill();
}

void thRoot_Write(thRootFileList *file)
{
  (file->tfile)->Flush();
}

void thRoot_Close(thRootFileList *file)
{
  TFile *hfile;
  if(--file->count <= 0){
    printf("Closing\n");
    hfile = file->tfile;
    hfile->Write();
    hfile->Close();
  } else {
    printf("Not Closing\n");
  }
}

