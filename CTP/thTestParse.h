/*
 * Revision History:
 * $Log: thTestParse.h,v $
 * Revision 1.1  2009/01/23 13:34:01  gaskelld
 * Initial revision
 *
 * Revision 1.1  1998/12/07 22:11:13  saw
 * Initial setup
 *
 *	  Revision 1.6  1995/08/03  14:37:13  saw
 *	  Add single argument functions
 *
 *	  Revision 1.5  1995/04/25  17:46:52  saw
 *	  Make compatible with OSF/Alpha (64 bit pointers)
 *
 *	  Revision 1.4  1994/11/07  14:38:43  saw
 *	  Add integer divide operator
 *
 *	  Revision 1.3  1994/07/21  20:52:33  saw
 *	  Add Revision history
 *
 */
#ifndef _TH_TESTPARSE_H
#define _TH_TESTPARSE_H

#ifndef _TH_UTILS_H
#include "thUtils.h"
#endif

/* thTestParse.h
 Header file for thTestParse.c and thTestExecute.c
 */

#define max(a,b) 		(a<b ? b : a)
#define min(a,b) 		(a>b ? b : a)

#ifndef NULL
#define NULL ((void *)0)
#endif

#define CODEGROWSIZE 250
#define CODESTARTSIZE 2500
typedef DAINT CODE;
typedef CODE * CODEPTR;



/* Operand types */

/* Hex charcter meanings
   H7 H6 H5 H4 H3 H2 H1 H0
   H7       Groups instructions together for execution convenience
   H7 H6 H5 Operator precidence ordering
   H4 H3    Operator code within precidence group
   H2       Type of Left operand
   H1       Type of Right operand
   H0       Type of the result

*/
/* Operand types */
#define OPRINT 0
#define OPRFLOAT 1
#define OPRDOUBLE 2
/* Operators */
#define OPPUSHGROUP  0x0000000
#define OPPUSHINT    0x0800000
#define OPPUSHFLOAT  0x0800001
#define OPPUSHDOUBLE 0x0800002
/* Next word is a pointer, push the value pointed to onto rpn stack */
#define OPPUSHINTP   0x0810000
#define OPPUSHFLOATP 0x0810001
#define OPPUSHDOUBLEP 0x0810002
/* Next word is a pointer, push the pointer onto the rpn stack */
#define OPPUSHPINT   0x0820000
#define OPPUSHPFLOAT 0x0820001
#define OPPUSHPDOUBLE 0x0820002
/* Next word is a index into a function table */
#define OPPUSHFUNCTION 0x0830000
/* Parenthesis operators */
#define OPLP      0x0100000
#define OPRP      0x0900000
#define OPRINDEX  0x0901000
#define OPRFARG   0x0B01000
#define OPLINDEXGROUP 0x1000000
#define OPLINDEX  0x1A00000
/* Fortran mode index */
#define OPLINDEXB 0x1A01000
/* These operators leave the pointer instead of the value on the stack */
/* They are used for the first indexing ( or [ before the = */
#define OPLINDEXP  0x1A10000
/* Fortran mode index */
#define OPLINDEXPB 0x1A11000
/* Function operator.  Is this the right place in precidence scheme? */
#define OPLFARG   0x1B00000
#define OPEOLGROUP 0x2000000
#define OPEOL     0x2E00000
#define OPCOMMA   0x2F00000

#define OPEQUAL   0x4000000

#define OPLOGGROUP 0x8000000
#define OPLOGOR   0x8300000
#define OPLOGXOR  0x8400000
#define OPLOGAND  0x8500000
#define OPBITOR   0x8600000
#define OPBITXOR  0x8700000
#define OPBITAND  0x8800000

#define OPCOMPGROUP 0xA000000
#define OPISEQUAL   0xA900000
#define OPISNOTEQUAL 0xA901000
#define OPISLT      0xAA00000
#define OPISLE      0xAA01000
#define OPISGT      0xAA02000
#define OPISGE      0xAA03000

#define OPSHIFTGROUP 0xB000000
#define OPSHL     0xBB00000
#define OPSHR     0xBB01000

#define OPADDGROUP 0xC000000
#define OPADD      0xC000000
#define OPSUB      0xC001000

#define OPMULGROUP 0xD000000
#define OPTIMES    0xD000000
#define OPDIV      0xD001000
#define OPIDIV     0xD002000
#define OPMOD      0xD003000

#define OPUNARY   0xE000000
#define OPNEG     0xE000000
#define OPNOT     0xE001000
#define OPCOMP    0xE002000

#define OPGROUPMASK 0xF000000
#define OPPRECMASK 0xFF00000
#define OPCODEMASK 0xFFFF000

#define OPLEFTTYPEMASK 0x0000F00
#define OPRIGHTTYPEMASK 0x00000F0
#define OPLRTYPEMASK 0x0000FF0
#define OPRESTYPEMASK  0x000000F

/* For Q like test package format */
typedef enum {
  tGATE, tPAT, tEQ, tBIT, tAND, tIOR, tXOR, tMAJ, tUSER, tBAD,
  PUSHI, PUSHITOFS, PUSHS, PUSHFTOIS, POPS
  } thTestType;

/* Operand types deduced from context */
typedef enum {
  otIMMED, otLOGIC, otVALUE, otRESULT} thOperandType;

typedef struct
{
  char *name;
  CODE result[3];
} INTRINSIC_FUNCTIONS;

char *thGetTok(char *linep,int *tokenid, char **tokstr, CODE *tokval,
	       void **tokptr, int expflag, daVarStructList **vlisthead);
thStatus thBookaTest(char *line, CODEPTR *codeheadp, CODEPTR *codenextp,
		     CODEPTR *codelimitp, CODEPTR *codelastop,
		     daVarStructList **vlisthead);
thStatus thExecuteCode(char *blockname, CODEPTR code, CODEPTR codelimit);
thOperandType thGetOperandType(char *soperand, char *rest, CODE lastop,
			       int expflag);
char **thGetClassList(thOperandType optype);

#endif
