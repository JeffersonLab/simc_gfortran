/*-----------------------------------------------------------------------------
 * Copyright (c) 1993-1995 Southeastern Universities Research Association,
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
 *  Test code executor
 *	
 * Author:  Stephen Wood, CEBAF Hall C
 *
 * Revision History:
 *   $Log: thTestExecute.c,v $
 *   Revision 1.2  2011/03/04 20:01:44  jones
 *   Add check for 64bit by looking for LP64
 *
 *   Revision 1.2.24.1  2007/09/10 21:32:47  pcarter
 *   Implemented changes to allow compilation on RHEL 3,4,5 and MacOSX
 *
 *   Revision 1.2  2003/02/21 20:55:25  saw
 *   Clean up some types and casts to reduce compiler warnings.
 *
 *   Revision 1.1  1998/12/07 22:11:13  saw
 *   Initial setup
 *
 *	  Revision 1.13  1996/08/01  01:33:31  saw
 *	  Add trig functions.  Print block name on errors.  Allow floating arguments
 *	  for mod (%).
 *
 *	  Revision 1.12  1995/08/03  13:56:00  saw
 *	  Add single argument functions
 *
 *	  Revision 1.11  1995/02/14  16:53:12  saw
 *	  Make compatible with OSF/Alpha (64 bit pointers)
 *
 *	  Revision 1.10  1994/10/03  12:39:45  saw
 *	  All "/" (division) has real result.  New op "//" has integerized result.
 *
 *	  Revision 1.9  1994/07/21  20:46:50  saw
 *	  Make some POP and PUSH stuff more portable
 *
 *	  Revision 1.8  1994/06/13  13:26:55  saw
 *	  Add some divide by zero checking
 *
 *	  Revision 1.7  1994/06/03  21:12:19  saw
 *	  Use memcpy for stack manipulation
 *
 *	  Revision 1.5  1994/02/10  18:40:25  saw
 *	  Typecasting of sp and pc pointer while incrementing or decrementing
 *	  doesn't work for ANSI C (Silicon Graphics).  Rewrote some POP and PUSH
 *	  macros and other codes for this (presumably) non-posix case.
 *
 *	  Revision 1.4  1993/12/02  21:33:17  saw
 *	  Fully allow the use of doubles in test expressions
 *
 *	  Revision 1.3  1993/11/22  20:42:12  saw
 *	  Add return of status codes on thExecuteCode
 *
 *	  Revision 1.2  1993/05/11  17:58:13  saw
 *	  Add copyright and header
 *
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "daVar.h"
#include "th.h"
#include "thUtils.h"
#include "thTestParse.h"
#include "thInternal.h"
#include "cfortran.h"


static DAINT stack[1000];			/* The stack */

/* Execute Test Package Pseudo Code
   code is first instruction, codelimit is pointer to location
   after last instruction.
   */

#ifdef __sgi
#define NOTPOSIX
#endif

#if defined(hpux) || defined(__sun)
#define USEMEMCPY
#endif

#ifdef __DECC
#define USEMEMCPY
#endif

#if (defined(__osf__) && defined(__alpha)) || defined(__LP64__)
#undef USEMEMCPY
#define NOTPOSIX
#define POINTER64
#endif

#ifdef USEMEMCPY
#define PUSHPOINTER(x) memcpy(((void **)sp)++, (void *)&x, sizeof(void *))
#else
//#define PUSHPOINTER(x) *((DAINT **)sp)++ = x
//#define PUSHPOINTER(x) *(DAINT **)sp = x; (DAINT **)(sp = (DAINT *)(DAINT **)((DAINT **)sp + 1))
#define PUSHPOINTER(x) *(DAINT **)sp = x; sp = (DAINT *) (DAINT **) ((DAINT **)sp + 1)
//#define PUSHPOINTER(x) *(DAINT **)sp = x; sp = (DAINT *) ((DAINT **)sp + 1)
#endif

#ifdef NOTPOSIX

/* Can't use --() constructs, must decrement stack pointer manually */ 

# ifdef POINTER64

#  define SAVEINT(x) **((DAINT **)(sp-2)) = x; sp--; sp--
#  define SAVEFLOAT(x) **((DAFLOAT **)(sp-2)) = x; sp--; sp--
#  define SAVEDOUBLE(x) **((DADOUBLE **)(sp-2)) = x; sp--; sp--
#  define FETCHIARRAY(x) x = *(*(((DAINT **)sp)-1) + index); sp--; sp--
#  define FETCHFARRAY(x) x = *(*(((DAFLOAT **)sp)-1) + index); sp--; sp--
#  define FETCHDARRAY(x) x = *(*(((DADOUBLE **)sp)-1) + index); sp--; sp--

# else /* 32 bit pointers */

#  define SAVEINT(x) **((DAINT **)(sp-1)) = x; sp--
#  define SAVEFLOAT(x) **((DAFLOAT **)(sp-1)) = x; sp--
#  define SAVEDOUBLE(x) **((DADOUBLE **)(sp-1)) = x; sp--
#  define FETCHIARRAY(x) x = *((DAINT *) *(sp-1) + index); sp--
#  define FETCHFARRAY(x) x = *((DAFLOAT *) *(sp-1) + index); sp--
#  define FETCHDARRAY(x) x = *((DADOUBLE *) *(sp-1) + index); sp--

# endif

#define POPDOUBLE(x) x = *((DADOUBLE *)(sp-2)); sp--; sp--
#define POPFLOAT(x) x = *((DAFLOAT *)(sp-1)); sp--
#define POPINT(x) x = *((DAINT *)(sp-1)); sp--
#define PUSHDOUBLE(x) *((DADOUBLE *)sp) = x; sp++; sp++
#define PUSHFLOAT(x) *((DAFLOAT *)sp++) = x;
#define PUSHINT(x) *((DAINT *)sp++) = x;

#else

//# define SAVEINT(x) **(--(DAINT **)sp) = x
//# define SAVEFLOAT(x) **(--(DAFLOAT **)sp) = x
//# define SAVEDOUBLE(x) **(--(DADOUBLE **)sp) = x
# define SAVEINT(x) sp = (DAINT *) (DAINT **) ((DAINT **)sp - 1); **(DAINT **)sp = x
# define SAVEFLOAT(x) sp = (DAINT *) (DAFLOAT **) ((DAFLOAT **)sp - 1); **(DAFLOAT **)sp = x
# define SAVEDOUBLE(x) sp = (DAINT *) (DADOUBLE **) ((DADOUBLE **)sp - 1); **(DADOUBLE **)sp = x
//# define FETCHIARRAY(x) x = (*(*(--(DAINT**)sp) + index));
//# define FETCHFARRAY(x) x = (*(*(--(DAFLOAT**)sp) + index));
//# define FETCHDARRAY(x) x = (*(*(--(DADOUBLE**)sp) + index));
# define FETCHIARRAY(x) sp = (DAINT *) (DAINT **) ((DAINT **)sp - 1); x = *(*(DAINT**)sp + index);
# define FETCHFARRAY(x) sp = (DAINT *) (DAFLOAT **) ((DAFLOAT **)sp - 1); x = *(*(DAFLOAT**)sp + index);
# define FETCHDARRAY(x) sp = (DAINT *) (DADOUBLE **) ((DADOUBLE **)sp - 1); x = *(*(DADOUBLE**)sp + index);

#ifdef USEMEMCPY

#define POPDOUBLE(x) ((DADOUBLE *)sp)--; memcpy((void *)&x,(void *)sp,sizeof(DADOUBLE))
#define POPFLOAT(x) x = *(--(DAFLOAT *)sp)
#define POPINT(x) x = *(--(DAINT *)sp)
#define PUSHDOUBLE(x) memcpy(((DADOUBLE *)sp)++,(void *)&x,sizeof(DADOUBLE));
#define PUSHFLOAT(x) *((DAFLOAT *)sp)++ = x;
#define PUSHINT(x) *((DAINT *)sp)++ = x;

#else

//#define POPDOUBLE(x) x = *(--(DADOUBLE *)sp)
#define POPDOUBLE(x) sp = (DAINT *) (DADOUBLE *) ((DADOUBLE *)sp - 1); x = *(DADOUBLE *)sp
//#define POPFLOAT(x) x = *(--(DAFLOAT *)sp)
#define POPFLOAT(x) sp = (DAINT *) (DAFLOAT *) ((DAFLOAT *)sp - 1); x = *(DAFLOAT *)sp
#define POPINT(x) x = *(--(DAINT *)sp)
//#define PUSHDOUBLE(x) *((DADOUBLE *)sp)++ = x;
#define PUSHDOUBLE(x) *(DADOUBLE *)sp = x; sp = (DAINT *) (DADOUBLE *) ((DADOUBLE *)sp + 1)
//#define PUSHFLOAT(x) *((DAFLOAT *)sp)++ = x
#define PUSHFLOAT(x) *(DAFLOAT *)sp = x; sp = (DAINT *) (DAFLOAT *) ((DAFLOAT *)sp + 1)
#define PUSHINT(x) *((DAINT *)sp)++ = x
#endif
#endif

#define GETNEXTINTP (((DAINT *)pc)++)
#define GETNEXTFLOATP (((DAFLOAT *)pc)++)
#define GETNEXTDOUBLEP (((DADOUBLE *)pc)++)
#define GETNEXTPOINTERP (((DAINT **)pc)++)

thStatus thExecuteCode(char *blockname,CODEPTR code, CODEPTR codelimit)
{
#ifdef PHILDEBUG
#ifdef NOTPOSIX
#warning Phil says NOTPOSIX!
#else
#warning Phil says not NOTPOSIX! i.e. POSIX!!
#endif
#ifdef POINTER64
#warning Phil says POINTER64!
#else
#warning Phil says not POINTER64!
#endif
#ifdef USEMEMCPY
#warning Phil says USEMEMCPY!
#else
#warning Phil says not USEMEMCPY!
#endif
#endif
  register CODEPTR pc;
  CODE rawopcode,opcode,ltype,rtype,lrtypes;
  DAINT nargs,result;
  register DAINT *sp;
  DAINT i,il,ir,*pi;
  DAFLOAT f,fl,fr,*pf;
  DADOUBLE d,dl,dr,*pd;
  DAINT index;

  sp = stack;
  pc = code;

  while(pc < codelimit){
/*    printf("PC=%x, Op code %x, Stack=%x, SP=%x\n",pc,*pc,stack,sp);*/
    rawopcode = *pc++;
    if(rawopcode >= OPLP){		/* New style */
      ltype = (rawopcode & OPLEFTTYPEMASK) >> 8;
      rtype = (rawopcode & OPRIGHTTYPEMASK) >> 4;
/*      lrtypes = opcode & OPLRTYPEMASK;*/
      opcode = rawopcode & OPCODEMASK;
      switch(opcode & OPGROUPMASK)
	{
	case OPPUSHGROUP:		/* Pushes */
	  switch(opcode)
	    {
#ifdef USEMEMCPY
	      void *tmpptr;
#endif
	    case OPPUSHINT:	/* Float included in pushes */
	      if((rawopcode & OPRESTYPEMASK) == OPRDOUBLE){
/*		printf("sp=%x, pc=%x\n",sp,pc);*/
#ifdef USEMEMCPY
		memcpy((void *)&d,((DADOUBLE *)pc)++,sizeof(DADOUBLE));
		PUSHDOUBLE(d);
#else
#ifdef __sgi
		PUSHDOUBLE(*((DADOUBLE *)pc)); pc++; pc++;
#else
		PUSHDOUBLE(*(DADOUBLE *)pc);/*phil*/
                pc = (CODEPTR) (DADOUBLE *) ((DADOUBLE *)pc + 1);
#endif
#endif
/*		printf("sp=%x, pc=%x\n",sp,pc);*/
	      } else {
		PUSHINT(*pc++);
	      }
	      break;
	    case OPPUSHPINT:	/*Push a pointer*/
#ifdef USEMEMCPY
	      PUSHPOINTER((memcpy(&tmpptr,(((DAINT **)pc)++),sizeof(void *))
			  ,tmpptr));
#else
              PUSHPOINTER(*(DAINT **)pc); /*phil*/
              pc = (CODEPTR)(DAINT **) ((DAINT **)pc + 1);
#endif
	      break;
	    case OPPUSHINTP:    /*Push what a pointer points to */
	      if((rawopcode & OPRESTYPEMASK) == OPRDOUBLE){
#ifdef USEMEMCPY
		memcpy(&tmpptr,(((DAINT **)pc)++),sizeof(void *));
		d = *(DADOUBLE *) tmpptr;
#else
		d = **(DADOUBLE **)pc;/*phil*/
                pc = (CODEPTR) (DADOUBLE **) ((DADOUBLE **)pc + 1);
#endif	      
		PUSHDOUBLE(d);/*phil*/
	      } else {
#ifdef USEMEMCPY
		memcpy(&tmpptr,(((DAINT **)pc)++),sizeof(void *));
		PUSHINT(*(DAINT *) tmpptr);
#else		
		PUSHINT(**(DAINT **)pc);/*phil*/
                pc = (CODEPTR) (DAINT **) ((DAINT **)pc + 1);
#endif
	      }
	      break;
	    case OPPUSHFUNCTION:	/*Push a intrinsic function code */
	      PUSHINT(*pc++);
	      break;
	    }
	  break;
	case OPEOLGROUP:
	  sp--;		/* Should empty the stack */
	  if(rtype == OPRDOUBLE) sp--; /* Double is two entries on stack */
	  break;
	case OPLINDEXGROUP:
	  if(opcode==OPLFARG) {
	    if(rtype==OPRINT) {POPINT(i);}
	    else if(rtype==OPRFLOAT) {POPFLOAT(f);}/*phil*/
	    else {POPDOUBLE(d);}/*phil*/
	    POPINT(index);	/* Pop the function code */
	    switch(index)
	      {
	      case 0:		/* abs */
		if(rtype==OPRINT) {
		  if(i<0) i = -i;
		  PUSHINT(i);
		} else if(rtype==OPRFLOAT) {
		  if(f<0.0) f = -f;
		  PUSHFLOAT(f);/*phil*/
		} else {
		  if(d<0.0) d = -d;
		  PUSHDOUBLE(d);/*phil*/
		}
		break;
	      case 1:		/* sqrt */
		if(rtype==OPRINT) d = i;
		else if(rtype==OPRFLOAT) d = f;
		if(d>=0) d = sqrt(d);
		else {
		  fprintf(STDERR,"Test block %s: sqrt(%f)\n",blockname,d);
		  d = 0;
		}
		PUSHDOUBLE(d);/*phil*/
		break;
	      case 2:		/* exp */
		if(rtype==OPRINT) d = i;
		else if(rtype==OPRFLOAT) d = f;
		d = exp(d);
		PUSHDOUBLE(d);/*phil*/
		break;
	      case 3:		/* sin */
		if(rtype==OPRINT) d = i;
		else if(rtype==OPRFLOAT) d = f;
		d = sin(d);
		PUSHDOUBLE(d);/*phil*/
		break;
	      case 4:		/* cos */
		if(rtype==OPRINT) d = i;
		else if(rtype==OPRFLOAT) d = f;
		d = cos(d);
		PUSHDOUBLE(d);/*phil*/
		break;
	      case 5:		/* tan */
		if(rtype==OPRINT) d = i;
		else if(rtype==OPRFLOAT) d = f;
		d = tan(d);
		PUSHDOUBLE(d);/*phil*/
		break;
	      }
	    break;
	  }
	  if(rtype==OPRFLOAT) {	/* Floating point index */
	    POPFLOAT(f);/*phil*/
	    index = floatToLong(f);
	  } else if(rtype==OPRDOUBLE) {	/* Double */
	    POPDOUBLE(d);/*phil*/
	    index = floatToLong(d);
	  } else {
	    POPINT(index);
	  }
	  index -= ((opcode & 0xF000) == 0x1000 ? 0 : 1);
	  /* ltype should always be == restype */
	  if(opcode == OPLINDEX || opcode == OPLINDEXB){
	    if(ltype == OPRDOUBLE) {
	      FETCHDARRAY(d);/*phil*/
	      PUSHDOUBLE(d);/*phil*/
	    } else if (ltype == OPRFLOAT) {
              FETCHFARRAY(f);/*phil*/
              PUSHFLOAT(f);/*phil*/
	    } else {
              FETCHIARRAY(i);/*phil*/
	      PUSHINT(i);
	    }
	  } else { /*pointer on stack*/
	    sp--;
#ifdef POINTER64
	    sp--;
#endif
	    if(ltype == OPRDOUBLE) {
	      /*	      *((DADOUBLE **)sp)++ =  (*((DADOUBLE **)sp)+index);*/
	      /* The following works better on the alpha */
	      pd = *((DADOUBLE **)sp);
	      pd += index;
	      PUSHPOINTER(pd);/*phil*/
	    } else {		/* Assume INT and FLOAT the same size */
	      /**((DAINT **)sp)++ =  (*((DAINT **)sp)+index);*/
	      /* The following works better on the alpha */
	      pi = *((DAINT **)sp);
	      pi += index;
	      PUSHPOINTER(pi);/*phil*/
	    }
	  }
	  break;
	case OPEQUAL:		/* Big ugly matrix of type conversions */
	  if(rtype==OPRINT) {
	    POPINT(i);
	    if(ltype==OPRINT) {
	      SAVEINT(i); /* Save result in result variable *//*phil*/
	      PUSHINT(i);	/* Put result back on stack */
	    } else if(ltype==OPRFLOAT) {
	      f = i;	/* Convert to floating */
	      SAVEFLOAT(f); /* Save variable *//*phil*/
	      PUSHFLOAT(f); /* Put back on stack *//*phil*/
	    } else {		/* if(ltype==OPRDOUBLE) */
	      d = i;
	      SAVEDOUBLE(d);/*phil*/
	      PUSHDOUBLE(d);/*phil*/
	    }
	  } else if(rtype==OPRFLOAT) {
	    POPFLOAT(f);/*phil*/
	    if(ltype==OPRINT) {
	      i = floatToLong(f);
	      SAVEINT(i); /* Save result in result variable *//*phil*/
	      *sp++ = i;
	    } else if(ltype==OPRFLOAT) {
	      SAVEFLOAT(f); /* Save variable *//*phil*/
	      *sp++ = *(DAINT *)&f;
	    } else {		/* if(ltype==OPRDOUBLE) */
	      d = f;
	      SAVEDOUBLE(d);/*phil*/
	      PUSHDOUBLE(d);/*phil*/
	    }
	  } else {		/* if(rtype==OPRDOUBLE) */
	    POPDOUBLE(d);/*phil*/
	    if(ltype==OPRINT) {
	      i = floatToLong(d);
	      SAVEINT(i); /* Save result in result variable *//*phil*/
	      *sp++ = i;
	    } else if(ltype==OPRFLOAT) {
	      f = d;
	      SAVEFLOAT(f); /* Save variable *//*phil*/
	      *sp++ = *(DAINT *)&f;
	    } else {		/* if(ltype==OPRDOUBLE) */
 	      SAVEDOUBLE(d);/*phil*/
	      PUSHDOUBLE(d);/*phil*/
	    }
	  }
	  break;
	case OPLOGGROUP:		/* Logic and Bit operations */
	case OPSHIFTGROUP:		/* Logic and Bit operations */
	  if(rtype==OPRINT) {
	    POPINT(ir);
	  } else if(rtype==OPRFLOAT) {
	    POPFLOAT(f);/*phil*/
	    ir = floatToLong(f);
	  } else {
	    POPDOUBLE(d);/*phil*/
	    ir = floatToLong(d);
	  }
	  if(ltype==OPRINT) {
	    POPINT(il);
	  } else if(ltype==OPRFLOAT) {
	    POPFLOAT(f);/*phil*/
	    il = floatToLong(f);
	  } else {
	    POPDOUBLE(d);/*phil*/
	    il = floatToLong(d);
	  }
	  switch(opcode)
	    {
	    case OPLOGOR:
	      *sp++ = il || ir;
	      break;
	    case OPLOGXOR:
	      *sp++ = (il != 0) ^ (ir != 0);
	      break;
	    case OPLOGAND:
	      *sp++ = il && ir;
	      break;
	    case OPBITOR:
	      *sp++ = il | ir;
	      break;
	    case OPBITXOR:
	      *sp++ = il ^ ir;
	      break;
	    case OPBITAND:
	      *sp++ = il & ir;
	      break;
	    case OPSHL:
	      *sp++ = il << ir;
	      break;
	    case OPSHR:
	      *sp++ = il >> ir;
	      break;
	    }
	  break;
	case OPCOMPGROUP:		/* Logic comparisons */
/* Result of Add amd MUL groups should now always be double */
	case OPADDGROUP:	/* Add and Subtract */
	case OPMULGROUP:	/* * / and % */
	  if(rtype==OPRINT) {
	    POPINT(ir);
	    dr = ir;
	  } else if (rtype==OPRFLOAT) {
	    POPFLOAT(fr);/*phil*/
	    dr = fr;
	  } else {
	    POPDOUBLE(dr);/*phil*/
	  }
	  if(ltype==OPRINT) {
	    POPINT(il);
	    dl = il;
	  } else if (ltype==OPRFLOAT) {
	    POPFLOAT(fl);/*phil*/
	    dl = fl;
	  } else {
	    POPDOUBLE(dl);/*phil*/
	  }
	  if(rtype!=OPRINT || ltype!=OPRINT){
	    switch(opcode)
	      {
	      case OPISEQUAL:
		*sp++ = dl == dr;
		break;
	      case OPISNOTEQUAL:
		*sp++ = dl != dr;
		break;
	      case OPISLT:
		*sp++ = dl < dr;
		break;
	      case OPISGT:
		*sp++ = dl > dr;
		break;
	      case OPISLE:
		*sp++ = dl <= dr;
		break;
	      case OPISGE:
		*sp++ = dl >= dr;
		break;
	      case OPADD:
		d = dl + dr;
		PUSHDOUBLE(d);/*phil*/
		break;
	      case OPSUB:
		d = dl - dr;
		PUSHDOUBLE(d);/*phil*/
		break;
	      case OPTIMES:
		d = dl * dr;	/* Need to deal with overflow */
		PUSHDOUBLE(d);/*phil*/
		break;
	      case OPIDIV:
/*		printf("OP=%x\n",rawopcode);*/
		if(dr == 0.0) {
		  fprintf(STDERR,"Test block %s: %f/0.0\n",blockname,dl);
		  d = 0.0;
		} else {
		  d = dl / dr;	/* Need to deal with overflow and div 0 */
		}
		*sp++ = floatToLong(d);
		break;
	      case OPDIV:
		if(dr == 0.0) {
		  fprintf(STDERR,"Test block %s: %f/0.0\n",blockname,dl);
		  d = 0.0;
		} else {
		  d = dl / dr;	/* Need to deal with overflow and div 0 */
		}
		PUSHDOUBLE(d);/*phil*/
		break;
	      case OPMOD:
		d = fmod(dl,dr);
		PUSHDOUBLE(d);/*phil*/
		break;
	      }
	  } else {		/* Both left and right are int */
	    switch(opcode)
	      {
	      case OPISEQUAL:
		*sp++ = il == ir;
		break;
	      case OPISNOTEQUAL:
		*sp++ = il != ir;
		break;
	      case OPISLT:
		*sp++ = il < ir;
		break;
	      case OPISGT:
		*sp++ = il > ir;
		break;
	      case OPISLE:
		*sp++ = il <= ir;
		break;
	      case OPISGE:
		*sp++ = il >= ir;
		break;
	      case OPADD:
		*sp++ = il + ir;
		break;
	      case OPSUB:
		*sp++ = il - ir;
		break;
	      case OPTIMES:
		*sp++ = il * ir; /* Need to deal with overflow */
		break;
	      case OPIDIV:
/*		printf("At OPIDIV all int branch\n");*/
		if(ir == 0) {
		  fprintf(STDERR,"Test block %s: %d/0.0\n",blockname,il);
		  *sp++ = 0;
		} else {
		  *sp++ = il / ir;
		}
		break;
	      case OPDIV:
		if(ir == 0) {
		  fprintf(STDERR,"Test block %s: %d/0.0\n",blockname,il);
		  d = 0.0;
		} else
		  d = dl / dr; /* Need to deal with overflow and div 0 */
		PUSHDOUBLE(d);/*phil*/
		break;
	      case OPMOD:
		*sp++ = il % ir; /* Need to deal with overflow and div 0 */
		break;
	      }
	  }
	  break;
	case OPUNARY:		/* Unary Operators */
	  switch(opcode)
	    {
	    case OPNEG:
	      if(rtype==OPRINT) {
		i = -(*--sp);
	        *sp++ = i;
	      } else if (rtype==OPRFLOAT) {
		f = *(DAFLOAT *)(--sp);
		f = -f;
		*sp++ = *(DAINT *)&f;
	      } else {
		POPDOUBLE(d);/*phil*/
		d = -d;
		PUSHDOUBLE(d);/*phil*/
	      }
	      break;
	    case OPNOT:
	    case OPCOMP:
	      if(rtype==OPRINT) {
		POPINT(i);
	      } else if(rtype==OPRFLOAT) {
		POPFLOAT(f);/*phil*/
		i = floatToLong(f);
	      } else {
		POPDOUBLE(d);/*phil*/
		i = floatToLong(d);
	      }
	      i = (opcode == OPNOT ? !i : ~i);
	      *sp++ = i;
	      break;
	    }
	  break;
	default:
	  fprintf(STDERR,"Test block %s: Operator %x not yet implimented\n",
		  blockname,opcode);
	  break;
	} /* Terminates switch */
    } else {	/* terminates if(rawopcode >=OPLP) *//* Old Style, May not work anymore */
      switch(*pc++)
	{
	case PUSHI:
	  *sp++ = *pc++;
	  break;
	case PUSHS:
	  *sp++ = *((DAINT *) *pc++);
	  pc++;			/* Skip variable name */
	  break;
	case PUSHFTOIS:
	  *sp++ = floatToLong(*((DAFLOAT *) *pc++));
	  pc++;			/* Skip variable name */
	  break;
	case PUSHITOFS:
	  *sp++ = *(DAINT *)&f;
	  pc++;			/* Skip variable name */
	  break;
	case POPS:
	  *((int *) *pc++) = *--sp;
/*	  printf("Putting result %d into %s\n",*sp,(char *) *pc);*/
	  pc++;			/* Skip variable name */
	  break;
	case tGATE:
	  nargs = *pc++;
/*	  printf("GATE: nargs=%d\n",nargs);*/
	  result = ((*((DAFLOAT *) sp-3) >= *((DAFLOAT *) sp-2)) 
		    && (*((DAFLOAT *) sp-1) > *((DAFLOAT *) sp-3)));
	  sp -= nargs;
	  *sp++ = result;
	  break;
	case tEQ:
	  nargs = *pc++;
	  result = (*(sp-1) == *(sp-2));
	  sp -= nargs;
	  *sp++ = result;
	  break;
	case tAND:
	  result = 1;
	  for(nargs = *pc++;(nargs > 0) && result; nargs--){
	    result = (*(--sp)) != 0;
	  }
	  sp -= nargs;
	  *sp++ = result;
	  break;
	case tIOR:
	  result = 0;
	  for(nargs = *pc++;(nargs > 0) && !result; nargs--){
	    result = (*(--sp)) != 0;
	  }
	  sp -= nargs;
	  *sp++ = result;
	  break;
	default:
	  fprintf(STDERR,"Test block %s: Opcode %d not defined\n",
		  blockname,*(pc-1));
	}
/*      printf("Stack depth %d\n",sp-stack);*/
    }
  }
  if(sp != stack){
    fprintf(STDERR,"\n");
    fprintf(STDERR,"Original SP %x\n",stack);
    fprintf(STDERR,"\n\n\n\n");
    fprintf(STDERR,"Final    SP %x\n",sp);
    fprintf(STDERR,"Items left on stack = %d\n",sp-stack);
    return(S_FAILURE);
  }
  return(S_SUCCESS);
}
