******************************* STRINGLIB.FOR **********************************
C+______________________________________________________________________________
!
!Facility: STRINGLIB
!Version:  2.1
!
!Purpose: A Library of string manipulation routines. See entry points for
!         details.
!
!Entry points:
!
!  last_char(string)
!  strip(cmd,p1)
!  rd_logical(cmd,l)
!  rd_int(cmd,i)
!  rd_real(cmd,x)
!  rd_hex(cmd,i)
!
!Author: David Potterveld - May 1985
!
!Modification History:
!
!  Jan-13-1990 (DHP) Added routine RD_HEX.
!  Mar-24-1992 (DHP) Conversion for UNIX F77 compatibility. STRIP, and RD_xxxx
!		routines are now logical functions. VAX RTL routines replaced.
!  Mar-25-1992 (DHP) Added routine RD_LOGICAL.
!  Aug-24-1993 (DHP) Fixed bug in STRIP.
C-______________________________________________________________________________

!###############################################################################

	integer*4 function last_char(string)
C+______________________________________________________________________________
!
! LAST_CHAR - Return the position of the last character in the string which
! is neither a space or a tab. Returns zero for null or empty strings.
C-______________________________________________________________________________

	implicit none
	integer*4 i
	character*(*) string
	character*1 sp/' '/
	character*1 tab/'	'/

	save

C ============================= Executable Code ================================

	last_char = 0
	do i = 1,len(string)
	   if (string(i:i).ne.sp.and.string(i:i).ne.tab) last_char = i
	enddo

	return
	end

!###############################################################################

	logical function strip(cmd,p1)
C+______________________________________________________________________________
!
! STRIP - strip space and tab separated substrings out of a string.
!
! Operation:
!   1 - Strip off leading separators (blanks and tabs) from CMD.
!   2 - If no non-separator characters encountered, return .FALSE.
!   3 - PUT all characters up to next separator into P1.
!   4 - Strip off all characters including separators from CMD until
!       either the next word is positioned at the beginning of CMD, or
!       the end of the orginal CMD is reached.
!   5 - Return .TRUE.
C-______________________________________________________________________________

	implicit none
	integer*4 char_cnt,len_cmd,len_word,len_p1
	integer*4 pos_of_sp,pos_of_sep,pos_of_tab
	character*(*) cmd,p1
	character*1 c1,sp,tab
	data	sp/' '/
	data	tab/'	'/
	logical more,stripped
	save

C ============================= Executable Code ================================

	char_cnt = 0				!count processed characters
	len_cmd = len(cmd)			!length of command string
	len_p1 = len(p1)			!length of substring
	more = .true.				!more work to do
	stripped = .false.			!haven't stripped word

	do while (more.and.char_cnt.lt.len_cmd)
	   c1 = cmd(1:1)			!get leading character
	   if (c1.eq.sp.or.c1.eq.tab) then	!strip away separators
	      cmd = cmd(2:)
	      char_cnt = char_cnt + 1
	   else					!Found non-separator char.
	      if (stripped) then
		more = .false.			!All done.
	      else	      			!Extract substring.
		stripped = .true.
		pos_of_sp = index(cmd,sp)		!find position of next sep.
		pos_of_tab = index(cmd,tab)
		if (pos_of_sp.eq.0) pos_of_sp = len_cmd + 1
		if (pos_of_tab.eq.0) pos_of_tab = len_cmd + 1
		pos_of_sep = min(pos_of_sp,pos_of_tab)
		len_word = pos_of_sep - 1
		p1 = cmd(1:min(len_word,len_p1))	!save the word
		cmd = cmd(min(pos_of_sep,len_cmd):len_cmd) !remove word from cmd
		char_cnt = char_cnt + len_word
	      endif
	   endif
	enddo

	strip = stripped			!return strip condition
	return
	end

!###############################################################################

	logical function rd_int(string,number)
C+______________________________________________________________________________
!
! Strip the leading substring out of STRING and try to convert it into
! An integer. If STRING is blank, or the conversion failed, the function
! returns .FALSE. Otherwise, .TRUE. is returned, and the integer value
! is returned in NUMBER.
C-______________________________________________________________________________

	implicit none
	integer*4 l,last_char

	integer*4 number
	character*(*) string
	character*32 str1
	logical strip
	save

C ============================= Executable Code ================================

	rd_int = .false.			!assume failure
	if (.not.strip(string,str1)) return	!Return if no words
	l = last_char(str1)
	read (str1(1:l),'(i)',err=999) number	!Try to read word as integer.
	rd_int = .true.				!Success.
999	return
	end

!###############################################################################

	logical function rd_logical(string,value)
C+______________________________________________________________________________
!
! Strip the leading substring out of STRING and try to convert it into
! a logical value. If STRING is blank, or the conversion failed, the function
! returns .FALSE. Otherwise, .TRUE. is returned, and the integer value
! is returned in value.
C-______________________________________________________________________________

	implicit none
	integer*4 l,last_char

	logical value
	character*(*) string
	character*32 str1
	logical strip
	save

C ============================= Executable Code ================================

	rd_logical = .false.			!assume failure
	if (.not.strip(string,str1)) return	!Return if no words
	l = last_char(str1)
	read (str1(1:l),'(L)',err=999) value	!Try to read word as logical.
	rd_logical = .true.			!Success.
999	return
	end

!###############################################################################

	logical function rd_hex(string,number)
C+______________________________________________________________________________
!
! Strip the leading substring out of STRING and try to convert it from a
! HEX string into an integer. If STRING is blank, or the conversion failed,
! the function returns .FALSE. Otherwise, .TRUE. is returned, and the
! integer value is returned in NUMBER.
C-______________________________________________________________________________

	implicit none
	integer*4 l,last_char

	integer*4 number
	character*(*) string
	character*32 str1
	logical strip
	save

C ============================= Executable Code ================================

	rd_hex = .false.			!assume failure
	if (.not.strip(string,str1)) return
	l = last_char(str1)
	read (str1(1:l),'(Z)',err=999) number
	rd_hex = .true.
999	return
	end

!###############################################################################

	logical function rd_real(string,number)
C+______________________________________________________________________________
!
! Strip the leading substring out of STRING and try to convert it into
! An R*4 number. If STRING is blank, or the conversion failed, the function
! returns .FALSE. Otherwise, .TRUE. is returned, and the R*4 value
! is returned in NUMBER.
C-______________________________________________________________________________

	implicit none
	integer*4 l,last_char

	real*8 number
	character*(*) string
	character*32 str1
	logical	strip
	save

C ============================= Executable Code ================================

	rd_real = .false.			!assume failure
	if (.not.strip(string,str1)) return
	l = last_char(str1)
	read (str1(1:l),*,err=999) number
	rd_real = .true.
999	return
	end
