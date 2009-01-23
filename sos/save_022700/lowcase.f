!--------------------------------------------------------------------
C       subroutine lowcase(str1,str2,length)
C       Which returns in str2 a lower cased copy of str1
C       Length is the length of the 2 strings (or max chars
C       to convert
C
C	pillaged (and slightly modified)
C	from the /src directory of T. Payerle, UMD.

        subroutine lowcase(str1,str2)
        implicit none
        character*(*) str1,str2

C       Local vars
        character*26 caps,lcase
        data caps /'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
        data lcase /'abcdefghijklmnopqrstuvwxyz' /
        integer i,j,l

	l = len(str1)

        do i=1,l
                j = index(caps,str1(i:i) )
C               index returns 0 if can't find substring
                if (j .gt. 0) then
                        str2(i:i) = lcase(j:j)
                else
                        str2(i:i) = str1(i:i)
                endif
        enddo

        end
