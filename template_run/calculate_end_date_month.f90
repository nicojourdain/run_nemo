program end_date

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain
!
! usage : ./calculate_end_date_month year month day duration isleap
!
!         isleap=1 if leap years are accounted for, =0 otherwise
!
! output : yearf, monthf, dayf (final if run length = duration)
!          durcorr (corrected duration to finish :
!                       - not after the end of current year
!                       - at the end of a finite month      )
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE

INTEGER :: an0, mois0, jour0, duration ! duration in nb days

INTEGER :: an, mois, jour, durleft, durcorr, nb_args, iargc, isleap

CHARACTER(LEN=8) :: an0_char, mois0_char, jour0_char, duration_char, isleap_char

INTEGER, DIMENSION(12) :: ndays_month

LOGICAL :: tocken

nb_args=iargc()

IF (nb_args .EQ. 5) THEN
  CALL getarg ( 1, an0_char )
  CALL getarg ( 2, mois0_char )
  CALL getarg ( 3, jour0_char )
  CALL getarg ( 4, duration_char )
  CALL getarg ( 5, isleap_char )
ELSE
  WRITE(*,*) "[calculate_end_date_month.f90] Error"
  STOP
ENDIF

READ(an0_char,*) an0
READ(mois0_char,*) mois0
READ(jour0_char,*) jour0
READ(duration_char,*) duration
READ(isleap_char,*) isleap

ndays_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)

tocken = .true.

durleft = duration
durcorr = duration
an      = an0
mois    = mois0
jour    = jour0

if ( ( isleap .eq. 1 ) .and. &
&    ( ( MOD(an,4) .eq. 0 .and. .not. ( MOD(an,100) .eq. 0 ) ) .or. MOD(an,400) .eq. 0 ) ) then
 ndays_month(2) = 29
endif

do while ( durleft .gt. 0 )

  if ( (jour+durleft) .gt. ndays_month(mois) ) then

    if ( mois .eq. 12 ) then

      durleft = durleft - ( ndays_month(mois) - jour + 1 )
      if ( tocken ) then
        durcorr = duration - durleft
        tocken = .false.
      endif
      an      = an + 1
      mois    = 1
      jour    = 1

    else

      durleft = durleft - ( ndays_month(mois) - jour + 1 )
      mois    = mois + 1
      jour    = 1

    endif

  else

    jour    = jour + durleft
    durleft = 0

  endif

enddo

!====================
if ( jour .ne. 1 .and. an .eq. an0 ) then
  durcorr = durcorr - jour + 1
  !jour = 1
endif
!====================

write(*,898) an, mois, jour, durcorr
898 FORMAT(i4.4,' ',i2.2,' ',i2.2,' ',i3)

end program end_date
