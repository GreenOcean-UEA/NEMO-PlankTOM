MODULE trcini_c14b
   !!======================================================================
   !!                         ***  MODULE trcini_c14b  ***
   !! TOP :   initialisation of the C14 bomb tracer
   !!======================================================================
   !! History :  1.0  ! 2005-10  (Z. Lachkar) Original code
   !!            2.0  ! 2007-12  (C. Ethe) 
   !!----------------------------------------------------------------------
#if defined key_c14b
   !!----------------------------------------------------------------------
   !!   'key_c14b'                                          C14 bomb tracer
   !!----------------------------------------------------------------------
   !! trc_ini_c14b      : C14 model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE trcsms_c14b     ! C14 sms trends

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_c14b   ! called by trcini.F90 module

   INTEGER  ::   nrec            ! number of year in CO2 Concentrations file
   INTEGER  ::   inum1, inum2    ! unit number

   REAL(wp) ::   ys40 = -40.     ! 40 degrees south
   REAL(wp) ::   ys20 = -20.     ! 20 degrees south
   REAL(wp) ::   yn20 =  20.     ! 20 degrees north
   REAL(wp) ::   yn40 =  40.     ! 40 degrees north

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcini_c14b.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_c14b
      !!-------------------------------------------------------------------
      !!                     ***  trc_ini_c14b  ***  
      !!
      !! ** Purpose :   initialization for C14 model
      !!
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl, jm
      !!----------------------------------------------------------------------

      !                     ! Allocate C14b arrays
      IF( trc_sms_c14b_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_c14b : unable to allocate C14b arrays' )

      CALL trc_ctl_c14b     !  Control consitency

      IF(lwp) WRITE(numout,*) ''
      IF(lwp) WRITE(numout,*) ' trc_ini_c14b: initialisation of Bomb C14 chemical model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~'
      IF(lwp) CALL FLUSH(numout)

      ! Initialization of boundaries conditions
      ! --------------------------------------- 
      qtr_c14(:,:) = 0._wp
      
      ! Initialization of qint in case of  no restart 
      !----------------------------------------------
      IF( .NOT. ln_rsttr ) THEN    
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'Initialization de qint ; No restart : qint equal zero '
         ENDIF
         trn     (:,:,:,jpc14) = 0._wp
         qint_c14(:,:        ) = 0._wp
      ENDIF


      ! Read CO2 atmospheric concentrations file...
      !------------------------------------------------

      nrec    = ( nyear_end - nyear_start + 1 )     ! number of year in CO2 Concentrations file
      nco2ti  = 1

      IF(lwp) WRITE(numout,*) 'Read nobomb CO2 atmospheric concentrations file '
      CALL ctl_opn( inum1, 'atmco2.dat', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      REWIND(inum1)
      DO jm = 1, jpmaxrec2
         READ(inum1, *)  co2year(jm), spco2(jm)
      END DO
      WRITE(numout,*)
      CLOSE(inum1)

      IF (lwp) WRITE(numout,*) 'Read C-14 atmospheric concentrations file '

      CALL ctl_opn( inum2, 'atmc14.dat', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      REWIND(inum2)

      ! Skip over 1st descriptor line
      READ(inum2, '(1x)')
      READ(inum2, '(1x)')

      ! READ FILE
      jl = nyear - nyear_start
      DO jm = 1, nrec
         READ(inum2,*) nobombyear(jm), nobomb(jm,1), nobomb(jm,2), nobomb(jm,3)
         IF (lwp .AND. (jm-jl.LT.3).AND. (jm-jl.GT.0)) WRITE(numout, '(f7.1, 3f9.4)') nobombyear(jm), nobomb(jm,1), nobomb(jm,2), nobomb(jm,3)
      END DO
      CLOSE(inum2)
      CALL ctl_opn( inum2, 'atmb14.dat', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      REWIND(inum2)

      ! Skip over 1st descriptor line
      READ(inum2, '(1x)')
      READ(inum2, '(1x)')

      ! READ FILE
      DO jm = 1, nrec
         READ(inum2,*) bombyear(jm), bomb(jm,1), bomb(jm,2), bomb(jm,3)
         IF (lwp .AND. (jm-jl.LT.3).AND. (jm-jl.GT.0)) WRITE(numout, '(f7.1, 3f9.4)') bombyear(jm), bomb(jm,1), bomb(jm,2), bomb(jm,3)
      END DO
      CLOSE(inum2)

      ! Conversion unit : Now atm units are in real C-14 [per mil]
      ! C-14(Orr) = C-14(per mil)/10.0
       DO jm = 1, nrec
         nobomb(jm,1) = ( nobomb(jm,1 ) + 17.40 ) * 0.1
         nobomb(jm,2) = ( nobomb(jm,2 ) + 10.40 ) * 0.1
         nobomb(jm,3) = ( nobomb(jm,3 ) + 14.65 ) * 0.1
         bomb(jm,1) = ( bomb(jm,1 ) + 17.40 ) * 0.1
         bomb(jm,2) = ( bomb(jm,2 ) + 10.40 ) * 0.1
         bomb(jm,3) = ( bomb(jm,3 ) + 14.65 ) * 0.1
       END DO

       ! Linear  interpolation of the C-14 source fonction
       ! in linear latitude band  (20N,40N) and (20S,40S)
       !------------------------------------------------------
       DO jj = 1 , jpj
          DO ji = 1 , jpi
            IF( gphit(ji,jj) >= yn40 ) THEN
                 fareaz(ji,jj,1) = 0.
                 fareaz(ji,jj,2) = 0.
                 fareaz(ji,jj,3) = 1.
            ELSE IF( gphit(ji,jj ) <= ys40) THEN
                 fareaz(ji,jj,1) = 1.
                 fareaz(ji,jj,2) = 0.
                 fareaz(ji,jj,3) = 0.
            ELSE IF( gphit(ji,jj) >= yn20 ) THEN
                 fareaz(ji,jj,1) = 0.
                 fareaz(ji,jj,2) = 2. * ( 1. - gphit(ji,jj) / yn40 )
                 fareaz(ji,jj,3) = 2. * gphit(ji,jj) / yn40 - 1.
            ELSE IF( gphit(ji,jj) <= ys20 ) THEN
                 fareaz(ji,jj,1) = 2. * gphit(ji,jj) / ys40 - 1.
                 fareaz(ji,jj,2) = 2. * ( 1. - gphit(ji,jj) / ys40 )
                 fareaz(ji,jj,3) = 0.
            ELSE
                 fareaz(ji,jj,1) = 0.
                 fareaz(ji,jj,2) = 1.
                 fareaz(ji,jj,3) = 0.
            ENDIF
         END DO
      END DO
      !
      IF(lwp) WRITE(numout,*) 'Initialization of C14 bomb tracer done'
      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) CALL FLUSH(numout)
   END SUBROUTINE trc_ini_c14b


   SUBROUTINE trc_ctl_c14b
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_ctl_c14b  ***
      !!
      !! ** Purpose :   control the cpp options, namelist and files 
      !!----------------------------------------------------------------------

      IF(lwp) THEN
          WRITE(numout,*) ' C14 bomb Model '
          WRITE(numout,*) ' '
      ENDIF
   END SUBROUTINE trc_ctl_c14b
   
#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                    No C14 bomb tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_c14b             ! Empty routine
   END SUBROUTINE trc_ini_c14b
#endif

   !!======================================================================
END MODULE trcini_c14b
