MODULE updtide
   !!======================================================================
   !!                       ***  MODULE  updtide  ***
   !! Initialization of tidal forcing
   !!======================================================================
   !! History :  9.0  !  07  (O. Le Galloudec)  Original code
   !!----------------------------------------------------------------------
#if defined key_tide
   !!----------------------------------------------------------------------
   !!   'key_tide' :                                        tidal potential
   !!----------------------------------------------------------------------
   !!   upd_tide       : update tidal potential
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O units
   USE phycst          ! physical constant
   USE sbctide         ! tide potential variable
   USE tideini, ONLY: ln_tide_ramp, rdttideramp

   IMPLICIT NONE
   PUBLIC

   PUBLIC   upd_tide   ! called in dynspg_... modules
  
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: updtide.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE upd_tide( kt, kit, kbaro, koffset )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE upd_tide  ***
      !!
      !! ** Purpose :   provide at each time step the astronomical potential
      !!
      !! ** Method  :   computed from pulsation and amplitude of all tide components
      !!
      !! ** Action  :   pot_astro   actronomical potential
      !!----------------------------------------------------------------------      
      INTEGER, INTENT(in)           ::   kt      ! ocean time-step index
      INTEGER, INTENT(in), OPTIONAL ::   kit     ! external mode sub-time-step index (lk_dynspg_ts=T only)
      INTEGER, INTENT(in), OPTIONAL ::   kbaro   ! number of sub-time-step           (lk_dynspg_ts=T only)
      INTEGER, INTENT(in), OPTIONAL ::   koffset ! time offset in number 
                                                 ! of sub-time-steps                 (lk_dynspg_ts=T only)
      !
      INTEGER  ::   joffset      ! local integer
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zt, zramp    ! local scalar
      REAL(wp), DIMENSION(nb_harmo) ::   zwt 
      !!----------------------------------------------------------------------      
      !
      !                               ! tide pulsation at model time step (or sub-time-step)
      zt = ( kt - kt_tide ) * rdt
      !
      joffset = 0
      IF( PRESENT( koffset ) )   joffset = koffset
      !
      IF( PRESENT( kit ) .AND. PRESENT( kbaro ) )   THEN
         zt = zt + ( kit + 0.5_wp * ( joffset - 1 ) ) * rdt / REAL( kbaro, wp )
      ELSE
         zt = zt + joffset * rdt
      ENDIF
      !
      zwt(:) = omega_tide(:) * zt

      pot_astro(:,:) = 0._wp          ! update tidal potential (sum of all harmonics)
      DO jk = 1, nb_harmo   
         pot_astro(:,:) = pot_astro(:,:) + amp_pot(:,:,jk) * COS( zwt(jk) + phi_pot(:,:,jk) )      
      END DO
      !
      IF( ln_tide_ramp ) THEN         ! linear increase if asked
         zt = ( kt - nit000 ) * rdt
         IF( PRESENT( kit ) .AND. PRESENT( kbaro ) )   zt = zt + kit * rdt / REAL( kbaro, wp )
         zramp = MIN(  MAX( zt / (rdttideramp*rday) , 0._wp ) , 1._wp  )
         pot_astro(:,:) = zramp * pot_astro(:,:)
      ENDIF
      !
   END SUBROUTINE upd_tide

#else
  !!----------------------------------------------------------------------
  !!   Dummy module :                                        NO TIDE
  !!----------------------------------------------------------------------
CONTAINS
  SUBROUTINE upd_tide( kt, kit, kbaro, koffset )          ! Empty routine
    INTEGER, INTENT(in)           ::   kt      !  integer  arg, dummy routine
    INTEGER, INTENT(in), OPTIONAL ::   kit     !  optional arg, dummy routine
    INTEGER, INTENT(in), OPTIONAL ::   kbaro   !  optional arg, dummy routine
    INTEGER, INTENT(in), OPTIONAL ::   koffset !  optional arg, dummy routine
    WRITE(*,*) 'upd_tide: You should not have seen this print! error?', kt
  END SUBROUTINE upd_tide

#endif

  !!======================================================================

END MODULE updtide
