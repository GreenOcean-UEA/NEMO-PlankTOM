MODULE agrif_opa_interp
   !!======================================================================
   !!                   ***  MODULE  agrif_opa_interp  ***
   !! AGRIF: interpolation package
   !!======================================================================
   !! History :  2.0  !  2002-06  (XXX)  Original cade
   !!             -   !  2005-11  (XXX) 
   !!            3.2  !  2009-04  (R. Benshila) 
   !!----------------------------------------------------------------------
#if defined key_agrif && ! defined key_offline
   !!----------------------------------------------------------------------
   !!   'key_agrif'                                              AGRIF zoom
   !!   NOT 'key_offline'                               NO off-line tracers
   !!----------------------------------------------------------------------
   !!   Agrif_tra     :
   !!   Agrif_dyn     : 
   !!   interpu       :
   !!   interpv       :
   !!----------------------------------------------------------------------
   USE par_oce
   USE oce
   USE dom_oce      
   USE sol_oce
   USE agrif_oce
   USE phycst
   USE in_out_manager
   USE agrif_opa_sponge
   USE lib_mpp
   USE wrk_nemo
   USE dynspg_oce

   IMPLICIT NONE
   PRIVATE

   ! Barotropic arrays used to store open boundary data during
   ! time-splitting loop:
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::  ubdy_w, vbdy_w, hbdy_w
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::  ubdy_e, vbdy_e, hbdy_e
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::  ubdy_n, vbdy_n, hbdy_n
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::  ubdy_s, vbdy_s, hbdy_s
    
   PUBLIC   Agrif_tra, Agrif_dyn, Agrif_ssh, Agrif_dyn_ts, Agrif_ssh_ts, Agrif_dta_ts
   PUBLIC   interpu, interpv, interpunb, interpvnb, interpsshn

#  include "domzgr_substitute.h90"  
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/NST 3.3 , NEMO Consortium (2010)
   !! $Id: agrif_opa_interp.F90 4486 2014-02-05 11:23:56Z jchanut $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   CONTAINS
   
   SUBROUTINE Agrif_tra
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_Tra  ***
      !!----------------------------------------------------------------------
      !!
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::   zrhox , alpha1, alpha2, alpha3
      REAL(wp) ::   alpha4, alpha5, alpha6, alpha7
      REAL(wp), POINTER, DIMENSION(:,:,:,:) :: ztsa
      !!----------------------------------------------------------------------
      !
      IF( Agrif_Root() )   RETURN

      CALL wrk_alloc( jpi, jpj, jpk, jpts, ztsa ) 

      Agrif_SpecialValue    = 0.e0
      Agrif_UseSpecialValue = .TRUE.
      ztsa(:,:,:,:) = 0.e0

      CALL Agrif_Bc_variable( ztsa, tsn_id, procname=interptsn )
      Agrif_UseSpecialValue = .FALSE.

      zrhox = Agrif_Rhox()

      alpha1 = ( zrhox - 1. ) * 0.5
      alpha2 = 1. - alpha1

      alpha3 = ( zrhox - 1. ) / ( zrhox + 1. )
      alpha4 = 1. - alpha3

      alpha6 = 2. * ( zrhox - 1. ) / ( zrhox + 1. )
      alpha7 =    - ( zrhox - 1. ) / ( zrhox + 3. )
      alpha5 = 1. - alpha6 - alpha7

      IF( nbondi == 1 .OR. nbondi == 2 ) THEN

         DO jn = 1, jpts
            tsa(nlci,:,:,jn) = alpha1 * ztsa(nlci,:,:,jn) + alpha2 * ztsa(nlci-1,:,:,jn)
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  IF( umask(nlci-2,jj,jk) == 0.e0 ) THEN
                     tsa(nlci-1,jj,jk,jn) = tsa(nlci,jj,jk,jn) * tmask(nlci-1,jj,jk)
                  ELSE
                     tsa(nlci-1,jj,jk,jn)=(alpha4*tsa(nlci,jj,jk,jn)+alpha3*tsa(nlci-2,jj,jk,jn))*tmask(nlci-1,jj,jk)
                     IF( un(nlci-2,jj,jk) > 0.e0 ) THEN
                        tsa(nlci-1,jj,jk,jn)=( alpha6*tsa(nlci-2,jj,jk,jn)+alpha5*tsa(nlci,jj,jk,jn)  &
                           &                 + alpha7*tsa(nlci-3,jj,jk,jn) ) * tmask(nlci-1,jj,jk)
                     ENDIF
                  ENDIF
               END DO
            END DO
         ENDDO
      ENDIF

      IF( nbondj == 1 .OR. nbondj == 2 ) THEN

         DO jn = 1, jpts
            tsa(:,nlcj,:,jn) = alpha1 * ztsa(:,nlcj,:,jn) + alpha2 * ztsa(:,nlcj-1,:,jn)
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  IF( vmask(ji,nlcj-2,jk) == 0.e0 ) THEN
                     tsa(ji,nlcj-1,jk,jn) = tsa(ji,nlcj,jk,jn) * tmask(ji,nlcj-1,jk)
                  ELSE
                     tsa(ji,nlcj-1,jk,jn)=(alpha4*tsa(ji,nlcj,jk,jn)+alpha3*tsa(ji,nlcj-2,jk,jn))*tmask(ji,nlcj-1,jk)        
                     IF (vn(ji,nlcj-2,jk) > 0.e0 ) THEN
                        tsa(ji,nlcj-1,jk,jn)=( alpha6*tsa(ji,nlcj-2,jk,jn)+alpha5*tsa(ji,nlcj,jk,jn)  &
                           &                 + alpha7*tsa(ji,nlcj-3,jk,jn) ) * tmask(ji,nlcj-1,jk)
                     ENDIF
                  ENDIF
               END DO
            END DO
         ENDDO 
      ENDIF

      IF( nbondi == -1 .OR. nbondi == 2 ) THEN
         DO jn = 1, jpts
            tsa(1,:,:,jn) = alpha1 * ztsa(1,:,:,jn) + alpha2 * ztsa(2,:,:,jn)
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  IF( umask(2,jj,jk) == 0.e0 ) THEN
                     tsa(2,jj,jk,jn) = tsa(1,jj,jk,jn) * tmask(2,jj,jk)
                  ELSE
                     tsa(2,jj,jk,jn)=(alpha4*tsa(1,jj,jk,jn)+alpha3*tsa(3,jj,jk,jn))*tmask(2,jj,jk)        
                     IF( un(2,jj,jk) < 0.e0 ) THEN
                        tsa(2,jj,jk,jn)=(alpha6*tsa(3,jj,jk,jn)+alpha5*tsa(1,jj,jk,jn)+alpha7*tsa(4,jj,jk,jn))*tmask(2,jj,jk)
                     ENDIF
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

      IF( nbondj == -1 .OR. nbondj == 2 ) THEN
         DO jn = 1, jpts
            tsa(:,1,:,jn) = alpha1 * ztsa(:,1,:,jn) + alpha2 * ztsa(:,2,:,jn)
            DO jk=1,jpk      
               DO ji=1,jpi
                  IF( vmask(ji,2,jk) == 0.e0 ) THEN
                     tsa(ji,2,jk,jn)=tsa(ji,1,jk,jn) * tmask(ji,2,jk)
                  ELSE
                     tsa(ji,2,jk,jn)=(alpha4*tsa(ji,1,jk,jn)+alpha3*tsa(ji,3,jk,jn))*tmask(ji,2,jk)
                     IF( vn(ji,2,jk) < 0.e0 ) THEN
                        tsa(ji,2,jk,jn)=(alpha6*tsa(ji,3,jk,jn)+alpha5*tsa(ji,1,jk,jn)+alpha7*tsa(ji,4,jk,jn))*tmask(ji,2,jk)
                     ENDIF
                  ENDIF
               END DO
            END DO
         ENDDO
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, jpk, jpts, ztsa ) 
      !
   END SUBROUTINE Agrif_tra


   SUBROUTINE Agrif_dyn( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_DYN  ***
      !!----------------------------------------------------------------------  
      !! 
      INTEGER, INTENT(in) ::   kt
      !!
      INTEGER :: ji,jj,jk
      REAL(wp) :: timeref
      REAL(wp) :: z2dt, znugdt
      REAL(wp) :: zrhox, zrhoy
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zua, zva
      REAL(wp), POINTER, DIMENSION(:,:) :: spgv1, spgu1, zua2d, zva2d
      !!----------------------------------------------------------------------  

      IF( Agrif_Root() )   RETURN

      CALL wrk_alloc( jpi, jpj, spgv1, spgu1, zua2d, zva2d )
      CALL wrk_alloc( jpi, jpj, jpk, zua, zva )

      zrhox = Agrif_Rhox()
      zrhoy = Agrif_Rhoy()

      timeref = 1.

      ! time step: leap-frog
      z2dt = 2. * rdt
      ! time step: Euler if restart from rest
      IF( neuler == 0 .AND. kt == nit000 ) z2dt = rdt
      ! coefficients
      znugdt =  grav * z2dt    

      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = ln_spc_dyn

      zua = 0.
      zva = 0.
      CALL Agrif_Bc_variable(zua,un_id,procname=interpu)
      CALL Agrif_Bc_variable(zva,vn_id,procname=interpv)
      zua2d = 0.
      zva2d = 0.

#if defined key_dynspg_flt
      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = ln_spc_dyn
      CALL Agrif_Bc_variable(zua2d,e1u_id,calledweight=1.,procname=interpu2d)
      CALL Agrif_Bc_variable(zva2d,e2v_id,calledweight=1.,procname=interpv2d)
#endif
      Agrif_UseSpecialValue = .FALSE.


      IF((nbondi == -1).OR.(nbondi == 2)) THEN

#if defined key_dynspg_flt
         DO jj=1,jpj
            laplacu(2,jj) = timeref * (zua2d(2,jj)/(zrhoy*e2u(2,jj)))*umask(2,jj,1)
         END DO
#endif

         DO jk=1,jpkm1
            DO jj=1,jpj
               ua(1:2,jj,jk) = (zua(1:2,jj,jk)/(zrhoy*e2u(1:2,jj)))
               ua(1:2,jj,jk) = ua(1:2,jj,jk) / fse3u_a(1:2,jj,jk)
            END DO
         END DO

#if defined key_dynspg_flt
         DO jk=1,jpkm1
            DO jj=1,jpj
               ua(2,jj,jk) = (ua(2,jj,jk) - z2dt * znugdt * laplacu(2,jj))*umask(2,jj,jk)
            END DO
         END DO

         spgu(2,:)=0.

         DO jk=1,jpkm1
            DO jj=1,jpj
               spgu(2,jj)=spgu(2,jj)+fse3u_a(2,jj,jk)*ua(2,jj,jk)
            END DO
         END DO

         DO jj=1,jpj
            IF (umask(2,jj,1).NE.0.) THEN
               spgu(2,jj)=spgu(2,jj)*hur_a(2,jj)
            ENDIF
         END DO
#else
         spgu(2,:) = ua_b(2,:)
#endif

         DO jk=1,jpkm1
            DO jj=1,jpj
               ua(2,jj,jk) = 0.25*(ua(1,jj,jk)+2.*ua(2,jj,jk)+ua(3,jj,jk))
               ua(2,jj,jk) = ua(2,jj,jk) * umask(2,jj,jk)
            END DO
         END DO

         spgu1(2,:)=0.

         DO jk=1,jpkm1
            DO jj=1,jpj
               spgu1(2,jj)=spgu1(2,jj)+fse3u_a(2,jj,jk)*ua(2,jj,jk)
            END DO
         END DO

         DO jj=1,jpj
            IF (umask(2,jj,1).NE.0.) THEN
               spgu1(2,jj)=spgu1(2,jj)*hur_a(2,jj)
            ENDIF
         END DO

         DO jk=1,jpkm1
            DO jj=1,jpj
               ua(2,jj,jk) = (ua(2,jj,jk)+spgu(2,jj)-spgu1(2,jj))*umask(2,jj,jk)
            END DO
         END DO

         DO jk=1,jpkm1
            DO jj=1,jpj
               va(2,jj,jk) = (zva(2,jj,jk)/(zrhox*e1v(2,jj)))*vmask(2,jj,jk)
               va(2,jj,jk) = va(2,jj,jk) / fse3v_a(2,jj,jk)
            END DO
         END DO

#if defined key_dynspg_ts
         ! Set tangential velocities to time splitting estimate
         spgv1(2,:)=0.
         DO jk=1,jpkm1
            DO jj=1,jpj
               spgv1(2,jj)=spgv1(2,jj)+fse3v_a(2,jj,jk)*va(2,jj,jk)
            END DO
         END DO

         DO jj=1,jpj
            spgv1(2,jj)=spgv1(2,jj)*hvr_a(2,jj)
         END DO

         DO jk=1,jpkm1
            DO jj=1,jpj
               va(2,jj,jk) = (va(2,jj,jk)+va_b(2,jj)-spgv1(2,jj))*vmask(2,jj,jk)
            END DO
         END DO
#endif

      ENDIF

      IF((nbondi == 1).OR.(nbondi == 2)) THEN
#if defined key_dynspg_flt
         DO jj=1,jpj
            laplacu(nlci-2,jj) = timeref * (zua2d(nlci-2,jj)/(zrhoy*e2u(nlci-2,jj)))
         END DO
#endif

         DO jk=1,jpkm1
            DO jj=1,jpj
               ua(nlci-2:nlci-1,jj,jk) = (zua(nlci-2:nlci-1,jj,jk)/(zrhoy*e2u(nlci-2:nlci-1,jj)))
               ua(nlci-2:nlci-1,jj,jk) = ua(nlci-2:nlci-1,jj,jk) / fse3u_a(nlci-2:nlci-1,jj,jk)
            END DO
         END DO

#if defined key_dynspg_flt
         DO jk=1,jpkm1
            DO jj=1,jpj
               ua(nlci-2,jj,jk) = (ua(nlci-2,jj,jk)- z2dt * znugdt * laplacu(nlci-2,jj))*umask(nlci-2,jj,jk)
            END DO
         END DO


         spgu(nlci-2,:)=0.

         do jk=1,jpkm1
            do jj=1,jpj
               spgu(nlci-2,jj)=spgu(nlci-2,jj)+fse3u_a(nlci-2,jj,jk)*ua(nlci-2,jj,jk)
            enddo
         enddo

         DO jj=1,jpj
            IF (umask(nlci-2,jj,1).NE.0.) THEN
               spgu(nlci-2,jj)=spgu(nlci-2,jj)*hur_a(nlci-2,jj)
            ENDIF
         END DO
#else
         spgu(nlci-2,:) = ua_b(nlci-2,:)
#endif

         DO jk=1,jpkm1
            DO jj=1,jpj
               ua(nlci-2,jj,jk) = 0.25*(ua(nlci-3,jj,jk)+2.*ua(nlci-2,jj,jk)+ua(nlci-1,jj,jk))

               ua(nlci-2,jj,jk) = ua(nlci-2,jj,jk) * umask(nlci-2,jj,jk)

            END DO
         END DO

         spgu1(nlci-2,:)=0.

         DO jk=1,jpkm1
            DO jj=1,jpj
               spgu1(nlci-2,jj)=spgu1(nlci-2,jj)+fse3u_a(nlci-2,jj,jk)*ua(nlci-2,jj,jk)*umask(nlci-2,jj,jk)
            END DO
         END DO

         DO jj=1,jpj
            IF (umask(nlci-2,jj,1).NE.0.) THEN
               spgu1(nlci-2,jj)=spgu1(nlci-2,jj)*hur_a(nlci-2,jj)
            ENDIF
         END DO

         DO jk=1,jpkm1
            DO jj=1,jpj
               ua(nlci-2,jj,jk) = (ua(nlci-2,jj,jk)+spgu(nlci-2,jj)-spgu1(nlci-2,jj))*umask(nlci-2,jj,jk)
            END DO
         END DO

         DO jk=1,jpkm1
            DO jj=1,jpj-1
               va(nlci-1,jj,jk) = (zva(nlci-1,jj,jk)/(zrhox*e1v(nlci-1,jj)))*vmask(nlci-1,jj,jk)
               va(nlci-1,jj,jk) = va(nlci-1,jj,jk) / fse3v_a(nlci-1,jj,jk)
            END DO
         END DO

#if defined key_dynspg_ts
         ! Set tangential velocities to time splitting estimate
         spgv1(nlci-1,:)=0._wp
         DO jk=1,jpkm1
            DO jj=1,jpj
               spgv1(nlci-1,jj)=spgv1(nlci-1,jj)+fse3v_a(nlci-1,jj,jk)*va(nlci-1,jj,jk)*vmask(nlci-1,jj,jk)
            END DO
         END DO

         DO jj=1,jpj
            spgv1(nlci-1,jj)=spgv1(nlci-1,jj)*hvr_a(nlci-1,jj)
         END DO

         DO jk=1,jpkm1
            DO jj=1,jpj
               va(nlci-1,jj,jk) = (va(nlci-1,jj,jk)+va_b(nlci-1,jj)-spgv1(nlci-1,jj))*vmask(nlci-1,jj,jk)
            END DO
         END DO
#endif

      ENDIF

      IF((nbondj == -1).OR.(nbondj == 2)) THEN

#if defined key_dynspg_flt
         DO ji=1,jpi
            laplacv(ji,2) = timeref * (zva2d(ji,2)/(zrhox*e1v(ji,2)))
         END DO
#endif

         DO jk=1,jpkm1
            DO ji=1,jpi
               va(ji,1:2,jk) = (zva(ji,1:2,jk)/(zrhox*e1v(ji,1:2)))
               va(ji,1:2,jk) = va(ji,1:2,jk) / fse3v_a(ji,1:2,jk)
            END DO
         END DO

#if defined key_dynspg_flt
         DO jk=1,jpkm1
            DO ji=1,jpi
               va(ji,2,jk) = (va(ji,2,jk) - z2dt * znugdt * laplacv(ji,2))*vmask(ji,2,jk)
            END DO
         END DO

         spgv(:,2)=0.

         DO jk=1,jpkm1
            DO ji=1,jpi
               spgv(ji,2)=spgv(ji,2)+fse3v_a(ji,2,jk)*va(ji,2,jk)
            END DO
         END DO

         DO ji=1,jpi
            IF (vmask(ji,2,1).NE.0.) THEN
               spgv(ji,2)=spgv(ji,2)*hvr_a(ji,2)
            ENDIF
         END DO
#else
         spgv(:,2)=va_b(:,2)
#endif

         DO jk=1,jpkm1
            DO ji=1,jpi
               va(ji,2,jk)=0.25*(va(ji,1,jk)+2.*va(ji,2,jk)+va(ji,3,jk))
               va(ji,2,jk)=va(ji,2,jk)*vmask(ji,2,jk)
            END DO
         END DO

         spgv1(:,2)=0.

         DO jk=1,jpkm1
            DO ji=1,jpi
               spgv1(ji,2)=spgv1(ji,2)+fse3v_a(ji,2,jk)*va(ji,2,jk)*vmask(ji,2,jk)
            END DO
         END DO

         DO ji=1,jpi
            IF (vmask(ji,2,1).NE.0.) THEN
               spgv1(ji,2)=spgv1(ji,2)*hvr_a(ji,2)
            ENDIF
         END DO

         DO jk=1,jpkm1
            DO ji=1,jpi
               va(ji,2,jk) = (va(ji,2,jk)+spgv(ji,2)-spgv1(ji,2))*vmask(ji,2,jk)
            END DO
         END DO

         DO jk=1,jpkm1
            DO ji=1,jpi
               ua(ji,2,jk) = (zua(ji,2,jk)/(zrhoy*e2u(ji,2)))*umask(ji,2,jk) 
               ua(ji,2,jk) = ua(ji,2,jk) / fse3u_a(ji,2,jk)
            END DO
         END DO

#if defined key_dynspg_ts
         ! Set tangential velocities to time splitting estimate
         spgu1(:,2)=0._wp
         DO jk=1,jpkm1
            DO ji=1,jpi
               spgu1(ji,2)=spgu1(ji,2)+fse3u_a(ji,2,jk)*ua(ji,2,jk)*umask(ji,2,jk)
            END DO
         END DO

         DO ji=1,jpi
            spgu1(ji,2)=spgu1(ji,2)*hur_a(ji,2)
         END DO

         DO jk=1,jpkm1
            DO ji=1,jpi
               ua(ji,2,jk) = (ua(ji,2,jk)+ua_b(ji,2)-spgu1(ji,2))*umask(ji,2,jk)
            END DO
         END DO
#endif
      ENDIF

      IF((nbondj == 1).OR.(nbondj == 2)) THEN

#if defined key_dynspg_flt
         DO ji=1,jpi
            laplacv(ji,nlcj-2) = timeref * (zva2d(ji,nlcj-2)/(zrhox*e1v(ji,nlcj-2)))
         END DO
#endif

         DO jk=1,jpkm1
            DO ji=1,jpi
               va(ji,nlcj-2:nlcj-1,jk) = (zva(ji,nlcj-2:nlcj-1,jk)/(zrhox*e1v(ji,nlcj-2:nlcj-1)))
               va(ji,nlcj-2:nlcj-1,jk) = va(ji,nlcj-2:nlcj-1,jk) / fse3v_a(ji,nlcj-2:nlcj-1,jk)
            END DO
         END DO

#if defined key_dynspg_flt
         DO jk=1,jpkm1
            DO ji=1,jpi
               va(ji,nlcj-2,jk) = (va(ji,nlcj-2,jk)-z2dt * znugdt * laplacv(ji,nlcj-2))*vmask(ji,nlcj-2,jk)
            END DO
         END DO

         spgv(:,nlcj-2)=0.

         DO jk=1,jpkm1
            DO ji=1,jpi
               spgv(ji,nlcj-2)=spgv(ji,nlcj-2)+fse3v_a(ji,nlcj-2,jk)*va(ji,nlcj-2,jk)
            END DO
         END DO

         DO ji=1,jpi
            IF (vmask(ji,nlcj-2,1).NE.0.) THEN
               spgv(ji,nlcj-2)=spgv(ji,nlcj-2)*hvr_a(ji,nlcj-2)
            ENDIF
         END DO
#else
         spgv(:,nlcj-2)=va_b(:,nlcj-2)
#endif

         DO jk=1,jpkm1
            DO ji=1,jpi
               va(ji,nlcj-2,jk)=0.25*(va(ji,nlcj-3,jk)+2.*va(ji,nlcj-2,jk)+va(ji,nlcj-1,jk))
               va(ji,nlcj-2,jk) = va(ji,nlcj-2,jk) * vmask(ji,nlcj-2,jk)
            END DO
         END DO

         spgv1(:,nlcj-2)=0.

         DO jk=1,jpkm1
            DO ji=1,jpi
               spgv1(ji,nlcj-2)=spgv1(ji,nlcj-2)+fse3v_a(ji,nlcj-2,jk)*va(ji,nlcj-2,jk)
            END DO
         END DO

         DO ji=1,jpi
            IF (vmask(ji,nlcj-2,1).NE.0.) THEN
               spgv1(ji,nlcj-2)=spgv1(ji,nlcj-2)*hvr_a(ji,nlcj-2)
            ENDIF
         END DO

         DO jk=1,jpkm1
            DO ji=1,jpi
               va(ji,nlcj-2,jk) = (va(ji,nlcj-2,jk)+spgv(ji,nlcj-2)-spgv1(ji,nlcj-2))*vmask(ji,nlcj-2,jk)
            END DO
         END DO

         DO jk=1,jpkm1
            DO ji=1,jpi
               ua(ji,nlcj-1,jk) = (zua(ji,nlcj-1,jk)/(zrhoy*e2u(ji,nlcj-1)))*umask(ji,nlcj-1,jk)
               ua(ji,nlcj-1,jk) = ua(ji,nlcj-1,jk) / fse3u_a(ji,nlcj-1,jk)
            END DO
         END DO

#if defined key_dynspg_ts
         ! Set tangential velocities to time splitting estimate
         spgu1(:,nlcj-1)=0._wp
         DO jk=1,jpkm1
            DO ji=1,jpi
               spgu1(ji,nlcj-1)=spgu1(ji,nlcj-1)+fse3u_a(ji,nlcj-1,jk)*ua(ji,nlcj-1,jk)
            END DO
         END DO

         DO ji=1,jpi
            spgu1(ji,nlcj-1)=spgu1(ji,nlcj-1)*hur_a(ji,nlcj-1)
         END DO

         DO jk=1,jpkm1
            DO ji=1,jpi
               ua(ji,nlcj-1,jk) = (ua(ji,nlcj-1,jk)+ua_b(ji,nlcj-1)-spgu1(ji,nlcj-1))*umask(ji,nlcj-1,jk)
            END DO
         END DO
#endif

      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, spgv1, spgu1, zua2d, zva2d )
      CALL wrk_dealloc( jpi, jpj, jpk, zua, zva )
      !
   END SUBROUTINE Agrif_dyn

   SUBROUTINE Agrif_dyn_ts( jn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_dyn_ts  ***
      !!----------------------------------------------------------------------  
      !! 
      INTEGER, INTENT(in) ::   jn
      !!
      INTEGER :: ji, jj
      !!----------------------------------------------------------------------  

      IF( Agrif_Root() )   RETURN

      IF((nbondi == -1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
            va_e(2,jj) = vbdy_w(jj) * hvr_e(2,jj)
! Specified fluxes:
            ua_e(2,jj) = ubdy_w(jj) * hur_e(2,jj)
! Characteristics method:
!alt            ua_e(2,jj) = 0.5_wp * ( ubdy_w(jj) * hur_e(2,jj) + ua_e(3,jj) &
!alt                       &           - sqrt(grav * hur_e(2,jj)) * (sshn_e(3,jj) - hbdy_w(jj)) )
         END DO
      ENDIF

      IF((nbondi == 1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
            va_e(nlci-1,jj) = vbdy_e(jj) * hvr_e(nlci-1,jj)
! Specified fluxes:
            ua_e(nlci-2,jj) = ubdy_e(jj) * hur_e(nlci-2,jj)
! Characteristics method:
!alt            ua_e(nlci-2,jj) = 0.5_wp * ( ubdy_e(jj) * hur_e(nlci-2,jj) + ua_e(nlci-3,jj) &
!alt                            &           + sqrt(grav * hur_e(nlci-2,jj)) * (sshn_e(nlci-2,jj) - hbdy_e(jj)) )
         END DO
      ENDIF

      IF((nbondj == -1).OR.(nbondj == 2)) THEN
         DO ji=1,jpi
            ua_e(ji,2) = ubdy_s(ji) * hur_e(ji,2)
! Specified fluxes:
            va_e(ji,2) = vbdy_s(ji) * hvr_e(ji,2)
! Characteristics method:
!alt            va_e(ji,2) = 0.5_wp * ( vbdy_s(ji) * hvr_e(ji,2) + va_e(ji,3) &
!alt                       &           - sqrt(grav * hvr_e(ji,2)) * (sshn_e(ji,3) - hbdy_s(ji)) )
         END DO
      ENDIF

      IF((nbondj == 1).OR.(nbondj == 2)) THEN
         DO ji=1,jpi
            ua_e(ji,nlcj-1) = ubdy_n(ji) * hur_e(ji,nlcj-1)
! Specified fluxes:
            va_e(ji,nlcj-2) = vbdy_n(ji) * hvr_e(ji,nlcj-2)
! Characteristics method:
!alt            va_e(ji,nlcj-2) = 0.5_wp * ( vbdy_n(ji) * hvr_e(ji,nlcj-2)  + va_e(ji,nlcj-3) &
!alt                            &           + sqrt(grav * hvr_e(ji,nlcj-2)) * (sshn_e(ji,nlcj-2) - hbdy_n(ji)) )
         END DO
      ENDIF
      !
   END SUBROUTINE Agrif_dyn_ts

   SUBROUTINE Agrif_dta_ts( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_dta_ts  ***
      !!----------------------------------------------------------------------  
      !! 
      INTEGER, INTENT(in) ::   kt
      !!
      INTEGER :: ji, jj
      LOGICAL :: ll_int_cons
      REAL(wp) :: zrhox, zrhoy, zrhot, zt
      REAL(wp) :: zaa, zab, zat
      REAL(wp) :: zt0, zt1
      REAL(wp), POINTER, DIMENSION(:,:) :: zunb, zvnb, zsshn
      REAL(wp), POINTER, DIMENSION(:,:) :: zuab, zvab, zubb, zvbb, zutn, zvtn
      !!----------------------------------------------------------------------  

      IF( Agrif_Root() )   RETURN

      ll_int_cons = ln_bt_fw ! Assume conservative temporal integration in
                             ! the forward case only

      zrhox = Agrif_Rhox()
      zrhoy = Agrif_Rhoy()
      zrhot = Agrif_rhot()

      IF ( kt==nit000 ) THEN ! Allocate boundary data arrays
         ALLOCATE( ubdy_w(jpj), vbdy_w(jpj), hbdy_w(jpj))
         ALLOCATE( ubdy_e(jpj), vbdy_e(jpj), hbdy_e(jpj))
         ALLOCATE( ubdy_n(jpi), vbdy_n(jpi), hbdy_n(jpi))
         ALLOCATE( ubdy_s(jpi), vbdy_s(jpi), hbdy_s(jpi))
      ENDIF

      CALL wrk_alloc( jpi, jpj, zunb, zvnb, zsshn )

      ! "Central" time index for interpolation:
      IF (ln_bt_fw) THEN
         zt = REAL(Agrif_NbStepint()+0.5_wp,wp) / zrhot
      ELSE
         zt = REAL(Agrif_NbStepint(),wp) / zrhot
      ENDIF

      ! Linear interpolation of sea level
      Agrif_SpecialValue    = 0.e0
      Agrif_UseSpecialValue = .TRUE.
      CALL Agrif_Bc_variable(zsshn, sshn_id,calledweight=zt, procname=interpsshn )
      Agrif_UseSpecialValue = .FALSE.

      ! Interpolate barotropic fluxes
      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = ln_spc_dyn

      IF (ll_int_cons) THEN ! Conservative interpolation
         CALL wrk_alloc( jpi, jpj, zuab, zvab, zubb, zvbb, zutn, zvtn )
         zuab(:,:) = 0._wp ; zvab(:,:) = 0._wp
         zubb(:,:) = 0._wp ; zvbb(:,:) = 0._wp
         zutn(:,:) = 0._wp ; zvtn(:,:) = 0._wp
         CALL Agrif_Bc_variable(zubb,unb_id ,calledweight=0._wp, procname=interpunb) ! Before
         CALL Agrif_Bc_variable(zvbb,vnb_id ,calledweight=0._wp, procname=interpvnb)
         CALL Agrif_Bc_variable(zuab,unb_id ,calledweight=1._wp, procname=interpunb) ! After
         CALL Agrif_Bc_variable(zvab,vnb_id ,calledweight=1._wp, procname=interpvnb)
         CALL Agrif_Bc_variable(zutn,ub2b_id,calledweight=1._wp, procname=interpub2b)! Time integrated
         CALL Agrif_Bc_variable(zvtn,vb2b_id,calledweight=1._wp, procname=interpvb2b)
         
         ! Time indexes bounds for integration
         zt0 = REAL(Agrif_NbStepint()  , wp) / zrhot
         zt1 = REAL(Agrif_NbStepint()+1, wp) / zrhot

         ! Polynomial interpolation coefficients:
         zaa = zrhot * (  zt1**2._wp * (       zt1 - 1._wp)        &
                 &      - zt0**2._wp * (       zt0 - 1._wp)        )
         zab = zrhot * (  zt1        * (       zt1 - 1._wp)**2._wp &
                 &      - zt0        * (       zt0 - 1._wp)**2._wp )
         zat = zrhot * (  zt1**2._wp * (-2._wp*zt1 + 3._wp)        &
                 &      - zt0**2._wp * (-2._wp*zt0 + 3._wp)        ) 

         ! Do time interpolation
         IF((nbondi == -1).OR.(nbondi == 2)) THEN
            DO jj=1,jpj
               zunb(2,jj) = zaa * zuab(2,jj) + zab * zubb(2,jj) + zat * zutn(2,jj)
               zvnb(2,jj) = zaa * zvab(2,jj) + zab * zvbb(2,jj) + zat * zvtn(2,jj)
            END DO
         ENDIF
         IF((nbondi == 1).OR.(nbondi == 2)) THEN
            DO jj=1,jpj
               zunb(nlci-2,jj) = zaa * zuab(nlci-2,jj) + zab * zubb(nlci-2,jj) + zat * zutn(nlci-2,jj)
               zvnb(nlci-1,jj) = zaa * zvab(nlci-1,jj) + zab * zvbb(nlci-1,jj) + zat * zvtn(nlci-1,jj)
            END DO
         ENDIF
         IF((nbondj == -1).OR.(nbondj == 2)) THEN
            DO ji=1,jpi
               zunb(ji,2) = zaa * zuab(ji,2) + zab * zubb(ji,2) + zat * zutn(ji,2)
               zvnb(ji,2) = zaa * zvab(ji,2) + zab * zvbb(ji,2) + zat * zvtn(ji,2)
            END DO
         ENDIF
         IF((nbondj == 1).OR.(nbondj == 2)) THEN
            DO ji=1,jpi
               zunb(ji,nlcj-1) = zaa * zuab(ji,nlcj-1) + zab * zubb(ji,nlcj-1) + zat * zutn(ji,nlcj-1)
               zvnb(ji,nlcj-2) = zaa * zvab(ji,nlcj-2) + zab * zvbb(ji,nlcj-2) + zat * zvtn(ji,nlcj-2)
            END DO
         ENDIF
         CALL wrk_dealloc( jpi, jpj, zuab, zvab, zubb, zvbb, zutn, zvtn )

      ELSE ! Linear interpolation
         zunb(:,:) = 0._wp ; zvnb(:,:) = 0._wp
         CALL Agrif_Bc_variable(zunb,unb_id,calledweight=zt, procname=interpunb)
         CALL Agrif_Bc_variable(zvnb,vnb_id,calledweight=zt, procname=interpvnb)
      ENDIF
      Agrif_UseSpecialValue = .FALSE.

      ! Fill boundary data arrays:
      IF((nbondi == -1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
               ubdy_w(jj) = (zunb(2,jj)/(zrhoy*e2u(2,jj))) * umask(2,jj,1)
               vbdy_w(jj) = (zvnb(2,jj)/(zrhox*e1v(2,jj))) * vmask(2,jj,1)
               hbdy_w(jj) = zsshn(2,jj) * tmask(2,jj,1)
         END DO
      ENDIF

      IF((nbondi == 1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
               ubdy_e(jj) = zunb(nlci-2,jj)/(zrhoy*e2u(nlci-2,jj)) * umask(nlci-2,jj,1)
               vbdy_e(jj) = zvnb(nlci-1,jj)/(zrhox*e1v(nlci-1,jj)) * vmask(nlci-1,jj,1)
               hbdy_e(jj) = zsshn(nlci-1,jj) * tmask(nlci-1,jj,1)
         END DO
      ENDIF

      IF((nbondj == -1).OR.(nbondj == 2)) THEN
         DO ji=1,jpi
               ubdy_s(ji) = zunb(ji,2)/(zrhoy*e2u(ji,2)) * umask(ji,2,1)
               vbdy_s(ji) = zvnb(ji,2)/(zrhox*e1v(ji,2)) * vmask(ji,2,1)
               hbdy_s(ji) = zsshn(ji,2) * tmask(ji,2,1)
         END DO
      ENDIF

      IF((nbondj == 1).OR.(nbondj == 2)) THEN
         DO ji=1,jpi
            ubdy_n(ji) = zunb(ji,nlcj-1)/(zrhoy*e2u(ji,nlcj-1)) * umask(ji,nlcj-1,1)
            vbdy_n(ji) = zvnb(ji,nlcj-2)/(zrhox*e1v(ji,nlcj-2)) * vmask(ji,nlcj-2,1)
            hbdy_n(ji) = zsshn(ji,nlcj-1) * tmask(ji,nlcj-1,1)
         END DO
      ENDIF

      CALL wrk_dealloc( jpi, jpj, zunb, zvnb, zsshn )

   END SUBROUTINE Agrif_dta_ts

   SUBROUTINE Agrif_ssh( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_DYN  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) ::   kt
      !!
      !!----------------------------------------------------------------------  

      IF( Agrif_Root() )   RETURN


      IF((nbondi == -1).OR.(nbondi == 2)) THEN
         ssha(2,:)=ssha(3,:)
         sshn(2,:)=sshn(3,:)
      ENDIF

      IF((nbondi == 1).OR.(nbondi == 2)) THEN
         ssha(nlci-1,:)=ssha(nlci-2,:)
         sshn(nlci-1,:)=sshn(nlci-2,:)        
      ENDIF

      IF((nbondj == -1).OR.(nbondj == 2)) THEN
         ssha(:,2)=ssha(:,3)
         sshn(:,2)=sshn(:,3)
      ENDIF

      IF((nbondj == 1).OR.(nbondj == 2)) THEN
         ssha(:,nlcj-1)=ssha(:,nlcj-2)
         sshn(:,nlcj-1)=sshn(:,nlcj-2)                
      ENDIF

   END SUBROUTINE Agrif_ssh

   SUBROUTINE Agrif_ssh_ts( jn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_ssh_ts  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) ::   jn
      !!
      INTEGER :: ji,jj
      !!----------------------------------------------------------------------  

      IF((nbondi == -1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
            ssha_e(2,jj) = hbdy_w(jj)
         END DO
      ENDIF

      IF((nbondi == 1).OR.(nbondi == 2)) THEN
         DO jj=1,jpj
            ssha_e(nlci-1,jj) = hbdy_e(jj)
         END DO
      ENDIF

      IF((nbondj == -1).OR.(nbondj == 2)) THEN
         DO ji=1,jpi
            ssha_e(ji,2) = hbdy_s(ji)
         END DO
      ENDIF

      IF((nbondj == 1).OR.(nbondj == 2)) THEN
         DO ji=1,jpi
            ssha_e(ji,nlcj-1) = hbdy_n(ji)
         END DO
      ENDIF

   END SUBROUTINE Agrif_ssh_ts

   SUBROUTINE interpsshn(tabres,i1,i2,j1,j2)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpsshn  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      !!
      INTEGER :: ji,jj
      !!----------------------------------------------------------------------  

      tabres(i1:i2,j1:j2) = sshn(i1:i2,j1:j2)

   END SUBROUTINE interpsshn

   SUBROUTINE interpu(tabres,i1,i2,j1,j2,k1,k2)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpu  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres
      !!
      INTEGER :: ji,jj,jk
      !!----------------------------------------------------------------------  

      DO jk=k1,k2
         DO jj=j1,j2
            DO ji=i1,i2
               tabres(ji,jj,jk) = e2u(ji,jj) * un(ji,jj,jk)
               tabres(ji,jj,jk) = tabres(ji,jj,jk) * fse3u_n(ji,jj,jk)
            END DO
         END DO
      END DO
   END SUBROUTINE interpu


   SUBROUTINE interpu2d(tabres,i1,i2,j1,j2)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpu2d  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      !!
      INTEGER :: ji,jj
      !!----------------------------------------------------------------------  

      DO jj=j1,j2
         DO ji=i1,i2
            tabres(ji,jj) = e2u(ji,jj) * ((gcx(ji+1,jj) - gcx(ji,jj))/e1u(ji,jj)) &
               * umask(ji,jj,1)
         END DO
      END DO

   END SUBROUTINE interpu2d


   SUBROUTINE interpv(tabres,i1,i2,j1,j2,k1,k2)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpv  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres
      !!
      INTEGER :: ji, jj, jk
      !!----------------------------------------------------------------------  

      DO jk=k1,k2
         DO jj=j1,j2
            DO ji=i1,i2
               tabres(ji,jj,jk) = e1v(ji,jj) * vn(ji,jj,jk)
               tabres(ji,jj,jk) = tabres(ji,jj,jk) * fse3v_n(ji,jj,jk)
            END DO
         END DO
      END DO

   END SUBROUTINE interpv


   SUBROUTINE interpv2d(tabres,i1,i2,j1,j2)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpu2d  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      !!
      INTEGER :: ji,jj
      !!----------------------------------------------------------------------  

      DO jj=j1,j2
         DO ji=i1,i2
            tabres(ji,jj) = e1v(ji,jj) * ((gcx(ji,jj+1) - gcx(ji,jj))/e2v(ji,jj)) &
               * vmask(ji,jj,1)
         END DO
      END DO

   END SUBROUTINE interpv2d

   SUBROUTINE interpunb(tabres,i1,i2,j1,j2)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpunb  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      !!
      INTEGER :: ji,jj
      !!----------------------------------------------------------------------  

      DO jj=j1,j2
         DO ji=i1,i2
            tabres(ji,jj) = un_b(ji,jj) * e2u(ji,jj) * hu(ji,jj) 
         END DO
      END DO

   END SUBROUTINE interpunb

   SUBROUTINE interpvnb(tabres,i1,i2,j1,j2)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpvnb  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      !!
      INTEGER :: ji,jj
      !!----------------------------------------------------------------------  

      DO jj=j1,j2
         DO ji=i1,i2
            tabres(ji,jj) = vn_b(ji,jj) * e1v(ji,jj) * hv(ji,jj)
         END DO
      END DO

   END SUBROUTINE interpvnb

   SUBROUTINE interpub2b(tabres,i1,i2,j1,j2)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpub2b  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      !!
      INTEGER :: ji,jj
      !!----------------------------------------------------------------------  

      DO jj=j1,j2
         DO ji=i1,i2
            tabres(ji,jj) = ub2_b(ji,jj) * e2u(ji,jj)
         END DO
      END DO

   END SUBROUTINE interpub2b

   SUBROUTINE interpvb2b(tabres,i1,i2,j1,j2)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpvb2b  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      !!
      INTEGER :: ji,jj
      !!----------------------------------------------------------------------  

      DO jj=j1,j2
         DO ji=i1,i2
            tabres(ji,jj) = vb2_b(ji,jj) * e1v(ji,jj)
         END DO
      END DO

   END SUBROUTINE interpvb2b

#else
   !!----------------------------------------------------------------------
   !!   Empty module                                          no AGRIF zoom
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE Agrif_OPA_Interp_empty
      WRITE(*,*)  'agrif_opa_interp : You should not have seen this print! error?'
   END SUBROUTINE Agrif_OPA_Interp_empty
#endif

   !!======================================================================
END MODULE agrif_opa_interp
