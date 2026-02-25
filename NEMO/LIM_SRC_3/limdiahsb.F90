MODULE limdiahsb
   !!======================================================================
   !!                       ***  MODULE limdia_hsb   ***
   !!  LIM-3 sea ice model :   diagnostics of ice model 
   !!======================================================================
   !! History :  3.4  ! 2012-10  (C. Rousset)  original code
   !!----------------------------------------------------------------------
#if defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3'                                       LIM3 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_dia_hsb        : computation and output of the time evolution of keys variables
   !!   lim_dia_hsb_init   : initialization and namelist read
   !!----------------------------------------------------------------------
   USE ice             ! LIM-3: sea-ice variable
   USE dom_ice         ! LIM-3: sea-ice domain
   USE dom_oce         ! ocean domain
   USE sbc_oce         ! surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: sea-ice fields
   USE daymod          ! model calendar
   USE phycst          ! physical constant
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE timing          ! preformance summary
   USE iom             ! I/O manager
   USE lib_fortran     ! glob_sum
   USE limrst          ! ice restart

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_diahsb        ! routine called by ice_step.F90

   real(wp) ::   frc_sal, frc_vol   ! global forcing trends
   real(wp) ::   bg_grme            ! global ice growth+melt trends

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.4 , NEMO Consortium (2012)
   !! $Id: limdiahsb.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE lim_diahsb
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE lim_diahsb  ***
      !!     
      !! ** Purpose: Compute the ice global heat content, salt content and volume conservation
      !!	
      !!---------------------------------------------------------------------------
      !!
      real(wp)   ::   zbg_ivo, zbg_svo, zbg_are, zbg_sal ,zbg_tem ,zbg_ihc ,zbg_shc
      real(wp)   ::   zbg_sfx, zbg_sfx_bri, zbg_sfx_bog, zbg_sfx_bom, zbg_sfx_sum, zbg_sfx_sni,   &
      &               zbg_sfx_opw, zbg_sfx_res, zbg_sfx_dyn 
      real(wp)   ::   zbg_vfx, zbg_vfx_bog, zbg_vfx_opw, zbg_vfx_sni, zbg_vfx_dyn
      real(wp)   ::   zbg_vfx_bom, zbg_vfx_sum, zbg_vfx_res, zbg_vfx_spr, zbg_vfx_snw, zbg_vfx_sub  
      real(wp)   ::   zbg_hfx_dhc, zbg_hfx_spr
      real(wp)   ::   zbg_hfx_res, zbg_hfx_sub, zbg_hfx_dyn, zbg_hfx_thd, zbg_hfx_snw, zbg_hfx_out, zbg_hfx_in   
      real(wp)   ::   zbg_hfx_sum, zbg_hfx_bom, zbg_hfx_bog, zbg_hfx_dif, zbg_hfx_opw
      real(wp)   ::   z_frc_vol, z_frc_sal, z_bg_grme 
      real(wp)   ::   z1_area                     !    -     -
      REAL(wp)   ::   ztmp
      !!---------------------------------------------------------------------------
      IF( nn_timing == 1 )   CALL timing_start('lim_diahsb')

      IF( numit == nstart ) CALL lim_diahsb_init 

      ! 1/area
      z1_area = 1._wp / MAX( glob_sum( e12t(:,:) * tmask(:,:,1) ), epsi06 )

      rswitch = MAX( 0._wp , SIGN( 1._wp , glob_sum( e12t(:,:) * tmask(:,:,1) ) - epsi06 ) )
      ! ----------------------- !
      ! 1 -  Content variations !
      ! ----------------------- !
      zbg_ivo = glob_sum( vt_i(:,:) * e12t(:,:) * tmask(:,:,1) ) ! volume ice 
      zbg_svo = glob_sum( vt_s(:,:) * e12t(:,:) * tmask(:,:,1) ) ! volume snow
      zbg_are = glob_sum( at_i(:,:) * e12t(:,:) * tmask(:,:,1) ) ! area
      zbg_sal = glob_sum( SUM( smv_i(:,:,:), dim=3 ) * e12t(:,:) * tmask(:,:,1) )       ! mean salt content
      zbg_tem = glob_sum( ( tm_i(:,:) - rt0 ) * vt_i(:,:) * e12t(:,:) * tmask(:,:,1) )  ! mean temp content

      !zbg_ihc = glob_sum( et_i(:,:) * e12t(:,:) * tmask(:,:,1) ) / MAX( zbg_ivo,epsi06 ) ! ice heat content
      !zbg_shc = glob_sum( et_s(:,:) * e12t(:,:) * tmask(:,:,1) ) / MAX( zbg_svo,epsi06 ) ! snow heat content

      ! Volume
      ztmp = rswitch * z1_area * r1_rau0 * rday
      zbg_vfx     = ztmp * glob_sum(     emp(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_vfx_bog = ztmp * glob_sum( wfx_bog(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_vfx_opw = ztmp * glob_sum( wfx_opw(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_vfx_sni = ztmp * glob_sum( wfx_sni(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_vfx_dyn = ztmp * glob_sum( wfx_dyn(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_vfx_bom = ztmp * glob_sum( wfx_bom(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_vfx_sum = ztmp * glob_sum( wfx_sum(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_vfx_res = ztmp * glob_sum( wfx_res(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_vfx_spr = ztmp * glob_sum( wfx_spr(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_vfx_snw = ztmp * glob_sum( wfx_snw(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_vfx_sub = ztmp * glob_sum( wfx_sub(:,:) * e12t(:,:) * tmask(:,:,1) )

      ! Salt
      zbg_sfx     = ztmp * glob_sum(     sfx(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_sfx_bri = ztmp * glob_sum( sfx_bri(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_sfx_res = ztmp * glob_sum( sfx_res(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_sfx_dyn = ztmp * glob_sum( sfx_dyn(:,:) * e12t(:,:) * tmask(:,:,1) )

      zbg_sfx_bog = ztmp * glob_sum( sfx_bog(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_sfx_opw = ztmp * glob_sum( sfx_opw(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_sfx_sni = ztmp * glob_sum( sfx_sni(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_sfx_bom = ztmp * glob_sum( sfx_bom(:,:) * e12t(:,:) * tmask(:,:,1) )
      zbg_sfx_sum = ztmp * glob_sum( sfx_sum(:,:) * e12t(:,:) * tmask(:,:,1) )

      ! Heat budget
      zbg_ihc      = glob_sum( et_i(:,:) * e12t(:,:) * 1.e-20 ) ! ice heat content  [1.e20 J]
      zbg_shc      = glob_sum( et_s(:,:) * e12t(:,:) * 1.e-20 ) ! snow heat content [1.e20 J]
      zbg_hfx_dhc  = glob_sum( diag_heat(:,:) * e12t(:,:) * tmask(:,:,1) ) ! [in W]
      zbg_hfx_spr  = glob_sum( hfx_spr(:,:) * e12t(:,:) * tmask(:,:,1) ) ! [in W]

      zbg_hfx_thd  = glob_sum( hfx_thd(:,:) * e12t(:,:) * tmask(:,:,1) ) ! [in W]
      zbg_hfx_dyn  = glob_sum( hfx_dyn(:,:) * e12t(:,:) * tmask(:,:,1) ) ! [in W]
      zbg_hfx_res  = glob_sum( hfx_res(:,:) * e12t(:,:) * tmask(:,:,1) ) ! [in W]
      zbg_hfx_sub  = glob_sum( hfx_sub(:,:) * e12t(:,:) * tmask(:,:,1) ) ! [in W]
      zbg_hfx_snw  = glob_sum( hfx_snw(:,:) * e12t(:,:) * tmask(:,:,1) ) ! [in W]
      zbg_hfx_sum  = glob_sum( hfx_sum(:,:) * e12t(:,:) * tmask(:,:,1) ) ! [in W]
      zbg_hfx_bom  = glob_sum( hfx_bom(:,:) * e12t(:,:) * tmask(:,:,1) ) ! [in W]
      zbg_hfx_bog  = glob_sum( hfx_bog(:,:) * e12t(:,:) * tmask(:,:,1) ) ! [in W]
      zbg_hfx_dif  = glob_sum( hfx_dif(:,:) * e12t(:,:) * tmask(:,:,1) ) ! [in W]
      zbg_hfx_opw  = glob_sum( hfx_opw(:,:) * e12t(:,:) * tmask(:,:,1) ) ! [in W]
      zbg_hfx_out  = glob_sum( hfx_out(:,:) * e12t(:,:) * tmask(:,:,1) ) ! [in W]
      zbg_hfx_in   = glob_sum(  hfx_in(:,:) * e12t(:,:) * tmask(:,:,1) ) ! [in W]
    
      ! --------------------------------------------- !
      ! 2 - Trends due to forcing and ice growth/melt !
      ! --------------------------------------------- !
      z_frc_vol = r1_rau0 * glob_sum( - emp(:,:) * e12t(:,:) * tmask(:,:,1) ) ! volume fluxes
      z_frc_sal = r1_rau0 * glob_sum(   sfx(:,:) * e12t(:,:) * tmask(:,:,1) ) ! salt fluxes
      z_bg_grme = glob_sum( - ( wfx_bog(:,:) + wfx_opw(:,:) + wfx_sni(:,:) + wfx_dyn(:,:) + &
                          &     wfx_bom(:,:) + wfx_sum(:,:) + wfx_res(:,:) + wfx_snw(:,:) + &
                          &     wfx_sub(:,:) ) * e12t(:,:) * tmask(:,:,1) ) ! volume fluxes
      !
      frc_vol  = frc_vol  + z_frc_vol  * rdt_ice
      frc_sal  = frc_sal  + z_frc_sal  * rdt_ice
      bg_grme  = bg_grme  + z_bg_grme  * rdt_ice
      
      ! difference
      !frc_vol = zbg_ivo - frc_vol
      !frc_sal = zbg_sal - frc_sal
      
      ! ----------------------- !
      ! 3 - Diagnostics writing !
      ! ----------------------- !
      rswitch = MAX( 0._wp , SIGN( 1._wp , zbg_ivo - epsi06 ) )
      !
      IF( iom_use('ibgvoltot') )   &
      CALL iom_put( 'ibgvoltot' , zbg_ivo * rhoic * r1_rau0 * 1.e-9        )   ! ice volume (km3 equivalent liquid)         
      IF( iom_use('sbgvoltot') )   &
      CALL iom_put( 'sbgvoltot' , zbg_svo * rhosn * r1_rau0 * 1.e-9        )   ! snw volume (km3 equivalent liquid)       
      IF( iom_use('ibgarea') )   &
      CALL iom_put( 'ibgarea'   , zbg_are * 1.e-6                          )   ! ice area   (km2)
      IF( iom_use('ibgsaline') )   &
      CALL iom_put( 'ibgsaline' , rswitch * zbg_sal / MAX( zbg_ivo, epsi06 ) )   ! ice saline (psu)
      IF( iom_use('ibgtemper') )   &
      CALL iom_put( 'ibgtemper' , rswitch * zbg_tem / MAX( zbg_ivo, epsi06 ) )   ! ice temper (C)
      CALL iom_put( 'ibgheatco' , zbg_ihc                                  )   ! ice heat content (1.e20 J)        
      CALL iom_put( 'sbgheatco' , zbg_shc                                  )   ! snw heat content (1.e20 J)
      IF( iom_use('ibgsaltco') )   &
      CALL iom_put( 'ibgsaltco' , zbg_sal * rhoic * r1_rau0 * 1.e-9        )   ! ice salt content (psu*km3 equivalent liquid)        

      CALL iom_put( 'ibgvfx'    , zbg_vfx                                  )   ! volume flux emp (m/day liquid)
      CALL iom_put( 'ibgvfxbog' , zbg_vfx_bog                              )   ! volume flux bottom growth     -(m/day equivalent liquid)
      CALL iom_put( 'ibgvfxopw' , zbg_vfx_opw                              )   ! volume flux open water growth -
      CALL iom_put( 'ibgvfxsni' , zbg_vfx_sni                              )   ! volume flux snow ice growth   -
      CALL iom_put( 'ibgvfxdyn' , zbg_vfx_dyn                              )   ! volume flux dynamic growth    -
      CALL iom_put( 'ibgvfxbom' , zbg_vfx_bom                              )   ! volume flux bottom melt       -
      CALL iom_put( 'ibgvfxsum' , zbg_vfx_sum                              )   ! volume flux surface melt      -
      CALL iom_put( 'ibgvfxres' , zbg_vfx_res                              )   ! volume flux resultant         -
      CALL iom_put( 'ibgvfxspr' , zbg_vfx_spr                              )   ! volume flux from snow precip         -
      CALL iom_put( 'ibgvfxsnw' , zbg_vfx_snw                              )   ! volume flux from snow melt         -
      CALL iom_put( 'ibgvfxsub' , zbg_vfx_sub                              )   ! volume flux from sublimation         -
          
      CALL iom_put( 'ibgsfx'    , zbg_sfx                                  )   ! salt flux         -(psu*m/day equivalent liquid)       
      CALL iom_put( 'ibgsfxbri' , zbg_sfx_bri                              )   ! salt flux brines  -      
      CALL iom_put( 'ibgsfxdyn' , zbg_sfx_dyn                              )   ! salt flux dynamic -    
      CALL iom_put( 'ibgsfxres' , zbg_sfx_res                              )   ! salt flux result  -    
      CALL iom_put( 'ibgsfxbog' , zbg_sfx_bog                              )   ! salt flux bottom growth   
      CALL iom_put( 'ibgsfxopw' , zbg_sfx_opw                              )   ! salt flux open water growth -
      CALL iom_put( 'ibgsfxsni' , zbg_sfx_sni                              )   ! salt flux snow ice growth   -
      CALL iom_put( 'ibgsfxbom' , zbg_sfx_bom                              )   ! salt flux bottom melt       -
      CALL iom_put( 'ibgsfxsum' , zbg_sfx_sum                              )   ! salt flux surface melt      -

      CALL iom_put( 'ibghfxdhc' , zbg_hfx_dhc                              )   ! Heat content variation in snow and ice [W]
      CALL iom_put( 'ibghfxspr' , zbg_hfx_spr                              )   ! Heat content of snow precip [W]

      CALL iom_put( 'ibghfxres' , zbg_hfx_res                              )   ! 
      CALL iom_put( 'ibghfxsub' , zbg_hfx_sub                              )   ! 
      CALL iom_put( 'ibghfxdyn' , zbg_hfx_dyn                              )   ! 
      CALL iom_put( 'ibghfxthd' , zbg_hfx_thd                              )   ! 
      CALL iom_put( 'ibghfxsnw' , zbg_hfx_snw                              )   ! 
      CALL iom_put( 'ibghfxsum' , zbg_hfx_sum                              )   ! 
      CALL iom_put( 'ibghfxbom' , zbg_hfx_bom                              )   ! 
      CALL iom_put( 'ibghfxbog' , zbg_hfx_bog                              )   ! 
      CALL iom_put( 'ibghfxdif' , zbg_hfx_dif                              )   ! 
      CALL iom_put( 'ibghfxopw' , zbg_hfx_opw                              )   ! 
      CALL iom_put( 'ibghfxout' , zbg_hfx_out                              )   ! 
      CALL iom_put( 'ibghfxin'  , zbg_hfx_in                               )   ! 

      CALL iom_put( 'ibgfrcvol' , frc_vol * 1.e-9                          )   ! vol - forcing     (km3 equivalent liquid) 
      CALL iom_put( 'ibgfrcsfx' , frc_sal * 1.e-9                          )   ! sal - forcing     (psu*km3 equivalent liquid)   
      IF( iom_use('ibgvolgrm') )   &
      CALL iom_put( 'ibgvolgrm' , bg_grme * r1_rau0 * 1.e-9                )   ! vol growth + melt (km3 equivalent liquid)         

      !
      IF( lrst_ice )   CALL lim_diahsb_rst( numit, 'WRITE' )
      !
      IF( nn_timing == 1 )   CALL timing_stop('lim_diahsb')
!
   END SUBROUTINE lim_diahsb


   SUBROUTINE lim_diahsb_init
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE lim_diahsb_init  ***
      !!     
      !! ** Purpose: Initialization for the heat salt volume budgets
      !!	
      !! ** Method : Compute initial heat content, salt content and volume
      !!
      !! ** Action : - Compute initial heat content, salt content and volume
      !!             - Initialize forcing trends
      !!             - Compute coefficients for conversion
      !!---------------------------------------------------------------------------
      INTEGER            ::   jk       ! dummy loop indice
      INTEGER            ::   ierror   ! local integer
      !!
      !!NAMELIST/namicehsb/ blabla
      !!----------------------------------------------------------------------
      !
      !!REWIND ( numnam_ice )              ! Read Namelist namicehsb 
      !!READ   ( numnam_ice, namicehsb )
      !
      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'lim_diahsb_init : check the heat and salt budgets'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      CALL lim_diahsb_rst( nstart, 'READ' )  !* read or initialize all required files
      !
   END SUBROUTINE lim_diahsb_init

   SUBROUTINE lim_diahsb_rst( kt, cdrw )
     !!---------------------------------------------------------------------
     !!                   ***  ROUTINE limdia_rst  ***
     !!                     
     !! ** Purpose :   Read or write DIA file in restart file
     !!
     !! ** Method  :   use of IOM library
     !!----------------------------------------------------------------------
     INTEGER         , INTENT(in) ::   kt     ! ice time-step
     CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
     !
     INTEGER ::   id1, id2, id3   ! local integers
     !!----------------------------------------------------------------------
     !
     IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise 
        IF( ln_rstart ) THEN                   !* Read the restart file
           !id1 = iom_varid( numrir, 'frc_vol'  , ldstop = .TRUE. )
           !
           IF(lwp) WRITE(numout,*) '~~~~~~~'
           IF(lwp) WRITE(numout,*) ' lim_diahsb_rst at it= ', kt,' date= ', ndastp
           IF(lwp) WRITE(numout,*) '~~~~~~~'
           CALL iom_get( numrir, 'frc_vol', frc_vol )
           CALL iom_get( numrir, 'frc_sal', frc_sal )
           CALL iom_get( numrir, 'bg_grme', bg_grme )
        ELSE
           IF(lwp) WRITE(numout,*) '~~~~~~~'
           IF(lwp) WRITE(numout,*) ' lim_diahsb at initial state '
           IF(lwp) WRITE(numout,*) '~~~~~~~'
           frc_vol  = 0._wp                                          
           frc_sal  = 0._wp                                                 
           bg_grme  = 0._wp                                       
       ENDIF

     ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
        !                                   ! -------------------
        IF(lwp) WRITE(numout,*) '~~~~~~~'
        IF(lwp) WRITE(numout,*) ' lim_diahsb_rst at it= ', kt,' date= ', ndastp
        IF(lwp) WRITE(numout,*) '~~~~~~~'
        CALL iom_rstput( kt, nitrst, numriw, 'frc_vol'   , frc_vol     )
        CALL iom_rstput( kt, nitrst, numriw, 'frc_sal'   , frc_sal     )
        CALL iom_rstput( kt, nitrst, numriw, 'bg_grme'   , bg_grme     )
        !
     ENDIF
     !
   END SUBROUTINE lim_diahsb_rst
 
#else
   !!----------------------------------------------------------------------
   !!   Default option :         Empty module          NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_diahsb          ! Empty routine
   END SUBROUTINE lim_diahsb
#endif
   !!======================================================================
END MODULE limdiahsb
