      SUBROUTINE bgcbio(kt)
#if defined key_planktom && defined key_top
!!======================================================================
!!                         ***  ROUTINE bgcbio  ***
!! TOP : PlankTOM ecosystem model 
!!======================================================================
!! History : PISCES       ! 2001 (O. Aumont) Original code 
!!           PlankTOM4    ! 2002 (E.T. Buitenhuis) improved zooplankton 
!!           PlankTOM5    ! 2006 (E.T. Buitenhuis, M. Vogt, M. Manizza)
!!                               5 Plankton Functional Types (PFTs)
!!           PlankTOM10   ! 2016 (E.T. Buitenhuis, C. Le Quéré) 10 PFTs
!!           PlankTOM12.1 ! 2024 (E.T. Buitenhuis, R. Wright, J. Guest) 
!!                               12 PFTs consolidated 
!!            -           ! 2025-01 (E.T. Buitenhuis, C. Le Quéré) 
!!                               Tidy up code and documentation
!!----------------------------------------------------------------------
!!    'key_c14b'          :       calculate 14C fluxes and concentration
!!    'key_trc_n2o'       :                        calculates N2O fluxes
!!    'key_trc_ch4'       :                        calculates CH4 fluxes
!!    'key_trc_dms'       :            includes the DMS cycle and fluxes
!!    'key_trc_foram'     :      includes zooplankton calcifiers (foram) 
!!    'key_trc_piic'      :    additional tracer with pre-industrial DIC 
!!    'key_trc_diaadd'    :       save tracer diagnostic files 2D and 3D
!!                                         (to use with key_iomput only)
!!    'key_iomput'        :                              use IOM library 
!!----------------------------------------------------------------------
!! ** Purpose : computes all ecosystem sources and sinks 
!! ** Action  : - calls : bgcpro
!!                        bgclos
!!                        bgcsnk
!!                        bgcnul
!!              - revises growth and loss variables after nul tests
!!              - add ecosystem sources and sinks to tracers  
!!
!! References : Buitenhuis et al. 2019 doi:10.1029/2018GB006110
!!              Le Quéré et al., Biogeosciences, 2016
!!              PlankTOM12 manual, Buitenhuis et al. 2023
!!                                 https://zenodo.org/records/8388158 
!! ---------------------------------------------------------------------
!
! parameters and commons
      USE trc
      USE trp_trc
      USE sms_planktom
      USE oce_trc
      USE lbclnk
      USE lib_mpp
      USE iom
      USE traqsr
#  if defined key_c14b
      USE trcsms_c14b , ONLY : rn_bfrd14
#  endif
      IMPLICIT NONE
! 
! local variables
      INTEGER ji, jj, jk, ikt, jl, jm, jn, iostat, kt, level
      REAL consum, consump, gongoc, dgongoc, mfecal, dumch4, temch4
      REAL omeara, omecal, remara, remco3
      REAL zgochoc(2), zgrdmp
#  if defined key_trc_n2o
      REAL(wp) dumn2s, lo2n2s
#  endif
#  include "domzgr_substitute.h90"
! 
! Initialisation of variables
! ----------------------------
! 
      grazing = 0.
      ppint   = 0.
      ppt     = 0.
      pptdoc  = 0.
      tchl    = 0.
      trophic = 0.
      out2d   = 0.
      out3d   = 0.
      zgochoc = 0.
!
      DO jk = 1, jpk-1
       DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
!
! Iron and Si deposition at the surface
!--------------------------------------
! The dust variable is in kgdust/m2/s (from Jickells et al. 2005)  
! We use 0.035gFe/gdust, 0.308 gSi/gdust (or 8.8gSi/gFe)  
! The solubility of Fe in dust is rn_fersol (usually 2%, set in namelist)  
! The solubility of Si in dust is 7.5%  
! Variables irondep and sidep are in mol/L/time_step for Fe and Si, respectively 
! kg/m3=g/L
! The molecular weight of Fe and Si is 55.85 and 28.01, respectively 
! 0.3314 is the correction so that the total iron added is the same as when
! iron was only added at the surface (0.80037 kmol/s), and is approximately
! the integral of
! (gdept(jk)/gdept(1))**(-0.858) over 30 layers
!
            irondep(ji,jj,jk) = dust(ji,jj)*rn_fersol*rfact             &
     &        *(gdept_0(ji,jj,jk)/gdept_0(ji,jj,1))**(-0.858)*0.3314    &
     &        /(fse3t(ji,jj,jk)*55.85) 
            sidep(ji,jj,jk) = irondep(ji,jj,jk)*                        &
     &                (0.308/0.035)*(rn_silsol/rn_fersol)* (55.85/28.01)
        END DO
       END DO
      END DO
!!
!! -----------------------------------------
!! Call biological production and loss terms 
!! -----------------------------------------
!!
! Call optical /production routine to compute phytoplankton growth rate
! over the global ocean  
! Growth rates for each element are computed (C, Si, Fe, Chl)
! -------------------------------------------------------------------
!
      IF ( kt .EQ. nit000 ) CALL tra_qsr(kt)
      CALL bgcpro(kt)
!
! Calculate loss rates
! -------------------------------------------------------------------
!     
      CALL bgclos
!
! Call subroutine for computation of the vertical flux 
! of particulate organic matter
! -------------------------------------------------------------------
!
      CALL bgcsnk
!
! Pre-compute concentration changes of the (rapidly
! varying) tracers for preventing them to fall below 0
! -------------------------------------------------------------------
!
      CALL bgcnul
!! 
!! --------------------------------------------------------------------
!! Recompute the sources ans sinks (SMS) related to zooplankton grazing
!! --------------------------------------------------------------------
!! 
      grazoc = 0.
      grazof = 0.
!
      DO jk = 1, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
!
! assign grazing to organic matter. For each food source, sum the ingestion of each zooplankton
! variable grazoc is for carbon, grazof is for iron  
            DO jn = 1, jpfoo
              jm = grizoo(1,jn)
              jl = grizoo(2,jn)
              grazoc(ji,jj,jk,jm) = grazoc(ji,jj,jk,jm)+grazoo(ji,jj,jk,jm,jl)
              grazof(ji,jj,jk,jm) = grazof(ji,jj,jk,jm)+grazoo(ji,jj,jk,jm,jl)*stofoo(ji,jj,jk,jl,2)
              level = min(max(int(float(jl-jpbac)/float(jpzft)+0.81),0),2)+1  !  1 = detriti, 2 = carni, 3 = herbivory
              trophic(ji,jj,jk,level) = trophic(ji,jj,jk,level)+grazoo(ji,jj,jk,jm,jl)
            END DO
!
! calculate GGE (gross gross efficiency, mgezoo) for each zooplankton
! calculate remineralisation. grarem is for carbon, grafer is for iron. 
! calculate fecal pellet production to sinking material grapoc 
            DO jm = 1, jpzft
              mgezoo(ji,jj,jk,jm) = min(1.-rn_unazoo(jm),rn_ggezoo(jm)     &
     &          +reszoo(ji,jj,jk,jm)/max(rtrn,grazoc(ji,jj,jk,jm)),        &
     &          grazof(ji,jj,jk,jm)*(1.-rn_unazoo(jm))                     &
     &          /max(grazoc(ji,jj,jk,jm)*ferat3,minfer))
              grafer(ji,jj,jk,jm) = grazof(ji,jj,jk,jm)*(1.-rn_unazoo(jm)) &
     &          -ferat3*mgezoo(ji,jj,jk,jm)*grazoc(ji,jj,jk,jm)
              grarem(ji,jj,jk,jm) = grazoc(ji,jj,jk,jm)                    &
     &          *(1.-mgezoo(ji,jj,jk,jm)-rn_unazoo(jm))
              grapoc(ji,jj,jk,jm) = grazoc(ji,jj,jk,jm)*rn_unazoo(jm)
            END DO
!
! calculate the bacterial growth efficiency
            bactge(ji,jj,jk) = min(bactge(ji,jj,jk),                    &
     &        (ubafer(ji,jj,jk)+remsfe(ji,jj,jk)+rembfe(ji,jj,jk)       &
     &        +remufe(ji,jj,jk))/max((remdoc(ji,jj,jk)+rempoc(ji,jj,jk) &
     &        +remgoc(ji,jj,jk)+remhoc(ji,jj,jk))*ferat3,minfer))
!
#  if defined key_trc_dms
! recalculate prodms, prodmd:
            DO jl = jpdia, jpdia+jppft-1
              zgrdmp = 0.
              DO jm = 1, jpzft
                zgrdmp = zgrdmp+grazoo(ji,jj,jk,jm,jl)                  &
     &            *(1.-mgezoo(ji,jj,jk,jm)-rn_assdms(jm)*rn_unazoo(jm))
              END DO
!
              prodmd(ji,jj,jk) = prodmd(ji,jj,jk)+((1.-rn_rdddms)*      &
     &          zgrdmp)*rphdmd(ji,jj,jk,jl)
!
              prodms(ji,jj,jk) = prodms(ji,jj,jk)+(rn_rdddms*zgrdmp+    &
     &          rn_xpldmd(jl)*rfact/rjjss*etot(ji,jj,jk)/rn_etomax*     &
     &          trn(ji,jj,jk,jl))*rphdmd(ji,jj,jk,jl)
            END DO
#  endif
!!     
!! -------------------------------------------------------
!! Determination of tracer concentrations as a function of 
!! biological sources and sinks
!! -------------------------------------------------------
!!
! calculate the sum of all the PFT respiration terms (bacteria, phytoplankton, zooplankton)
            consum = (1.-bactge(ji,jj,jk))*(remdoc(ji,jj,jk)            &
     &        +rempoc(ji,jj,jk)+remhoc(ji,jj,jk))+resbac(ji,jj,jk)/3.0
            DO jl = jpdia, jpdia+jppft-1
              consum = consum +resphy(ji,jj,jk,jl,1)
            END DO
            DO jm = 1, jpzft
              consum = consum+grarem(ji,jj,jk,jm)*rn_sigzoo(jm)         &
     &          +reszoo(ji,jj,jk,jm)
            END DO
!
! calculate the regeneration of PO4 (including all fixed N:C ratio terms)
            consump = consum+remgon(ji,jj,jk)*ratc2n-bactge(ji,jj,jk)*remgoc(ji,jj,jk)
!
! add variable N:C ratio to respiration term
            consum = consum+(1.-bactge(ji,jj,jk))*remgoc(ji,jj,jk)
            delo2(ji,jj,jk) = min(rato2c*consum,delo2(ji,jj,jk))
!
! calculate denitrification
            denitr(ji,jj,jk) = 0.8*(rato2c*consum-delo2(ji,jj,jk))
!
! calculate the N2O production
!
!            prbn2s = denitr(ji,jj,jk)*rn_denn2o
!            trn(ji,jj,jk,jpn2o)= trn(ji,jj,jk,jpn2o)+prbn2s
!            trc3d(ji,jj,jk,17) = prbn2s*rfactr*1e3
!            denno3=max(-0.5,(rn_doxn2o-trno2(ji,jj,jk))     &
!       &      /(2.*rn_doxn2o+trno2(ji,jj,jk)))
!            denno3=(sin(denno3*rpi)+1.)/2.
!            degn2o=min(trn(ji,jj,jk,jpn2o), &
!       &      0.8*rato2c*consum*denno3*rn_degn2o)
!            trn(ji,jj,jk,jpn2o)= trn(ji,jj,jk,jpn2o)-degn2o
!            trc3d(ji,jj,jk,18) = degn2o*rfactr*1e3
!!
!! ---------------------------------------------------------
!! Apply the sources and sinks terms to the tracer evolution
!! ---------------------------------------------------------
!!
! Evolution of O2
! ---------------
!
! add to oxygen the consumption by nitrogen fixation, which follows the reaction: 
! 0.5 N2 + 1.25 O2 + 0.5 H2O -> HNO3
            trn(ji,jj,jk,jpoxy) = trn(ji,jj,jk,jpoxy)-delo2(ji,jj,jk)    &
     &        -(prophy(ji,jj,jk,jpfix,1)+docphy(ji,jj,jk,jpfix))*ratn2c* &
     &        (1.-dinpft(ji,jj,jk,jpfix))*1.25
!
            prodt(ji,jj,jk) = 0.
!
! Evolution of phytoplankton
! --------------------------
!
! loop over phytoplankton types
            DO jl = jpdia, jpdia+jppft-1
!
              trn(ji,jj,jk,jppo4) = trn(ji,jj,jk,jppo4)                &
     &          -prophy(ji,jj,jk,jl,1)-docphy(ji,jj,jk,jl)
!
              trn(ji,jj,jk,jpdin) = trn(ji,jj,jk,jpdin)                &
     &          -(prophy(ji,jj,jk,jl,1)+docphy(ji,jj,jk,jl))           &
     &          *ratn2c*dinpft(ji,jj,jk,jl)
!
! calculate productivity (prodt) to add to all phytoplankton
              prodt(ji,jj,jk) = prodt(ji,jj,jk)+prophy(ji,jj,jk,jl,1)  &
     &          +docphy(ji,jj,jk,jl)
! create ppint and pptdoc for diagnostics
              ppint(ji,jj) = ppint(ji,jj)+prophy(ji,jj,jk,jl,1)*fse3t(ji,jj,jk)*1e3*rfactr
              pptdoc(ji,jj,jk) = pptdoc(ji,jj,jk)+docphy(ji,jj,jk,jl)*1e3*rfactr
!
! loop over elements (C, Fe, Chl)
              DO jm = 0, 2
!
                trn(ji,jj,jk,jl+jm*jppft) = trn(ji,jj,jk,jl+jm*jppft)  &
     &            +prophy(ji,jj,jk,jl,jm+1)                            &
     &            -resphy(ji,jj,jk,jl,jm+1)
                DO jn = 1, jpzft
                trn(ji,jj,jk,jl+jm*jppft) = trn(ji,jj,jk,jl+jm*jppft)  &
     &            -grazoo(ji,jj,jk,jn,jl)*stofoo(ji,jj,jk,jl,jm+1)
                END DO
              END DO
!
              trn(ji,jj,jk,jpdoc) = trn(ji,jj,jk,jpdoc)+docphy(ji,jj,jk,jl)
!
              trn(ji,jj,jk,jpoxy) = trn(ji,jj,jk,jpoxy)                &
     &          +rato2c*(prophy(ji,jj,jk,jl,1)+docphy(ji,jj,jk,jl))
!
              trn(ji,jj,jk,jpfer) = trn(ji,jj,jk,jpfer)                &
     &          -prophy(ji,jj,jk,jl,2)+resphy(ji,jj,jk,jl,2)
            END DO
!
! Evolution of Calcite
! --------------------
!
! production of attached CaCO3 (no tracer, but subtracted from DIC and Alk)
! cal is for calcite, ara is for aragonite 
            prcaca(ji,jj,jk) =                                          &
#    if defined key_trc_foram
     &        rn_forcal*mgezoo(ji,jj,jk,2)*grazoc(ji,jj,jk,2)+          &
#    endif
     &        rn_coccal*prophy(ji,jj,jk,jpcoc,1)
#    if defined key_trc_foram
            proara(ji,jj,jk) = rn_pteara*mgezoo(ji,jj,jk,3)             &
     &         *grazoc(ji,jj,jk,3)
!
! transfer of attached CaCO3 to sinking CaCO3 by loss processes
            losara(ji,jj,jk) = rn_pteara*reszoo(ji,jj,jk,3)
            losfor(ji,jj,jk) = rn_forcal*reszoo(ji,jj,jk,2)
            DO jm = 3, jpzft
              losfor(ji,jj,jk) = losfor(ji,jj,jk)+rn_forcal*grazoo(ji,jj,jk,jm,jpfor)
            END DO
            trn(ji,jj,jk,jpcal) = trn(ji,jj,jk,jpcal)+(1.-rn_disfor)*losfor(ji,jj,jk)
            DO jm = 4, jpzft
              losara(ji,jj,jk) = losara(ji,jj,jk)+rn_pteara*grazoo(ji,jj,jk,jm,jppte)
            END DO
#    else
            proara(ji,jj,jk) = rn_pteara*mgezoo(ji,jj,jk,2)             &
     &         *grazoc(ji,jj,jk,2)
!
! transfer of attached CaCO3 to sinking CaCO3 by loss processes
            losara(ji,jj,jk) = rn_pteara*reszoo(ji,jj,jk,2)
            DO jm = 3, jpzft
              losara(ji,jj,jk) = losara(ji,jj,jk)+rn_pteara*grazoo(ji,jj,jk,jm,jppte)
            END DO
#    endif
            trn(ji,jj,jk,jpara) = trn(ji,jj,jk,jpara)+(1.-rn_disara)*losara(ji,jj,jk)  &
     &        +(snkara(ji,jj,jk)-snkara(ji,jj,jk+1))/fse3t(ji,jj,jk)
            trn(ji,jj,jk,jpdic) = trn(ji,jj,jk,jpdic)-proara(ji,jj,jk)  &
#    if defined key_trc_foram
     &        +rn_disfor*losfor(ji,jj,jk)                               &
#    endif
     &        +rn_disara*losara(ji,jj,jk)
#    if defined key_trc_piic
            trn(ji,jj,jk,jppiic) = trn(ji,jj,jk,jppiic)-proara(ji,jj,jk)  &
#      if defined key_trc_foram
     &        +rn_disfor*losfor(ji,jj,jk)                               &
#      endif
     &        +rn_disara*losara(ji,jj,jk)
#    endif
            trn(ji,jj,jk,jptal) = trn(ji,jj,jk,jptal)-2.*proara(ji,jj,jk) &
#    if defined key_trc_foram
     &        +2.*rn_disfor*losfor(ji,jj,jk)                            &
#    endif
     &        +2.*rn_disara*losara(ji,jj,jk)
!
! transfer of attached CaCO3 to sinking CaCO3 by loss processes
            loscal(ji,jj,jk) = rn_coccal*resphy(ji,jj,jk,jpcoc,1)
            DO jm = 1, jpzft
              loscal(ji,jj,jk) = loscal(ji,jj,jk)+rn_coccal*grazoo(ji,jj,jk,jm,jpcoc)
            END DO
!
            trn(ji,jj,jk,jpcal) = trn(ji,jj,jk,jpcal)+(1.-rn_discal)*loscal(ji,jj,jk) &
     &        +(snkcal(ji,jj,jk)-snkcal(ji,jj,jk+1))/fse3t(ji,jj,jk)
#  if defined key_c14b
            d14pro(ji,jj,jk) = rn_bfrd14/10.*prodt(ji,jj,jk)/max(trn(ji,jj,jk,jpdic),rtrn)
            d14res(ji,jj,jk) = (trn(ji,jj,jk,jpd14)-(trn(ji,jj,1,jpd14)-rn_bfrd14/10.)) &
     &        *consum/max(trn(ji,jj,jk,jpdic),rtrn)
            trn(ji,jj,jk,jpd14) = trn(ji,jj,jk,jpd14)+(d14pro(ji,jj,jk)-d14res(ji,jj,jk))*tmask(ji,jj,jk)
#  endif
!
! Evolution of PO4
! ----------------
!
            trn(ji,jj,jk,jppo4) = trn(ji,jj,jk,jppo4)+consump           &
     &        +deppo4(ji,jj,jk)
!
            trn(ji,jj,jk,jpdin) = trn(ji,jj,jk,jpdin)-denitr(ji,jj,jk) &
     &        +consump*ratn2c+depnit(ji,jj,jk)+atmdin(ji,jj,jk)
!
            prodt(ji,jj,jk) = prodt(ji,jj,jk)-consum
!    
! Evolution of zooplankton
! ------------------------
!     
            DO jm = 1, jpzft
              trn(ji,jj,jk,jpbac+jm) = trn(ji,jj,jk,jpbac+jm)           &
     &      -reszoo(ji,jj,jk,jm)+mgezoo(ji,jj,jk,jm)*grazoc(ji,jj,jk,jm)
            END DO
#  if defined key_trc_foram
            DO jm = 0, jpzft-1
              DO jn = jm+1, jpzft
#  else
            DO jm = 0, jpzft
              DO jn = 1, jpzft
#  endif
                trn(ji,jj,jk,jpbac+jm) = trn(ji,jj,jk,jpbac+jm)         &
     &            -grazoo(ji,jj,jk,jn,jpbac+jm)
              END DO
            END DO
!
! Evolution of Macrozooplankton (crustaceans and jellyfish)
! ---------------------------------------------------------
!
            trn(ji,jj,jk,jpcru) = trn(ji,jj,jk,jpcru)-torcru(ji,jj,jk)
#  if ! defined key_trc_foram
            trn(ji,jj,jk,jpgel) = trn(ji,jj,jk,jpgel)-torgel(ji,jj,jk)
#  endif
!
! Evolution of Bacteria
! ---------------------
!
            trn(ji,jj,jk,jpbac) = trn(ji,jj,jk,jpbac)+bactge(ji,jj,jk)*&
     &        (remdoc(ji,jj,jk)+rempoc(ji,jj,jk)+remgoc(ji,jj,jk)       &
     &        +remhoc(ji,jj,jk))-resbac(ji,jj,jk)
!    
! Evolution of DOC and particles 
! ------------------------------
!     
            trn(ji,jj,jk,jpdoc) = trn(ji,jj,jk,jpdoc)-remdoc(ji,jj,jk)  &
     &        +resbac(ji,jj,jk)/3.0  &
     &        -xaggdoc(ji,jj,jk)-xaggdoc2(ji,jj,jk)+depdoc(ji,jj,jk)    
            gongoc = trn(ji,jj,jk,jpgon)/max(trn(ji,jj,jk,jpgoc),rtrn)
            DO jm = 1, jpzft
              trn(ji,jj,jk,jpdoc) = trn(ji,jj,jk,jpdoc)                 &
     &          +grarem(ji,jj,jk,jm)*(1.-rn_sigzoo(jm))
              trn(ji,jj,jk,jppoc) = trn(ji,jj,jk,jppoc)                 &
     &          -grazoo(ji,jj,jk,jm,jppoc)
              trn(ji,jj,jk,jpgoc) = trn(ji,jj,jk,jpgoc)                 &
     &          -grazoo(ji,jj,jk,jm,jpgoc)
              trn(ji,jj,jk,jpgon) = trn(ji,jj,jk,jpgon)                 &
     &          -grazoo(ji,jj,jk,jm,jpgoc)*gongoc
              trn(ji,jj,jk,jphoc) = trn(ji,jj,jk,jphoc)                 &
     &          -grazoo(ji,jj,jk,jm,jphoc)
              trn(ji,jj,jk,jppoc+nn_sizzoo(jm))=                        &
     &          trn(ji,jj,jk,jppoc+nn_sizzoo(jm))+grapoc(ji,jj,jk,jm)
              trn(ji,jj,jk,jpfer) = trn(ji,jj,jk,jpfer)                 &
     &          +grafer(ji,jj,jk,jm)+ferat3*reszoo(ji,jj,jk,jm)
              trn(ji,jj,jk,jpsfe) = trn(ji,jj,jk,jpsfe)                 &
     &          -grazoo(ji,jj,jk,jm,jppoc)*stofoo(ji,jj,jk,jppoc,2)
              trn(ji,jj,jk,jpbfe) = trn(ji,jj,jk,jpbfe)                 &
     &          -grazoo(ji,jj,jk,jm,jpgoc)*stofoo(ji,jj,jk,jpgoc,2)
              trn(ji,jj,jk,jpufe) = trn(ji,jj,jk,jpufe)                 &
     &          -grazoo(ji,jj,jk,jm,jphoc)*stofoo(ji,jj,jk,jphoc,2)
              trn(ji,jj,jk,jpsfe+nn_sizzoo(jm))=                        &
     &          trn(ji,jj,jk,jpsfe+nn_sizzoo(jm))                       &
     &          +rn_unazoo(jm)*grazof(ji,jj,jk,jm)
            END DO
            zgochoc(ngochoc) =                                          &
#  if ! defined key_trc_foram
     &        +torgel(ji,jj,jk)                                         &
#  endif        
     &        +torcru(ji,jj,jk)
!
! Evolution of nitrogen in large particles
! ----------------------------------------
!
            dgongoc=gongoc-ratn2c
            DO jm = nn_gonmin, nn_gonmax
              trn(ji,jj,jk,jpgon) = trn(ji,jj,jk,jpgon)                 &
     &          +float(nn_sizzoo(jm))*(grapoc(ji,jj,jk,jm)*ratn2c       &
     &          +grazoo(ji,jj,jk,jm,jpgoc)*rn_unazoo(jm)*dgongoc)
            END DO
!
! Evolution of small POC
! ----------------------
!     
            trn(ji,jj,jk,jppoc) = trn(ji,jj,jk,jppoc)-rempoc(ji,jj,jk)  &
     &        +resbac(ji,jj,jk)/3.0                                     &
     &        +(snkpoc(ji,jj,jk)-snkpoc(ji,jj,jk+1))/fse3t(ji,jj,jk)    &
     &        -xagg(ji,jj,jk)+xaggdoc(ji,jj,jk)+deppoc(ji,jj,jk)
!    
! Evolution of big POC
! --------------------
!
            zgrdmp = zgochoc(1)+xagg(ji,jj,jk)+xaggdoc2(ji,jj,jk)
            trn(ji,jj,jk,jpgoc) = trn(ji,jj,jk,jpgoc)-remgoc(ji,jj,jk)  &
     &        +zgrdmp                                                   &
     &        +(snkgoc(ji,jj,jk)-snkgoc(ji,jj,jk+1))/fse3t(ji,jj,jk)
            trn(ji,jj,jk,jpgon) = trn(ji,jj,jk,jpgon)-remgon(ji,jj,jk)  &
     &        +zgrdmp*ratn2c                                            &
     &        +(snkgon(ji,jj,jk)-snkgon(ji,jj,jk+1))/fse3t(ji,jj,jk)
!
! Evolution of huge POC
! --------------------
!
            trn(ji,jj,jk,jphoc) = trn(ji,jj,jk,jphoc)-remhoc(ji,jj,jk)  &
     &        +zgochoc(2)                                               &
     &        +(snkhoc(ji,jj,jk)-snkhoc(ji,jj,jk+1))/fse3t(ji,jj,jk)
!
! Evolution of dissolved Iron
! ---------------------------
!
            trn(ji,jj,jk,jpfer) = trn(ji,jj,jk,jpfer)                   &
     &        +rbafer(ji,jj,jk)-ubafer(ji,jj,jk)                        &
     &        -xscave(ji,jj,jk)+irondep(ji,jj,jk)+depfer(ji,jj,jk)
!
! Evolution of small biogenic Iron
! --------------------------------
!
            trn(ji,jj,jk,jpsfe) = trn(ji,jj,jk,jpsfe)-remsfe(ji,jj,jk)  &
     &        -xaggfe(ji,jj,jk)                                         &
     &        +(snksfe(ji,jj,jk)-snksfe(ji,jj,jk+1))/fse3t(ji,jj,jk)    &
     &        +deppoc(ji,jj,jk)*ferat3
!
! Evolution of big biogenic Iron
! ------------------------------
!
            trn(ji,jj,jk,jpbfe) = trn(ji,jj,jk,jpbfe)+ferat3*zgochoc(1) &
     &        +xaggfe(ji,jj,jk)-rembfe(ji,jj,jk)                        &
     &        +(snkbfe(ji,jj,jk)-snkbfe(ji,jj,jk+1))/fse3t(ji,jj,jk)
!
! Evolution of huge biogenic Iron
! ------------------------------
!
            trn(ji,jj,jk,jpufe) = trn(ji,jj,jk,jpufe)+ferat3*zgochoc(2) &
     &        -remufe(ji,jj,jk)                                         &
     &        +(snkufe(ji,jj,jk)-snkufe(ji,jj,jk+1))/fse3t(ji,jj,jk)
!
! Evolution of biogenic Silica
! ----------------------------
!
            trn(ji,jj,jk,jpbsi) = trn(ji,jj,jk,jpbsi)                   &
     &        +prorca3(ji,jj,jk)-losbsi(ji,jj,jk)
!
! Evolution of sinking biogenic silica
! ------------------------------------
!
            trn(ji,jj,jk,jpdsi) = trn(ji,jj,jk,jpdsi)                   &
     &        +losbsi(ji,jj,jk)-remdsi(ji,jj,jk)                        &
     &        +(snkdsi(ji,jj,jk)-snkdsi(ji,jj,jk+1))/fse3t(ji,jj,jk)
!
! Evolution of dissolved Silica
! -----------------------------
!
            trn(ji,jj,jk,jpsil) = trn(ji,jj,jk,jpsil)                   &
     &        -prorca3(ji,jj,jk)+remdsi(ji,jj,jk)                       &
     &        +depsil(ji,jj,jk)+sidep(ji,jj,jk)
#  if defined key_trc_dms
!
! Evolution of DMS,DMD
! --------------------
!
            trn(ji,jj,jk,jpdmd) = trn(ji,jj,jk,jpdmd)+                  &
     &        prodmd(ji,jj,jk) - degdmd(ji,jj,jk)
!
            trn(ji,jj,jk,jpdms) = trn(ji,jj,jk,jpdms)+                  &
     &    prodms(ji,jj,jk)+0.5*dmddms(ji,jj,jk)-degdms(ji,jj,jk)
#  endif
!
! Consumption of Total (12C)O2
! ----------------------------
!     
            trn(ji,jj,jk,jpdic) = trn(ji,jj,jk,jpdic)                      &
     &        -prodt(ji,jj,jk)-prcaca(ji,jj,jk)+rn_discal*loscal(ji,jj,jk) &
     &        +depdic(ji,jj,jk)
#  if defined key_trc_piic
            trn(ji,jj,jk,jppiic) = trn(ji,jj,jk,jppiic)                    &
     &        -prodt(ji,jj,jk)-prcaca(ji,jj,jk)+rn_discal*loscal(ji,jj,jk) &
     &        +depdic(ji,jj,jk)
#  endif
!     
! Consumption of alkalinity due to ca++ uptake and increase 
! of alkalinity due to nitrate consumption during organic 
! soft tissue production
! ---------------------------------------------------------
!     
            trn(ji,jj,jk,jptal) = trn(ji,jj,jk,jptal)                  &
     &        +alknut*prodt(ji,jj,jk)-2.*prcaca(ji,jj,jk)+             &
     &        2.*rn_discal*loscal(ji,jj,jk)+depdic(ji,jj,jk)+denitr(ji,jj,jk)
!            deln2o = nitrif*rn_aoun2o
!            trn(ji,jj,jk,jpn2o)= trn(ji,jj,jk,jpn2o)+deln2o
!            trc3d(ji,jj,jk,16) = deln2o*rfactr*1e3
          END DO
        END DO
      END DO
!!
!! --------------------------------------------------------------------
!! Apply the biological sources and sinks to the other greenhouse gases
!! --------------------------------------------------------------------
!!
! CH4 model 
! ---------
!
#  if defined key_trc_ch4
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
          ikt=max(mbathy(ji,jj)-1,1)
          IF (ikt .LE. 24) THEN
            DO jk = 1, 24
            STOP 'implement CH4 in TOM12'
            trn(ji,jj,jk,jpch1)  = trn(ji,jj,jk,jpch1)+mfecal*rn_proch2*0.5
            proch1(ji,jj,jk)     = mfecal*rn_proch1*1e3*rfactr
            trn(ji,jj,jk,jpch2)  = trn(ji,jj,jk,jpch2)+mfecal*rn_proch2
            proch2(ji,jj,jk)     = mfecal*rn_proch2*1e3*rfactr
            trn(ji,jj,jk,jpch3)  = trn(ji,jj,jk,jpch3)+mfecal*rn_proch2*1.5
            proch3(ji,jj,jk)     = mfecal*rn_proch3*1e3*rfactr
            trn(ji,jj,jk,jpch4)  = trn(ji,jj,jk,jpch4)+mfecal*rn_proch3
            proch4(ji,jj,jk)     = mfecal*rn_proch4*1e3*rfactr
            trn(ji,jj,jk,jpch5)  = trn(ji,jj,jk,jpch5)+mfecal*rn_proch2*0.5
            trn(ji,jj,jk,jpch6)  = trn(ji,jj,jk,jpch6)+mfecal*rn_proch2*1.5
            trn(ji,jj,jk,jpch7)  = trn(ji,jj,jk,jpch7)+mfecal*rn_proch2*0.5
            trn(ji,jj,jk,jpch8)  = trn(ji,jj,jk,jpch8)+mfecal*rn_proch2*1.5
            trn(ji,jj,jk,jpch9)  = trn(ji,jj,jk,jpch9)+mfecal*rn_proch2*0.5
            trn(ji,jj,jk,jpch10) = trn(ji,jj,jk,jpch10)+mfecal*rn_proch2
            trn(ji,jj,jk,jpch11) = trn(ji,jj,jk,jpch11)+mfecal*rn_proch2*1.5
            trn(ji,jj,jk,jpch12) = trn(ji,jj,jk,jpch12)+mfecal*rn_proch3
            trn(ji,jj,jk,jpch13) = trn(ji,jj,jk,jpch13)+mfecal*rn_proch4
            trn(ji,jj,jk,jpch14) = trn(ji,jj,jk,jpch14)+mfecal*rn_proch2*0.5
            trn(ji,jj,jk,jpch15) = trn(ji,jj,jk,jpch15)+mfecal*rn_proch2
            trn(ji,jj,jk,jpch16) = trn(ji,jj,jk,jpch16)+mfecal*rn_proch2*1.5
            trn(ji,jj,jk,jpch17) = trn(ji,jj,jk,jpch17)+mfecal*rn_proch2*0.75
            trn(ji,jj,jk,jpch18) = trn(ji,jj,jk,jpch18)+mfecal*rn_proch2*1.25
            trn(ji,jj,jk,jpch19) = trn(ji,jj,jk,jpch19)+mfecal*rn_proch2
            trn(ji,jj,jk,jpch20) = trn(ji,jj,jk,jpch20)+mfecal*rn_proch2
            trn(ji,jj,jk,jpch21) = trn(ji,jj,jk,jpch21)+mfecal*rn_proch2*1.25
            trn(ji,jj,jk,jpch22) = trn(ji,jj,jk,jpch22)+mfecal*rn_proch2*0.5
            trn(ji,jj,jk,jpch23) = trn(ji,jj,jk,jpch23)+mfecal*rn_proch2*1.5
            trn(ji,jj,jk,jpch24) = trn(ji,jj,jk,jpch24)+mfecal*rn_proch2*0.5
            trn(ji,jj,jk,jpch25) = trn(ji,jj,jk,jpch25)+mfecal*rn_proch2*1.5
            END DO
          ELSE
            DO jk = 1, 24
            mfecal = (grapoc(ji,jj,jk,2)+grapoc(ji,jj,jk,3))
            temch4 = min(trn(ji,jj,jk,jpch1),mfecal*rn_botch2*-5.)
            trn(ji,jj,jk,jpch5) = trn(ji,jj,jk,jpch1)-temch4
            temch4 = min(trn(ji,jj,jk,jpch2),mfecal*rn_botch2*-5.)
            trn(ji,jj,jk,jpch6) = trn(ji,jj,jk,jpch2)-temch4
            temch4 = min(trn(ji,jj,jk,jpch3),mfecal*rn_botch2*-5.)
            trn(ji,jj,jk,jpch7) = trn(ji,jj,jk,jpch3)-temch4
            temch4 = min(trn(ji,jj,jk,jpch4),mfecal*rn_botch2*-5.)
            trn(ji,jj,jk,jpch8) = trn(ji,jj,jk,jpch4)-temch4
            trn(ji,jj,jk,jpch1) = trn(ji,jj,jk,jpch1)+mfecal*rn_botch1
            trn(ji,jj,jk,jpch2) = trn(ji,jj,jk,jpch2)+mfecal*rn_botch1
            trn(ji,jj,jk,jpch3) = trn(ji,jj,jk,jpch3)+mfecal*rn_botch1
            trn(ji,jj,jk,jpch4) = trn(ji,jj,jk,jpch4)+mfecal*rn_botch1
            temch4 = min(trn(ji,jj,jk,jpch5),mfecal*rn_botch2*-1.)
            trn(ji,jj,jk,jpch5) = trn(ji,jj,jk,jpch5)-temch4
            proch5(ji,jj,jk)    = -temch4*1e3*rfactr
            temch4 = min(trn(ji,jj,jk,jpch6),mfecal*rn_botch2*-1.)
            trn(ji,jj,jk,jpch6) = trn(ji,jj,jk,jpch6)-temch4
            proch6(ji,jj,jk)    = -temch4*1e3*rfactr
            temch4 = min(trn(ji,jj,jk,jpch7),mfecal*rn_botch3*-1.)
            trn(ji,jj,jk,jpch7) = trn(ji,jj,jk,jpch7)-temch4
            temch4 = min(trn(ji,jj,jk,jpch8),mfecal*rn_botch3*-1.)
            trn(ji,jj,jk,jpch8) = trn(ji,jj,jk,jpch8)-temch4
            proch8(ji,jj,jk)    = -temch4*1e3*rfactr
            temch4 = min(trn(ji,jj,jk,jpch9),mfecal*rn_botch3*-10.)
            trn(ji,jj,jk,jpch9) = trn(ji,jj,jk,jpch9)-temch4
            proch9(ji,jj,jk)    = -temch4*1e3*rfactr
            temch4 = min(trn(ji,jj,jk,jpch10),mfecal*rn_botch3*-10.)
            trn(ji,jj,jk,jpch10) = trn(ji,jj,jk,jpch10)-temch4
            proch10(ji,jj,jk)    = -temch4*1e3*rfactr
            temch4 = min(trn(ji,jj,jk,jpch11),mfecal*rn_botch3*-10.)
            trn(ji,jj,jk,jpch11) = trn(ji,jj,jk,jpch11)-temch4
            proch11(ji,jj,jk)    = -temch4*1e3*rfactr
            temch4 = min(trn(ji,jj,jk,jpch12),mfecal*rn_botch3*-10.)
            trn(ji,jj,jk,jpch12) = trn(ji,jj,jk,jpch12)-temch4
            proch12(ji,jj,jk)    = -temch4*1e3*rfactr
            temch4 = min(trn(ji,jj,jk,jpch13),mfecal*rn_botch3*-10.)
            trn(ji,jj,jk,jpch13) = trn(ji,jj,jk,jpch13)-temch4
            proch13(ji,jj,jk)    = -temch4*1e3*rfactr
            temch4 = min(trn(ji,jj,jk,jpch14),mfecal*rn_botch2*-2.)
            trn(ji,jj,jk,jpch14) = trn(ji,jj,jk,jpch14)-temch4
            proch14(ji,jj,jk)    = -temch4*1e3*rfactr
            temch4 = min(trn(ji,jj,jk,jpch15),mfecal*rn_botch2*-2.)
            trn(ji,jj,jk,jpch15) = trn(ji,jj,jk,jpch15)-temch4
            proch7(ji,jj,jk)     = -temch4*1e3*rfactr
            temch4 = min(trn(ji,jj,jk,jpch16),mfecal*rn_botch2*-2.)
            trn(ji,jj,jk,jpch16) = trn(ji,jj,jk,jpch16)-temch4
            temch4 = min(trn(ji,jj,jk,jpch17),mfecal*rn_botch2*-1.)
            trn(ji,jj,jk,jpch17) = trn(ji,jj,jk,jpch17)-temch4
            temch4 = min(trn(ji,jj,jk,jpch18),mfecal*rn_botch2*-1.)
            trn(ji,jj,jk,jpch18) = trn(ji,jj,jk,jpch18)-temch4
            temch4 = min(trn(ji,jj,jk,jpch19),mfecal*rn_botch2*-1.5)
            trn(ji,jj,jk,jpch19) = trn(ji,jj,jk,jpch19)-temch4
            temch4 = min(trn(ji,jj,jk,jpch20),mfecal*rn_botch2*-0.5)
            trn(ji,jj,jk,jpch20) = trn(ji,jj,jk,jpch20)-temch4
            temch4 = min(trn(ji,jj,jk,jpch21),mfecal*rn_botch2*-2.)
            trn(ji,jj,jk,jpch21) = trn(ji,jj,jk,jpch21)-temch4
            trn(ji,jj,jk,jpch22) = trn(ji,jj,jk,jpch22)+mfecal*rn_botch1*5./6.
            trn(ji,jj,jk,jpch23) = trn(ji,jj,jk,jpch23)+mfecal*rn_botch1*5./6.
!           trn(ji,jj,jk,jpch24) = trn(ji,jj,jk,jpch24)+mfecal*rn_botch5*0.8379
!           trn(ji,jj,jk,jpch25) = trn(ji,jj,jk,jpch25)+mfecal*rn_botch5*0.8379
            END DO
          ENDIF
        END DO
      END DO
#  endif
!
! N2O model
! ---------
! 
#  if defined key_trc_n2o
      DO jk = 1, nn_deun2s-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
            out3d(ji,jj,jk) = 0.
          END DO
        END DO
      END DO
      DO jk = nn_deun2s, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
            dumn2s = delo2(ji,jj,jk)
            ppt(ji,jj,jk) = dumn2s*rn_aoun2s
            lo2n2s = dumn2s*rn_betn2s                                   &
     &        *exp(rn_decn2s*(trn(ji,jj,jk,jpoxy)-rn_omxn2s)/rn_omxn2s)
            trn(ji,jj,jk,jpn2s) = trn(ji,jj,jk,jpn2s)+ppt(ji,jj,jk)+lo2n2s
            out3d(ji,jj,jk)     = lo2n2s*1e3*rfactr
          END DO
        END DO
      END DO
#    if defined key_trc_diaadd && key_iomput
        CALL iom_put( "N2SprodB", out3d )
          out3d = ppt*1e3*rfactr
        CALL iom_put( "N2SprodA", out3d )
#    endif
      ppt=0.
#  endif
#  if defined key_trc_diaadd
!!
!! ---------------------------
!! Save biological diagnostics
!! ---------------------------
!!
#    if defined key_iomput
        out3d(:,:,:) = (snkhoc(:,:,:)+snkgoc(:,:,:)+snkpoc(:,:,:))*1e3*rfactr
        CALL iom_put("EXP", out3d )

          out3d = grazoc(:,:,:,1)*1e3*rfactr
        CALL iom_put("GRAPRO", out3d(:,:,:) )
!        
        DO jl = jpdia, jpdia+jppft-1   
            ppt(2:nlci-1,2:nlcj-1,1:jpk-1) = &
     &           ppt(2:nlci-1,2:nlcj-1,1:jpk-1)+ &
     &           prophy(2:nlci-1,2:nlcj-1,1:jpk-1,jl,1)*rfactr*1e3
            tchl(2:nlci-1,2:nlcj-1,1:jpk-1) = &
     &           tchl(2:nlci-1,2:nlcj-1,1:jpk-1)+ &
     &           trn(2:nlci-1,2:nlcj-1,1:jpk-1,jl+2*jppft)
            grazing(2:nlci-1,2:nlcj-1,1:jpk-1,1) = & 
     &           grazing(2:nlci-1,2:nlcj-1,1:jpk-1,1)+ &
     &           grazoo(2:nlci-1,2:nlcj-1,1:jpk-1,1,jl)*rfactr*1e3
#        if defined key_trc_foram
            grazing(2:nlci-1,2:nlcj-1,1:jpk-1,2) = &
     &           grazing(2:nlci-1,2:nlcj-1,1:jpk-1,2) +  &
     &           grazoo(2:nlci-1,2:nlcj-1,1:jpk-1,4,jl)*rfactr*1e3

            grazing(2:nlci-1,2:nlcj-1,1:jpk-1,3) = &
     &           grazing(2:nlci-1,2:nlcj-1,1:jpk-1,3) + &
     &           grazoo(2:nlci-1,2:nlcj-1,1:jpk-1,5,jl)*rfactr*1e3
 
#        else
            grazing(2:nlci-1,2:nlcj-1,1:jpk-1,2) = &
     &           grazing(2:nlci-1,2:nlcj-1,1:jpk-1,2) +  &
     &           grazoo(2:nlci-1,2:nlcj-1,1:jpk-1,2,jl)*rfactr*1e3
            grazing(2:nlci-1,2:nlcj-1,1:jpk-1,3) = &
     &           grazing(2:nlci-1,2:nlcj-1,1:jpk-1,3) + &
     &           grazoo(2:nlci-1,2:nlcj-1,1:jpk-1,3,jl)*rfactr*1e3
            grazing(2:nlci-1,2:nlcj-1,1:jpk-1,4) = &
     &           grazing(2:nlci-1,2:nlcj-1,1:jpk-1,4) + &
     &           grazoo(2:nlci-1,2:nlcj-1,1:jpk-1,4,jl)*rfactr*1e3
            grazing(2:nlci-1,2:nlcj-1,1:jpk-1,5) = &
     &           grazing(2:nlci-1,2:nlcj-1,1:jpk-1,5) + &
     &           grazoo(2:nlci-1,2:nlcj-1,1:jpk-1,5,jl)*rfactr*1e3
#         endif
        END DO
!
        CALL iom_put("PPT", ppt(:,:,:) )
        CALL iom_put("PPTDOC",pptdoc)
        CALL iom_put("TChl",tchl)
!
        IF( iom_use('PPT_DIA') ) CALL iom_put("PPT_DIA", prophy(:,:,:,jpdia,1)*rfactr*1e3   )
        IF( iom_use('PPT_MIX') ) CALL iom_put("PPT_MIX", prophy(:,:,:,jpdia+1,1)*rfactr*1e3 )
        IF( iom_use('PPT_COC') ) CALL iom_put("PPT_COC", prophy(:,:,:,jpdia+2,1)*rfactr*1e3 )
        IF( iom_use('PPT_PIC') ) CALL iom_put("PPT_PIC", prophy(:,:,:,jpdia+3,1)*rfactr*1e3 )
        IF( iom_use('PPT_PHA') ) CALL iom_put("PPT_PHA", prophy(:,:,:,jpdia+4,1)*rfactr*1e3 )
        IF( iom_use('PPT_FIX') ) CALL iom_put("PPT_FIX", prophy(:,:,:,jpdia+5,1)*rfactr*1e3 )
!
        IF( iom_use('Detrit') ) CALL iom_put("Detrit",trophic(:,:,:,1)*rfactr*1.e3 )
        IF( iom_use('Carniv') ) CALL iom_put("Carniv",trophic(:,:,:,2)*rfactr*1.e3 )
        IF( iom_use('Herbiv') ) CALL iom_put("Herbiv",trophic(:,:,:,3)*rfactr*1.e3 )
!
        IF( iom_use('GRAPROPHY') ) CALL iom_put("GRAPROPHY", grazing(:,:,:,1) )
        IF( iom_use('GRAPTEPHY') ) CALL iom_put("GRAPTEPHY", grazing(:,:,:,2) )
        IF( iom_use('GRAMESPHY') ) CALL iom_put("GRAMESPHY", grazing(:,:,:,3) )
        IF( iom_use('GRAGELPHY') ) CALL iom_put("GRAGELPHY", grazing(:,:,:,4) )
        IF( iom_use('GRACRUPHY') ) CALL iom_put("GRACRUPHY", grazing(:,:,:,5) )
!
          out3d = (prophy(:,:,:,jpfix,1)             &
     &        +docphy(:,:,:,jpfix))*ratn2c*(1.-dinpft(:,:,:,jpfix))*rfactr*1.e3
        CALL iom_put("nitrfix", out3d) 
!
          out3d = denitr(:,:,:)*rfactr*1e3 
        CALL iom_put("denitr", out3d )
!
          out3d = delo2(:,:,:)*1e3*rfactr
        CALL iom_put("DELO2", out3d)
          out3d = snkdsi(:,:,:)*1e3*rfactr
        CALL iom_put("sinksil", out3d)
          out3d = (snkcal(:,:,:)+snkara(:,:,:))*1e3*rfactr
        CALL iom_put("ExpCO3", out3d )
          out3d = snkara(:,:,:)*1e3*rfactr
        CALL iom_put("ExpARA", out3d )
          out3d = grazoc(:,:,:,5)*1e3*rfactr
        CALL iom_put("GRACRU", out3d(:,:,:) )
          out3d = proara(:,:,:)*rfactr*1e3
        CALL iom_put("proara", out3d )
#        if defined key_trc_foram
          out3d = grazoc(:,:,:,2)*1e3*rfactr
        CALL iom_put("GRAFOR", out3d )
          out3d = grazoc(:,:,:,3)*1e3*rfactr
        CALL iom_put("GRAPTE", out3d )
          out3d = grazoc(:,:,:,4)*1e3*rfactr
        CALL iom_put("GRAMES", out3d(:,:,:) )
          out3d = rn_forcal*mgezoo(:,:,:,2)*grazoc(:,:,:,2)*1e3*rfactr
        CALL iom_put("profor", out3d )
#        else
          out3d = grazoc(:,:,:,2)*1e3*rfactr
        CALL iom_put("GRAPTE", out3d )
          out3d = grazoc(:,:,:,3)*1e3*rfactr
        CALL iom_put("GRAMES", out3d(:,:,:) )
          out3d = grazoc(:,:,:,4)*1e3*rfactr
        CALL iom_put( "GRAGEL", out3d(:,:,:) )
#        endif
          out3d = rn_coccal*prophy(:,:,:,jpcoc,1)*rfactr*1e3
        CALL iom_put("prococ", out3d )
        CALL iom_put("probsi",prorca3*rfactr*1e3)
        CALL iom_put("losbsi",losbsi*rfactr*1e3)
        CALL iom_put("remdsi",remdsi*rfactr*1e3)
#      if defined key_trc_dms
        CALL iom_put("dmspp", dmspp(:,:,:) )
        CALL iom_put("degdms", degdms(:,:,:) )
        CALL iom_put("prodms", prodms(:,:,:) )
#      endif
#      if defined key_trc_ch4
        CALL iom_put("PROCH1",proch1)
        CALL iom_put("PROCH2",proch2)
        CALL iom_put("PROCH3",proch3)
        CALL iom_put("PROCH4",proch4)
        CALL iom_put("PROCH5",proch5)
        CALL iom_put("PROCH6",proch6)
        CALL iom_put("PROCH7",proch7)
        CALL iom_put("PROCH8",proch8)
        CALL iom_put("PROCH9",proch9)
        CALL iom_put("PROCH10",proch10)
        CALL iom_put("PROCH11",proch11)
        CALL iom_put("PROCH12",proch12)
        CALL iom_put("PROCH13",proch13)
        CALL iom_put("PROCH14",proch14)
#      endif
#      if defined key_c14b
        CALL iom_put("D14PRO",d14pro)
        CALL iom_put("D14RES",d14res)
#      endif
        CALL iom_put("PPINT",ppint)
#    else
      STOP 'key_trc_diaadd without key_iomput is no longer supported'
#    endif
#  endif
#endif
      RETURN
      END
