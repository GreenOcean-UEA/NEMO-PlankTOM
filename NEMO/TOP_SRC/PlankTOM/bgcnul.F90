      SUBROUTINE bgcnul
#if defined key_planktom && defined key_top
!!======================================================================
!!                         ***  ROUTINE bgcnul  ***
!! TOP : PlankTOM ecosystem model 
!!======================================================================
!! History : original pisces       ! 2001 (O. Aumont) 
!!           original p4znull : 2002    E. Buitenhuis
!!           PlankTOM see bgcbio
!!----------------------------------------------------------------------
!!    'key_trc_dms'       :            includes the DMS cycle and fluxes
!!    'key_trc_foram'     :      includes zooplankton calcifiers (foram) 
!!----------------------------------------------------------------------
!! ** Purpose : Sets some processes to zero to protect passive tracers 
!!              from going negative.
!!              variables using xcond are based on 0 or 1 multiplying factor
!!              variables using zprotec are based on a gradual multiplying factor
!! 
!! ** Action  : - 
!! ** External : NONE
!!
!! References : PlankTOM12 manual, Buitenhuis et al. 2023
!!                                 https://zenodo.org/records/8388158 
!!----------------------------------------------------------------------
!
! parameters and commons
      USE trp_trc
      USE sms_planktom
      USE oce_trc
      IMPLICIT NONE
! 
! local variables
      INTEGER ji, jj, jk, jl, jm, jn
      REAL xcond,phosph2,zoo2,mes2,dsil2,doc212
      REAL bac212,bsil2,din2,iron2,phyto2,poc2(0:3),piron2(0:2)
      REAL consum,consump,nitrfac,zlosss,zprotec,zminc,zloszoo(jpzft)
      REAL zlosfe(0:2),zprzoo(jpzft),gongoc,dgongoc
#    if defined key_trc_dms
      REAL dmd212, dms212
#    endif
#  include "domzgr_substitute.h90"
!
!  minimum concentration for PFT biomass
      consum = 0.
      zminc  = 1.e-10
!
!     Vertical loop to pre-compute concentration changes of the (rapidly
!     varying) tracers for preventing them to fall below 0
!     ------------------------------------------------------------------
!
      DO jk = 1, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
!
! sum the remineralisation rates for P (consump) and N:C:O (consum)
            consum = (1.-bactge(ji,jj,jk))*(remdoc(ji,jj,jk)               &
     &        +rempoc(ji,jj,jk)+remhoc(ji,jj,jk))+resbac(ji,jj,jk)/3.0
            DO jl = jpdia, jpdia+jppft-1
              consum = consum +resphy(ji,jj,jk,jl,1)
            END DO
            DO jm = 1, jpzft
              consum = consum+grarem(ji,jj,jk,jm)*rn_sigzoo(jm)         &
     &          +reszoo(ji,jj,jk,jm)
            END DO
            consump = consum+(1.-bactge(ji,jj,jk))*remgon(ji,jj,jk)*ratc2n
            consum  = consum+(1.-bactge(ji,jj,jk))*remgoc(ji,jj,jk)
!
! calculate the fraction of remineralisation that is supported by nitrate
            nitrfac = max(-0.5,(rn_denno3-trn(ji,jj,jk,jpoxy))/           &
     &        (2.*rn_denno3+trn(ji,jj,jk,jpoxy)))
            nitrfac = (sin(nitrfac*rpi)+1.)/2.
!
! C122N16H120 + 172 O2 -> 122 CO2 + 52 H2O + 16 HNO3
! 0.8 HNO3 -> 0.4 N2 + O2 + 0.4 H2O  (0.8 = 1/1.25)
            delo2(ji,jj,jk)  = rato2c*consum*(1.-nitrfac)
            denitr(ji,jj,jk) = 0.8*(rato2c*consum-delo2(ji,jj,jk))
!     
!     Evolution of PO4
!     ----------------
!     
            zlosss = max(bactge(ji,jj,jk)*(remgoc(ji,jj,jk)-remgon(ji,jj,jk)*ratc2n),0.)
            DO jl = jpdia, jpdia+jppft-1
              zlosss = zlosss +prophy(ji,jj,jk,jl,1)+docphy(ji,jj,jk,jl)
            END DO 
            phosph2 = trn(ji,jj,jk,jppo4)+consump &
     &       +deppo4(ji,jj,jk)
!     
!     Nullity test for PO4
!     --------------------
!     
            zprotec = max(0.,min(1.,(phosph2-1e-10)/(zlosss+rtrn)))
            DO jl = jpdia, jpdia+jppft-1
              prophy(ji,jj,jk,jl,1) = prophy(ji,jj,jk,jl,1)*zprotec
              docphy(ji,jj,jk,jl)   = docphy(ji,jj,jk,jl)  *zprotec
            END DO
            remgoc(ji,jj,jk) = remgoc(ji,jj,jk)*zprotec
!
!     Evolution of DIN
!     ----------------
!
            din2 = trn(ji,jj,jk,jpdin) &
     &        +consump*ratn2c+depnit(ji,jj,jk)+atmdin(ji,jj,jk)-1e-10
!
            zlosss = max(bactge(ji,jj,jk)*(remgoc(ji,jj,jk)*ratn2c-remgon(ji,jj,jk)),0.)
            DO jl = jpdia, jpdia+jppft-1
              zlosss = zlosss+(prophy(ji,jj,jk,jl,1)+docphy(ji,jj,jk,jl))&
     &          *ratn2c*dinpft(ji,jj,jk,jl)
            END DO
!
!     Nullity test for DIN
!     --------------------
!
            zprotec = max(0.,min(1.,din2/(zlosss+rtrn)))
            DO jl = jpdia, jpdia+jppft-1
              prophy(ji,jj,jk,jl,1) = prophy(ji,jj,jk,jl,1)*zprotec
              docphy(ji,jj,jk,jl)   = docphy(ji,jj,jk,jl)  *zprotec
            END DO
            remgoc(ji,jj,jk) = remgoc(ji,jj,jk)*zprotec
!
!     if denitr > din2+consum decrease consum of (O2+macronutrient)
            xcond = (0.5+sign(0.5,din2-zlosss*zprotec-denitr(ji,jj,jk)))
            din2  = trn(ji,jj,jk,jpdin)+depnit(ji,jj,jk)+atmdin(ji,jj,jk) &
     &        -zlosss*zprotec-1e-10
            zprotec = max(xcond,min(1.,din2*(0.8*rato2c)/ &
     &        ((0.8*rato2c*ratc2n-1.)*denitr(ji,jj,jk)+rtrn)))
            remdoc(ji,jj,jk) = remdoc(ji,jj,jk)*zprotec
            rempoc(ji,jj,jk) = rempoc(ji,jj,jk)*zprotec
            remgoc(ji,jj,jk) = remgoc(ji,jj,jk)*zprotec
            remhoc(ji,jj,jk) = remhoc(ji,jj,jk)*zprotec
            DO jn = 1, jpfoo
              jm = grizoo(1,jn)
              jl = grizoo(2,jn)
              grazoo(ji,jj,jk,jm,jl) = grazoo(ji,jj,jk,jm,jl)*zprotec
            END DO
            DO jl = jpdia, jpdia+jppft-1
              DO jm = 1, 3
                resphy(ji,jj,jk,jl,jm) = resphy(ji,jj,jk,jl,jm)*zprotec
              END DO
            END DO
            DO jm = 1, jpzft
              reszoo(ji,jj,jk,jm) = reszoo(ji,jj,jk,jm)*zprotec
            END DO
            resbac(ji,jj,jk) = resbac(ji,jj,jk)*zprotec
!
!     Evolution of dissolved iron
!     ---------------------------
!
            iron2 = trn(ji,jj,jk,jpfer)                               &
     &        +remsfe(ji,jj,jk)+rembfe(ji,jj,jk)                           &
     &       +rbafer(ji,jj,jk)+irondep(ji,jj,jk)+depfer(ji,jj,jk)
            DO jl = jpdia, jpdia+jppft-1
              iron2 = iron2 + resphy(ji,jj,jk,jl,2)
            END DO
            DO jm = 1, jpzft
              iron2 = iron2+grafer(ji,jj,jk,jm)+ferat3*reszoo(ji,jj,jk,jm)
            END DO
            zlosss = ubafer(ji,jj,jk)+xscave(ji,jj,jk)
            DO jl = jpdia, jpdia+jppft-1
              zlosss = zlosss +prophy(ji,jj,jk,jl,2)
            END DO
!
!     Nullity test for iron
!     ---------------------
!
            zprotec = max(0.,min(1.,(iron2-1e-13)/(zlosss+rtrn)))
            DO jl = jpdia, jpdia+jppft-1
              prophy(ji,jj,jk,jl,1) = prophy(ji,jj,jk,jl,1)*zprotec
              prophy(ji,jj,jk,jl,2) = prophy(ji,jj,jk,jl,2)*zprotec
              prophy(ji,jj,jk,jl,3) = prophy(ji,jj,jk,jl,3)*zprotec
            END DO
            ubafer(ji,jj,jk) = ubafer(ji,jj,jk)*zprotec
!    
!     Evolution of phytoplankton
!     --------------------------
!     
              DO jl = jpdia, jpdia+jppft-1
                DO jm = 1, 3
                  phyto2 = trn(ji,jj,jk,jl+(jm-1)*jppft)                &
     &              +prophy(ji,jj,jk,jl,jm)-resphy(ji,jj,jk,jl,jm)
                  DO jn = 1, jpzft
                    phyto2 = phyto2                                &
     &                -grazoo(ji,jj,jk,jn,jl)*stofoo(ji,jj,jk,jl,jm)
                  END DO
!     
!     Nullity test for phytoplankton
!     ------------------------------
!     
                  xcond=(0.5+sign(0.5,phyto2))
                  DO jn = 1, jpzft
                    grazoo(ji,jj,jk,jn,jl) = grazoo(ji,jj,jk,jn,jl)*xcond
                  END DO
                  resphy(ji,jj,jk,jl,jm) = resphy(ji,jj,jk,jl,jm)*xcond
                END DO
              END DO
!
!     Nullity test for SiO3
!     ---------------------
!
              dsil2 = trn(ji,jj,jk,jpsil)+remdsi(ji,jj,jk)                &
     &        +depsil(ji,jj,jk)+sidep(ji,jj,jk)
              zlosss = prorca3(ji,jj,jk)
              zprotec = max(0.,min(1.,(dsil2-1e-10)/(zlosss+rtrn)))
              prorca3(ji,jj,jk) = prorca3(ji,jj,jk)*zprotec
!
!     Evolution of biogenic silica in diatoms
!     ---------------------------------------
!
              bsil2 = trn(ji,jj,jk,jpbsi)                               &
     &        +prorca3(ji,jj,jk)-losbsi(ji,jj,jk)
!
!     Nullity test for biogenic silica in diatoms
!     -------------------------------------------
!
              xcond = (0.5+sign(0.5,bsil2))
              losbsi(ji,jj,jk) = losbsi(ji,jj,jk)*xcond
!    
!     Evolution of microzooplankton
!     -----------------------------
!    
            zloszoo = 0.
#  if ! defined key_trc_foram
            zloszoo(jpzft-1) = torgel(ji,jj,jk)
#  endif
            zloszoo(jpzft) = torcru(ji,jj,jk)
            DO jm = 1, jpzft
              zoo2 = trn(ji,jj,jk,jpbac+jm)                             &
     &          +mgezoo(ji,jj,jk,jm)*grazoc(ji,jj,jk,jm)
              zloszoo(jm) = zloszoo(jm)+reszoo(ji,jj,jk,jm)
#  if defined key_trc_foram
              DO jn = jm+1, jpzft
#  else
              DO jn = 1, jpzft
#  endif
                zloszoo(jm) = zloszoo(jm)+grazoo(ji,jj,jk,jn,jpbac+jm)
              END DO
!    
!     Nullity test for zooplankton
!     ----------------------------
!    
              zprzoo(jm) = max(0.,min(1.,(zoo2-zminc)/(zloszoo(jm)+rtrn)))
              reszoo(ji,jj,jk,jm) = reszoo(ji,jj,jk,jm)*zprzoo(jm)
#  if defined key_trc_foram
              DO jn = jm+1, jpzft
#  else
              DO jn = 1, jpzft
#  endif
                grazoo(ji,jj,jk,jn,jpbac+jm) = grazoo(ji,jj,jk,jn,jpbac+jm)*zprzoo(jm)
              END DO
            END DO
#  if ! defined key_trc_foram
              torgel(ji,jj,jk) = torgel(ji,jj,jk)*zprzoo(jpzft-1)
#  endif
              torcru(ji,jj,jk) = torcru(ji,jj,jk)*zprzoo(jpzft)
!
!     Nullity test for DOC
!     --------------------
!
            doc212 = trn(ji,jj,jk,jpdoc)                                &
     &        +depdoc(ji,jj,jk)+resbac(ji,jj,jk)/3.0 
            DO jm = 1, jpzft
              doc212 = doc212+grarem(ji,jj,jk,jm)*(1.-rn_sigzoo(jm))
            END DO
          zlosss  = remdoc(ji,jj,jk)+xaggdoc(ji,jj,jk)+xaggdoc2(ji,jj,jk)
          zprotec = max(0.,min(1.,(doc212-1e-10)/(zlosss+rtrn)))
          remdoc(ji,jj,jk)   = remdoc(ji,jj,jk)*zprotec
          xaggdoc(ji,jj,jk)  = xaggdoc(ji,jj,jk)*zprotec
          xaggdoc2(ji,jj,jk) = xaggdoc2(ji,jj,jk)*zprotec
!
!     Nullity test for bacteria
!     -------------------------
          bac212 = trn(ji,jj,jk,jpbac)+bactge(ji,jj,jk)* &
     &      (remdoc(ji,jj,jk)+rempoc(ji,jj,jk)+remgoc(ji,jj,jk)         &
     &      +remhoc(ji,jj,jk))
          zlosss = resbac(ji,jj,jk)
            DO jm = 1, jpzft
              zlosss = zlosss+grazoo(ji,jj,jk,jm,jpbac)
            END DO
          zprotec = max(0.,min(1.,(bac212-1e-10)/(zlosss+rtrn)))
          resbac(ji,jj,jk) = resbac(ji,jj,jk)*zprotec
            DO jm = 1, jpzft
              grazoo(ji,jj,jk,jm,jpbac) = grazoo(ji,jj,jk,jm,jpbac)*zprotec
            END DO
!
!     Evolution of all POC
!     --------------------
!
              poc2(0) = trn(ji,jj,jk,jppoc)-rempoc(ji,jj,jk)              &
     &       -xagg(ji,jj,jk)+xaggdoc(ji,jj,jk)+deppoc(ji,jj,jk)         &
     &       +(snkpoc(ji,jj,jk)-snkpoc(ji,jj,jk+1))/fse3t(ji,jj,jk)   &
     &        +resbac(ji,jj,jk)/3.0
              xcond = xagg(ji,jj,jk)+xaggdoc2(ji,jj,jk)
              poc2(1) = trn(ji,jj,jk,jpgoc)-remgoc(ji,jj,jk)             &
     &       +xcond                                                     &
     &       +(snkgoc(ji,jj,jk)-snkgoc(ji,jj,jk+1))/fse3t(ji,jj,jk)
              poc2(2) = trn(ji,jj,jk,jphoc)+snkhoc(ji,jj,jk)/fse3t(ji,jj,jk)
              poc2(ngochoc) = poc2(ngochoc)                               &
#  if ! defined key_trc_foram
     &          +torgel(ji,jj,jk)                                       &
#  endif 
     &          +torcru(ji,jj,jk)
            DO jm = 1, jpzft
              poc2(nn_sizzoo(jm)) = poc2(nn_sizzoo(jm))+grapoc(ji,jj,jk,jm)
              poc2(0) = poc2(0)-grazoo(ji,jj,jk,jm,jppoc)
              poc2(1) = poc2(1)-grazoo(ji,jj,jk,jm,jpgoc)
            END DO
!
!     Nullity test for GOC
!     --------------------
!
              xcond = (0.5+sign(0.5,poc2(1)))
              snkgoc(ji,jj,jk+1) = snkgoc(ji,jj,jk+1)*xcond
              snkgon(ji,jj,jk+1) = snkgon(ji,jj,jk+1)*xcond
              remgoc(ji,jj,jk)   = remgoc(ji,jj,jk)*xcond
              remgon(ji,jj,jk)   = remgon(ji,jj,jk)*xcond
            DO jm = 1, jpzft
              grazoo(ji,jj,jk,jm,jpgoc) = grazoo(ji,jj,jk,jm,jpgoc)*xcond
            END DO
!
!     Evolution of HOC
!     ----------------
!
              zlosss = remhoc(ji,jj,jk)+snkhoc(ji,jj,jk+1)/fse3t(ji,jj,jk)
            DO jm = 1, jpzft
              zlosss = zlosss+grazoo(ji,jj,jk,jm,jphoc)
            END DO
!
!     Nullity test for HOC
!     --------------------
!
              zprotec = max(0.,min(1.,(poc2(2)-zminc)/(zlosss+rtrn)))
              snkhoc(ji,jj,jk+1) = snkhoc(ji,jj,jk+1)*zprotec
              remhoc(ji,jj,jk)   = remhoc(ji,jj,jk)*zprotec
            DO jm = 1, jpzft
              grazoo(ji,jj,jk,jm,jphoc) = grazoo(ji,jj,jk,jm,jphoc)*zprotec
            END DO
!
!     Evolution of GON
!     ----------------
!
              poc2(3) = trn(ji,jj,jk,jpgon)-remgon(ji,jj,jk)             &
     &       +xcond*ratn2c                                              &
     &       +(snkgon(ji,jj,jk)-snkgon(ji,jj,jk+1))/fse3t(ji,jj,jk)
            gongoc  = trn(ji,jj,jk,jpgon)/max(trn(ji,jj,jk,jpgoc),rtrn)
            dgongoc = gongoc-ratn2c
            DO jm = nn_gonmin, nn_gonmax
              poc2(3) = poc2(3)+float(nn_sizzoo(jm))*(grapoc(ji,jj,jk,jm)&
     &          *ratn2c+grazoo(ji,jj,jk,jm,jpgoc)*rn_unazoo(jm)*dgongoc)
            END DO
            DO jm = 1, jpzft
              poc2(3) = poc2(3)-grazoo(ji,jj,jk,jm,jpgoc)*gongoc
            END DO
              xcond=(0.5+sign(0.5,poc2(3)))
              snkgon(ji,jj,jk+1) = snkgon(ji,jj,jk+1)*xcond
              remgon(ji,jj,jk)   = remgon(ji,jj,jk)*xcond
            DO jm = 1, jpzft
              grazoo(ji,jj,jk,jm,jpgoc) = grazoo(ji,jj,jk,jm,jpgoc)*xcond
            END DO
!     
!     Nullity test for POC
!     --------------------
!     
              xcond = (0.5+sign(0.5,poc2(0)))
              snkpoc(ji,jj,jk+1) = snkpoc(ji,jj,jk+1)*xcond
              rempoc(ji,jj,jk)   = rempoc(ji,jj,jk)*xcond
              xagg(ji,jj,jk)     = xagg(ji,jj,jk)*xcond
            DO jm = 1, jpzft
              grazoo(ji,jj,jk,jm,jppoc) = grazoo(ji,jj,jk,jm,jppoc)*xcond
            END DO
!
!     Evolution of small biogenic iron
!     --------------------------------
!
              piron2(0) = trn(ji,jj,jk,jpsfe)+ferat3*deppoc(ji,jj,jk)     &
     &        -remsfe(ji,jj,jk)-xaggfe(ji,jj,jk)           &
     &        +(snksfe(ji,jj,jk)-snksfe(ji,jj,jk+1))/fse3t(ji,jj,jk)
              piron2(1) = trn(ji,jj,jk,jpbfe)-rembfe(ji,jj,jk)             &
     &        +xaggfe(ji,jj,jk)                                         &
     &        +(snkbfe(ji,jj,jk)-snkbfe(ji,jj,jk+1))/fse3t(ji,jj,jk)
              piron2(2) = trn(ji,jj,jk,jpufe)+snkufe(ji,jj,jk)/fse3t(ji,jj,jk)
              piron2(ngochoc) = piron2(ngochoc)                           &
#  if ! defined key_trc_foram
     &          +ferat3*torgel(ji,jj,jk)                                &
#  endif
     &          +ferat3*torcru(ji,jj,jk)
              zlosfe(2) = remufe(ji,jj,jk)+snkufe(ji,jj,jk+1)/fse3t(ji,jj,jk)
            DO jm = 1, jpzft
              piron2(0) = piron2(0)-grazoo(ji,jj,jk,jm,jppoc)*stofoo(ji,jj,jk,jppoc,2)
              zlosfe(2) = zlosfe(2)+grazoo(ji,jj,jk,jm,jphoc)*stofoo(ji,jj,jk,jphoc,2)
              piron2(nn_sizzoo) = piron2(nn_sizzoo)+rn_unazoo(jm)*grazof(ji,jj,jk,jm)
            END DO
!
!     Nullity test for small biogenic iron
!     ------------------------------------
!
              xcond = (0.5+sign(0.5,piron2(0)))
              snksfe(ji,jj,jk+1) = snksfe(ji,jj,jk+1)*xcond
              remsfe(ji,jj,jk)   = remsfe(ji,jj,jk)*xcond
              xaggfe(ji,jj,jk)   = xaggfe(ji,jj,jk)*xcond
!
!     Evolution of big biogenic iron
!     ------------------------------
!
              DO jl = jpdia, jpdia+jppft-1
                piron2(1) = piron2(1)+rn_unazoo(2)*grazoo(ji,jj,jk,2,jl)*stofoo(ji,jj,jk,jl,2)
              END DO
            DO jm = 1, jpzft 
              grazoo(ji,jj,jk,jm,jppoc) = grazoo(ji,jj,jk,jm,jppoc)*xcond
              piron2(1) = piron2(1)-grazoo(ji,jj,jk,jm,jpgoc)*stofoo(ji,jj,jk,jpgoc,2)
            END DO
!
!     Nullity test for big biogenic iron
!     ----------------------------------
!
              xcond = (0.5+sign(0.5,piron2(1)))               
              snkbfe(ji,jj,jk+1) = snkbfe(ji,jj,jk+1)*xcond             
              rembfe(ji,jj,jk)   = rembfe(ji,jj,jk)*xcond
            DO jm = 1, jpzft
              grazoo(ji,jj,jk,jm,jpgoc) = grazoo(ji,jj,jk,jm,jpgoc)*xcond
            END DO
!
!     Nullity test for huge biogenic iron
!     -----------------------------------
!
              zprotec = max(0.,min(1.,(piron2(2)-rtrn)/(zlosfe(2)+rtrn)))
              remufe(ji,jj,jk)   = remufe(ji,jj,jk)*zprotec
              snkufe(ji,jj,jk+1) = snkufe(ji,jj,jk+1)*zprotec
            DO jm = 1, jpzft
              grazoo(ji,jj,jk,jm,jphoc) = grazoo(ji,jj,jk,jm,jphoc)*zprotec
            END DO
!
!     Evolution of sinking biogenic silica
!     ------------------------------------
!
            dsil2 = trn(ji,jj,jk,jpdsi)+losbsi(ji,jj,jk)                  &
     &        +snkdsi(ji,jj,jk)/fse3t(ji,jj,jk)
            zlosss = remdsi(ji,jj,jk)+snkdsi(ji,jj,jk+1)/fse3t(ji,jj,jk)
!
!     Nullity test for sinking biogenic silica
!     ----------------------------------------
!
             zprotec = max(0.,min(1.,(dsil2-1e-10)/(zlosss+rtrn)))
             snkdsi(ji,jj,jk+1) = snkdsi(ji,jj,jk+1)*zprotec            
             remdsi(ji,jj,jk)   = remdsi(ji,jj,jk)*zprotec
!
!     Protection of CaCO3 concentration
!     ---------------------------------
!
               snkcal(ji,jj,jk+1) = min(trn(ji,jj,jk,jpcal), &
     &           snkcal(ji,jj,jk+1)/fse3t(ji,jj,jk+1))*fse3t(ji,jj,jk+1)
               snkara(ji,jj,jk+1) = min(trn(ji,jj,jk,jpara), &
     &           snkara(ji,jj,jk+1)/fse3t(ji,jj,jk+1))*fse3t(ji,jj,jk+1)
#    if defined key_trc_dms
!
!     Protection of DMS(Pd)
!     ---------------------
!
            dmd212  = trn(ji,jj,jk,jpdmd)+prodmd(ji,jj,jk)
            zlosss  = degdmd(ji,jj,jk)
            zprotec = max(0.,min(1.,(dmd212-1.e-12)/(zlosss+rtrn)))
            if ( zprotec .gt. 1.-1.e-12 .or. &
     &           zprotec .lt. 1.+1.e-12 ) dmd_snkcount=dmd_snkcount+1   
            degdmd(ji,jj,jk) = degdmd(ji,jj,jk)*zprotec
            dmddms(ji,jj,jk) = dmddms(ji,jj,jk)*zprotec
            dms212  = trn(ji,jj,jk,jpdms)+prodms(ji,jj,jk)+dmddms(ji,jj,jk)
            zlosss  = degdms(ji,jj,jk)
            zprotec = max(0.,min(1.,(dms212-1.e-12)/(zlosss+rtrn)))
            degdms(ji,jj,jk) = degdms(ji,jj,jk)*zprotec
            if ( zprotec .gt. 1.-1.e-12 .or. &
     &           zprotec .lt. 1.+1.e-12 ) dms_snkcount=dms_snkcount+1   
#    endif
          END DO
        END DO
      END DO
#endif
      RETURN
      END
