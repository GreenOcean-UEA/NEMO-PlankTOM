      SUBROUTINE bgcsnk
#if defined key_planktom  && key_top
!!======================================================================
!!                         ***  ROUTINE bgcsnk  ***
!! TOP : PlankTOM ecosystem model 
!!======================================================================
!! History : PISCES       ! 2001 (O. Aumont) Original code 
!!           PlankTOM see bgcbio
!!----------------------------------------------------------------------
!!    'key_trc_dms'       :            includes the DMS cycle and fluxes
!!    'key_trc_piic'      :    additional tracer with pre-industrial DIC 
!!    'key_iomput'        :                              use IOM library 
!!----------------------------------------------------------------------
!! ** Purpose : Compute vertical flux of particulate matter due to
!!              gravitational sinking
!! ** Action  : Pelagic model
!!              Sediment model 
!!
!! References : Buitenhuis et al. 2013 doi:10.1002/gbc.20074
!!              PlankTOM12 manual, Buitenhuis et al. 2023
!!                                 https://zenodo.org/records/8388158 
!! ---------------------------------------------------------------------
!
! parameters and commons
      USE trc
      USE trp_trc
      USE sms_planktom
      USE oce_trc
#if defined key_iomput
      USE iom
#endif
      IMPLICIT NONE
! 
! local variables
      INTEGER ji, jj, jk,jl
      REAL bacfer,maxlig,omeara,omecal
      REAL remik,sedflx,siremin,snkspd(jpdsi:jphoc),totprt
      REAL xagg1,xagg2,xagg3,xagg4,xdens,xfeequi,xkeq,xlam1b,xlamfc
      REAL zdenom,zrfe2
      REAL xlibad,xlibas,xphdms,zbldmd,xcldmd
#  include "domzgr_substitute.h90"
!!
!! -------------
!! Pelagic model
!! -------------
!!
!
! Computation of the vertical sinking speed
! -----------------------------------------
!
      DO jk = 1, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
!    
!     Sinking flux of pFe, POC, GOC, BSi, and Cal
!     ------------------------------------------------------------------
!
! these water column equations are also used in the sediment model, so
! WHEN YOU CHANGE THE EQUATIONS THAT FOLLOW, CHANGE THE SEDIMENT
! MODEL BELOW AS WELL
          snkpoc(ji,jj,jk+1) =                       &
     &           rn_snkpoc/rjjss*trn(ji,jj,jk,jppoc) &
     &           *rfact*tmask(ji,jj,jk+1)
          snksfe(ji,jj,jk+1) =                       &
     &           rn_snkpoc/rjjss*trn(ji,jj,jk,jpsfe) &
     &           *rfact*tmask(ji,jj,jk+1)
!
! rn_siegoc and rn_sisgoc are derived by fitting sinking speeds 
! calculated with the sedimentation function of Buitenhuis et al. 2001
! to an exponential function 
! snkmax is calculated in trcini_planktom.F90 to achieve a CFL condition
!
! calculate the density of large particles 
          xdens = (trn(ji,jj,jk,jpgoc)*moworg+trn(ji,jj,jk,jpcal)*mowco3 &
     &     +trn(ji,jj,jk,jpara)*mowco3                                  &
     &     +trn(ji,jj,jk,jpdsi)*mwsio2)*1e6/                            &
     &     MAX(trn(ji,jj,jk,jpgoc)*moworg/rhoorg                        &
     &     +trn(ji,jj,jk,jpcal)*mowco3/rhoco3                           &
     &     +trn(ji,jj,jk,jpara)*mowco3/rhoco3                           &
     &     +trn(ji,jj,jk,jpdsi)*mwsio2/rhsio2,rtrn)-rhop(ji,jj,jk)*1000.
!
! calculate the sinking speed of big particles 
          xvsink(ji,jj,jk) = MIN(rn_sisgoc*MAX(xdens,dnsmin)**rn_siegoc,snkmax(jk))
!
! calculate the sedimentation of big particles and mineral content 
          snkgoc(ji,jj,jk+1) =                    &
     &            xvsink(ji,jj,jk)/rjjss*trn(ji,jj,jk,jpgoc) &
     &            *rfact*tmask(ji,jj,jk+1)
          snkgon(ji,jj,jk+1) =                    &
     &            xvsink(ji,jj,jk)/rjjss*trn(ji,jj,jk,jpgon) &
     &            *rfact*tmask(ji,jj,jk+1)
!
          snkbfe(ji,jj,jk+1) =                      &
     &              xvsink(ji,jj,jk)/rjjss*trn(ji,jj,jk,jpbfe) &
     &              *rfact*tmask(ji,jj,jk+1)
 
!
          snkdsi(ji,jj,jk+1) = xvsink(ji,jj,jk)*trn(ji,jj,jk,jpdsi) &
     &           /rjjss*rfact*tmask(ji,jj,jk+1)
!
          snkcal(ji,jj,jk+1) = xvsink(ji,jj,jk)*trn(ji,jj,jk,jpcal) &
     &           /rjjss*rfact*tmask(ji,jj,jk+1)
          snkara(ji,jj,jk+1) = xvsink(ji,jj,jk)*trn(ji,jj,jk,jpara)  &
     &           /rjjss*rfact*tmask(ji,jj,jk+1)
          snkhoc(ji,jj,jk+1) = rn_snkhoc/rjjss*trn(ji,jj,jk,jphoc) &
     &            *rfact*tmask(ji,jj,jk+1)
          snkufe(ji,jj,jk+1) = rn_snkhoc/rjjss*trn(ji,jj,jk,jpufe) &
     &            *rfact*tmask(ji,jj,jk+1)
          END DO
        END DO
      END DO
!
!  Exchange between organic matter compartments due to
!  coagulation/disaggregation
!  ---------------------------------------------------
!
      DO jk = 1, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
!
!    Part I : Coagulation dependent on turbulence
!    --------------------------------------------
!
       xagg1 = rn_ag1poc/rjjss*rfact*min(avt(ji,jj,jk)/5.E-4,1.)* &
     &       trn(ji,jj,jk,jppoc)**2
       xagg2 = rn_ag2poc/rjjss*rfact*min(avt(ji,jj,jk)/5.E-4,1.)* &
     &       trn(ji,jj,jk,jppoc)*trn(ji,jj,jk,jpgoc)
!
!    Aggregation of small into large particles
!    Part II : Differential settling
!    -----------------------------------------
!
       xagg3 = rn_ag3poc/rjjss*rfact*                  &
     &       trn(ji,jj,jk,jppoc)*trn(ji,jj,jk,jpgoc)
       xagg4 = rn_ag4poc/rjjss*rfact* &
     &       trn(ji,jj,jk,jppoc)**2
       xagg(ji,jj,jk)   = xagg1+xagg2+xagg3+xagg4
       xaggfe(ji,jj,jk) = xagg(ji,jj,jk)*trn(ji,jj,jk,jpsfe)/        &
     &                  (trn(ji,jj,jk,jppoc)+rtrn)*tmask(ji,jj,jk)
!
!    Aggregation of DOC to small and large particles
!    -----------------------------------------------
!
        xaggdoc(ji,jj,jk) = (rn_ag5doc*trn(ji,jj,jk,jpdoc)              &
     & +rn_ag7doc*trn(ji,jj,jk,jppoc))/rjjss*rfact                      &
     & *min(avt(ji,jj,jk)/5.E-4,1.)*trn(ji,jj,jk,jpdoc)*tmask(ji,jj,jk)
        xaggdoc2(ji,jj,jk)=rn_ag6doc*trn(ji,jj,jk,jpgoc)*rfact          &
     & /rjjss*min(avt(ji,jj,jk)/5.E-4,1.)*trn(ji,jj,jk,jpdoc)           &
     & *tmask(ji,jj,jk)
          END DO
        END DO
      END DO
      DO jk = 1, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
!
!    TOC remineralization. Depends on temperature and nutrients
!    ----------------------------------------------------------
!
            zdenom = rn_kmobac+gbadoc*trn(ji,jj,jk,jpdoc) &
     &        +gbapoc*trn(ji,jj,jk,jppoc) &
     &        +gbagoc*trn(ji,jj,jk,jpgoc)                            &
     &        +gbahoc*trn(ji,jj,jk,jphoc)
            remik = rn_grabac/rjjss*rfact*tmask(ji,jj,jk) &
     &        *tgfunc(ji,jj,jk,jpbac) &
     &        *(trn(ji,jj,jk,jpoxy)+3.E-6)/(10.E-6+trn(ji,jj,jk,jpoxy)) &
     &        *trn(ji,jj,jk,jpbac)
!
! remineralization of DOC
            remdoc(ji,jj,jk) = remik*gbadoc*trn(ji,jj,jk,jpdoc)/zdenom
!
! remineralization of POC
! These water column equations are also used in the sediment model, so
! WHEN YOU CHANGE THE EQUATIONS THAT FOLLOW, CHANGE THE SEDIMENT
! MODEL BELOW AS WELL
            remik = rn_grabac/rjjss*rfact*tmask(ji,jj,jk) &
     &        *tgfunc(ji,jj,jk,jpbac) &
     &        *(trn(ji,jj,jk,jpoxy)+3.E-6)/(10.E-6+trn(ji,jj,jk,jpoxy)) &
     &        *max(trn(ji,jj,jk,jpbac),rn_rembac)
            rempoc(ji,jj,jk) = remik*gbapoc*trn(ji,jj,jk,jppoc)/zdenom
!
! remineralization of GOC and HOC
            remgoc(ji,jj,jk) = remik*gbagoc*trn(ji,jj,jk,jpgoc)/zdenom
            remgon(ji,jj,jk) = remik*gbagon*trn(ji,jj,jk,jpgon)/zdenom
            remhoc(ji,jj,jk) = remik*gbahoc*trn(ji,jj,jk,jphoc)/zdenom
!
! remineralization of iron in POC, GOC and HOC
            remsfe(ji,jj,jk) = remik*gbapoc*trn(ji,jj,jk,jpsfe)/zdenom
            rembfe(ji,jj,jk) = remik*gbagoc*trn(ji,jj,jk,jpbfe)/zdenom
            remufe(ji,jj,jk) = remik*gbahoc*trn(ji,jj,jk,jpufe)/zdenom
!
            bactge(ji,jj,jk) = rn_ggebac-rn_ggtbac*tsnbio(ji,jj,jk,1)
!
! Fe that is taken up by bacteria minus what is available in
! DOFe, sPOFe and bPOFe
            bacfer = bactge(ji,jj,jk)*ferat3*(remdoc(ji,jj,jk)             &
     &        +rempoc(ji,jj,jk)+remgoc(ji,jj,jk)+remhoc(ji,jj,jk))      &
     &        -remsfe(ji,jj,jk)-rembfe(ji,jj,jk)-remufe(ji,jj,jk)
!
! If there is not enough Fe in DOFe, sPOFe and bPOFe, then take up dissolved Fe
            ubafer(ji,jj,jk) = max(bacfer*trn(ji,jj,jk,jpfer)/            &
     &        (trn(ji,jj,jk,jpfer)+rn_kmfbac),0.)
!
! If there is not enough dissolved Fe, then decrease ggebac
            bactge(ji,jj,jk) = min(bactge(ji,jj,jk),                      &
     &        (ubafer(ji,jj,jk)+remsfe(ji,jj,jk)+rembfe(ji,jj,jk)       &
     &        +remufe(ji,jj,jk))/max((remdoc(ji,jj,jk)+rempoc(ji,jj,jk) &
     &        +remgoc(ji,jj,jk)+remhoc(ji,jj,jk))                       &
     &        *ferat3,minfer))
!
! Fe in excess of that taken up by bacteria
            rbafer(ji,jj,jk) = max(-bacfer,0.)
            resbac(ji,jj,jk) = rn_resbac*(rn_retbac**tsnbio(ji,jj,jk,1)) &
     &        /rjjss*rfact*max(trn(ji,jj,jk,jpbac)-1e-10,0.)*tmask(ji,jj,jk)
!
!     Remineralisation rate of BSi dependent on T and O2
!     --------------------------------------------------
!
            siremin = min(rn_remdsi*exp(rn_retdsi/(273.15+tsnbio(ji,jj,jk,1))) &
     &        ,rn_readsi)/rjjss*rfact*tmask(ji,jj,jk) &
     &        *(trn(ji,jj,jk,jpoxy)+3.E-6)/(10.E-6+trn(ji,jj,jk,jpoxy))
            remdsi(ji,jj,jk) = siremin*trn(ji,jj,jk,jpdsi)
!
!     Scavenging rate of iron based on Parekh et al., GBC 2005 as implemented by 
!     Aumont & Bopp GBC 2006
!     -------------------------------------------------------------------------
!
         xkeq    = 10**(17.27 - 1565.7 / ( 273.15 + tsnbio(ji,jj,jk,1) +         &
     &              (1.-tmask(ji,jj,jk))*20. ) )

         xfeequi = (-(1.+ligfer(jj,jk)*xkeq-xkeq*trn(ji,jj,jk,jpfer))+            &
     &             ((1.+ligfer(jj,jk)*xkeq-xkeq*trn(ji,jj,jk,jpfer))**2           &
     &             +4.*trn(ji,jj,jk,jpfer)*xkeq)**0.5)/(2.*xkeq)
         totprt  = trn(ji,jj,jk,jppoc)+trn(ji,jj,jk,jpgoc)              &
     &     +trn(ji,jj,jk,jphoc)+trn(ji,jj,jk,jpcal)+trn(ji,jj,jk,jpdsi)
         xlam1b  = rn_scmfer+rn_scofer*totprt*1E6

         xscave(ji,jj,jk) = xfeequi*xlam1b/rjjss*rfact*tmask(ji,jj,jk)
#    if defined key_trc_dms
!
!     Bacterial production and degradation of DMS(Pd)
!     -----------------------------------------------
!
         zbldmd = min(1., &
     &            max(0.66,1.-(etot(ji,jj,jk)/rn_etomax)**6.+0.18))
         xlibad = min(xlimbac(ji,jj,jk),trn(ji,jj,jk,jpdmd)/      &
     &            (trn(ji,jj,jk,jpdmd)+rn_xkdmd))                 &
     &            * zbldmd
!
!     Cleavage of dmspd
!     -----------------
!
         xcldmd = rn_xcldmd*rfact/rjjss*trn(ji,jj,jk,jpdmd)

         degdmd(ji,jj,jk) = 2*rn_grabac*xlibad  &
     &    *1.12**(tsnbio(ji,jj,jk,1))*2.*trn(ji,jj,jk,jpbac) &
     &    * rfact/rjjss *trn(ji,jj,jk,jpdmd)

         dmddms(ji,jj,jk) = rn_dmsyld*degdmd(ji,jj,jk)+xcldmd
         degdmd(ji,jj,jk) = degdmd(ji,jj,jk)+xcldmd
         xlibas = min(xlimbac(ji,jj,jk), trn(ji,jj,jk,jpdms)/ &
     &            (trn(ji,jj,jk,jpdms)+rn_xkdms)) &
     &            * zbldmd
!
!     Photolysis of dms
!     -----------------
!
         degdms(ji,jj,jk) = 2*(rn_grabac*xlibas  &
     &        *1.12**(tsnbio(ji,jj,jk,1))*0.6*trn(ji,jj,jk,jpbac) &
     &        +rn_xpodms*etot(ji,jj,jk))*rfact/rjjss*trn(ji,jj,jk,jpdms)
#    endif
          END DO
        END DO
      END DO
!
!! -----------------------------------------------------------------
!!     Apply the sediment model described in the PlankTOM12.2 manual
!!                           (September 2023)
!! -----------------------------------------------------------------
! These sediment equations are also used in the water column model, so
! WHEN YOU CHANGE THE EQUATIONS THAT FOLLOW, CHANGE THE WATER COLUMN
! MODEL ABOVE AS WELL
! Calculate the sediment remineralisation as a function of T, BAC, and
! OXY, using values in overlying water (because we do not calculate T,
! BAC and OXY in the sediment). This uses the same remineralisation
! formulation and rates as for particles in the open water.
!
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
!
!     Sinking fluxes out of water layer above sediment
!
          jk = max(2,mbathy(ji,jj))-1
!
! First calculate the remineralisation rate for carbon and iron (both
! remid) and silica (siremin)
! -------------------------------------------------------------------
!
          remik = rn_grabac/rjjss*tmask(ji,jj,jk)*tgfunc(ji,jj,jk,jpbac) &
     &      *(trn(ji,jj,jk,jpoxy)+3.E-6)/(10.E-6+trn(ji,jj,jk,jpoxy)) &
     &      *max(trn(ji,jj,jk,jpbac),rn_rembac)
          siremin = min(rn_remdsi*exp(rn_retdsi/(273.15+tsnbio(ji,jj,jk,1))) &
     &      ,rn_readsi)/rjjss*tmask(ji,jj,jk) &
     &      *(trn(ji,jj,jk,jpoxy)+3.E-6)/(10.E-6+trn(ji,jj,jk,jpoxy))
!
! Then apply remineralisation rate to all POM pools in the sediment.
! There is no DOC in the sediment
! ------------------------------------------------------------------
!
          zdenom = rn_kmobac+gbapoc*trnsed(ji,jj,jppoc)                   &
     &      +gbagoc*trnsed(ji,jj,jpgoc)+gbahoc*trnsed(ji,jj,jphoc)
!
          remsed(ji,jj,jppoc) = remik*gbapoc*trnsed(ji,jj,jppoc)/zdenom
          remsed(ji,jj,jpgoc) = remik*gbagoc*trnsed(ji,jj,jpgoc)/zdenom
          remsed(ji,jj,jpgon) = remik*gbagon*trnsed(ji,jj,jpgon)/zdenom
          remsed(ji,jj,jphoc) = remik*gbahoc*trnsed(ji,jj,jphoc)/zdenom
          remsed(ji,jj,jpsfe) = remik*gbapoc*trnsed(ji,jj,jpsfe)/zdenom
          remsed(ji,jj,jpbfe) = remik*gbagoc*trnsed(ji,jj,jpbfe)/zdenom
          remsed(ji,jj,jpufe) = remik*gbahoc*trnsed(ji,jj,jpufe)/zdenom
          remsed(ji,jj,jpdsi) = siremin*trnsed(ji,jj,jpdsi)
!
! Calculate dissolution rates for calcite and aragonite, and apply this
! to the calcite and aragonite particulate pools
! ---------------------------------------------------------------------
!
          omecal = co3(ji,jj,jk)*calcon/aksp(ji,jj,jk)
          remsed(ji,jj,jpcal) = max(trnsed(ji,jj,jpcal)*rn_lyscal/rjjss*(1.-omecal),0.)**rn_lyoco3
          omeara = co3(ji,jj,jk)*calcon/aksara(ji,jj,jk)
          remsed(ji,jj,jpara) = max(trnsed(ji,jj,jpara)*rn_lysara/rjjss*(1.-omeara),0.)**rn_lyoco3
!
! Calculate the sinking rate of POM pools
! ---------------------------------------
!
          xdens = (trn(ji,jj,jk,jpgoc)*moworg+trn(ji,jj,jk,jpcal)*mowco3 &
     &     +trn(ji,jj,jk,jpara)*mowco3                                  &
     &     +trn(ji,jj,jk,jpdsi)*mwsio2)*1e6/                            &
     &     MAX(trn(ji,jj,jk,jpgoc)*moworg/rhoorg                        &
     &     +trn(ji,jj,jk,jpcal)*mowco3/rhoco3                           &
     &     +trn(ji,jj,jk,jpara)*mowco3/rhoco3                           &
     &     +trn(ji,jj,jk,jpdsi)*mwsio2/rhsio2,rtrn)-rhop(ji,jj,jk)*1000.
          xvsink(ji,jj,jk) = MIN(rn_sisgoc*MAX(xdens,dnsmin)**rn_siegoc,snkmax(jk))
          snkspd = xvsink(ji,jj,jk)/rjjss
          snkspd(jppoc) = rn_snkpoc/rjjss
          snkspd(jpsfe) = rn_snkpoc/rjjss
          snkspd(jphoc) = rn_snkhoc/rjjss                                 
          snkspd(jpufe) = rn_snkhoc/rjjss
!
! Calculate the concentration of the POM pools in the sediment (trnsed),
! which is the old pool minus remineralisation plus sinking, and takes
! it out of the above water. Decrease the remnineralisation (remsed) if
! there is not enough mass in the sediment
! ----------------------------------------------------------------------
!
          DO jl = jpdsi, jphoc
            sedflx = snkspd(jl)*trn(ji,jj,jk,jl)*tmask(ji,jj,jk)/fse3t(ji,jj,jk)
            remsed(ji,jj,jl) = max(min((trnsed(ji,jj,jl)-rn_gramin)*rfactr  &
     &        +sedflx,remsed(ji,jj,jl)),0.)
            trn(ji,jj,jk,jl) = trn(ji,jj,jk,jl)-sedflx*rfact
            trnsed(ji,jj,jl) = trnsed(ji,jj,jl)+(sedflx-remsed(ji,jj,jl))*rfact
          END DO
!
! Nutrients are instantaneously added to the overlying water
! ----------------------------------------------------------
!
          trn(ji,jj,jk,jpsil) = trn(ji,jj,jk,jpsil)+remsed(ji,jj,jpdsi)*rfact
          trn(ji,jj,jk,jppo4) = trn(ji,jj,jk,jppo4)+(remsed(ji,jj,jppoc) &
     &      +remsed(ji,jj,jpgon)*ratc2n+remsed(ji,jj,jphoc))*rfact
          remik = remsed(ji,jj,jppoc)+remsed(ji,jj,jpgoc)+remsed(ji,jj,jphoc)
          sedflx = min(remik*rato2c,trn(ji,jj,jk,jpoxy)*rfactr)
          trn(ji,jj,jk,jpoxy) = trn(ji,jj,jk,jpoxy)-sedflx*rfact
!
! in contrast to the water column, denitrification starts when O2=0
          trn(ji,jj,jk,jpdin) = trn(ji,jj,jk,jpdin)+(0.8*(sedflx-remik*rato2c) &
     &      +(remsed(ji,jj,jppoc)+remsed(ji,jj,jphoc))*ratn2c           &
     &      +remsed(ji,jj,jpgon))*rfact
          trn(ji,jj,jk,jpfer) = trn(ji,jj,jk,jpfer)+(remsed(ji,jj,jpsfe) &
     &      +remsed(ji,jj,jpbfe)+remsed(ji,jj,jpufe))*rfact
          trn(ji,jj,jk,jpdic) = trn(ji,jj,jk,jpdic)+(remik              &
     &      +remsed(ji,jj,jpcal)+remsed(ji,jj,jpara))*rfact
#  if defined key_trc_piic
!
! apply sediment model to pre-industrial DIC
          trn(ji,jj,jk,jppiic) = trn(ji,jj,jk,jppiic)+(remik            &
     &      +remsed(ji,jj,jpcal)+remsed(ji,jj,jpara))*rfact
#  endif
          trn(ji,jj,jk,jptal) = trn(ji,jj,jk,jptal)+(-alknut*remik      &
     &      +2.*remsed(ji,jj,jpcal)+2.*remsed(ji,jj,jpara))*rfact
        END DO
      END DO
#if defined key_iomput
      WHERE (tmask(:,:,:) .EQ. 0 ) 
        out3d = ncf_fill
      ELSEWHERE
        out3d = xvsink/rjjss
      END WHERE
!
      CALL iom_put("vsink",    out3d )
      CALL iom_put("remdsised",remsed(:,:,jpdsi)*1e3)
#endif
 
#endif
      RETURN
      END SUBROUTINE bgcsnk
