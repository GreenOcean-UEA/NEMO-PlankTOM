MODULE trcini_planktom
#if defined key_planktom
!!======================================================================
!!                         ***  MODULE trcini_planktom  ***
!! TOP : PlankTOM ecosystem model 
!!======================================================================
!! History : original  : 1988-07  E. MAIER-REIMER      MPI HAMBURG
!!           additions : 1999-10  O. Aumont and C. Le Quere
!!           additions : 2001-03  0. Aumont and E. Kastenare : Changes in 
!!             computation of the export profiles for silicate and calcite
!!           deletions : 2002     Erik Buitenhuis
!!           see bgcbio
!!----------------------------------------------------------------------
!!    'key_planktom'      :             use the PlankTOM ecosystem model
!!    'key_trc_n2o'       :                        calculates N2O fluxes
!!    'key_trc_ch4'       :                        calculates CH4 fluxes
!!    'key_trc_piic'      :    additional tracer with pre-industrial DIC 
!!----------------------------------------------------------------------
!! ** Purpose : Initialisation of biological and chemical variables
!! ** Action  : 
!! 
!! ** External: 
!!              SEAICE
!!              WSJCD
!!              WSEK
!!              RHO
!!
!! ** Method  :
!!         1) SET CONSTANTS FOR CARBONATE CHEMISTRY AS DESCRIBED IN
!!            IN BROECKER ET AL. (1982, GEOSECS) AND EDMOND A. GIESKES
!!            (1970)
!!         2) INITIATE [CO3--] AND PH-VALUE BY ITERATION
!!            (NEWTON-RAPHSON METHOD FOR SOLVING NONLINEAR SIMULTANEOUS
!!             EQUATIONS, SEE E.G. SCARBOROUGH, J. (1958))
!!
!! ** References for biology:
!!
!!         DEGENS, E.T, S. KEMPE, AND A. SPITZY (1984)
!!         CARBON DIOXIDE: A BIOGEOCHEMICAL PORTRAIT.
!!         IN: THE HANDBOOK OF ENVIRONMENTAL CHEMISTRY, VOLUME 1/
!!         PART C, O. HUTZINGER, ED., SPRINGER-VERLAG, BERLIN,
!!         HEIDELBERG, PP. 127-215.
!!
!!         DUGDALE. R.C. (1967)
!!         NUTRIENT LIMITATION IN THE SEA: DYNAMICS, IDENTIFICATION
!!         AND SIGNIFICANCE.
!!         LIMNOLOGY AND OCEANOGRAPHY, VOL.12, 685-695.
!!
!!         PARSONS, T.R., AND M. TAKAHASHI (1973)
!!         BIOLOGICAL OCEANOGRAPHIC PROCESSES.
!!         PERGAMON PRESS, 186 PP.
!!
!!         TAKAHASHI, T., W.S. BROECKER, AND S. LANGER (1985)
!!         REDFIELD RATIO BASED ON CHEMICAL DATA FROM ISOPYCNAL
!!         SURFACES.
!!         JOURNAL OF GEOPHYSICAL RESEARCH, 90(C4), 6907-6924.
!!
!! ** References for chemistry:
!!
!!         BERNER, R. A. (1976)
!!         THE SOLUBILITY OF CALCITE AND ARAGONITE IN SEA WATER
!!         AT ATMOSPHERIC PRESSURE AND 34.5 O/OO SALINITY.
!!         AMERICAN JOURNAL OF SCIENCE, VOL. 276, 713-730.
!!         (K'SP(ARAGONITE)=1.45 K'SP(CALCITE))
!!
!!         BROECKER, W.S., D.W. SPENCER, AND H. CRAIG (1982)
!!         GEOSECS PACIFIC EXPEDITION. VOL. 3.. HYDROGRAPHIC DATA
!!         1973-1974, SUPERINTENDANT OF DOCUMENTS, U.S. GOVERNMENT
!!         PRINTING OFFICE, WASHINGTON, D.C., 137 PP..
!!
!!         CULBERSON, C.H., AND R.M. PYTKOWICZ (1968)
!!         EFFECT ON PRESSURE ON CARBONIC ACID, BORIC ACID AND THE PH
!!         IN SEA WATER.
!!         LIMNOLOGY AND OCEANOGRAPHY, VOL. 13, 403-417.
!!
!!         DICKSON, A.G., AND J.P. RILEY (1979)
!!         THE ESTIMATION OF ACID DISSOCIATION CONSTANTS IN SEAWATER
!!         MEDIA FROM POTENTIOMETRIC TITRATIONS WITH STRONG BASE.
!!         I. THE IONIC PRODUCT OF WATER - KW.
!!         MARINE CHEMISTRY, VOL. 7, 89-99.
!!
!!         EDMOND, J.M., AND J.M.T.M. GIESKES (1970)
!!         ON THE CALCULATION OF THE DEGREE OF SATURATION OF SEA WATER
!!         WITH RESPECT TO CALCIUM CARBONATE UNDER IN SITU CONDITIONS.
!!         GEOCHIM. ET COSMOCHIM. ACTA, 34, 1261-1291.
!!
!!         INGLE, S.E. (1800)
!!         SOLUBILITY OF CALCITE IN THE OCEAN.
!!         MARINE CHEMISTRY, VOL. 3, 301-319.
!!
!!         INGLE, S.E., C.H. CULBERSON, J.E. HAWLEY, AND R.M. PYTKOWICZ
!!         (1973) THE SOLUBILITY OF CALCITE IN SEAWATER AT ATMOSPHERIC
!!         PRESSURE AND 35 O/OO SALINITY.
!!         MARINE CHEMISTRY, VOL. 1, 295-307.
!!
!!         RILEY, J. P., AND G. SKIRROW, EDS. (1965)
!!         CHEMICAL OCEANOGRAPHY. VOL. 1, 712 PP., ACADEMIC PRESS,
!!         LONDON A. NEW YORK.
!!
!!         SCARBOROUGH, J. (1958) NUMERICAL MATHEMATICAL ANALYSIS.
!!         OXFORD UNIVERSITY PRESS, LONDON, 4TH ED., 576 PP..
!!
!!         WEISS, R. F. (1970) THE SOLUBILITY OF NITROGEN
!!         OXYGEN AND ARGON IN WATER AND SEAWATER.
!!         DEEP-SEA RESEARCH, VOL. 17, 721-735.
!!
!!         WEISS, R. F. (1974)
!!         CARBON DIOXIDE IN WATER AND SEAWATER: THE SOLUBILITY OF A
!!         NON IDEAL GAS. MARINE CHEMISTRY, VOL. 2, 203-215.
!!
!!         WOOSTER, W.S., A.J. LEE, AND G. DIETRICH (1969)
!!         REDEFINITION OF SALINITY. Z. GEOPHYS., VOL.35, 611-613.
!!
!!         BROECKER, W.S., D.W. SPENCER, AND H. CRAIG (1982)
!!         GEOSECS PACIFIC EXPEDITION. VOL. 3.. HYDROGRAPHIC DATA
!!         1973-1974, SUPERINTENDANT OF DOCUMENTS, U.S. GOVERNMENT
!!         PRINTING OFFICE, WASHINGTON, D.C., 137 PP..
!!
!!!---------------------------------------------------------------------
   CONTAINS
      subroutine trc_ini_planktom
!
! parameters and commons
   USE iom
   USE oce_trc
   USE phycst
   USE sms_planktom
   USE trc
   USE trp_trc
   USE trcrst
   USE trcdta    
   USE zpshde          ! partial step: hor. derivative 
   USE in_out_manager  ! I/O manager
   USE prtctl_trc      ! Print control passive tracers (prt_ctl_trc_init routine)
   USE lib_mpp         ! distributed memory computing library
   
   IMPLICIT NONE
#  include "domzgr_substitute.h90" 
! 
! local variables
      INTEGER mo, numdust, ji, jj, jk, jl, jm, ierr, ik, yearrc
      REAL ztest, sum, sum1, dumco2
      REAL denitide, expide
      REAL(wp) :: zsecond, zdate0
      INTEGER , DIMENSION (1) :: istep0
      INTEGER , PARAMETER :: jpmois = 12, jpan = 1
      INTEGER , DIMENSION (jpmois) :: istep
      INTEGER :: numriv, numbath, ipi, ipj, ipk, itime
      INTEGER :: numatm
      REAL(wp), DIMENSION (jpk) :: zlev
      REAL(wp) , DIMENSION (jpi,jpj) :: zlon, zlat
      CHARACTER (len=34) :: clname
      REAL sumdic, sumdoc, sumpoc, sumnit, sumpo4, sumsil, sumfer, sumfes, sumani
      REAL riverdic(jpi,jpj), riverdoc(jpi,jpj), riverpoc(jpi,jpj)
      REAL rivernit(jpi,jpj), riverpo4(jpi,jpj), riversil(jpi,jpj)
      REAL riverfer(jpi,jpj)
      REAL atmosdin(jpi,jpj)
      REAL cmask(jpi,jpj,jpk)
      REAL sumdepsi, sumdepfe
      REAL maxlig, ligz, bot, h3o, dic, calk, ah2
      REAL totbio,biowprf
      INTEGER numtemp, numsal  !  files for reading in observed temp and salinity
!
! allocate PlankTOM arrays
      ierr =         sms_planktom_alloc()          
!
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'planktom_alloc: unable to allocate PlankTOM arrays' )
!
! Initialization
! --------------
!
! initialisation of local variables
      ferat3 = 2.0e-6
      rfact  = rdt !* float(ndttrc)
      rfactr = 1./rfact
      IF(lwp) WRITE(numout,*) ' Tracer time step=',rfact,' rdt=',rdt
      IF(lwp) WRITE(numout,*) ' ferat3 =',ferat3 
      rjjss  = rjjhh * rhhmm * rmmss                ! number of seconds in day
      raajj  = REAL(nyear_len(1),wp)              ! number of days / year
      raass  = raajj * rjjss
      rmoss  = raass/12.                            ! number of seconds in year

!
! initialise dust input from the atmosphere
      CALL iom_open('dust.orca.nc',numdust)
      DO mo = 1, jpmois
        CALL iom_get(numdust,jpdom_data,'DEP',dustmo(:,:,mo),mo)
      END DO
      CALL iom_close(numdust)
!!
!!     -------------------------------------------------
!!     Computation of the total atmospheric supply of Si
!!     -------------------------------------------------
!!
!     Iron and Si deposition at the surface
!     the dust variable is in kgdust/m2/s (from Jickells et al. 2005).
!     We use 0.035gFe/gdust, 0.308 gSi/gdust (or 8.8gSi/gFe).
!     The solubility of Fe in dust is rn_fersol (usually 2%, set in namelist).
!     The solubility of Si in dust is 7.5%.
!
!     for Mahowald et al. 2009:
!     the dust variable is in ngFe(soluble)/m2/s 
!     We use 8.8gSi/gFe (0.308/0.035).
!     rn_fersol (set in namelist) should be = 1 and is used as a scaling factor.
!     The solubility of Si in dust is 7.5%. We use a mean solubility of Fe of 
!                                     1.25% to scale the Si flux.
!
!     Variables irondep and sidep are in mol/L/time_step for Fe and Si, respectively.
!     The molecular weight of Fe and Si is 55.85 and 28.01, respectively.
!
!     -------------------------------------------------------------------
!
      DO jk = 1, jpk
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
            irondep(ji,jj,jk) = 0.
            sidep(ji,jj,jk)   = 0.
          END DO
        END DO
      END DO
!
      sumdepfe = 0.
      sumdepsi = 0.
!
! use this if dustmo is dust depostion
!       DO mo = 1, 12
!         DO jj = 2, nlcj-1
!           DO ji = 2, nlci-1
!             sumdepfe=sumdepfe+dustmo(ji,jj,mo)*0.035*rn_fersol/55.85* &
!      &                           e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,1)
!             sumdepsi=sumdepsi+dustmo(ji,jj,mo)*0.308*rn_silsol/28.01 *    &
!      &                           e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,1)
!           END DO
!         END DO
!       END DO
!
! 'use this if dustmo is iron in dust deposition
      DO mo = 1, 12
       DO jk = 1, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
            sumdepfe = sumdepfe+dustmo(ji,jj,mo)*rn_fersol/55.85         &
     &        *(gdept_0(ji,jj,jk)/gdept_0(ji,jj,1))**(-0.858)*0.3314     &
     &        *e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,1)
          END DO
        END DO
       END DO
      END DO
      sumdepsi = sumdepfe*(0.308/0.035)*(rn_silsol/rn_fersol)*(55.85/28.01)
!!
!!    -------------------------------------------------------
!!    Computation of the coastal mask.
!!    Computation of an island mask to enhance coastal supply
!!    of iron
!!    -------------------------------------------------------
!!
      CALL iom_open('bathy.orca.nc',numbath)
      CALL iom_get(numbath,jpdom_data,'bathy',cmask(:,:,:),jpan)
!     volumt = 0.
      DO jk = 1, jpk
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
            expide   = min(8.,(gdept_0(ji,jj,jk)/500.)**(-1.5))
            denitide = -0.9543+0.7662*log(expide)-0.235*log(expide)**2
            cmask(ji,jj,jk)  = cmask(ji,jj,jk)*exp(denitide)/0.6858
            volumt(ji,jj,jk) = e1t(ji,jj) * e2t(ji,jj) &
     &        * fse3t(ji,jj,jk) * tmask(ji,jj,jk) *1000.
          END DO
        END DO
      END DO
!!
!!    ---------------------------------------
!!    Initialise the nutrient input by rivers
!!    ---------------------------------------
!!
      CALL iom_open('river.nc',numriv)
      CALL iom_get(numriv,jpdom_data,'RIVERDIC',riverdic(:,:),1)
      CALL iom_get(numriv,jpdom_data,'RIVERDOC',riverdoc(:,:),1)
      CALL iom_get(numriv,jpdom_data,'RIVERPOC',riverpoc(:,:),1)
      CALL iom_get(numriv,jpdom_data,'RIVERNIT',rivernit(:,:),1)
      CALL iom_get(numriv,jpdom_data,'RIVERPHOS',riverpo4(:,:),1)
      CALL iom_get(numriv,jpdom_data,'RIVERSIL',riversil(:,:),1)
      CALL iom_get(numriv,jpdom_data,'RIVERFe80',riverfer(:,:),1)
      CALL iom_close(numriv)
      IF(lwp) WRITE(numout,*) ' rn_rivdic=',rn_rivdic,' rfact=',rfact
      IF(lwp) WRITE(numout,*) nlcj, nlci, jpi, jpj
!
!     convert river file [Tg/gridcell/y] to [mol/L/timestep]
!     *1e12 g/Tg /1000. L/m3 * rfact s/timestep /
!     (12 g/mol * raass s/y * e1t*e2t*e3t m3/gridcell)
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
          depdic(ji,jj,1) = rn_rivdic*riverdic(ji,jj)*       &
     &      1E9*rfact*tmask(ji,jj,1)/                        &
     &      (12.*raass*e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,1))
          depdoc(ji,jj,1) = rn_rivdoc*riverdoc(ji,jj)*       &
     &      1E9*rfact*tmask(ji,jj,1)/                        &
     &      (12.*raass*e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,1))
          deppoc(ji,jj,1) = rn_rivpoc*riverpoc(ji,jj)*       &
     &      1E9*rfact*tmask(ji,jj,1)/                        &
     &      (12.*raass*e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,1))
          depnit(ji,jj,1) = rn_rivnit*rivernit(ji,jj)*       &
     &      1E9*rfact*tmask(ji,jj,1)/                        &
     &      (14.*raass*e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,1))
          deppo4(ji,jj,1) = rn_rivpo4*riverpo4(ji,jj)*       &
     &      122.*1E9*rfact*tmask(ji,jj,1)/                   &
     &      (31.*raass*e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,1))
          depsil(ji,jj,1) = rn_rivsil*riversil(ji,jj)*       &
     &      1E9*rfact*tmask(ji,jj,1)/                        &
     &      (28.1*raass*e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,1))
          depfer(ji,jj,1) = (rn_rivfer*riverfer(ji,jj)       &
     &      +rn_rivdoc*riverdoc(ji,jj)*ferat3)*              &
     &      1E9*rfact*tmask(ji,jj,1)/                        &
     &      (55.8*raass*e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,1))
        END DO
      END DO
!!
!!    ---------------------------------------------------------------------
!!    Initialise NO3 deposition from the atmosphere (from Duce et al. 2008)
!!    original data in mgN/m2/y
!!    ---------------------------------------------------------------------
!!
      CALL iom_open('ndeposition.nc',numatm)
      CALL iom_get(numatm,jpdom_data,'ndep',atmosdin(:,:),1)
      CALL iom_close(numatm)
! 
!     atmospheric deposition in molN/L (/1000 mgN/g /1000 L/m3)
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
          DO jk = 1, jpk
            atmdin(ji,jj,jk) = 0.
          END DO
        END DO
      END DO
!
      sumdic = 0.
      sumdoc = 0.
      sumpoc = 0.
      sumnit = 0.
      sumpo4 = 0.
      sumsil = 0.
      sumfer = 0.
      sumani = 0.
!
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
          sumdic = sumdic+rn_rivdic*riverdic(ji,jj)*tmask(ji,jj,1)
          sumdoc = sumdoc+rn_rivdoc*riverdoc(ji,jj)*tmask(ji,jj,1)
          sumpoc = sumpoc+rn_rivpoc*riverpoc(ji,jj)*tmask(ji,jj,1)
          sumnit = sumnit+rn_rivnit*rivernit(ji,jj)*tmask(ji,jj,1)
          sumpo4 = sumpo4+rn_rivpo4*riverpo4(ji,jj)*tmask(ji,jj,1)
          sumsil = sumsil+rn_rivsil*riversil(ji,jj)*tmask(ji,jj,1)
          sumfer = sumfer+rn_rivfer*riverfer(ji,jj)*tmask(ji,jj,1)
          sumani = sumani+atmosdin(ji,jj)*tmask(ji,jj,1)
        END DO
      END DO
      sumfes = 0.
!
!     rn_sedfer is in mol/m2/d 
      DO jk = 1, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
          depfer(ji,jj,jk) = depfer(ji,jj,jk)+rn_sedfer*cmask(ji,jj,jk) &
     &      /(fse3t(ji,jj,jk)*rjjss)*rfact
          sumfes = sumfes+rn_sedfer*cmask(ji,jj,jk) &
     &      *raajj*e1t(ji,jj)*e2t(ji,jj)
          END DO
        END DO
      END DO
!
      IF( lk_mpp ) THEN
        CALL mpp_sum( sumdic )
        CALL mpp_sum( sumdoc )
        CALL mpp_sum( sumpoc )
        CALL mpp_sum( sumnit )
        CALL mpp_sum( sumpo4 )
        CALL mpp_sum( sumsil )
        CALL mpp_sum( sumdepsi )
        CALL mpp_sum( sumdepfe )
        CALL mpp_sum( sumfer )
        CALL mpp_sum( sumfes )
        CALL mpp_sum( sumani )
      ENDIF
!
! convert from Tg/y to mol/timestep
      extinp(1) = (sumsil/28.1*1e12/raass+sumdepsi*1000.)*rfact
      extinp(2) = (sumpoc/12.+sumpo4/31.*122.)*1e12/raass*rfact
      extinp(3) = (((sumfer/55.8+(sumpoc+sumdoc)/12.*ferat3)*1e12+sumfes)/raass+sumdepfe*1000.)*rfact
      extinp(4) = sumdoc*1e12/12./raass*rfact
      extinp(5) = sumdic/12.*1e12/raass*rfact/2.
      extinp(6) = sumnit/14.*1e12/raass*rfact
      extinp(7) = (sumdic/12.-sumpo4/31.*122.)*1e12/raass*rfact
!
      IF (lwp) THEN
        write(numout,*) ' river inputs DIC, DOC, POC, NO3, PO4, Si, Fe [Tg/y]' &
     &                   ,sumdic,sumdoc,sumpoc,sumnit,sumpo4,sumsil,sumfer
        write(numout,*) ' dust inputs Si, Fe [kmol/s] ',sumdepsi,sumdepfe
        write(numout,*) ' sediment input Fe [mol/y]',sumfes, rn_sedfer
        write(numout,*) ' atmospheric input N [mol N/y]',sumani
      ENDIF

      IF (lwp) WRITE(numout,*) ' corrections of Si, POC, Fe, DOC, CO3, PON, DIC [mol/timestep] ' &
     &,extinp(1:7)
      IF (lwp) write(numout,31) "corrections in Gmol/timestep:", (extinp(jl)*1.e-9,jl=1,7)
31    FORMAT(a12,7f10.3)
      sedcor = 0.
!
! Parameters for ballasting of sinking material
! 24, 240, 100, 60.1 are the molar weight of dry organic material, wet org mat, CaCO3, SiO2
! 1.08, 1.335, 1.2 are the density of organic material,CaCO3,SiO2
! rn_siegoc and rn_sisgoc are derived by fitting sinking speeds 
! calculated with the sedimentation function of Buitenhuis et al. 2001
! to an exponential function 
!------------------------------------------------------------------------------------------
!
      moworg  = 240.
!     moworg  = 24.
      mowco3  = 100.
      mwsio2  = 60.1
      rhoorg  = 1.08
      rhoco3  = 1.335
      rhsio2  = 1.2
!
      DO jk = 1, jpk
         snkmax(jk) = min(0.9*e3t_1d(jk)*rday/rfact,rn_snkmax)
      END DO
!!
!!    --------------------------------------------------
!!    Read in 4D observed temperature and salinity files
!!    --------------------------------------------------
!!
      IF(lwp) THEN
        WRITE(numout,*)
        WRITE(numout,*) 'biodts : temperature diagnostic for ecosystem and chemical processes'
        WRITE(numout,*) '~~~~~~~~~~~~~'

        SELECT CASE ( nn_biodts )  !  this value is specified in namelist.trc.sms
          CASE ( 0 )  ;  WRITE(numout,*) ' nn_biodts = 0 : using default modelled T and S'
          CASE ( 1 )  ;  WRITE(numout,*) ' nn_biodts = 1 : biological subroutines use obs T and S'
          CASE ( 2 )  ;  WRITE(numout,*) ' nn_biodts = 2 : chemical subroutines use obs T and S'
          CASE ( 3 )  ;  WRITE(numout,*) ' nn_biodts = 3 : bio. and chem. subroutines use obs T and S'
        END SELECT

        WRITE(numout,*)
      ENDIF

      CALL iom_open ( 'data_1m_potential_temperature_nomask.nc', numtemp )
      DO mo = 1, jpmois
        CALL iom_get  ( numtemp, jpdom_data, 'votemper', obstempmo(:,:,:,mo), mo )
      END DO
      CALL iom_close( numtemp )

      CALL iom_open ( 'data_1m_salinity_nomask.nc', numsal )
      DO mo = 1, jpmois
        CALL iom_get  ( numsal, jpdom_data, 'vosaline', obssalmo(:,:,:,mo), mo )
      END DO
      CALL iom_close( numsal )

      IF(lwp) WRITE(numout,*)
!!
!!--------------------------------
!! Initialize biological variables 
!!--------------------------------
!
! Set biological ratios
! ---------------------
!
      alknut = (16.+6.+1.)/122. ! N+S+P Wolf-Gladrow et al. 2007
      minfer = 1e-15
      ratn2c = 16./122.
      ratc2n = 122./16.
      rato2c = 172./122.
!
      IF(lwp) WRITE(numout,*) ' minfer =',minfer 
!
! Set fractionation factors of 13C, 14C
! -------------------------------------
!
      pdb     = 0.011112
      plafr13 = 0.980
!
! Assign additional variables based on namelist options
! -------------------------------------------------------------
!
      dnsmin = (rn_snkpoc/rn_sisgoc)**(1/rn_siegoc)
      IF(lwp) WRITE(numout,*) ' dnsmin =', dnsmin
! 
! Assign grazing preferences (prfzoo) based on namelist options
! -------------------------------------------------------------
! 
! totbio   :  sum of biomasses rn_bmspft(:) where rn_rpfzoo(:,:) =/= 0
! biowprf  :  sum of ( rn_bmspft(:) * rn_rpfzoo(:,:) ) for each PFT
!
      jpfoo = 0
!
      DO jm = 1, jpzft
!
! initialise variables
        totbio  = 0.
        biowprf = 0.
!
! calculate total biomass and total available biomass for each PFT
        DO jl = jppoc,jpdia+jppft-1
          IF (rn_prfzoo(jm,jl) .GT. rtrn) THEN
            totbio  = totbio+rn_bmspft(jl)
            biowprf = biowprf+rn_bmspft(jl)*rn_prfzoo(jm,jl)

            jpfoo = jpfoo+1
            IF (jpfoo .GT. jpfmx) STOP 'increase jpfmx'
            grizoo(1,jpfoo) = jm
            grizoo(2,jpfoo) = jl
          ENDIF
        END DO
!
! calculate preference of each zPFT for each food source
        DO jl = jppoc,jpdia+jppft-1
          prfzoo(jm,jl) = totbio/biowprf*rn_prfzoo(jm,jl)
        END DO
      END DO
! 
! Assign the preferences of bacteria remineralisation (gbaXXX) for each detrital pool 
! -----------------------------------------------------------------------------------
! 
! calculate total biomass and total available biomass for bacteria
      totbio  = rn_bmsdet(1)+rn_bmsdet(2)+rn_bmsdet(3)+rn_bmsdet(4)
      biowprf = rn_bmsdet(1)*rn_gbadoc+rn_bmsdet(2)*rn_gbapoc+rn_bmsdet(3)*rn_gbagoc+rn_bmsdet(4)*rn_gbahoc
!
! calculate preference of bacteria for each organic matter source
      gbadoc = totbio/biowprf*rn_gbadoc
      gbapoc = totbio/biowprf*rn_gbapoc
      gbagoc = totbio/biowprf*rn_gbagoc
      gbahoc = totbio/biowprf*rn_gbahoc
      gbagon = totbio/biowprf*rn_gbagon
!!
!!------------------------------
!! Initialize chemical variables 
!!------------------------------
!!
! Set pre-industrial o2/n2 ratio
! ------------------------------
!
      atcox = 0.20946
!
! Set half precision constants
! ----------------------------
!
      smicr  = 1.E-6
      thousi = 1./1000.
      perc   = 0.01
      third  = 1./3.
!
! Set lower/upper limits for temperature and salinity
! ---------------------------------------------------
!
      salchl = 1./1.80655
      temzer = 273.16
      calcon = 1.03E-2
!
! Set coefficients for apparent solubility equilibrium
!   of calcite (Ingle, 1800, eq. 6)
! ----------------------------------------------------
!
      akcc1 = -34.452
      akcc2 = -39.866
      akcc3 = 110.21
      akcc4 = -7.5752E-6
!
      arafra = 0.
      calfra = 1.-arafra
      aracal = arafra*1.45+calfra
!
! Set coefficients for seawater pressure correction
! -------------------------------------------------
!
      devk1  = 24.2
      devk2  = 16.4
      devkb  = 27.5
      devk1t = 0.085
      devk2t = 0.04
      devkbt = 0.095
!
      devkst = 0.23
      devks  = 32.8*arafra+35.4*calfra
!
! Set universal gas constants
! ---------------------------
!
      rgas = 83.143
      oxyco = 1./22.4144
!
! Set boron constants
! -------------------
!
      bor1 = 0.00023
      bor2 = 1./10.82
!
! Set volumetric solubility constants for co2 in ml/l (Weiss, 1974)
! -----------------------------------------------------------------
!
      c00 = -58.0931
      c01 = 90.5069
      c02 = 22.2940
      c03 = 0.027766
      c04 = -0.025888
      c05 = 0.0050578
!
! Set coeff. for 1. dissoc. of carbonic acid (Edmond and Gieskes, 1970)
! ---------------------------------------------------------------------
!
      c10 = -2307.1266
      c11 = 2.83655
      c12 = -1.5529413
      c13 = -4.0484
      c14 = -0.20760841
      c15 = 0.08468345
      c16 = -0.00654208
      c17 = -0.001005
!
! Set coeff. for 2. dissoc. of carbonic acid (Edmond and Gieskes, 1970)
! ---------------------------------------------------------------------
!
      c20 = -3351.6106
      c21 = -9.226508
      c22 = -0.2005743
      c23 = -23.9722
      c24 = -0.106901773
      c25 = 0.1130822
      c26 = -0.00846934
      c27 = -0.001005
!
! Set coeff. for 1. dissoc. of boric acid (Edmond and Gieskes, 1970)
! ------------------------------------------------------------------
!
      cb0  = -8966.90
      cb1  = -2890.53
      cb2  = -77.942
      cb3  = 1.728
      cb4  = -0.0996
      cb5  = 148.0248
      cb6  = 137.1942
      cb7  = 1.62142
      cb8  = -24.4344
      cb9  = -25.085
      cb10 = -0.2474
      cb11 = 0.053105
!
! Set coeff. for dissoc. of water (Dickson and Riley, 1979, 
! eq. 7, coefficient cw2 corrected from 0.9415 to 0.09415 
! after pers. commun. to B. Bacastow, 1988)
! ---------------------------------------------------------
!
      cw0 = -13847.26
      cw1 = 148.9652
      cw2 = -23.6521
      cw3 = 118.67
      cw4 = -5.977
      cw5 = 1.0495
      cw6 = -0.01615
!
! Set volumetric solubility constants for o2 in ml/l (Weiss, 1970)
! ----------------------------------------------------------------
!
      ox0 = -58.3877
      ox1 = 85.8079
      ox2 = 23.8439
      ox3 = -0.034892
      ox4 = 0.015568
      ox5 = -0.0019387
!      
! Calculate chemical tracers
! --------------------------
!
      CALL bgcche
!
#if defined key_trc_atmco2
!
! Read atmospheric trace gas concentrations
! -----------------------------------------
!
!
! read atmospheric CO2 concentration
      nutatm = 73
      IF (lwp) WRITE (numout,*) ' Reading real atmospheric CO2 option '
      OPEN(nutatm,FILE='atmco2.dat')
!
      DO ji = 1, nmaxrec
#  if defined key_trc_ch4
!
! read atmospheric CH4 concentration
        read(nutatm, *,END=152) yrco2(ji),sipco2(ji),sipn2o(ji),sipch4(ji)
#  else
#    if defined key_trc_n2o
!
! read atmospheric N2O concentration
        read(nutatm, *,END=152) yrco2(ji),sipco2(ji),sipn2o(ji)
#    else
        read(nutatm, *,END=152) yrco2(ji),sipco2(ji)
#    endif
#  endif
      END DO
!
! read atmospheric CFC concentration
152   nutatm = 1173
      IF (lwp) WRITE (numout,*) ' Reading atmospheric CFC '
      OPEN(nutatm,FILE='atmco2cfc.dat')
      DO ji = 1, 336
        read(nutatm, *) yrcfc(ji),dumco2,sipcfc(ji,1,:),sipcfc(ji,2,:),sipcfc(ji,3,:)
      END DO
      yearrc = ndastp/10000-int(yrcfc(1))+1
      IF(lwp) WRITE(numout,*) ' atm CFC11',ndastp,yrcfc(yearrc),sipcfc(yearrc,1,:)
!
! Interpolate atmospheric CO2 concentrations
! ------------------------------------------
!
      CALL bgcatm
!
#endif
!
! Calculate the temperature-dependence of PFT growth rates
! --------------------------------------------------------
!
      DO jl = jpdia, jpdia+jppft-1
        pcmax(jl) = rn_mumpft(jl)*(1.+rn_resphy(jl))/rjjss
      END DO
!
! Fix the sub-surface iron ligand concentration
! ---------------------------------------------
!
      DO jj = 2, nlcj-1
!
! between rn_liglat-10 degrees South and rn_liglat degrees South, 
! the subsurface ligand concentration changes from rn_ligfrn to rn_ligfrs
        maxlig = rn_ligfrs+min(max(gphit(2,jj)+rn_liglat,0.),10.)/10.*(rn_ligfrn-rn_ligfrs)
!
! ligz increases from rn_ligdpn (depth in the north) to rn_ligdps (depth in the south)
        ligz  = rn_ligdpn+min(max(-gphit(2,jj)-rn_liglat,0.),10.)/10.*(rn_ligdps-rn_ligdpn)
!
        DO jk = 1, jpk
!
! ligand concentration (ligfer) is set to 0 at depths shallower than ligz, 
! and to its (latitude-dependent) maximum below 2*ligz
          ligfer(jj,jk) = max(min(gdept_1d(jk)-ligz,ligz),0.)/ligz*maxlig
        END DO
      END DO
!
! Initialise the fraction of phyto growth that is supported by NH4 to 1
! ---------------------------------------------------------------------
!
      DO jl = jpdia, jpdia+jppft-1
        DO jk = 1, jpk
          DO jj = 1, jpj
            DO ji = 1, jpi
              dinpft(ji,jj,jk,jl) = tmask(ji,jj,jk)
            END DO
          END DO
        END DO
      END DO
!
! Identify the bottom depth
! -------------------------
!
      DO jj = 1, jpj
        DO ji = 1, jpi
          ik = mbathy(ji,jj)
          IF( ik > 0 ) mdept(ji,jj) = gdept_0(ji,jj,ik)
        END DO
      END DO
!
! Assign the initial values for the carbonate variables 
! -----------------------------------------------------
!
      hi = 1e-8
#  if defined key_trc_piic
      pihi = 1e-8
#  endif
      DO jk = 1,jpkm1
         DO jj = 1,jpj
            DO ji = 1,jpi
              DO jl = 1, 10
                bot  = borat(ji,jj,jk)
                h3o  = amax1(hi(ji,jj,jk),1.E-10)
                dic  = trn(ji,jj,jk,jpdic)
                calk = trn(ji,jj,jk,jptal) -                                      &
     &                 (akw3(ji,jj,jk)/h3o-h3o+bot/(1.+h3o/akb3(ji,jj,jk)))
                ah2  = sqrt((dic-calk)**2+4.*(calk*ak23(ji,jj,jk)/ak13(ji,jj,jk)) &
     &                *(2*dic-calk))
                hi(ji,jj,jk)=0.5*ak13(ji,jj,jk)/calk*((dic-calk)+ah2)
#  if defined key_trc_piic
                h3o  = amax1(pihi(ji,jj,jk),1.E-10)
                dic  = trn(ji,jj,jk,jppiic)
                calk = trn(ji,jj,jk,jptal)- &
     &                 (akw3(ji,jj,jk)/h3o-h3o+bot/(1.+h3o/akb3(ji,jj,jk)))
                ah2  = sqrt((dic-calk)**2+4.*(calk*ak23(ji,jj,jk)/ak13(ji,jj,jk)) &
     &                *(2*dic-calk))
                pihi(ji,jj,jk) = 0.5*ak13(ji,jj,jk)/calk*((dic-calk)+ah2)
#  endif
              END DO
            END DO
         END DO
      END DO
!
! Initialise further ecosystem and output variables
! -------------------------------------------------
!
      denitr = 0.
      out2d  = 0.
      out3d  = 0.
      xvsink = 0.
!
! Don't calculate grazing terms if prfzoo==0.
      grazoo = 0.
      delo2  = 0.
!
! Set sinking into surface layer to 0.
      snkara  = 0.
      snkbfe  = 0.
      snkcal  = 0.
      snkdsi  = 0.
      snkgoc  = 0.
      snkgon  = 0.
      snkhoc  = 0.
      snkpoc  = 0.
      snksfe  = 0.
      snkufe  = 0.
      ngochoc = maxval(nn_sizzoo)
      remsed  = 0.
      stofoo(:,:,:,:,1) = 1.     ! C:C ratio
      stofoo(:,:,:,:,2) = ferat3 ! Fe:C ratio of tracers that don't have explicit Fe
!
#  if defined key_trc_ch4
      proch1 = 0.
      proch2 = 0.
      proch3 = 0.
      proch4 = 0.
      proch5 = 0.
#  endif
      ncf_fill = -1e35
      IF(lwp) WRITE(numout,*) ' Initialisation of PlankTOM done',jpi,jpj,jpk,numatm
      return
      end subroutine
#else
   !!----------------------------------------------------------------------   
   !!  Empty module :                                     NO PlankTOM model
   !!----------------------------------------------------------------------
#endif
   
   !!======================================================================  
END MODULE trcini_planktom
