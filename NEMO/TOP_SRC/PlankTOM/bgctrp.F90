      SUBROUTINE bgctrp(kt,kindic)
#    if defined key_planktom && defined key_trc_diaadd && defined key_nonexist
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE bgctrp
!!!                     ******************
!!!
!!!  Purpose :
!!!  ---------
!!!     output of sinking rates where sediment trap data are available
!!!
!!   Method :
!!   -------
!!      At the beginning of the first time step (nit000) read
!!	 coordinates of trap data
!!      At each time step calculate the mean fluxes
!!	Each nwritesed time steps, output the mean fields
!!
!!   Input :
!!   -----
!!             sedim.nc
!!   Output :
!!   ------
!!             *sedtrp.nc
!!   Workspace :
!!   ---------
!!   EXTERNAL :
!!   --------
!!   MODIFICATIONS:
!!   --------------
!!      original sedtrp : 2002  (Erik Buitenhuis)
!!----------------------------------------------------------------------
! 
! local variables
      INTEGER kt,kindic,ji,jj,jk,jn,nitb,lonid,latid,depid,status,numid
      INTEGER start(2),edge(2),nwritesed,nwritecal,dimsurf(3),blid
      INTEGER start2(3),edge2(3)
      PARAMETER (nwritesed=15,nwritecal=120)
      REAL sedmnt(3,dimlen)
      REAL step
      CHARACTER*40 clhstnam
#  include "domzgr_substitute.h90"
!!!---------------------------------------------------------------------
!!!  OPA8, LODYC (15/11/96)
!!!---------------------------------------------------------------------
      IF (kt.eq.nit000.and.kindic.eq.1) THEN
!
! 0.1 Initialisation
! -----------------
        status=nf_open('sedim.nc',NF_NOWRITE,nitb)
!          status=nf_inq_dimid(nitb,'NUMBER',inpid)
!          status=nf_inq_dimlen(nitb,inpid,dimlen)
!          ALLOCATE(sedimi(dimlen))
!          ALLOCATE(sedimj(dimlen))
!          ALLOCATE(sedimk(dimlen))
!          ALLOCATE(sedmnt(3,dimlen))
        status=nf_inq_varid(nitb,'SEDI',lonid)
        status=nf_inq_varid(nitb,'SEDJ',latid)
        status=nf_inq_varid(nitb,'SEDK',depid)
        status = nf_get_var_int(nitb,lonid,sedimi)
        status = nf_get_var_int(nitb,latid,sedimj)
        status = nf_get_var_int(nitb,depid,sedimk)
        status = nf_close(nitb)
        status=nf_open('climppin.nc',NF_NOWRITE,nitb)
        IF (status .NE. 0) STOP 'climppin.nc'
        status=nf_inq_varid(nitb,'ORCAI',lonid)
        IF (status .NE. 0) STOP 'ORCAI'
        status=nf_inq_varid(nitb,'ORCAJ',latid)
        IF (status .NE. 0) STOP 'ORCAJ'
        status=nf_inq_varid(nitb,'DATE',depid)
        IF (status .NE. 0) STOP 'DATE'
        status = nf_get_var_int(nitb,lonid,climi)
        IF (status .NE. 0) STOP 'climi'
        status = nf_get_var_int(nitb,latid,climj)
        status = nf_get_var_int(nitb,depid,climd)
        status = nf_close(nitb)
!
! 0.2 Define NETCDF files and fields at beginning of first time step
! ----------------------------------------------------------------------
!
        CALL dianam(clhstnam,nwrite,'sedtrp.nc')
        IF(lwp)WRITE(numout,*)" Name of NETCDF file ",clhstnam
        status = nf_create(clhstnam,NF_CLOBBER,nitc)
        status = nf_def_dim(nitc,'NUMBER',dimlen,dimid(1))
        status = nf_def_dim(nitc,'TIME',NF_UNLIMITED,dimid(2))
        status = nf_def_var(nitc,'TIME',nf_double,1,dimid(2),timeid)
!        status = nf_def_var(nitc,'TIME',nf_int,1,dimid(2),timeid)
        status = nf_def_var(nitc,'Export',nf_double,2,dimid,orgid)
        status = nf_def_var(nitc,'Expcal',nf_double,2,dimid,calid)
        status = nf_def_var(nitc,'Expsil',nf_double,2,dimid,silid)
        status = nf_enddef(nitc)
        sedoav = 0.
        sedcav = 0.
        sedsav = 0.
        climpar = 0.
        climsst = 0.
        climmld = 0.
        climchl = 0.
        climppt = 0.
        CALL dianam(clhstnam,nwrite,'calc.nc')
        IF(lwp)WRITE(numout,*)" Name of NETCDF file ",clhstnam
        status = nf_create(clhstnam,NF_CLOBBER,nitf)
        IF (status .NE. 0) WRITE(numout,*) 'calc cannot create',clhstnam
        status = nf_def_dim(nitf,'LONGITUDE',jpi,dimsurf(1))
        IF (status .NE. 0) WRITE(numout,*) 'calc cannot def LON'
        status = nf_def_dim(nitf,'LATITUDE',jpj,dimsurf(2))
        IF (status .NE. 0) WRITE(numout,*) 'calc cannot def LAT'
        status = nf_def_dim(nitf,'TIME',NF_UNLIMITED,dimsurf(3))
        IF (status .NE. 0) WRITE(numout,*) 'calc cannot def TIM'
        status = nf_def_var(nitf,'TIME',nf_double,1,dimsurf(3),timeid2)
        IF (status .NE. 0) WRITE(numout,*) 'calc cannot var TIM'
        status = nf_def_var(nitf,'CaCO3',nf_double,3,dimsurf,blmid)
        IF (status .NE. 0) WRITE(numout,*) 'calc cannot var CaCO3'
        status = nf_enddef(nitf)
        IF (status .NE. 0) WRITE(numout,*) 'calc cannot enddef'
        caco3av = 0.
        IF (mod( kt-1, nwritesed ).NE.0) THEN
          status=nf_open('sedtrp_old.nc',NF_NOWRITE,nitb)
          IF (status.NE.NF_NOERR) THEN
            WRITE(numout,*)" ERROR: sedtrp_old.nc needed, program continues"
          ELSE
            WRITE(numout,*)"sedtrp_old.nc used"
            status=nf_inq_varid(nitb,'Export',lonid)
            status=nf_inq_varid(nitb,'Expcal',latid)
            status=nf_inq_varid(nitb,'Expsil',depid)
            status=NF_GET_VAR_DOUBLE(nitb,lonid,sedoav)
            status=NF_GET_VAR_DOUBLE(nitb,latid,sedcav)
            status=NF_GET_VAR_DOUBLE(nitb,depid,sedsav)
            status=NF_CLOSE(nitb)
          ENDIF
        ENDIF
        IF(mod(mod(kt-1,5475),nwritecal).NE.0) THEN
          status=nf_open('calc_old.nc',NF_NOWRITE,nitb)
          IF (status.NE.NF_NOERR) THEN
            WRITE(numout,*)" ERROR: calc_old.nc needed, program continues"
          ELSE
            WRITE(numout,*)"calc_old.nc used"
            status=nf_inq_varid(nitb,'CaCO3',blid)
            status=NF_GET_VAR_DOUBLE(nitb,blid,caco3av)
            status=NF_CLOSE(nitb)
          ENDIF
        ENDIF
      ELSE
!
! 1. Every time step: calculate sediment trap averages
! ----------------------------------------------------------------------
!
        jn=1
        step = float(mod(kt-1,nwritesed)+1)
        DO jk = 1, jpk
          DO jj = 2, nlcj-1
            DO ji = 2, nlci-1
              IF ((jk.EQ.sedimk(jn)).AND.(jj.EQ.sedimj(jn)).AND. &
     &          (ji.EQ.sedimi(jn))) THEN
                sedmnt(1,jn)=(snkpoc(ji,jj,jk)+snkgoc(ji,jj,jk)) &
     &            *rjjss/rfact*1000
                sedmnt(2,jn)=snkcal(ji,jj,jk)*rjjss/rfact*1000
                sedmnt(3,jn)=snkdsi(ji,jj,jk)*rjjss/rfact*1000
                sedoav(jn)=sedoav(jn)+(sedmnt(1,jn)-sedoav(jn))/step
                sedcav(jn)=sedcav(jn)+(sedmnt(2,jn)-sedcav(jn))/step
                sedsav(jn)=sedsav(jn)+(sedmnt(3,jn)-sedsav(jn))/step
                jn=jn+1
              ENDIF
            END DO
          END DO
        END DO
!        WRITE(0,*)'sedmnt',(sedmnt(1,jn),jn=1,10)
!        WRITE(0,*)'sedoav',(sedoav(jn),jn=1,10)
!  2. Every day: write sediment trap data
!
        IF(mod( kt, nwritesed ).EQ.0) THEN
          IF(lwp) WRITE(numout,*) '**** sedtrp : write NetCDF array',kt
          start(1) = 1
          start(2) = (kt-nit000)/nwritesed+1
          edge(1) = dimlen
          edge(2) = 1
          status = NF_PUT_VARA_DOUBLE(nitc,timeid,start(2),1, &
     &      float(ndastp))
!          status = NF_PUT_VARA_INT(nitc,timeid,start(2),1,njulian)
          status=NF_PUT_VARA_DOUBLE(nitc,orgid,start,edge,sedoav)
          status=NF_PUT_VARA_DOUBLE(nitc,calid,start,edge,sedcav)
          status=NF_PUT_VARA_DOUBLE(nitc,silid,start,edge,sedsav)
          sedoav = 0.
          sedcav = 0.
          sedsav = 0.
        ELSE IF (kt.EQ.nitend) THEN
          status=nf_create('sedtrp.nc',NF_CLOBBER,nitb)
          status = nf_def_dim(nitb,'NUMBER',dimlen,numid)
          status = nf_def_var(nitb,'Export',nf_double,1,numid,lonid)
          status = nf_def_var(nitb,'Expcal',nf_double,1,numid,latid)
          status = nf_def_var(nitb,'Expsil',nf_double,1,numid,depid)
          status = nf_enddef(nitb)
          status=NF_PUT_VAR_DOUBLE(nitb,lonid,sedoav)
          status=NF_PUT_VAR_DOUBLE(nitb,latid,sedcav)
          status=NF_PUT_VAR_DOUBLE(nitb,depid,sedsav)
          status=NF_CLOSE(nitb)
        ENDIF
!
! On the correct day save climPP data
!
        DO jn = 1, clilen
!          IF (ndastp .EQ. climd(jn)) THEN
          IF ((ndastp-10000*(ndastp/10000)) .EQ.  &
     &      (climd(jn)-10000*(climd(jn)/10000))) THEN
            climpar(jn)=climpar(jn)+qsr(climi(jn),climj(jn))/15.
            climsst(jn)=climsst(jn)+tn(climi(jn),climj(jn),1)/15.
            DO ji=jpdia,jpdia+jppft-1
              climchl(jn)=climchl(jn)+trn(climi(jn),climj(jn),1, &
     &          ji+2*jppft)/15.
              DO jk = 1, jpk-1
                climppt(jn)=climppt(jn)+fse3t(ji,jj,jk)*1000.*rfactr*rjjss/15.* &
     &            pprphy(climi(jn),climj(jn),jk,ji,1)*12.
              END DO 
            END DO
            climmld(jn)=climmld(jn)+hmld(climi(jn),climj(jn))/15.
            IF (mod(kt,15) .EQ. 0) WRITE(numout,*) '**** clim', &
     &        climpar(jn),climsst(jn),climmld(jn) &
     &        ,climchl(jn),climppt(jn)
            IF ((mod(kt,15) .EQ. 0) .AND. (ndastp .EQ. climd(jn))) &
     &        WRITE(numout,*) 'climPPARR ',jn,ndastp,climpar(jn), &
     &          climsst(jn),climmld(jn),climchl(jn),climppt(jn)
          ENDIF
        END DO
        IF (kt.EQ.nitend) THEN
        CALL dianam(clhstnam,nwrite,'climppout.nc')
        IF(lwp)WRITE(numout,*)" Name of NETCDF file ",clhstnam
          status = nf_create(clhstnam,NF_CLOBBER,nitg)
          status = nf_def_dim(nitg,'NUMBER',clilen,climid)
          status = nf_def_var(nitg,'PAR',nf_double,1,climid,parid)
          status = nf_def_var(nitg,'SST',nf_double,1,climid,sstid)
          status = nf_def_var(nitg,'CHL',nf_double,1,climid,chlid)
          status = nf_def_var(nitg,'MLD',nf_double,1,climid,mldid)
          status = nf_def_var(nitg,'PPT',nf_double,1,climid,pptid)
          status = nf_enddef(nitg)
          status=NF_PUT_VAR_DOUBLE(nitg,parid,climpar)
          status=NF_PUT_VAR_DOUBLE(nitg,sstid,climsst)
          status=NF_PUT_VAR_DOUBLE(nitg,chlid,climchl)
          status=NF_PUT_VAR_DOUBLE(nitg,mldid,climmld)
          status=NF_PUT_VAR_DOUBLE(nitg,pptid,climppt)
          status=NF_CLOSE(nitg)
        ENDIF
! 3. Every time step: calculate detached CaCO3 averages
! ----------------------------------------------------------------------
!
        step = float(mod(mod(kt-1,5475),nwritecal)+1)
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
            caco3av(ji,jj)=caco3av(ji,jj)+ &
     &        (trn(ji,jj,1,jpcal)-caco3av(ji,jj))/step
          END DO
        END DO
!  2. Every 8 days: write CaCO3 data
        IF(mod( mod(kt,5475), nwritecal).EQ.0) THEN
          IF(lwp) WRITE(numout,*) '**** calc : write NetCDF array'
          start2(1) = 1
          start2(2) = 1
          start2(3) = (kt-nit000)/nwritecal+1
          IF (mod(kt,5475).EQ.0) start2(3)=start2(3)+1
          edge2(1) = jpi
          edge2(2) = jpj
          edge2(3) = 1
          status = NF_PUT_VARA_DOUBLE(nitf,timeid2,start2(3),1, &
     &      float(ndastp))
        IF(lwp) WRITE(numout,*) '**** calc', status,start2,edge2,ndastp
          status=NF_PUT_VARA_DOUBLE(nitf,blmid,start2,edge2,caco3av)
        IF(lwp) WRITE(numout,*) '**** calc', status,caco3av(75,75)
          caco3av = 0.
        ELSE IF (kt.EQ.nitend) THEN
          status=nf_create('calc.nc',NF_CLOBBER,nitb)
          status = nf_def_dim(nitb,'LONGITUDE',jpi,dimsurf(1))
          status = nf_def_dim(nitb,'LATITUDE',jpj,dimsurf(2))
          status = nf_def_var(nitb,'caco3',nf_double,2,dimsurf,blid)
          status = NF_PUT_VAR_DOUBLE(nitb,blid,caco3av)
          status = NF_CLOSE(nitb)
        ENDIF
      ENDIF
!
! 3. Closing file
! --------------------
      IF (kt.EQ.nitend.OR.kindic.LT.0) THEN
        status=NF_CLOSE(nitc)
        status=NF_CLOSE(nitf)
        IF(lwp) WRITE(numout,*) '**** calc close', status
      ENDIF
#    endif
      RETURN
      END
