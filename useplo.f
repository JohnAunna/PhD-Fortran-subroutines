c=======================================================================
      SUBROUTINE useplo(mat,mph,ifile)
c=======================================================================
c
c     calculation of viscous term indirectly
c-----------------------------------------------------------------------
c-----
      USE comm0
      USE comm1
      USE prec_mod, ONLY: prec     !...precision in FIRE
      IMPLICIT NONE
!      include 'com007.inc'
      INCLUDE 'SwiftIO_FortranFunctions.inc'
      INTEGER  mat   !< material index
      INTEGER, INTENT(IN) :: mph   !multiphase flag
      INTEGER, INTENT(IN) :: ifile !file number
      INTEGER          ::isel,i,ic,iic,ityp,nelem,nc,j,jc,no,PBs
       Real          ::sumnode,sumPB !cell selection vari
      REAL(prec), DIMENSION(ncell) ::conv_x,conv_y,conv_z !convection t
       REAL(prec), DIMENSION(ncell) ::p_x,p_y,p_z   !pressure gradients
       REAL(prec), DIMENSION(ncell) ::del2u,del2v,del2w,del2
      CHARACTER*256 string
c--------------------------------------------------------------------------local variables
c-----ARGUMENTS
      REAL(prec), ALLOCATABLE, DIMENSION(:) :: floatcell, floatbd
      REAL(prec), ALLOCATABLE, DIMENSION(:) :: nodevsc,nodevscbd
      REAL(prec), ALLOCATABLE, DIMENSION(:) :: PBvsc,PBvscbd
c
c-----------------------------------------------------------------------allocating the viscousity term in all cells
c-----
       IF(I_USEPLO == 1) THEN 
c---------------------------------------------------------------------------------------------------------------------------------------------------------viscous term calculation
         call gradfi(u,ub,g1,3,1,1)
	   call gradfi(u,ub,g2,3,2,1)
	   call gradfi(u,ub,g3,3,3,1)
c---------------------------------------------------------------------------
          Do ic = 1, ncell  
c------------------------------------------------calculating the convection term for all cells
      conv_x(ic)=den(ic)*vol(ic)*(u(1,ic)*g1(1,ic)+u(2,ic)*g1(2,ic)
     X +u(3,ic)*g1(3,ic))
	conv_y(ic)=den(ic)*vol(ic)*(u(1,ic)*g2(1,ic)+u(2,ic)*g2(2,ic)
     X  +u(3,ic)*g2(3,ic))
	conv_z(ic)=den(ic)*vol(ic)*(u(1,ic)*g3(1,ic)+u(2,ic)*g3(2,ic)
     X +u(3,ic)*g3(3,ic))
	    End do
c---------------------------------------------------------------------------
        call gradfi(p,pb,g1,1,1,1)
          Do ic = 1, ncell  
c------------------------------------------------calculating the pressure term for all cells
      p_x(ic)=g1(1,ic)
      p_y(ic)=g1(2,ic)
      p_z(ic)=g1(3,ic)
           End do
		 do ic = 1, ncell  
        del2u(ic)=conv_x(ic)+p_x(ic)+0.0  !Convection term+Pressure trem
        del2v(ic)=conv_y(ic)+p_y(ic)+0.0  !Convection term+Pressure trem
        del2w(ic)=conv_z(ic)+p_z(ic)+9.81*den(ic)*vol(ic)  !.+.-Gravity
c---------------------------------------------------------------------non-dimensionalized (magnitude of viscous term/ gravity term)
        del2(ic)=sqrt((del2u(ic))**2+(del2v(ic))**2+(del2w(ic))**2) 
     X              /(9.81*den(ic)) ! rhou.g
		 END DO
c---------------------------------------------------------------------getting the average of viscous term in node		 
		 string = 'just_node'
        isel = SWIFTIO_INDEX_OF_SELECTION(string)
        IF(isel.GT.0) THEN
          ityp= SWIFTIO_TYPE_OF_SELECTION(ISEL)
          nelem = SWIFTIO_NELEMENTS_OF_SELECTION(ISEL)
          IF(ityp.EQ.3) THEN
c---------------------------------------------------------------
        DO j = 1,nelem
         jc = swiftio_element_of_selection(isel,j)
c       
         sumnode= sumnode+del2(jc)    !sum of visc term in node
        End do
        no=nelem
          write(6,*) no,sumnode   !average of visc term
          End If
        End If
c---------------------------------------------------------------------getting the average of viscous term in top PB		 
		 string = 'just_PB'
        isel = SWIFTIO_INDEX_OF_SELECTION(string)
        IF(isel.GT.0) THEN
          ityp= SWIFTIO_TYPE_OF_SELECTION(ISEL)
          nelem = SWIFTIO_NELEMENTS_OF_SELECTION(ISEL)
          IF(ityp.EQ.3) THEN
c---------------------------------------------------------------
        DO i = 1,nelem
         iic = swiftio_element_of_selection(isel,i)
c       
         sumPB= sumPB+del2(iic)    !sum of visc term in PB
        End do
        PBs=nelem
          write(6,*) PBs,sumPB   !average of visc term in PB
          End If
        End If
c----------------------------------------------------------------------------End of the calculation of the average viscosity term in node
cc-------------------------------------------------------------------------------------------------------allocating the non-deimensionalized viscous term
       ALLOCATE (floatcell(ncell)); floatcell = zero
       ALLOCATE (floatbd(nbfac));  floatbd = zero
       do ic = 1, ncell
         floatcell(ic) = del2(ic) 
       end do
       call Write_User_Fl3 (mat,mph,ifile,'viscous term/rho.g',
     &                       '0',ncell,nbfac,floatcell,floatbd)

       DEALLOCATE (floatcell)
       DEALLOCATE (floatbd) 
c------------------------------------------------------------------------------start calculation of the average viscosity term in node
c----------------------------------------------------------------------------allocating the average viscosity term in node
       ALLOCATE (nodevsc(ncell));      nodevsc = zero
       ALLOCATE (nodevscbd(nbfac));  nodevscbd = zero
        do ic = 1, ncell
        nodevsc(ic)=sumnode/no !sumavg 
        End do       
       call Write_User_Fl3 (mat,mph,ifile,'node viscous average',
     &                       '0',ncell,nbfac,nodevsc,nodevscbd)

       DEALLOCATE (nodevsc)
       DEALLOCATE (nodevscbd) 
c----------------------------------------------------------------------------allocating the average viscosity term in node
       ALLOCATE (PBvsc(ncell));      PBvsc = zero
       ALLOCATE (PBvscbd(nbfac));  PBvscbd = zero
       do ic = 1, ncell
         PBvsc(ic) = sumPB/PBs 
       end do
       call Write_User_Fl3 (mat,mph,ifile,'PB viscous average',
     &                       '0',ncell,nbfac,PBvsc,PBvscbd)

       DEALLOCATE (PBvsc)
       DEALLOCATE (PBvscbd)       
      END IF
c-----
      RETURN
      END
