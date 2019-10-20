c=======================================================================
      SUBROUTINE usebnd(k,mph)
c=======================================================================
c     author:Hamid Anazad         
c     date:30 October 2015             
c     project:Simulation of Foam          
c     project manager:Jamal Naser
c     description: user subroutines for modification of Boundary Condition in interface of air and PB in Film and PB plus Marangoni area
c------------------------------------------------------------------------------------
      USE comm0
      USE comm1
      USE prec_mod, ONLY: prec     !...precision in FIRE
!      IMPLICIT NONE
      include 'com007.inc'
      INCLUDE 'SwiftIO_FortranFunctions.inc'
      INTEGER  mat   !< material index
      INTEGER k !< value which shows the boundary region, see .log file
      INTEGER mph !< phase index
c-----------------------------------------------------------------------Boundary variables 
      CHARACTER*256 string
      INTEGER            :: ip       !< cell index 
      REAL(prec)         :: Bo1,Bo2,N,E,g_eq, g_ex        ! < constants
      REAL(prec)         :: Mar_x,Mar_y,Mar_z
      INTEGER         :: isel, ityp ,nelem,i,nc,ib,idir,nfc,ic
      REAL(prec), DIMENSION(ncell) :: p_face_area,esb_x,esb_y,esb_z
      REAL(prec), DIMENSION(ncell) :: Uc_pbx,Uc_pby,Uc_pbz,p_zc
      REAL(prec), DIMENSION(nbfac) :: dt2_I,B1,B2,B3,force
      REAL(prec), DIMENSION(nbfac) :: ddelu_I,ddelv_I,ddelw_I
c--------------------------------------------------------------------variables for viscous term
      REAL(prec), DIMENSION(:)    :: delu(3,ncell),delv(3,ncell)
      REAL(prec), DIMENSION(:)    :: delw(3,ncell)
      REAL(prec), DIMENSION(:)    :: delub(3,ncell),delvb(3,ncell)
      REAL(prec), DIMENSION(:)    :: delwb(3,ncell)
      REAL(prec), DIMENSION(ncell)    ::del2u,del2v,del2w
      REAL(prec), DIMENSION(ncell) :: ddudxdx,ddudxdy,ddudxdz
	REAL(prec), DIMENSION(ncell) :: ddudydx,ddudydy,ddudydz   
	REAL(prec), DIMENSION(ncell) :: ddudzdx,ddudzdy,ddudzdz   
	REAL(prec), DIMENSION(ncell) :: ddvdxdx,ddvdxdy,ddvdxdz   
	REAL(prec), DIMENSION(ncell) :: ddvdydx,ddvdydy,ddvdydz   
	REAL(prec), DIMENSION(ncell) :: ddvdzdx,ddvdzdy,ddvdzdz
	REAL(prec), DIMENSION(ncell) :: ddwdxdx,ddwdxdy,ddwdxdz   
	REAL(prec), DIMENSION(ncell) :: ddwdydx,ddwdydy,ddwdydz   
	REAL(prec), DIMENSION(ncell) :: ddwdzdx,ddwdzdy,ddwdzdz 
c-----------------------------------------------------------------------Boussinesque number
	     Bo1=0.00000035/1.        !constant for SDS 0.5
           Bo2=0.00000035/1.        !constant for 
c---------------------------------------------------------------------Marangoni variables
!           N=0.02            !constant for adapting the delta_n
!           E=0.01           !Gibbs elasticity
!           g_ex=0.00000001  !instead of 10e-06, 110 has been used
!           g_eq=0.000001  !equilibrium concentration gama
c----------------------------------------------------------------------------------
      IF(I_USEBND == 1) THEN
c---------------------------------------------------------------------------------------------------------------------------------------------------------viscous term calculation
         call gradfi(u,ub,g1,3,1,1)
	   call gradfi(u,ub,g2,3,2,1)
	   call gradfi(u,ub,g3,3,3,1)
c---------------------------------------------------------------------------
          Do ic = 1, ncell  !defining the second gradient for all cells
c------------------------------------------------putting velocity first gradient components into defined variables
         delu(1,ic)=g1(1,ic)    !delu/delx
		 delu(2,ic)=g1(2,ic)    !delu/dely
		 delu(3,ic)=g1(3,ic)    !delu/delz
		 delv(1,ic)=g2(1,ic)    !delv/delx
		 delv(2,ic)=g2(2,ic)    !delv/dely
		 delv(3,ic)=g2(3,ic)    !delv/delz
		 delw(1,ic)=g3(1,ic)    !delw/delx
		 delw(2,ic)=g3(2,ic)    !delw/dely
		 delw(3,ic)=g3(3,ic)    !delw/delz
c                                   boundary values
         delub(1,ic)=g1(1,ic)
		 delub(2,ic)=g1(2,ic)
		 delub(3,ic)=g1(3,ic)
		 delvb(1,ic)=g2(1,ic)
		 delvb(2,ic)=g2(2,ic)
		 delvb(3,ic)=g2(3,ic)
		 delwb(1,ic)=g3(1,ic)
		 delwb(2,ic)=g3(2,ic)
		 delwb(3,ic)=g3(3,ic)
		 END DO
c-----------------------------------------------calling second gradient of velocity
	   call gradfi(delu,delub,g1,3,1,1)
	   call gradfi(delu,delub,g2,3,2,1)
	   call gradfi(delu,delub,g3,3,3,1)
c----------------------------------------------
	    do ic = 1, ncell  
C---------------------------------------------------------------------------------------
	   ddudxdx(ic)=g1(1,ic)   !del(delu/delx)/delx
		 ddudxdy(ic)=g1(2,ic)   !del(delu/delx)/dely
		 ddudxdz(ic)=g1(3,ic)   !del(delu/delx)/delz
		 ddudydx(ic)=g2(1,ic)   !del(delu/dely)/delx
		 ddudydy(ic)=g2(2,ic)   !del(delu/dely)/dely
		 ddudydz(ic)=g2(3,ic)   !del(delu/dely)/delz
		 ddudzdx(ic)=g3(1,ic)   !del(delu/delz)/delx
		 ddudzdy(ic)=g3(2,ic)   !del(delu/delz)/dely
		 ddudzdz(ic)=g3(3,ic)   !del(delu/delz)/delz
		 END DO
c-----------------
         call gradfi(delv,delvb,g1,3,1,1)
	   call gradfi(delv,delvb,g2,3,2,1)
	   call gradfi(delv,delvb,g3,3,3,1)
c------------------------------------------------------------------------------------
         do ic = 1, ncell  
c----------------------------------------------------------
         ddvdxdx(ic)=g1(1,ic)   !del(delv/delx)/delx
		 ddvdxdy(ic)=g1(2,ic)   !del(delv/delx)/dely
		 ddvdxdz(ic)=g1(3,ic)   !del(delv/delx)/delz
		 ddvdydx(ic)=g2(1,ic)   !del(delv/dely)/delx
		 ddvdydy(ic)=g2(2,ic)   !del(delv/dely)/dely
		 ddvdydz(ic)=g2(3,ic)   !del(delv/dely)/delz
		 ddvdzdx(ic)=g3(1,ic)   !del(delv/delz)/delx
		 ddvdzdy(ic)=g3(2,ic)   !del(delv/delz)/dely
		 ddvdzdz(ic)=g3(3,ic)   !del(delv/delz)/delz
		  END DO
c-----------------
	   call gradfi(delw,delwb,g1,3,1,1)
	   call gradfi(delw,delwb,g2,3,2,1)
	   call gradfi(delw,delwb,g3,3,3,1)
c-------------------------------------------------------------------------
        do ic = 1, ncell  
         ddwdxdx(ic)=g1(1,ic)   !del(delw/delx)/delx
		 ddwdxdy(ic)=g1(2,ic)   !del(delw/delx)/dely
		 ddwdxdz(ic)=g1(3,ic)   !del(delw/delx)/delz
		 ddwdydx(ic)=g2(1,ic)   !del(delw/dely)/delx
		 ddwdydy(ic)=g2(2,ic)   !del(delw/dely)/dely
		 ddwdydz(ic)=g2(3,ic)   !del(delw/dely)/delz
		 ddwdzdx(ic)=g3(1,ic)   !del(delw/delz)/delx
		 ddwdzdy(ic)=g3(2,ic)   !del(delw/delz)/dely
		 ddwdzdz(ic)=g3(3,ic)   !del(delw/delz)/delz
		 END DO
		 do ic = 1, ncell  
        del2u(ic)=2.*ddudxdx(ic)
     X            +ddudydy(ic)+ddvdydx(ic)
     X            +ddudzdz(ic)+ddwdxdz(ic)
        del2v(ic)=ddudydx(ic)+ddvdxdx(ic)
     X           +2.*ddvdydy(ic)
     X           +ddvdzdz(ic)+ddwdydz(ic)
        del2w(ic)=ddudzdx(ic)+ddwdxdx(ic)
     X           +ddvdzdy(ic)+ddwdydy(ic)
     X           +2.*ddwdzdz(ic)
        del2(ic)=sqrt((del2u(ic))**2+(del2v(ic))**2+(del2w(ic))**2)
		 END DO
c---------------------------------------------------------------------------------------------------------------------------------------------------------boundary definition
            string = 'intfc'
        isel = SWIFTIO_INDEX_OF_SELECTION(string)
        IF(isel.GT.0) THEN
          ityp  = SWIFTIO_TYPE_OF_SELECTION(ISEL)
          nelem = SWIFTIO_NELEMENTS_OF_SELECTION(ISEL)
          IF(ityp.EQ.2) THEN
c---------------------------------------------------------------
        DO i = 1,nelem,2
         nc = swiftio_element_of_selection(isel,i)
         idir = swiftio_element_of_selection(isel,i+1)
         ib = swiftio_Index_of_Sel_Bnd_Face(nc,idir)
         nfc = lb(ib)
              ip=nfc
c
      p_face_area(ib)=sqrt(SB(1,ib)**2 +SB(2,ib)**2+SB(3,ib)**2)
c---------------------------------------------------------------------- Boundary face unit normal vector components
        esb_x(ib) = sb(1,ib)/p_face_area(ib)
        esb_y(ib) = sb(2,ib)/p_face_area(ib)
        esb_z(ib) = sb(3,ib)/p_face_area(ib)
c-----------------------------------------------------Boundary face connected cell velocity component parallel to boundary face
       Uc_pbx(ip) = u(1,ip)-u(1,ip)*esb_x(ib)
       Uc_pby(ip) = u(2,ip)-u(2,ip)*esb_y(ib)
       Uc_pbz(ip) = u(3,ip)-u(3,ip)*esb_z(ib)
C---------------------------------------------------------- PERPENDICULAR distance from boundary face centre to cell centre-Delta n
          p_zc(ib) = sqrt((DB(1,ib)*esb_x(ib))**2 
     x         + (DB(2,ib)*esb_y(ib))**2 
     x         + (DB(3,ib)*esb_z(ib))**2)  
c------------------------------------------------------------distance between two boundary centres in tangential direction_{(i-j)+(j-k)+(i-k)}/3.0
        dt2_I(ib)=(XB(1,ib)-XB(1,ib_1(ib)))**2 
     x         + (XB(2,ib)-XB(2,ib_1(ib)))**2 
     x         + (XB(3,ib)-XB(3,ib_1(ib)))**2
c--------------------------------------------------------------------tri-mesh new approach for calculation of second gradient of boundary suface velocity
        ddelu_I(ib)=sqrt((ub(1,ib_1(ib))-2*ub(1,ib)+ub(1,ib_2(ib)))**2
     x         +(ub(1,ib_2(ib))-2*ub(1,ib)+ub(1,ib_3(ib)))**2
     x         +(ub(1,ib_1(ib))-2*ub(1,ib)+ub(1,ib_3(ib)))**2)/dt2_I(ib)
     
         ddelv_I(ib)=sqrt((ub(2,ib_1(ib))-2*ub(2,ib)+ub(2,ib_2(ib)))**2
     x        +(ub(2,ib_2(ib))-2*ub(2,ib)+ub(2,ib_3(ib)))**2
     x        +(ub(2,ib_1(ib))-2*ub(2,ib)+ub(2,ib_3(ib)))**2)/dt2_I(ib)
     
          ddelw_I(ib)=sqrt((ub(1,ib_1(ib))-2*ub(1,ib)+ub(1,ib_2(ib)))**2
     x         +(ub(3,ib_2(ib))-2*ub(3,ib)+ub(3,ib_3(ib)))**2
     x         +(ub(3,ib_1(ib))-2*ub(3,ib)+ub(3,ib_3(ib)))**2)/dt2_I(ib)
!c------------------------------------------------------------------------ Surface velocity at boundary interface- Ito ensure deduction of magnitude of velocity, IF condition is used
!      IF(Uc_pbx(ip)*(M*ddelu_I(ib)*p_zc(ib))<0.0) THEN
!      ub(1,ib)=Uc_pbx(ip)+M*ddelu_I(ib)*p_zc(ib) 
!      ELSE 
!      ub(1,ib)=Uc_pbx(ip)-M*ddelu_I(ib)*p_zc(ib)
!      END IF
      IF(Uc_pby(ip)>0.0) THEN
      ub(2,ib)=Uc_pby(ip)-Bo1*ddelv_I(ib)*p_zc(ib)
      ELSE
      Uc_pby(ip)=1e-30
      END IF
      ub(3,ib)=Uc_pbz(ip)+Bo1*ddelw_I(ib)*p_zc(ib)
      ENDDO
	write(6,*)'ddelu,v,w_I=second gradients for intfc'
	write(6,*)ddelv_I(15909)
	write(6,*)ddelw_I(15909)
	write(6,*)'                                       '
	write(6,*)ddelv_I(14135)
	write(6,*)ddelw_I(14135)
	write(6,*)'                                       '	
	write(6,*)ddelv_I(18465)
	write(6,*)ddelw_I(18465)
	write(6,*)'                                       '
	write(6,*)ddelv_I(13780)
	write(6,*)ddelw_I(13780)
	write(6,*)'                                       '
	write(6,*)ddelv_I(17620)
	write(6,*)ddelw_I(17620)
	write(6,*)'                                       '
	write(6,*)'Debit term'
	write(6,*) Bo1*ddelv_I(15909)*p_zc(15909)
	write(6,*) Bo1*ddelw_I(15909)*p_zc(15909)
	write(6,*)'                                       '
	write(6,*) Bo1*ddelv_I(14135)*p_zc(14135)
	write(6,*) Bo1*ddelw_I(14135)*p_zc(14135)
	write(6,*)'                                       '
	write(6,*) Bo1*ddelv_I(18465)*p_zc(18465)
	write(6,*) Bo1*ddelw_I(18465)*p_zc(18465)
	write(6,*)'                                       '
	write(6,*) Bo1*ddelv_I(13780)*p_zc(13780)
	write(6,*) Bo1*ddelw_I(13780)*p_zc(13780)
	write(6,*)'                                       '
	write(6,*) Bo1*ddelv_I(17620)*p_zc(17620)
	write(6,*) Bo1*ddelw_I(17620)*p_zc(17620)
	write(6,*)' Ub(1,2,3)                                      '
      write(6,*)ub(2,15909)
      write(6,*)ub(3,15909)
      	write(6,*)'                                       '
      	 write(6,*)ub(2,14135)
      write(6,*)ub(3,14135)
      	write(6,*)'                                       '
      	 write(6,*)ub(2,18465)
      write(6,*)ub(3,18465)
      	write(6,*)'                                       '
      	 write(6,*)ub(2,13780)
      write(6,*)ub(3,13780)
      	write(6,*)'                                       '
      	 write(6,*)ub(2,17620)
      write(6,*)ub(3,17620)
      	write(6,*)'                                       '
              ENDIF
c-----
        END IF
c---------------------------------------------------------------------------------------       
         string = 'intfc_small'
        isel = SWIFTIO_INDEX_OF_SELECTION(string)
        IF(isel.GT.0) THEN
          ityp  = SWIFTIO_TYPE_OF_SELECTION(ISEL)
          nelem = SWIFTIO_NELEMENTS_OF_SELECTION(ISEL)
          IF(ityp.EQ.2) THEN
c---------------------------------------------------------------
        DO i = 1,nelem,2
         nc = swiftio_element_of_selection(isel,i)
         idir = swiftio_element_of_selection(isel,i+1)
         ib = swiftio_Index_of_Sel_Bnd_Face(nc,idir)
         nfc = lb(ib)
              ip=nfc
c
      p_face_area(ib)=sqrt(SB(1,ib)**2 +SB(2,ib)**2+SB(3,ib)**2)
c---------------------------------------------------------------------- Boundary face unit normal vector components
        esb_x(ib) = sb(1,ib)/p_face_area(ib)
        esb_y(ib) = sb(2,ib)/p_face_area(ib)
        esb_z(ib) = sb(3,ib)/p_face_area(ib)
c-----------------------------------------------------Boundary face connected cell velocity component parallel to boundary face
       Uc_pbx(ip) = u(1,ip)-u(1,ip)*esb_x(ib)
       Uc_pby(ip) = u(2,ip)-u(2,ip)*esb_y(ib)
       Uc_pbz(ip) = u(3,ip)-u(3,ip)*esb_z(ib)
C---------------------------------------------------------- PERPENDICULAR distance from boundary face centre to cell centre-Delta n
          p_zc(ib) = sqrt((DB(1,ib)*esb_x(ib))**2 
     x         + (DB(2,ib)*esb_y(ib))**2 
     x         + (DB(3,ib)*esb_z(ib))**2)  
c------------------------------------------------------------distance between two boundary centres in tangential direction_{(i-j)+(j-k)+(i-k)}/3.0
        dt2_I(ib)=(XB(1,ib)-XB(1,ib_1(ib)))**2 
     x         + (XB(2,ib)-XB(2,ib_1(ib)))**2 
     x         + (XB(3,ib)-XB(3,ib_1(ib)))**2
c--------------------------------------------------------------------tri-mesh new approach for calculation of second gradient of boundary suface velocity
        ddelu_I(ib)=sqrt((ub(1,ib_1(ib))-2*ub(1,ib)+ub(1,ib_2(ib)))**2
     x         +(ub(1,ib_2(ib))-2*ub(1,ib)+ub(1,ib_3(ib)))**2
     x         +(ub(1,ib_1(ib))-2*ub(1,ib)+ub(1,ib_3(ib)))**2)/dt2_I(ib)
     
         ddelv_I(ib)=sqrt((ub(2,ib_1(ib))-2*ub(2,ib)+ub(2,ib_2(ib)))**2
     x        +(ub(2,ib_2(ib))-2*ub(2,ib)+ub(2,ib_3(ib)))**2
     x        +(ub(2,ib_1(ib))-2*ub(2,ib)+ub(2,ib_3(ib)))**2)/dt2_I(ib)
     
          ddelw_I(ib)=sqrt((ub(1,ib_1(ib))-2*ub(1,ib)+ub(1,ib_2(ib)))**2
     x         +(ub(3,ib_2(ib))-2*ub(3,ib)+ub(3,ib_3(ib)))**2
     x         +(ub(3,ib_1(ib))-2*ub(3,ib)+ub(3,ib_3(ib)))**2)/dt2_I(ib)
c---------------------------------------------------------------------Surface velocity at boundary interface-to ensure deduction of magnitude of velocity, IF condition is used
      ub(2,ib)=Uc_pby(ip)-Bo2*ddelv_I(ib)*p_zc(ib)
      ENDDO
	write(6,*)'ddelv=second gradients for intfc_small'
	write(6,*)ddelv_I(17428)
	write(6,*)ddelv_I(18115)
	write(6,*)'Debit term'
	write(6,*) Bo2*ddelv_I(17428)*p_zc(17428)
	write(6,*) Bo2*ddelv_I(18115)*p_zc(18115)
	write(6,*)'Boundary velocity'
      write(6,*)ub(2,17428)
      write(6,*)ub(2,18115)
              ENDIF
c-----------------------------------------------------------For checking the PB
        END IF
          string = 'intfc_PB'
        isel = SWIFTIO_INDEX_OF_SELECTION(string)
        IF(isel.GT.0) THEN
          ityp  = SWIFTIO_TYPE_OF_SELECTION(ISEL)
          nelem = SWIFTIO_NELEMENTS_OF_SELECTION(ISEL)
          IF(ityp.EQ.2) THEN
c---------------------------------------------------------------
        DO i = 1,nelem,2
         nc = swiftio_element_of_selection(isel,i)
         idir = swiftio_element_of_selection(isel,i+1)
         ib = swiftio_Index_of_Sel_Bnd_Face(nc,idir)
         nfc = lb(ib)
              ip=nfc
c
      p_face_area(ib)=sqrt(SB(1,ib)**2 +SB(2,ib)**2+SB(3,ib)**2)
c---------------------------------------------------------------------- Boundary face unit normal vector components
        esb_x(ib) = sb(1,ib)/p_face_area(ib)
        esb_y(ib) = sb(2,ib)/p_face_area(ib)
        esb_z(ib) = sb(3,ib)/p_face_area(ib)
c-----------------------------------------------------Boundary face connected cell velocity component parallel to boundary face
       Uc_pbx(ip) = u(1,ip)-u(1,ip)*esb_x(ib)
       Uc_pby(ip) = u(2,ip)-u(2,ip)*esb_y(ib)
       Uc_pbz(ip) = u(3,ip)-u(3,ip)*esb_z(ib)
C---------------------------------------------------------- PERPENDICULAR distance from boundary face centre to cell centre-Delta n
          p_zc(ib) = sqrt((DB(1,ib)*esb_x(ib))**2 
     x         + (DB(2,ib)*esb_y(ib))**2 
     x         + (DB(3,ib)*esb_z(ib))**2)  
c------------------------------------------------------------distance between two boundary centres in tangential direction_{(i-j)+(j-k)+(i-k)}/3.0
        dt2_I(ib)=(XB(1,ib)-XB(1,ib_1(ib)))**2 
     x         + (XB(2,ib)-XB(2,ib_1(ib)))**2 
     x         + (XB(3,ib)-XB(3,ib_1(ib)))**2
c--------------------------------------------------------------------tri-mesh new approach for calculation of second gradient of boundary suface velocity
        ddelu_I(ib)=sqrt((ub(1,ib_1(ib))-2*ub(1,ib)+ub(1,ib_2(ib)))**2
     x         +(ub(1,ib_2(ib))-2*ub(1,ib)+ub(1,ib_3(ib)))**2
     x         +(ub(1,ib_1(ib))-2*ub(1,ib)+ub(1,ib_3(ib)))**2)/dt2_I(ib)
     
         ddelv_I(ib)=sqrt((ub(2,ib_1(ib))-2*ub(2,ib)+ub(2,ib_2(ib)))**2
     x        +(ub(2,ib_2(ib))-2*ub(2,ib)+ub(2,ib_3(ib)))**2
     x        +(ub(2,ib_1(ib))-2*ub(2,ib)+ub(2,ib_3(ib)))**2)/dt2_I(ib)
     
          ddelw_I(ib)=sqrt((ub(1,ib_1(ib))-2*ub(1,ib)+ub(1,ib_2(ib)))**2
     x         +(ub(3,ib_2(ib))-2*ub(3,ib)+ub(3,ib_3(ib)))**2
     x         +(ub(3,ib_1(ib))-2*ub(3,ib)+ub(3,ib_3(ib)))**2)/dt2_I(ib)
c------------------------------------------------------------------------ Surface velocity at boundary interface- Ito ensure deduction of magnitude of velocity, IF condition is used
      ub(3,ib)=Uc_pbz(ip)+Bo1*ddelw_I(ib)*p_zc(ib)
       ENDDO
!      write(6,*)'PB PB PB PB'
!      write(6,*)ib_1 (4272),ib_2(4272),ib_3(4272)
!      write(6,*) p_zc(4272)
!      write(6,*) dt2_I(4272)
!	write(6,*)'ddelu,v,w_I=second gradients'
!	write(6,*)ddelw_I(4272)
!	write(6,*)'                                       '
!	write(6,*)'Debit term'
!	write(6,*) Bo1*ddelw_I(4272)*p_zc(4272)
!	write(6,*)' Ub(1,2,3)                                      '
!	write(6,*)ub(1,4272)
!      write(6,*)ub(2,4272)
!      write(6,*)ub(3,4272)
              ENDIF
c-----
        END IF
!      string = 'Mar_Intfc'
!        isel = SWIFTIO_INDEX_OF_SELECTION(string)
!        IF(isel.GT.0) THEN
!          ityp  = SWIFTIO_TYPE_OF_SELECTION(ISEL)
!          nelem = SWIFTIO_NELEMENTS_OF_SELECTION(ISEL)
!          IF(ityp.EQ.2) THEN
!c---------------------------------------------------------------calling first gradint of gamma 
!         DO i = 1,nelem,2
!         nc = swiftio_element_of_selection(isel,i)
!         idir = swiftio_element_of_selection(isel,i+1)
!         ib = swiftio_Index_of_Sel_Bnd_Face(nc,idir)
!         nfc = lb(ib)
!              ip=nfc
!          phib(ib)=phi(ip)    
!         ENDDO    
!         call gradfi(phi,phib,g1,1,1,1)
!c---------------------------------------------------------------loop over the boundary-adjacent cells
!        DO i = 1,nelem,2
!         nc = swiftio_element_of_selection(isel,i)
!         idir = swiftio_element_of_selection(isel,i+1)
!         ib = swiftio_Index_of_Sel_Bnd_Face(nc,idir)
!         nfc = lb(ib)
!              ip=nfc
!c
!      p_face_area(ib)=sqrt(SB(1,ib)**2 +SB(2,ib)**2+SB(3,ib)**2)
!c---------------------------------------------------------------------- Boundary face unit normal vector components
!        esb_x(ib) = sb(1,ib)/p_face_area(ib)
!        esb_y(ib) = sb(2,ib)/p_face_area(ib)
!        esb_z(ib) = sb(3,ib)/p_face_area(ib)
!c-----------------------------------------------------Boundary face connected cell velocity component parallel to boundary face
!       Uc_pbx(ip) = u(1,ip)-u(1,ip)*esb_x(ib)
!       Uc_pby(ip) = u(2,ip)-u(2,ip)*esb_y(ib)
!       Uc_pbz(ip) = u(3,ip)-u(3,ip)*esb_z(ib)
!C---------------------------------------------------------- PERPENDICULAR distance from boundary face centre to cell centre-Delta n
!          p_zc(ib) = sqrt((DB(1,ib)*esb_x(ib))**2 
!     x         + (DB(2,ib)*esb_y(ib))**2 
!     x         + (DB(3,ib)*esb_z(ib))**2)  
!c----------------------------------------------------------------Marangoni force direction and magnitude
!      Mar_x=0.0
!      Mar_y=-0.8480   !-sin(58)
!      Mar_z=0.5299    !cos(58)
!      force(ib)=abs(g1(3,ip))*p_zc(ib)*E*N*g_ex/vim(ip)/g_eq
!c--------------------------------------------------------------------Surface velocity at boundary interface---ub-uc=*deltan*(dgama/dz)*E/(viscosity*Equilibrium surface concentration)
!      ub(1,ib)=Uc_pbx(ip)+Mar_x*force(ib)
!      ub(2,ib)=Uc_pby(ip)+Mar_y*force(ib)
!      ub(3,ib)=Uc_pbz(ip)+Mar_z*force(ib)
!c
!          ENDDO
!!       write(6,*)'Mar_Intfc is working'
!!      write(6,*)'cell velocity'
!!      write(6,*)u(1,23845)
!!      write(6,*)u(2,23845)
!!      write(6,*)u(3,23845)
!!      write(6,*) 'gama gradient'
!!      write(6,*) g1(1,23845)     
!!      write(6,*) g1(2,23845)
!!      write(6,*) abs(g1(3,23845))
!!      write(6,*) force(15079)
!!      write(6,*) Mar_x,Mar_y,Mar_z
!!	write(6,*)'                                       '
!!	write(6,*)'Boundary velocity'
!!      write(6,*)ub(1,15079)
!!      write(6,*)ub(2,15079)
!!      write(6,*)ub(3,15079)
!!      write(6,*)'------------------------------------------------------'
!        END IF 
!        ENDIF
        END IF    
      RETURN
      END
      
  