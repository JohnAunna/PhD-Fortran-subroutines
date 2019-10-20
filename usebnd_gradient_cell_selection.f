c=======================================================================
      SUBROUTINE usebnd(k,mph)
c=======================================================================
c     author:John Aunna           
c     date:30 October 2015             
c     project:Simulation of Foam          
c     project manager:Jamal Naser
c     description: user subroutines for modification of Boundary Condition in interface of air and PB
c------------------------------------------------------------------------------------
      USE comm0
      USE comm1
      USE prec_mod, ONLY: prec     !...precision in FIRE
      IMPLICIT NONE
      INCLUDE 'SwiftIO_FortranFunctions.inc'
      INTEGER  mat   !< material index
      INTEGER k !< value which shows the boundary region, see .log file
      INTEGER mph !< phase index
c--------------------------------------------------------------------------local variables
      CHARACTER*256 string
      INTEGER            :: ip       !< cell index 
      REAL(prec)         :: M  ! < constant dependent on the surfactant
      INTEGER         :: strin, isel, ityp ,nelem,i,nc,ib1,idir,ip1,nfc 
      REAL(prec), DIMENSION(ncell) :: p_face_area,esb_x,esb_y,esb_z
      REAL(prec), DIMENSION(ncell) :: Uc_pbx,Uc_pby,Uc_pbz,p_zc
      REAL(prec), DIMENSION(ncell) :: ddelu,ddelv,ddelw
      REAL(prec), DIMENSION(ncell) :: ddelu_t,ddelv_t,ddelw_t
      REAL(prec), DIMENSION(:)     :: delu(3,ncell),delv(3,ncell)
      REAL(prec), DIMENSION(:)    :: delw(3,ncell),B(3,ncell)
      REAL(prec), DIMENSION(:)   :: delub(3,nbfac),delvb(3,nbfac)
      REAL(prec), DIMENSION(:)   :: delwb(3,nbfac),B_b(3,nbfac)
!          M=0.018               !constant For BSA
           M=0.000035            !constant for SDS
c----------------------------------------------------------------------------------
      IF(I_USEBND == 1) THEN
c-----------------------------------------------------------------------face selection with the name of 'fac'
            string = 'fac'
        isel = SWIFTIO_INDEX_OF_SELECTION(string)
        IF(isel.GT.0) THEN
          ityp  = SWIFTIO_TYPE_OF_SELECTION(ISEL)
          nelem = SWIFTIO_NELEMENTS_OF_SELECTION(ISEL)
          IF(ityp.EQ.2) THEN
c---------------------------------------------------------------
        DO i = 1,nelem,2
         nc = swiftio_element_of_selection(isel,i)
         idir = swiftio_element_of_selection(isel,i+1)
         ib1 = swiftio_Index_of_Sel_Bnd_Face(nc,idir)
         nfc = lb(ib1)
              ip1=nfc
          B(1,ip1)=ub(1,ib1)
          B(2,ip1)=ub(2,ib1)
          B(3,ip1)=ub(3,ib1)
          B_b(1,ib1)=ub(1,ib1)
          B_b(2,ib1)=ub(1,ib1)
          B_b(3,ib1)=ub(1,ib1)       
          End Do   
c---------------------------------------------------------------calling first gradint of velocity 
         call gradfi(B,B_b,g1,3,1,1)
	   call gradfi(B,B_b,g2,3,2,1)
	   call gradfi(B,B_b,g3,3,3,1)
!c---------------------------------------------------------------------------
!         call gradfi(u,ub,g1,3,1,1)
!	   call gradfi(u,ub,g2,3,2,1)
!	   call gradfi(u,ub,g3,3,3,1)
c---------------------------------------------------------------------------
         DO i = 1,nelem,2
         nc = swiftio_element_of_selection(isel,i)
         idir = swiftio_element_of_selection(isel,i+1)
         ib1 = swiftio_Index_of_Sel_Bnd_Face(nc,idir)
         nfc = lb(ib1)
              ip1=nfc
c------------------------------------------------putting velocity first gradient components into defined variables
         delu(1,ip1)=g1(1,ip1)
		 delu(2,ip1)=g1(2,ip1)
		 delu(3,ip1)=g1(3,ip1)
		 delv(1,ip1)=g2(1,ip1)
		 delv(2,ip1)=g2(2,ip1)
		 delv(3,ip1)=g2(3,ip1)
		 delw(1,ip1)=g3(1,ip1)
		 delw(2,ip1)=g3(2,ip1)
		 delw(3,ip1)=g3(3,ip1)
c                                   boundary values
         delub(1,ib1)=g1(1,ib1)
		 delub(2,ib1)=g1(2,ib1)
		 delub(3,ib1)=g1(3,ip1)
		 delvb(1,ib1)=g2(1,ip1)
		 delvb(2,ib1)=g2(2,ib1)
		 delvb(3,ib1)=g2(3,ib1)
		 delwb(1,ib1)=g3(1,ib1)
		 delwb(2,ib1)=g3(2,ib1)
		 delwb(3,ib1)=g3(3,ib1)
		 END DO
c-----------------------------------------------calling second gradient of velocity
	   call gradfi(delu,delub,g1,3,1,1)
	   call gradfi(delu,delub,g2,3,2,1)
	   call gradfi(delu,delub,g3,3,3,1)
c----------------------------------------------putting these three gradient in one vector(ddel_u)-[d2u/dx2, du2/dy2, d2u/dz2]
	   DO i = 1,nelem,2
         nc = swiftio_element_of_selection(isel,i)
         idir = swiftio_element_of_selection(isel,i+1)
         ib1 = swiftio_Index_of_Sel_Bnd_Face(nc,idir)
            nfc = lb(ib1)
              ip1=nfc
C
	   ddelu(ip1)=g1(1,ip1)+g2(2,ip1)+g3(3,ip1)
		 END DO
c-----------------
         call gradfi(delv,delvb,g1,3,1,1)
	   call gradfi(delv,delvb,g2,3,2,1)
	   call gradfi(delv,delvb,g3,3,3,1)
c------------------------------putting these three gradient in one vector(ddel_v)-[d2v/dx2, dv2/dy2, d2v/dz2]
        DO i = 1,nelem,2
         nc = swiftio_element_of_selection(isel,i)
         idir = swiftio_element_of_selection(isel,i+1)
         ib1 = swiftio_Index_of_Sel_Bnd_Face(nc,idir)
          nfc = lb(ib1)
              ip1=nfc
C
	   ddelv(ip1)=g1(1,ip1)+g2(2,ip1)+g3(3,ip1)
		  END DO
c-----------------
	   call gradfi(delw,delwb,g1,3,1,1)
	   call gradfi(delw,delwb,g2,3,2,1)
	   call gradfi(delw,delwb,g3,3,3,1)
c------------------------------putting these three gradient in one vector(ddel_W)-[d2w/dx2, dw2/dy2, d2w/dz2]
        DO i = 1,nelem,2
         nc = swiftio_element_of_selection(isel,i)
         idir = swiftio_element_of_selection(isel,i+1)
         ib1 = swiftio_Index_of_Sel_Bnd_Face(nc,idir)
             nfc = lb(ib1)
              ip1=nfc
C
	   ddelw(ip1)=g1(1,ip1)+g2(2,ip1)+g3(3,ip1)
		 END DO
c-------------------------------sum of the second gradient of u,v,w[(d2u/dx2+d2v/dx2+d2w/dx2),(du2/dy2+dv2/dy2+dw2/dy2),(d2u/dz2+d2v/dz2+d2w/dz2)]
            DO i = 1,nelem,2
         nc = swiftio_element_of_selection(isel,i)
         idir = swiftio_element_of_selection(isel,i+1)
         ib1 = swiftio_Index_of_Sel_Bnd_Face(nc,idir)
         nfc = lb(ib1)
              ip1=nfc
!c-----------------------------------------------------Tangential components of second gradients
       ddelu_t(ip1) = ddelu(ip1)-ddelu(ip1)*esb_x(ib1)
       ddelv_t(ip1) = ddelv(ip1)-ddelv(ip1)*esb_y(ib1)
       ddelw_t(ip1) = ddelw(ip1)-ddelw(ip1)*esb_z(ib1)
c
      p_face_area(ib1)=sqrt(SB(1,ib1)**2 +SB(2,ib1)**2+SB(3,ib1)**2)
c---------------------------------------------------------------------- Boundary face unit normal vector components
        esb_x(ib1) = sb(1,ib1)/p_face_area(ib1)
        esb_y(ib1) = sb(2,ib1)/p_face_area(ib1)
        esb_z(ib1) = sb(3,ib1)/p_face_area(ib1)
c----- -------------------------------------------------Distance from boundary face centre to cell centre
c       p_zc = sqrt(DB(1,ib1)**2 + DB(2,ib1)**2 + DB(3,ib1)**2)
C---------------------------------------------------------- PERPENDICULAR distance from boundary face centre to cell centre
          p_zc(ib1) = sqrt((DB(1,ib1)*esb_x(ib1))**2 
     x         + (DB(2,ib1)*esb_y(ib1))**2 
     x         + (DB(3,ib1)*esb_z(ib1))**2)
c-----------------------------------------------------Boundary face connected cell velocity component parallel to boundary face
       Uc_pbx(ip1) = u(1,ip1)-u(1,ip1)*esb_x(ib1)
       Uc_pby(ip1) = u(2,ip1)-u(2,ip1)*esb_y(ib1)
       Uc_pbz(ip1) = u(3,ip1)-u(3,ip1)*esb_z(ib1)
c------------------------------------------------------------------------ Surface velocity at boundary interface
      ub(1,ib1)=Uc_pbx(ip1)-M*ddelu_t(ip1)*p_zc(ib1) 
      ub(2,ib1)=Uc_pby(ip1)-M*ddelv_t(ip1)*p_zc(ib1)
      ub(3,ib1)=Uc_pbz(ip1)-M*ddelw_t(ip1)*p_zc(ib1)
c
          ENDDO
      write(6,*)ib1,ip1
      write(6,*)'cell velocity'
      write(6,*)u(1,117720)
      write(6,*)u(2,117720)
      write(6,*)u(3,117720)
      write(6,*)'                                       '
      write(6,*)'sb(1,2,3)',sb(1,11800),sb(2,11800),sb(3,11800)
      write(6,*)'                                               '
      write(6,*)'db(1,2,3)', db(1,11800),db(2,11800),db(3,11800)
      write(6,*)'                                               '
      write(6,*)'p_face_area',p_face_area(11800)
       write(6,*)'                                               '
      write(6,*)'  p_zc',p_zc(11800)
      write(6,*)'                                               '
      write(6,*)'esb_x,y,z',esb_x(11800),esb_y(11800),esb_z(11800)
      write(6,*)'                                       ' 
      Write(6,*) 'B',B(1,33655),B(2,33655),B(3,33655)  
       write(6,*)'                                       ' 
       Write(6,*) 'B_b',B_b(1,33655),B_b(2,33655),B_b(3,33655)  
      write(6,*)'                                       '     
      write(6,*)'delu,v,w=values of gradients'
      write(6,*)delu(1,117720),delu(2,117720),delu(3,117720)
      write(6,*)delv(1,117720),delv(2,117720),delv(3,117720)
	write(6,*)delw(1,117720),delw(2,117720),delw(3,117720)
	write(6,*)'                                       '
	write(6,*)'delb_u,v,w=boundary values of gradients'
	write(6,*)delub(1,11800),delub(2,11800),delub(3,11800)
      write(6,*)delvb(1,11800),delvb(2,11800),delvb(3,11800)
	write(6,*)delwb(1,11800),delwb(2,11800),delwb(3,11800)
	write(6,*)'                                       '
	write(6,*)'ddelu,v,w=second gradients'
	write(6,*)ddelu(117720)
	write(6,*)ddelv(117720)
	write(6,*)ddelw(117720)
	write(6,*)'                                       '
	write(6,*)'ddelu_t,v,w=second gradients'
	write(6,*)ddelu_t(117720)
	write(6,*)ddelv_t(117720)
	write(6,*)ddelw_t(117720)
	write(6,*)'                                       '
      write(6,*)'ddelu*p_zc', ddelu_t(117720)*p_zc(11800)
      write(6,*)'ddelv*p_zc',ddelv_t(117720)*p_zc(11800)
      write(6,*)'ddelw*p_zc', ddelw_t(117720)*p_zc(11800)
      write(6,*)'                                       '
      write(6,*)'M*ddelu_t*p_zc', M*ddelu_t(117720)*p_zc(11800)
      write(6,*)'M*ddelv_t*p_zc',M*ddelv_t(117720)*p_zc(11800)
      write(6,*)'M*ddelw_t*p_zc', M*ddelw_t(117720)*p_zc(11800)
      write(6,*)'                                       '
      write(6,*)'Uc_pb=cell component parallel to boundary face'
      write(6,*)Uc_pbx(117720)
      write(6,*)Uc_pby(117720)
      write(6,*)Uc_pbz(117720)
      write(6,*)'-----                                     '
      write(6,*)'ub=Boundary velocity'
      write(6,*)ub(1,11800)
      write(6,*)ub(2,11800)
      write(6,*)ub(3,11800)
      write(6,*)'------------------------------------------------------'
      write(6,*)'                                       '
              ENDIF
c-----
        END IF
        END IF     
      RETURN
      END
      
      
