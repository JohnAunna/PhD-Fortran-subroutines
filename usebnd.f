c=======================================================================
      SUBROUTINE usebnd(k,mph)
c=======================================================================
c     author:Hamid Anazadeh           
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
      INTEGER         :: isel, ityp ,nelem,i,nc,ib,idir,nfc 
      REAL(prec), DIMENSION(ncell) :: p_face_area,esb_x,esb_y,esb_z
      REAL(prec), DIMENSION(ncell) :: Uc_pbx,Uc_pby,Uc_pbz,p_zc
      REAL(prec), DIMENSION(nbfac) :: dt2_I,dt2_II
      REAL(prec), DIMENSION(nbfac) :: ddelu_I,ddelv_I,ddelw_I
      REAL(prec), DIMENSION(nbfac) :: ddelu_II,ddelv_II,ddelw_II
          M=0.00000018               !constant For BSA
 !          M=0.000000035            !constant for SDS
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
c-------------------------------------------------------Distance from boundary face centre to cell centre
c       p_zc = sqrt(DB(1,ib)**2 + DB(2,ib)**2 + DB(3,ib)**2)
C---------------------------------------------------------- PERPENDICULAR distance from boundary face centre to cell centre-Delta n
          p_zc(ib) = sqrt((DB(1,ib)*esb_x(ib))**2 
     x         + (DB(2,ib)*esb_y(ib))**2 
     x         + (DB(3,ib)*esb_z(ib))**2)     
c------------------------------------------------------------distance between two boundary centres in tangential direction_I_delta t1**2
        dt2_I(ib)=(XB(1,ib)-XB(1,(ib-1)))**2 
     x         + (XB(2,ib)-XB(2,(ib-1)))**2 
     x         + (XB(3,ib)-XB(3,(ib-1)))**2
c------------------------------------------------------------distance between two boundary centres in tangential direction_II_deltat2**2
        dt2_II(ib)=(XB(1,ib)-XB(1,(ib-120)))**2 
     x         + (XB(2,ib)-XB(2,(ib-120)))**2 
     x         + (XB(3,ib)-XB(3,(ib-120)))**2     
c-----------------------------------------------------second gradient of surface velocity_consecutive faces
      ddelu_I(ib)= ub(1,(ib+1))-2*ub(1,ib)+ub(1,(ib-1))/dt2_I(ib)
      ddelv_I(ib)= ub(2,(ib+1))-2*ub(2,ib)+ub(2,(ib-1))/dt2_I(ib)
      ddelw_I(ib)= ub(3,(ib+1))-2*ub(3,ib)+ub(3,(ib-1))/dt2_I(ib)
c-----------------------------------------------------second gradient of surface velocity_upper and lower faces
      ddelu_II(ib)= ub(1,ib+120)-2*ub(1,ib)+ub(1,ib-120)/dt2_II(ib)
      ddelv_II(ib)= ub(2,ib+120)-2*ub(2,ib)+ub(2,ib-120)/dt2_II(ib)
      ddelw_II(ib)= ub(3,ib+120)-2*ub(3,ib)+ub(3,ib-120)/dt2_II(ib)
c------------------------------------------------------------------------ Surface velocity at boundary interface
!      ub(1,ib)=Uc_pbx(ip)-M*(ddelu_I(ib)+ddelu_II(ib))*p_zc(ib) 
!      ub(2,ib)=Uc_pby(ip)-M*(ddelv_I(ib)+ddelv_II(ib))*p_zc(ib)
      ub(3,ib)=Uc_pbz(ip)-M*(ddelw_I(ib)+ddelw_II(ib))*p_zc(ib)
c
!	write(6,*)'ddelu,v,w_II=second gradients'
!	write(6,*)ddelu_II(ib)
!	write(6,*)ddelv_II(ib)
!	write(6,*)ddelw_II(ib)
!      write(6,*)'                   '
!      write(6,*)'boundary velocity',ib,ub(1,ib),ub(2,ib),ub(3,ib)
!      write(6,*)'                   '
!      write(6,*)'cell velocity',ip,u(1,ip),u(2,ip),u(3,ip)
!      write(6,*)'                   '
!      write(6,*)'ddelu,v,w_I=second gradients'
!	write(6,*)ddelu_I(ib)
!	write(6,*)ddelv_I(ib)
!	write(6,*)ddelw_I(ib)
!	write(6,*)'               '
!      write(6,*) ib,ip
          ENDDO
!      write(6,*)
      write(6,*)'cell velocity'
      write(6,*)u(1,72117),u(2,72117),u(3,72117)
      write(6,*)'                                       '
!      write(6,*)'sb(1,2,3)',sb(1,16517),sb(2,16517),sb(3,16517)
!      write(6,*)'                                               '
!      write(6,*)'db(1,2,3)', db(1,16517),db(2,16517),db(3,16517)
!      write(6,*)'                                               '
!      write(6,*)'p_face_area',p_face_area(16517)
!       write(6,*)'                                               '
!      write(6,*)'p_zc',p_zc(16517)
!      write(6,*)'                                               '
!      write(6,*)'esb_x,y,z',esb_x(16517),esb_y(16517),esb_z(16517)
!      write(6,*)'                                       ' 
!      Write(6,*) 'dt2_I,dt2_II',dt2_I(16517),dt2_II(16517)
!      Write(6,*) sqrt(dt2_I(16517)),sqrt(dt2_II(16517))
!       write(6,*)'                                       ' 
	write(6,*)'ddelu,v,w_I=second gradients'
	write(6,*)ddelu_I(16517)
	write(6,*)ddelv_I(16517)
	write(6,*)ddelw_I(16517)
	write(6,*)'                                       '
	write(6,*)'ddelu,v,w_II=second gradients'
	write(6,*)ddelu_II(16517)
	write(6,*)ddelv_II(16517)
	write(6,*)ddelw_II(16517)
	write(6,*)'                                       '
      write(6,*)'delu*p_zc',(ddelu_I(16517)+ddelu_II(16517))*p_zc(16517)
      write(6,*)'delv*p_zc',(ddelv_I(16517)+ddelv_II(16517))*p_zc(16517)
      write(6,*)'delw*p_zc',(ddelw_I(16517)+ddelw_II(16517))*p_zc(16517)
      write(6,*)'                                       '
      write(6,*)'M*du*pz',M*(ddelu_I(16517)+ddelu_II(16517))*p_zc(16517)
      write(6,*)'M*du*pz',M*(ddelv_I(16517)+ddelv_II(16517))*p_zc(16517)
      write(6,*)'M*du*pz',M*(ddelw_I(16517)+ddelw_II(16517))*p_zc(16517)
      write(6,*)'                                       '
      write(6,*)'Uc_pb=cell component parallel to boundary face'
      write(6,*)Uc_pbx(72117)
      write(6,*)Uc_pby(72117)
      write(6,*)Uc_pbz(72117)
      write(6,*)'-----                                     '
      write(6,*)'ub=Boundary velocity'
      write(6,*)ub(1,16517),ub(2,16517),ub(3,16517)
!      write(6,*)ub(3,5920)
!      write(6,*)ub(3,5921)
!      write(6,*)ub(3,5922)
!      write(6,*)ub(3,6000)
!      write(6,*)'   '
!      write(6,*)ub(3,6040)
!      write(6,*)ub(3,6041)
!      write(6,*)ub(3,6042)
!      write(6,*)'   '
!      write(6,*)ub(3,6160)
!      write(6,*)ub(3,6161)
!      write(6,*)ub(3,6162)
!      write(6,*)'   '
!      write(6,*)ub(3,6280)
!      write(6,*)ub(3,6281)
!      write(6,*)ub(3,6282)
!      write(6,*)'   '
!      write(6,*)ub(3,9283)
!      write(6,*)ub(3,10000)
!      write(6,*)ub(3,11779)
!      write(6,*)ub(3,16517)
!       write(6,*)ub(3,11801)
      write(6,*)'------------------------------------------------------'
      write(6,*)'                                       '
              ENDIF
c-----
        END IF
        END IF     
      RETURN
      END
      
  