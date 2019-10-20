c------Convection term dropping
c      author:John Aunna    date:27 October 2015      
c      SUBROUTINE usebod(iconv,mph)       
c
      SUBROUTINE usebod(mat,mph)                
      USE prec_mod, ONLY : prec
      USE comm0
      USE comm1
      IMPLICIT NONE
c
c      INCLUDE 'comdp.inc'
c      INCLUDE 'SwiftIO_FortranFunctions.inc'
c      
      INTEGER, INTENT(IN) :: mat   !< material index
      INTEGER, INTENT(IN) :: mph   !< phase index (multiphase)
c-----
c-----local variables
c-----
      CHARACTER*256 string
      INTEGER            :: ip       !< cell index 
      REAL(prec)         :: mass ! < specific fluid mass in cell(kg/m3) 
      REAL(prec), DIMENSION(ncell) :: Sum_conu,Sum_conv,Sum_conw
c-------------------------------------------------------------
c      mat=1
c
	IF(i_usebod > 0) THEN
      call gradfi(u,ub,g1,3,1,mat)
	call gradfi(u,ub,g2,3,2,mat)
	call gradfi(u,ub,g3,3,3,mat)
c	
	write(*,*) 'g1(1,68000)',g1(1,68000) 
c	
	Do ip=1,ncell
      Sum_conu(ip)=den(ip)*vol(ip)*(u(1,ip)*g1(1,ip)+u(2,ip)*g1(2,ip)
     X +u(3,ip)*g1(3,ip))
	Sum_conv(ip)=den(ip)*vol(ip)*(u(1,ip)*g2(1,ip)+u(2,ip)*g2(2,ip)
     X  +u(3,ip)*g2(3,ip))
	Sum_conw(ip)=den(ip)*vol(ip)*(u(1,ip)*g3(1,ip)+u(2,ip)*g3(2,ip)
     X +u(3,ip)*g3(3,ip))
	End do
c	
c
c      
	Do ip=1,ncell
	su1(ip)=su1(ip)+Sum_conu(ip)
	su2(ip)=su2(ip)+Sum_conv(ip)
	su3(ip)=su3(ip)+Sum_conw(ip)
	End do
c
      write(*,*) 'Sumw',Sum_conw(68000)
c      
	END IF
	 RETURN
      END SUBROUTINE usebod
