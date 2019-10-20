c=======================================================================
      SUBROUTINE useini(k,mph)
c=======================================================================
c
c     USEINI IS A SPECIAL PURPOSE ROUTINE FOR INITIALISAION 
c     BASED ON USER CODING
C=======================================================================
c     author:Hamid Anazad         
c     date:7 November 2016             
c     project:Simulation of Foam          
c     project manager:Jamal Naser
c     description: user subroutines for defining the indexes of neighbour boundary faces for each boundary faces of air-liquid interface (intfc)
c------------------------------------------------------------------------------------
c-----
      USE comm0
      USE comm1
      USE prec_mod, ONLY: prec     !...precision in FIRE
!      IMPLICIT NONE
      include 'com007.inc'
      INCLUDE 'SwiftIO_FortranFunctions.inc'
      INTEGER  mat   !< material index
      INTEGER k !< value which shows the boundary region, see .log file
      INTEGER mph !< phase index
c--------------------------------------------------------------------------local variables 
      CHARACTER*256 string
      INTEGER         :: ip       !< cell index 
      INTEGER         :: isel, ityp ,nelem,i,nc,ib,idir,nfc,nbface
      INTEGER         :: j, jb, n
      REAL(prec), DIMENSION(nbfac)     :: dis !,ib_1,ib_2,ib_3
      REAL(prec), DIMENSION(3,nbfac)     :: nghbr
      n=1
c----------------------------------------------------------------------------------
      IF(I_USEINI == 1) THEN
c-----------------------------------------------------------------------face selection with the name of 'fac'
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
C---------------------------------------------------------- selecting 3 neighbour boundary face 
            DO j = 1,nelem,2
                   nc = swiftio_element_of_selection(isel,j)
                  idir = swiftio_element_of_selection(isel,j+1)
                   jb = swiftio_Index_of_Sel_Bnd_Face(nc,idir)
          dis(jb) = sqrt((XB(1,jb)-XB(1,ib))**2 
     x         +(XB(2,jb)-XB(2,ib))**2 
     x         +(XB(3,jb)-XB(3,ib))**2) 
!         write(6,*) ip,ib,jb
!         write(6,*) XB(1,ib),XB(2,ib),XB(3,ib)
!         write(6,*) XB(1,jb),XB(2,jb),XB(3,jb)
!        write(6,*) dis(jb)
           IF (dis(jb)<0.000038.AND.dis(jb)>0.00000000000001) THEN 
        nghbr(n,ib)=jb
          n=n+1
!      write(6,*) 'click'
           END IF
         END DO
            IF (nghbr(1,ib)<10000000.AND.nghbr(1,ib)>0) THEN
             ib_1(ib)=nghbr(1,ib)      !define this in the com007
            ELSE 
            ib_1(ib)=0
            END IF
           IF (nghbr(2,ib)<10000000.AND.nghbr(2,ib)>0)THEN 
             ib_2(ib)=nghbr(2,ib)      !define this in the com007
             ELSE
              ib_2(ib)=0
              END IF
           IF (nghbr(3,ib)<10000000.AND.nghbr(3,ib)>0.0) THEN
             ib_3(ib)=nghbr(3,ib)      !define this in the com007
             ELSE
              ib_3(ib)=0
              END IF
             n=1        !resetting the n number
        write(6,*)'ngbr------------------------------------------------'
!          write(6,*) ib,ip
!          write(6,*) ib_1(ib),ib_2(ib),ib_3(ib)
        END DO
        ENDIF
c-----
        END IF
c-----
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
C---------------------------------------------------------- selecting 3 neighbour boundary face 
            DO j = 1,nelem,2
                   nc = swiftio_element_of_selection(isel,j)
                  idir = swiftio_element_of_selection(isel,j+1)
                   jb = swiftio_Index_of_Sel_Bnd_Face(nc,idir)
          dis(jb) = sqrt((XB(1,jb)-XB(1,ib))**2 
     x         +(XB(2,jb)-XB(2,ib))**2 
     x         +(XB(3,jb)-XB(3,ib))**2) 
           IF (dis(jb)<0.0000045.AND.dis(jb)>0.00000000000001) THEN 
        nghbr(n,ib)=jb
          n=n+1
!      write(6,*) 'click'
           END IF
         END DO
            IF (nghbr(1,ib)<10000000.AND.nghbr(1,ib)>0) THEN
             ib_1(ib)=nghbr(1,ib)      !define this in the com007
            ELSE 
            ib_1(ib)=0
            END IF
           IF (nghbr(2,ib)<10000000.AND.nghbr(2,ib)>0)THEN 
             ib_2(ib)=nghbr(2,ib)      !define this in the com007
             ELSE
              ib_2(ib)=0
              END IF
           IF (nghbr(3,ib)<10000000.AND.nghbr(3,ib)>0.0) THEN
             ib_3(ib)=nghbr(3,ib)      !define this in the com007
             ELSE
              ib_3(ib)=0
              END IF
             n=1        !resetting the n number
        write(6,*)'PB neighbours---------------------------------------'
!          write(6,*) ib,ip
!          write(6,*) ib_1(ib),ib_2(ib),ib_3(ib)
        END DO
        ENDIF
c-----
        END IF
c-----------------------------------------------------------------------face selection with the name of 'intfc_small for lower face of geometry'
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
C---------------------------------------------------------- selecting 3 neighbour boundary face 
            DO j = 1,nelem,2
                   nc = swiftio_element_of_selection(isel,j)
                  idir = swiftio_element_of_selection(isel,j+1)
                   jb = swiftio_Index_of_Sel_Bnd_Face(nc,idir)
          dis(jb) = sqrt((XB(1,jb)-XB(1,ib))**2 
     x         +(XB(2,jb)-XB(2,ib))**2 
     x         +(XB(3,jb)-XB(3,ib))**2) 
!         write(6,*) ip,ib,jb
!         write(6,*) XB(1,ib),XB(2,ib),XB(3,ib)
!         write(6,*) XB(1,jb),XB(2,jb),XB(3,jb)
!        write(6,*) dis(jb)
           IF (dis(jb)<0.0000024.AND.dis(jb)>0.00000000000001) THEN 
        nghbr(n,ib)=jb
          n=n+1
!      write(6,*) 'click'
           END IF
         END DO
            IF (nghbr(1,ib)<10000000.AND.nghbr(1,ib)>0) THEN
             ib_1(ib)=nghbr(1,ib)      !define this in the com007
            ELSE 
            ib_1(ib)=0
            END IF
           IF (nghbr(2,ib)<10000000.AND.nghbr(2,ib)>0)THEN 
             ib_2(ib)=nghbr(2,ib)      !define this in the com007
             ELSE
              ib_2(ib)=0
              END IF
           IF (nghbr(3,ib)<10000000.AND.nghbr(3,ib)>0.0) THEN
             ib_3(ib)=nghbr(3,ib)      !define this in the com007
             ELSE
              ib_3(ib)=0
              END IF
             n=1        !resetting the n number
        write(6,*)'ngbr_small lower--------------------------------'
!          write(6,*) ib,ip
!          write(6,*) ib_1(ib),ib_2(ib),ib_3(ib)
        END DO
        ENDIF
c-----
        END IF
        END IF    
c-----
      RETURN
      END