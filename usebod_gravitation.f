c=======================================================================
!> @brief Special purpose routine for adding the body force
!> @details All user-specified expressions must be written in standard
!! Fortran 90 or in the Fortran version available on your machine. <br>
!! To activate the statements, please remove the 'c' from the
!! first column. <br>
!! This routine will be executed after each iteration. <br>
!! See the list of variables for the use of this routine.
!=======================================================================
      SUBROUTINE usebod(mat,mph)
!=======================================================================
!
!.....contact fire@avl.com
!
!     TO ACTIVATE THE EXAMPLE CODE SET VARIABLE "IS_EXAMPLE_CODE" TO
!     "FALSE".
!
!     Note: "I_USEBOD" is the integer value specified in the GUI at
!           "User-functions / Activation".
!
!-----------------------------------------------------------------------
!
      USE comm0
      USE comm1
      USE comm2
      USE prec_mod, ONLY: prec     !...precision in FIRE
!
      IMPLICIT NONE
!
!
!     arguments
      INTEGER, INTENT(IN) :: mat   !< material index
      INTEGER, INTENT(IN) :: mph   !< phase index (multiphase)
!
!     local variables
      INTEGER            :: ip                 !< cell index
      REAL(prec)         :: mass  !< specific fluid mass in cell (kg/m3)
      REAL(prec), DIMENSION(:)  :: bforvec(3)!< body force vector (m/s2)
!
!
!-----------------------------------------------------------------------
      IF(I_USEBOD == 1) THEN
!
!
!-----------------------------------------------------------------------
!     Example: Apply fix body force (for multiphase and single phase)
!-----------------------------------------------------------------------
!      ELSE IF(I_USEBOD == 2) THEN
!!
!       IF (.NOT. IS_EXAMPLE_CODE) THEN
!
!       set body-force vector
        bforvec(1) = 0.0
        bforvec(2) = 0.0
        bforvec(3) = -5.19
!
!       loop over all cells in the current domain (material)
        DO ip = nsp(mat), nep(mat)
!
!          calculate specific fluid mass in cell
           mass = den(ip) * vf(ip)
!
!          add source
           su1(ip) = su1(ip) + mass * bforvec(1) * vol(ip)
           su2(ip) = su2(ip) + mass * bforvec(2) * vol(ip)
           su3(ip) = su3(ip) + mass * bforvec(3) * vol(ip)
!
        END DO
!
       END IF
      write(6,*) 'g is working'
      RETURN
      END SUBROUTINE usebod
!