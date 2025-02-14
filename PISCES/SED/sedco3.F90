#include "cppdefs.h"

MODULE sedco3
   !!======================================================================
   !!              ***  MODULE  sedco3  ***
   !!    Sediment : carbonate in sediment pore water
   !!=====================================================================
   !! * Modules used
#if defined key_sediment
   USE sed     ! sediment global variable
   USE sedchem
   USE lib_mpp         ! distribued memory computing library


   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC sed_co3     

   !!* Substitution
#  include "ocean2pisces.h90"

   !!----------------------------------------------------------------------
   !!   OPA 9.0   !   LODYC-IPSL   (2003)
   !!----------------------------------------------------------------------

   !! $Id: sedco3.F90 10222 2018-10-25 09:42:23Z aumont $
CONTAINS


   SUBROUTINE sed_co3( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_co3  ***
      !!
      !! ** Purpose :  carbonate ion and proton concentration 
      !!               in sediment pore water
      !!
      !! ** Methode :  - solving nonlinear equation for [H+] with given alkalinity
      !!               and total co2 
      !!               - one dimensional newton-raphson algorithm for [H+])
      !!
      !!   History :
      !!        !  98-08 (E. Maier-Reimer, Christoph Heinze )  Original code
      !!        !  04-10 (N. Emprin, M. Gehlen ) coupled with PISCES
      !!        !  06-04 (C. Ethe)  Re-organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(in)  :: kt   ! time step
      !
      !---Local variables
      INTEGER  :: ji, jk           ! dummy loop indices

      REAL(wp), DIMENSION(jpoce,jpksed) :: zhinit, zhi
     !!----------------------------------------------------------------------

      IF( ln_timing )  CALL timing_start('sed_co3')

      DO jk = 1, jpksed
         DO ji = 1, jpoce
            zhinit(ji,jk)   = hipor(ji,jk) / densSW(ji)
         END DO
      END DO

      !     -------------------------------------------
      !     COMPUTE [CO3--] and [H+] CONCENTRATIONS
      !     -------------------------------------------

      CALL solve_at_general_sed(zhinit, zhi)

      DO jk = 1, jpksed
         DO ji = 1, jpoce
            co3por(ji,jk) = pwcp(ji,jk,jwdic) * ak1s(ji) * ak2s(ji) / (zhi(ji,jk)**2   &
            &               + ak1s(ji) * zhi(ji,jk) + ak1s(ji) * ak2s(ji) + rtrn )
            hipor(ji,jk)  = zhi(ji,jk) * densSW(ji)
         END DO
      END DO

     IF( ln_timing )  CALL timing_stop('sed_co3')

   END SUBROUTINE sed_co3

#endif

END MODULE sedco3
