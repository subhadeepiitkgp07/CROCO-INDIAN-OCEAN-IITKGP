#include "cppdefs.h"

#if defined MUSTANG

      module plug_MUSTANG_CROCO

      USE module_MUSTANG
      USE initMUSTANG, ONLY : MUSTANG_init
      USE sed_MUSTANG, ONLY : MUSTANG_update
      USE sed_MUSTANG, ONLY : MUSTANG_deposition
# ifdef MORPHODYN
      USE sed_MUSTANG, ONLY : MUSTANG_morpho
# endif

      IMPLICIT NONE

      PRIVATE

      PUBLIC   mustang_update_main
      PUBLIC   mustang_deposition_main
      PUBLIC   mustang_init_main
# ifdef MORPHODYN
      PUBLIC   mustang_morpho_main
# endif

CONTAINS
!
!-----------------------------------------------------------------------
!
      subroutine mustang_update_main (tile)

      REAL    :: TEMPREF_LIN
      REAL    :: SALREF_LIN
      INTEGER :: tile
# include "ocean2d.h"
# include "compute_tile_bounds.h"

      TEMPREF_LIN = 10.0 
      SALREF_LIN  = 35.0
      CALL MUSTANG_update (Istr, Iend, Jstr, Jend,  & 
                   t, zob, zeta,                    &
# if defined key_MUSTANG_lateralerosion || defined key_MUSTANG_bedload
                   ubar, vbar,                      &
# endif
                   SALREF_LIN, TEMPREF_LIN, dt)
      end subroutine
!
!-----------------------------------------------------------------------
!
      subroutine mustang_deposition_main (tile)

      integer :: tile
# include "ocean2d.h"
# include "compute_tile_bounds.h"
      CALL MUSTANG_deposition (Istr, Iend, Jstr, Jend, zeta, t)
      end subroutine
!
!-----------------------------------------------------------------------
!
      subroutine mustang_init_main (tile)

      REAL :: h0fond
      integer :: tile
# include "ocean2d.h"
# include "compute_tile_bounds.h"

# ifdef WET_DRY
            h0fond = D_wetdry
# else
            h0fond = 0.
# endif

      CALL MUSTANG_init (Istr, Iend, Jstr, Jend,  &
                    zeta,                         &
# if defined MORPHODYN
                    dh,                           &
# endif
                    h0fond, zob, t)
      end subroutine
!
!-----------------------------------------------------------------------
# ifdef MORPHODYN
      subroutine mustang_morpho_main (tile)

      integer :: tile
# include "ocean2d.h"
# include "compute_tile_bounds.h"

      CALL MUSTANG_morpho (Istr, Iend, Jstr, Jend, dh)

      end subroutine
# endif
!-----------------------------------------------------------------------
!
      end module plug_MUSTANG_CROCO

#else

      module plug_MUSTANG_CROCO_empty
      end module plug_MUSTANG_CROCO_empty
      
#endif /* MUSTANG */
