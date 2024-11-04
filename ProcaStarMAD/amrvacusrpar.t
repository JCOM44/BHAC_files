!##############################################################################
! include amrvacusrpar - grmhd_torus

! This file should contain the number of PROBLEM dependent equation parameters,
! the index names for them with values neqpar+1..neqpar+nspecialpar,
! and a string giving the names for the file header. For example:
!
! INTEGER,PARAMETER:: mass_=neqpar+1, nspecialpar=1
! CHARACTER*4,PARAMETER:: specialparname='mass'
!
! By default there are no special parameters
 INTEGER,PARAMETER:: rin_=neqpar+1,l_=neqpar+2,rhomin_=neqpar+3,pmin_=neqpar+4,kappa_=neqpar+5,beta_=neqpar+6

INTEGER,PARAMETER:: nspecialpar=6
CHARACTER*1,PARAMETER:: specialparname=' '
CHARACTER*20,PARAMETER:: typeuser='grmhd_torus'
! end include amrvacusrpar - grmhd_torus
!##############################################################################
