module dynsoil_lithium
implicit none

contains

subroutine lithium(reg_P_diss, temp, runoff, h_soil , Li_Friv, Li_Fsp, Li_driv, list_cont_pxl, ncontpxl)

  use dynsoil_empirical_laws, only: SP_fraction, Li_fractionation
  use dynsoil_physical_parameters, only: nlon, nlat, nlitho, rho_UCC, Li_UCC, Li_enrichment, dLi_UCC
  integer, parameter:: npxl = nlon*nlat

  double precision, dimension(nlitho,npxl), intent(in)::  reg_P_diss, h_soil
  double precision, dimension(npxl),        intent(in)::  temp, runoff
  double precision, dimension(nlitho,npxl), intent(out):: Li_Friv, Li_Fsp, Li_driv
  integer, dimension(npxl), intent(in):: list_cont_pxl
  integer, intent(in):: ncontpxl
  double precision:: Li_Fdiss, Deltaland
  integer:: i, j, j0


  do j0 = 1,ncontpxl
    j = list_cont_pxl(j0)

    if (runoff(j)>0) then !

      do i = 1,nlitho-1 ! Skip last lithology class => carbonates

        !===============================================================================================================!
        !==============     Secondary phases fluxes and lithium partitioning and isotope fractionation      ============!
        !===============================================================================================================!
        !                                                                                                               !
        !         compute riverine/sec phases partioning with water residence time. Lithium equations come from         !
        !                               Vigier and Godderis (2015) Clim Past, 11, 635-645                               !
        !                                                                                                               !
        !===============================================================================================================!

        ! WARNING: runoff conversion: cm/y -> m/y

        ! Lithium fluxes:
        Li_Fdiss      = Li_UCC * rho_UCC * reg_P_diss(i,j)
        Li_Fsp(i,j)   = Li_Fdiss * Li_enrichment * SP_fraction(temp(j), runoff(j)/100, h_soil(i,j))
        Li_Friv(i,j)  = Li_Fdiss - Li_Fsp(i,j)
        ! WARNING: SP_fraction must be lower or equal to 1/Li_enrichment to avoid negative Li Friv

        ! isotopic factionation:
        Deltaland = Li_fractionation(temp(j))

        ! Lithium isotopic budget:
        if (Li_Fdiss>0) then
          Li_driv(i,j) = dLi_UCC  +  DeltaLand*(Li_Fsp(i,j)/Li_Fdiss)
        else
          Li_driv(i,j) = 0
        end if

      end do
    end if
  end do


  end subroutine

end module
