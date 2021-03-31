module dynsoil_empirical_laws
implicit none

contains


!----------------------------------------------------------------------------------------------------------------------------------!

  function erosion(temp,runoff,slope,klith)
    use dynsoil_physical_parameters, only: ke, a, b
    double precision, intent(in):: temp, runoff, slope ! runoff (cm/y) ; temperature (°C) !
    integer, intent(in):: klith
    double precision:: erosion
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !erosion = ke(klith) * ((runoff/100)**a) * (slope**b) * max( 2. , temp ) ! BQART exponents
    erosion = ke(klith) * ((runoff/100)**a) * (slope**b)
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function reg_prod_opt(temp,runoff,klith)
    use dynsoil_physical_parameters, only: krp, Ea_rp, T0_rp, Rgas
    double precision, intent(in):: temp, runoff ! runoff (cm/y) ; temperature (°C) !
    integer, intent(in):: klith
    double precision:: reg_prod_opt
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    reg_prod_opt = krp(klith) * (runoff/100) * exp(-(Ea_rp(klith)/Rgas)*(1./(temp+273.15)-1./T0_rp(klith)))
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   Optimal soil production rate. Should be multiplied by the soil production function to get the actual soil production rate
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function soil_prod_func(h_soil)
    use dynsoil_physical_parameters, only: h0 !k1, d1, d2
    double precision, intent(in):: h_soil
    double precision:: soil_prod_func
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !soil_prod_func = ( exp(-1*h_soil/d1) - k1*exp(-1*h_soil/d2) ) / knorm ! normalized HUMPED FUNCTION
    !soil_prod_func = h0/h_soil                                            ! inverse function
    soil_prod_func = exp(-1*h_soil/h0)                                    ! exponential SPF
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function eq_reg_thick(RPopt,E)
    use dynsoil_physical_parameters, only: h0 !d1
    double precision, intent(in):: RPopt,E
    double precision:: eq_reg_thick
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !eq_reg_thick = d1*log(RPopt/E) ! for humped SPF
    !eq_reg_thick = h0*RPopt/E      ! for inverse SPF
    eq_reg_thick = h0*log(RPopt/E) ! for exponential SPF
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   Estimated soil equilibrium thickness by approximating the SPF.
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function dissolution_constant(temp,runoff,klith)
    use dynsoil_physical_parameters, only: Kwest, kw, Ea, T0, Rgas
    double precision, intent(in):: temp, runoff ! runoff (cm/y) ; temperature (°C) !
    integer, intent(in):: klith
    double precision:: dissolution_constant
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    dissolution_constant = Kwest(klith) * (1-exp(-kw(klith)*(runoff/100))) * exp( (Ea(klith)/Rgas) * (1/T0(klith)-1/(temp+273.15)) )
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   rock dissolution constant, computed from climatic variables
  end function

!----------------------------------------------------------------------------------------------------------------------------------!

  function SP_fraction(temp,runoff,h_soil)
    use dynsoil_physical_parameters, only: xSP_max, xSP_min, reg_porosity, TauSP0
    double precision, intent(in):: temp, runoff, h_soil ! temperature (°C) ; runoff (m/y) ; h_soil (m)
    double precision:: TauWat, TauSP ! water residence time (y) ; SP charcteristic time (y)
    double precision:: SP_fraction
    if (runoff == 0) then
      SP_fraction = xSP_max
    elseif (h_soil == 0) then
      SP_fraction = xSP_min
    else
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      TauWat = reg_porosity * h_soil / runoff
      TauSP = TauSP0 !* exp((Ea_SP/Rgas)*(1./(temp+273.15)-1./T0))
!     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      SP_fraction = xSP_min  +  (xSP_max-xSP_min) / ( 1 + (TauSP/TauWat) )
!     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    end if
  end function SP_fraction

!----------------------------------------------------------------------------------------------------------------------------------!

  function Li_fractionation(temp)
   use dynsoil_physical_parameters, only: aDland, bDland
   double precision, intent(in):: temp ! temperature (°C)
   double precision:: Li_fractionation
   Li_fractionation = aDland/(temp+273.15)**2 + bDland
  end function


end module
