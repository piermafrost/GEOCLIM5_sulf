! Define couplings with GEOCLIM modules
logical, parameter :: coupling_dynsoil = .true.
logical, parameter :: use_dynsoil_steady_state = .false.
logical, parameter :: coupling_ecogeo = .false.
logical, parameter :: lock_oxygen_cycle = .false.
logical, parameter :: lock_sulfur_cycle = .false.
character(len=*), parameter :: interpolation_mode = 'log'
