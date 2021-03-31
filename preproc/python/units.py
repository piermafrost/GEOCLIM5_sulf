'''
Define a bunch a accepted units for several physical quantities (temperature, ...)
'''



class units():
    '''
    "units" object contains a reference unit name, and a list of accepted units name,
    and corresponding list factors and offsets to convert the accepted units in the
    reference unit.
    '''

    def __init__(self, reference_units, accepted_units, conversion_factor, conversion_offset=None):
        '''
        Create a "units" object.
        Mandatory arguments:
          * reference_units: string. Name of the reference units
          * accepted_units: list. Names of the accepted units. The element of the list
                            must be strings, of lists of string (if several accepted
                            units have the same conversion parameters).
          * conversion_factor: list of numbers (real of integer): factor to convert the accepted
                               unit into the reference unit
        Optional arguments:
          * conversion_offset: list of numbers (real of integer): offset to convert the accepted
                               unit into the reference unit. By default, it is 0.
        
        NOTE: the unit conversion is 'factor*value + offset'
        '''

        # Input arguments
        # ---------------


        # reference units

        if isinstance(reference_units, str):
            self.reference = reference_units
        else:
            raise TypeError('"reference_units" argument must be a string')


        # list of accepted units

        if isinstance(accepted_units, list):
            passed = True
            for element in accepted_units:
                if isinstance(element, str):
                    loc_passed = True
                elif isinstance(element, list):
                    loc_passed = True
                    for el in element:
                        if not isinstance(el, str):
                            loc_passed = False

                else:
                    loc_passed = False

                if not loc_passed:
                    passed = False

        else:
            passed = False
                
        if passed:
            self.accepted = accepted_units
        else:
            raise TypeError('"accepted_units" argument must be a list of strings, or of lists of strings')


        # List of conversion factors

        if isinstance(conversion_factor, list):
            if len(conversion_factor) == len(accepted_units):
                self.factor = conversion_factor
            else:
                raise ValueError('"accepted_units" and "conversion_factor" arguments must have the same length')
        
        else:
            raise TypeError('"conversion_factor" must be a list of numbers')


        # List of conversion offsets

        if conversion_offset is None:

            self.offset = len(accepted_units)*[0]

        else:

            if isinstance(conversion_offset, list):
                if len(conversion_offset) == len(accepted_units):
                    self.offset = conversion_offset
                else:
                    raise ValueError('"accepted_units" and "conversion_offset" arguments must have the same length')
            
            else:
                raise TypeError('"conversion_offset" must be a list of numbers')




    # Class methods
    # -------------

    def convert(self, value, unit):
        '''
        convert the value "value" of a given unit "unit" into the reference unit
        of the object, if it belongs tothe list of accepted units.
        "value" must be number of array of number
        "unit" must be a string
        '''
        if unit == self.reference:
            return value

        for i,known in enumerate(self.accepted):
            if isinstance(known, list):
                for loc_known in known:
                    if unit == loc_known:
                        return self.factor[i]*value + self.offset[i]

            else:
                if unit == known:
                    return self.factor[i]*value + self.offset[i]

        # If unit does not match any accepted unit:
        raise ValueError('Not recognized unit "'+str(unit)+'"')




#################################################



temperature_units = units('degrees_celsius',
                          [['C', '°', '°C', 'deg C', 'deg Celsius', 'deg celsius', 'degree Celsius', 'degree celsius', 'degrees Celsius',
                            'degrees celsius', 'degree_Celsius', 'degree_celsius', 'degrees_Celsius'],
                           ['K', 'Kelvin', 'kelvin']],
                          [1, 1],
                          [0, -273.15])

length_units = units('m',
                     [['meter', 'meters'],
                      ['cm', 'centimeter', 'centimeters'],
                      ['mm', 'millimeter', 'millimeters'],
                      ['km', 'kilometer', 'kilometers'],
                      ['1e6m', '1e6 m', 'Mm', 'megameter', 'megameters'],
                      ['1e6km', '1e6 km', 'Tm', 'terameter', 'terameters']],
                     [1, 1e-2, 1e-3, 1e3, 1e6, 1e12])

latitude_units = units('degrees_north',
                       [['degree_north', 'degrees_North', 'degree_North', 'degrees north', 'degree north', 'degrees North', 'degree North',
                         'deg N', 'deg North', 'deg north', 'deg', 'degree', 'degrees', '°'],
                        ['rad', 'radian', 'radians']],
                       [1, 3.14159265358979/180])


