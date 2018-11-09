"""
Parses chemkin input files and sets up the system
"""

from monty.json import MSONable
import cantera.ck2cti

class GasCK(MSONable):
    """
    Parses the text in gas.inp input file of Chemkin to read phases and 
    reactions Alternatively generates the gas.inp file from the supplied 
    elements, phases and reactions

    Uses cantera's ck2cti functionality to read and setup
    """
    def __init__(self, energy_units='cal/mol', quantity_units='mol', 
                 elements=[], element_weights={}, species=[], reactions=[]):
        self.energy_units = energy_units,
        self.quantity_units = quantity_units
        self.elements = elements
        self.element_weights = element_weights
        self.species = species
        self.reactions = reactions

    @classmethod
    def from_file(cls, filename, skip_undeclared_species=True):
        parser = cantera.ck2cti.Parser()
        parser.loadChemkinFile(filename, 
                               skipUndeclaredSpecies=skip_undeclared_species)
        return cls(parser.energy_units, parser.quantity_units, 
                   parser.elements, parser.element_weights, parser.speciesList,
                   parser.reactions)

    def from_string(self, gasck_string):
        pass

    def write(self, filename):
        """Generate Chemkin readable input file"""
        pass

    def as_dict(self):
        """Generate MontyJSON serializable dictionary"""
        pass

    @classmethod
    def from_dict(cls, d):
        """Initialize from dictionary"""
        pass


class SurfaceCK(MSONable):
    """
    Parses the surface.inp input file of Chemkin to read surface phases 
    and surface reactions

    Uses cantera's ck2cti functionality to read and setup
    """
    def __init__(self, energy_units='cal/mol', quantity_units='mol', 
                 elements=[], element_weights={}, species=[], reactions=[]):
        self.energy_units = energy_units,
        self.quantity_units = quantity_units
        self.elements = elements
        self.element_weights = element_weights
        self.species = species
        self.reactions = reactions

    @classmethod
    def from_file(cls, filename, skip_undeclared_species=True):
        parser = cantera.ck2cti.Parser()
        parser.loadChemkinFile(filename, 
                               skipUndeclaredSpecies=skip_undeclared_species)
        print('processed units', parser.processed_units)
        print('energy units', parser.energy_units)
        print('quantity units', parser.quantity_units)
        print('motz_wise', parser.motz_wise)
        print('warnings_as_error', parser.warning_as_error)
        print('elements', parser.elements)
        print('element weights', parser.element_weights)
        print('species', parser.speciesList)
        print('species dict', parser.speciesDict)
        print('surfaces', parser.surfaces)
        print('reactions', parser.reactions)
        print('finalReactionComment', parser.finalReactionComment)
        print('headerLines', parser.headerLines)
        return cls(parser.energy_units, parser.quantity_units, 
                   parser.elements, parser.element_weights, parser.speciesList,
                   parser.reactions)

    def from_string(self, surfck_string):
        pass

    def write(self, filename):
        """Generate Chemkin readable input file"""
        pass

    def as_dict(self):
        """Generate MontyJSON serializable dictionary"""
        pass

    @classmethod
    def from_dict(cls, d):
        """Initialize from dictionary"""
        pass

ass

class TubeCK(MSONable):
    """
    Parses the tube.inp input file to setup reactor conditons

    """
    pass

class ThermDat(MSONable):
    pass
