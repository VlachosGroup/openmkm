"""
Redefine the classes defined in cantera to add as_dict method
"""

import cantera as ct
from monty.json import MSONable

class Species(ct.ck2cti.Species, MSONable):
    def __init__(self, *args, **kwargs):
        ck.ck2cti.Species.__init__(self, *args, **kwargs)

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "label": self.label,
                "thermo": self.thermo,
                "transport": self.transport,
                "note": self.note,
                "sites": self.sites,
                "composiiton": self.composition
                }

    @classmethod
    def from_dict(cls, d):
        species = cls(d['label'], d['sites'])
        if d['thermo']:
            species.thermo = ThermoModel.from_dict(d['thermo'])
        if d['transport']:
            species.transport = TransportModel.from_dict(d['transport'])
        species.note = d.get('note', None)
        species.composition = d.get('composition', None)


class ThermoModel(ct.ck2cti.ThermoModel, MSONable):
    def __init__(self, *args, **kwargs):
        ct.ck2cti.ThermoModel.__init__(self, *args, **kwargs)

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "Tmin": self.Tmin,
                "Tmax": self.Tmax,
                "comment": self.comment
                }

    @classmethod
    def from_dict(cls, d):
        return cls(d['Tmin'], d['Tmax'], d.get('comment', ''))


        
