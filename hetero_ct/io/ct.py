"""
Redefine the classes defined in cantera to add as_dict method
"""

import cantera as ct
from monty.json import MSONable


class Species(ct.ck2cti.Species, MSONable):
    def __init__(self, *args, **kwargs):
        ct.ck2cti.Species.__init__(self, *args, **kwargs)

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


class NASA(ThermoModel, MSONable):
    """
    TODO: Coping the NASA code from ct. Violates the DRY principle.
    TODO: Find alternative approach.

    """
    def __init__(self, coeffs, **kwargs):
        ThermoModel.__init__(self, **kwargs)
        if len(coeffs) not in (7, 9):
            raise ct.ck2cti.InputParseError('Invalid number of NASA polynomial coefficients; '
                                  'should be 7 or 9.')
        self.coeffs = coeffs

    def to_cti(self, indent=0):
        prefix = ' ' * indent
        if len(self.coeffs) == 7:
            vals = ['{0: 15.8E}'.format(i) for i in self.coeffs]
            lines = ['NASA([{0:.2f}, {1:.2f}],'.format(self.Tmin[0], self.Tmax[0]),
                     prefix + '     [{0}, {1}, {2},'.format(*vals[0:3]),
                     prefix + '      {0}, {1}, {2},'.format(*vals[3:6]),
                     prefix + '      {0}]),'.format(vals[6])]
        else:
            vals = ['{0: 15.9E}'.format(i) for i in self.coeffs]
            lines = ['NASA9([{0:.2f}, {1:.2f}],'.format(self.Tmin[0], self.Tmax[0]),
                     prefix + '      [{0}, {1}, {2},'.format(*vals[0:3]),
                     prefix + '       {0}, {1}, {2},'.format(*vals[3:6]),
                     prefix + '       {0}, {1}, {2}]),'.format(*vals[6:9])]

        return '\n'.join(lines)

    def as_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "coeffs": self.coeffs
             }
        thermo_d = ThermoModel.as_dict(self)
        del thermo_d["@module"]
        del thermo_d["@class"]
        d.update(thermo_d)
        return d

    @classmethod
    def from_dict(cls, d):
        d.pop("@module", None)
        d.pop("@class", None)
        coeff = d.pop("coeffs") # This line may not be required.
        return cls(coeff, **d)


class MultiNASA(ThermoModel, MSONable):
    """
    TODO: Coping the MultiNASA code from ct. Violates the DRY principle.
    TODO: Find alternative approach.

    """

    def __init__(self, polynomials, **kwargs):
        ThermoModel.__init__(self, **kwargs)
        self.polynomials = polynomials or []

    def to_cti(self, indent=0):
        prefix = ' ' * indent
        lines = []
        for i, p in enumerate(self.polynomials):
            if i == 0:
                if len(self.polynomials) == 1:
                    lines.append('({0})'.format(p.to_cti(indent + 1)))
                else:
                    lines.append('({0}'.format(p.to_cti(indent + 1)))
            elif i != len(self.polynomials) - 1:
                lines.append(prefix + ' {0}'.format(p.to_cti(indent + 1)))
            else:
                lines.append(prefix + ' {0})'.format(p.to_cti(indent + 1)[:-1]))

        return '\n'.join(lines)

    def as_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "polynomials": self.polynomials
             }
        thermo_d = ThermoModel.as_dict(self)
        del thermo_d["@module"]
        del thermo_d["@class"]
        d.update(thermo_d)
        return d

    @classmethod
    def from_dict(cls, d):
        d.pop("@module", None)
        d.pop("@class", None)
        polynomials = d.pop("polynomials")  # This line may not be required.
        return cls(polynomials, **d)


class Reaction(ct.ck2cti.Reaction, MSONable):
    def __init__(self, *args, **kwargs):
        ct.ck2cti.Reaction.__init__(self, *args, **kwargs)

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