import re

"""

Contains dictionary of units conversions to and from atomic units

"""


units = {
    "au_per_ang": 1.8897261328856432e00,  # au of length per Angstrom, PSI4
    "au_per_amu": 1.8228884840986366e03,  # au of mass per atomic mass unit, PSI4
    "au_per_debeye": 3.9343030716327643e-01,  # au of dipole per Debeye, PSI4
    "au_per_ev": 3.6749330610942925e-02,  # au of energy per eV, PSI4
    "au_per_kcal": 1.5936013717720609e-03,  # au of energy per kcal mol^-1, PSI4
    "au_per_cminv": 4.5563359040180500e-06,  # au of energy per cm^{-1} ("wavenumber"), PSI4
    "au_per_K": 1.0
    / 3.1577464e5,  # temperature (TODO: This is probably OK, maybe a bit low precision)
    "au_per_fs": 1.0 / 2.418884326505e-2,  # time
}
for k in list(units.keys()):
    v = units[k]
    mobj = re.match("(\S+)_per_(\S+)", k)
    units["%s_per_%s" % (mobj.group(2), mobj.group(1))] = 1.0 / v
