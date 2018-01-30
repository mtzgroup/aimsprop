import re

"""

Contains dictionary of units conversions to and from atomic units

"""


units = {
    'au_per_amu'   : 1.8228884855409500E+03,        # mass
    'au_per_cminv' : 1.0 / 219474.6305,             # ???
    'au_per_ang'   : 1.0 / 0.5291772109217,         # length
    'au_per_K'     : 1.0 / 3.1577464E5,             # temperature
    'au_per_fs'    : 1.0 / 2.418884326505E-2,       # time
}
for k in units.keys():
    v = units[k]
    mobj = re.match('(\S+)_per_(\S+)',k)
    units['%s_per_%s' % (mobj.group(2),mobj.group(1))] = 1.0 / v
    

# charge
# energy
    # kJ
    # kcal
    # nm
    # eV
# time
    # fs
# temperature
    # K
# electric field
# electric dipole moment
# Constants (yes!)
    # hbar
    # kb
    # c (137)
    
