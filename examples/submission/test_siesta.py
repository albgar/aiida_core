#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.5.0"
__contributors__ = "Andrea Cepellotti, Victor Garcia-Suarez, Alberto Garcia"

import sys
import os

from aiida.common.example_helpers import test_and_get_code
from aiida.common.exceptions import NotExistent

################################################################

PsfData = DataFactory('psf')
ParameterData = DataFactory('parameter')
KpointsData = DataFactory('array.kpoints')
StructureData = DataFactory('structure')

try:
    dontsend = sys.argv[1]
    if dontsend == "--dont-send":
        submit_test = True
    elif dontsend == "--send":
        submit_test = False
    else:
        raise IndexError
except IndexError:
    print >> sys.stderr, ("The first parameter can only be either "
                          "--send or --dont-send")
    sys.exit(1)

try:
    codename = sys.argv[2]
except IndexError:
    codename = 'Siesta-4.0@rinaldo'


queue = None
settings = None
#####

code = test_and_get_code(codename, expected_code_type='siesta')

alat = 15. # angstrom
cell = [[alat, 0., 0.,],
        [0., alat, 0.,],
        [0., 0., alat,],
       ]

# Benzene molecule
# Note an atom tagged (for convenience) with a different label
#
s = StructureData(cell=cell)
s.append_atom(position=(0.000,0.000,0.468),symbols=['H'])
s.append_atom(position=(0.000,0.000,1.620),symbols=['C'])
s.append_atom(position=(0.000,-2.233,1.754),symbols=['H'])
s.append_atom(position=(0.000,2.233,1.754),symbols=['H'])
s.append_atom(position=(0.000,-1.225,2.327),symbols='C',name="Cred")
s.append_atom(position=(0.000,1.225,2.327),symbols=['C'])
s.append_atom(position=(0.000,-1.225,3.737),symbols=['C'])
s.append_atom(position=(0.000,1.225,3.737),symbols=['C'])
s.append_atom(position=(0.000,-2.233,4.311),symbols=['H'])
s.append_atom(position=(0.000,2.233,4.311),symbols=['H'])
s.append_atom(position=(0.000,0.000,4.442),symbols=['C'])
s.append_atom(position=(0.000,0.000,5.604),symbols=['H'])

elements = list(s.get_symbols_set())

parameters = ParameterData(dict={
'pao:basistype': 'split',
'pao:splitnorm': 0.150,
'pao:energyshift': '0.020 Ry',
'xc:functional': 'LDA',
'xc:authors': 'CA',
'spinpolarized': True,
'noncollinearspin': False,
'meshcutoff': '200.000 Ry',
'maxscfiterations': 1000,
'dm:numberpulay': 5,
'dm:mixingweight': 0.050,
'dm:tolerance': 1.e-4,
'dm:mixscf1': True,
'neglnonoverlapint': False,
'solutionmethod': 'diagon',
'electronictemperature': '100.000 k',
'# md:typeofrun': 'cg',
'# md:numcgsteps': 1000,
'# md:maxcgdispl': '0.200 bohr',
'# md:maxforcetol': '0.050 ev/ang',
'writeforces': True,
'writecoorstep': True,
'xml:write': True,
'writemullikenpop': 1,
'md:usesavexv': True,
'dm:usesavedm': True
})

#
#  This entry is handled in a special way, interpreted as
#  a pao-basis-sizes block
#
basis = ParameterData(dict={
                'C': 'SZP',
                'Cred': 'SZ',
                'H': 'SZP',
                })

kpoints = KpointsData()

# method mesh
kpoints_mesh = 1
kpoints.set_kpoints_mesh([kpoints_mesh,kpoints_mesh,kpoints_mesh])

# to retrieve the bands
# (the object settings is optional)
#settings_dict={'also_bands': True}
#settings = ParameterData(dict=settings_dict)

## For remote codes, it is not necessary to manually set the computer,
## since it is set automatically by new_calc
#computer = code.get_remote_computer()
#calc = code.new_calc(computer=computer)

calc = code.new_calc()
calc.label = "Test Siesta. Benzene molecule"
calc.description = "Test calculation with the Siesta code. Benzene molecule"
calc.set_max_wallclock_seconds(30*60) # 30 min

calc.set_resources({"num_machines": 1})

if queue is not None:
    calc.set_queue_name(queue)

calc.use_structure(s)
calc.use_parameters(parameters)
calc.use_basis(basis)

#
# This exemplifies the handling of pseudos for different species
# Those sharing the same pseudo should be indicated.
# Families support is not yet available for this.
#
raw_pseudos = [ ("C.psf", ['C', 'Cred']),
                ("H.psf", 'H')]

for fname, kinds, in raw_pseudos:
    absname = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                            "data",fname))
    pseudo, created = PsfData.get_or_create(absname,use_first=True)
    if created:
        print "Created the pseudo for {}".format(kinds)
    else:
        print "Using the pseudo for {} from DB: {}".format(kinds,pseudo.pk)
        
    # Attach pseudo node to the calculation
    calc.use_pseudo(pseudo,kind=kinds)

calc.use_kpoints(kpoints)

if settings is not None:
    calc.use_settings(settings)
#from aiida.orm.data.remote import RemoteData
#calc.set_outdir(remotedata)

if submit_test:
    subfolder, script_filename = calc.submit_test()
    print "Test_submit for calculation (uuid='{}')".format(
        calc.uuid)
    print "Submit file in {}".format(os.path.join(
        os.path.relpath(subfolder.abspath),
        script_filename
        ))
else:
    calc.store_all()
    print "created calculation; calc=Calculation(uuid='{}') # ID={}".format(
        calc.uuid,calc.dbnode.pk)
    calc.submit()
    print "submitted calculation; calc=Calculation(uuid='{}') # ID={}".format(
        calc.uuid,calc.dbnode.pk)

