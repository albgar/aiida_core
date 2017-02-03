#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.5.0"
__contributors__ = "Andrea Cepellotti, Victor Garcia-Suarez, Alberto Garcia"

import sys
import os
import pymatgen as mg

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


#####    SET CODE!!!    TO BE MODIFIED   ######
try:
    codename = sys.argv[2]
except IndexError:
    codename = 'SiestaLast@parsons'
code = test_and_get_code(codename, expected_code_type='siesta')


#####     SETTINGS             ####
settings = None


####       STRUCTURE      #####
structure = mg.Structure.from_file("O2_ICSD_173933.cif", primitive=False)
from_struct = StructureData(pymatgen_structure=structure)
elements = list(from_struct.get_symbols_set())


######   PSEUDOS. TO BE MODIFIED   ######
# If auto_pseudo True, load the pseudos from the family specified below
# Otherwise, use static files provided
auto_pseudos = True
if auto_pseudos:
    valid_pseudo_groups = PsfData.get_psf_groups(filter_elements=elements)

    try:
        #pseudo_family = sys.argv[3]
        pseudo_family = '1stGenCorrCorrect'
    except IndexError:
        print >> sys.stderr, "Error, auto_pseudos set to True. You therefore need to pass as second parameter"
        print >> sys.stderr, "the pseudo family name."
        print >> sys.stderr, "Valid PSF families are:"
        print >> sys.stderr, "\n".join("* {}".format(i.name) for i in valid_pseudo_groups)
        sys.exit(1)
    try:
        PsfData.get_psf_group(pseudo_family)
    except NotExistent:
        print >> sys.stderr, "auto_pseudos is set to True and pseudo_family='{}',".format(pseudo_family)
        print >> sys.stderr, "but no group with such a name found in the DB."
        print >> sys.stderr, "Valid PSF groups are:"
        print >> sys.stderr, ",".join(i.name for i in valid_pseudo_groups)
        sys.exit(1)


######   PARAMETERS   ##########
parameters = ParameterData(dict={
                'xc:functional': 'LDA',
                'xc:authors': 'CA',
                'spinpolarized': True,
                'meshcutoff': '40.000 Ry',
                'maxscfiterations': 50,
                'dm:numberpulay': 4,
                'dm:mixingweight': 0.3,
                'dm:tolerance': 1.e-3,
                'Solution-method': 'diagon',
                'electronic-temperature': '25 meV',
                'md:typeofrun': 'cg',
                'md:numcgsteps': 3,
                'md:maxcgdispl': '0.1 Ang',
                'md:maxforcetol': '0.04 eV/Ang',
                'writeforces': True,
                'writecoorstep': True,
                'xml:write': True,
                'dm:usesavedm': True
                })


########   BASIS SET   #########
basis = ParameterData(dict={
'pao-energy-shift': '300 meV',
'%block pao-basis-sizes': """
O SZP                    """,
})


####            K-POINTS for bands, uncomment your favourite                          ###
# NOTE: bandskpoints.set_cell(from_struct.cell, from_struct.pbc) HAS TO BE SET ALWAYS ###
bandskpoints = KpointsData()

##..Set a path, label needed, 40 is number of kp between W-L and between L-G..##
#kpp = [('W',  (0.500,  0.250, 0.750), 'L', (0.500,  0.500, 0.500), 40),
#        ('L', (0.500,  0.500, 0.500), 'G', (0., 0., 0.), 40)]
#bandskpoints.set_cell(from_struct.cell, from_struct.pbc)
#bandskpoints.set_kpoints(kpp)


##..........................Only points, no labels............................##
#kpp = [(0.500,  0.250, 0.750), (0.500,  0.500, 0.500), (0., 0., 0.)]
#bandskpoints.set_cell(from_struct.cell, from_struct.pbc)
#bandskpoints.set_kpoints(kpp)

##..kp path automatically generated from structure (all high-simmetry point)..##
##.....labels automatically included, 0.05 is the distance between kpoints....##
bandskpoints.set_cell(from_struct.cell, from_struct.pbc)
bandskpoints.set_kpoints_path(kpoint_distance = 0.05)


##########            KPOINTS   ################
kts = KpointsData()
kpoints_mesh = 4
kts.set_kpoints_mesh([kpoints_mesh,kpoints_mesh,kpoints_mesh])


#########    SET CALCULATION   #########
calc = code.new_calc()
calc.label = "O_el_cell_spin"
calc.set_max_wallclock_seconds(30*60) # 30 min
calc.set_resources({"num_machines": 1, "num_mpiprocs_per_machine": 8})
queue = "debug"
if queue is not None:
    calc.set_queue_name(queue)
calc.use_structure(from_struct)
calc.use_parameters(parameters)
calc.use_basis(basis)
if auto_pseudos:
    try:
        calc.use_pseudos_from_family(pseudo_family)
        print "Pseudos successfully loaded from family {}".format(pseudo_family)
    except NotExistent:
        print ("Pseudo or pseudo family not found. You may want to load the "
               "pseudo family, or set auto_pseudos to False.")
        raise
else:
    raw_pseudos = [("Si.psf", 'Si')]
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
calc.use_kpoints(kts)
calc.use_bandskpoints(bandskpoints)
if settings is not None:
    calc.use_settings(settings)
calc.description = "Test calculation with the Siesta code. O elementary cell spin polarized"
print calc.description



###########  SEND or JUST PRINT    ##########
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
