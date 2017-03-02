Standard Siesta plugin
++++++++++++++++++++++

Description
-----------

A plugin for Siesta's basic functionality. There remain some details to address.

These docs are for version: *aiida-0.7--plugin-0.6.0* of the plugin,
which is compatible with AiiDA v.0.7.0 (modified to include a kpoints.py
module from v0.7.1)

Supported Siesta versions
-------------------------

In principle, all modern Siesta versions (>= 3.X) should work, but
some extra functionality (i.e., warnings handling and consistency of
the CML file in the face of early termination) might only be available
for recent versions. We are preparing a special set of patches for the
already-released 4.0 version that will make it compatible with the
features of the parser. It is work in progress in the '4.0-aiida'
branch of the SIESTA Launchpad platform (http://launchpad.net/siesta/)
and should be finished soon.

This will be documented more fully when v1.0 of the plugin is finalized.

Inputs
------

* **structure**, class :py:class:`StructureData <aiida.orm.data.structure.StructureData>`
A structure. Siesta employs "species labels" to implement special
conditions (such as basis set characteristics) for specific atoms
(e.g., surface atoms might have a richer basis set). This is
implemented through the 'name' attribute of the Site objects. For example::

  alat = 15. # angstrom
  cell = [[alat, 0., 0.,],
    [0., alat, 0.,],
    [0., 0., alat,],
   ]

   # Benzene molecule with a special carbon atom
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


This feature works well in normal use, but the association of
pseudos to multiple species is not currently working for restarts.
    
* **parameters**, class :py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>`
A dictionary with scalar fdf variables and blocks, which are the
basic elements of any Siesta input file. A given Siesta fdf file
can be cast almost directly into this dictionary form, except that
some items (for structure data) are blocked. Any units are
specified for now as part of the value string. Blocks are entered
by using an appropriate key and Python's multiline string
constructor. For example::

    {
      "mesh-cutoff": "200 Ry",
      "dm-tolerance": "0.0001",
	  "%block example-block": """
	  first line
	  second line             """,
    }

Note that Siesta fdf keywords allow '.', '-', or nothing as
internal separators. AiiDA does not allow the use of '.' in
nodes to be inserted in the database, so it should not be used
in the input script (or removed before assigning the dictionary to
the ParameterData instance).

* **pseudo**, class :py:class:`PsfData <aiida.orm.data.psf.PsfData>`

The PsfData class has been implemented along the lines of the Upf class for QE.

One pseudopotential file per atomic element. Several species in the
Siesta sense can share the same pseudopotential. For the example
above::

  pseudo_file_to_species_map = [ ("C.psf", ['C', 'Cred']),
                            ("H.psf", 'H')
			    ]


Alternatively, a pseudo for every atomic species can be set with the
**use_pseudos_from_family**  method, if a family of pseudopotentials
has been installed. (But the family approach will not support
multiple species sharing the same pseudopotential.)

* **basis**, class :py:class:`ParameterData  <aiida.orm.data.parameter.ParameterData>`
  
A dictionary specifically intended for basis set information. It
follows the same structure as the **parameters** element, including
the allowed use of fdf-block items. This raw interface allows a
direct translation of the myriad basis-set options supported by the
program. In future we might have a more structured input for
basis-set information.

* **kpoints**, class :py:class:`KpointsData <aiida.orm.data.array.kpoints.KpointsData>`
  
Reciprocal space points for the full sampling of the BZ during the
self-consistent-field iteration. It must be given in mesh form. There is no support
yet for Siesta's kgrid-cutoff keyword.

If this node is not present, only the Gamma point is used for sampling.

* **bandskpoints**, class :py:class:`KpointsData
  <aiida.orm.data.array.kpoints.KpointsData>`
  
Reciprocal space points for the calculation of bands.  They can be
given as a simple list of k-points, as segments with start and end
point and number of points, or as a complete automatic path, using the
functionality of modern versions of the class.

If this node is not present, no band structure is computed.

* **settings**, class
  :py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>`
      
An optional dictionary that activates non-default operations. For a list of possible
values to pass, see the section on :ref:`advanced features <siesta-advanced-features>`.

Outputs
-------

There are several output nodes that can be created by the plugin,
according to the calculation details.  All output nodes can be
accessed with the ``calculation.out`` method.

The output parser takes advantage of the structured output available
in Siesta as a Chemical Markup Language (CML) file. The CML-writer
functionality should be compiled in and active in the run!

* **output_parameters** :py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>` 
  (accessed by ``calculation.res``)

A dictionary with metadata, scalar result values, and a warnings
list.  Units are specified by means of an extra item with '_units'
appended to the key::

    {
      "siesta:Version": "siesta-4.0-540",
      "E_fermi": "-3.24",
	  "E_fermi_units": "eV",
      "Free_EK": "-6656.2343"
	  "Free_EK_units": "eV",
	  "warnings": [ "INFO: Job Completed"]
	}

The 'warnings' list contains program messages, labeled as INFO,
WARNING, or FATAL, read directly from a MESSAGES file produced by
Siesta, which include items from the execution of the program and
also a possible 'out of time' condition. This is implemented by
passing to the program the wallclock time specified in the script,
and checking at each scf step for the walltime consumed. This
'warnings' list can be examined by the parser itself to raise an
exception in the FATAL case.

* **output_array** :py:class:`ArrayData <aiida.orm.data.array.ArrayData>`

Contains the final forces (eV/Angstrom) and stresses (GPa) in array form.
  

* **output_structure** :py:class:`StructureData
  <aiida.orm.data.structure.StructureData>`
  
Present only if the calculation is moving the ions.  Cell and ionic
positions refer to the last configuration.

* **bands_array**, :py:class:`BandsData
  <aiida.orm.data.array.bands.BandsData>`
  
Present only if a band calculation is requested (signaled by the
presence of a 'bandskpoints' input node of class KpointsData)
Contains the list of electronic energies for every kpoint. For
spin-polarized calculations, the 'bands' array has an extra dimension
for spin.
  
No trajectories have been implemented yet.

Errors
------

Errors of the parsing are reported in the log of the calculation (accessible 
with the ``verdi calculation logshow`` command). 
Moreover, they are stored in the ParameterData under the key ``warnings``, and are
accessible with ``Calculation.res.warnings``.

.. _siesta-advanced-features:

Additional advanced features
----------------------------

While the input link with name **parameters** is used for the main
Siesta options (as would be given in an fdf file), additional settings
can be specified in the **settings** input, also of type ParameterData.

Below we summarise some of the options that you can specify, and their effect.
In each case, after having defined the content of ``settings_dict``, you can use
it as input of a calculation ``calc`` by doing::

  calc.use_settings(ParameterData(dict=settings_dict))

The keys of the settings dictionary are internally converted to
uppercase by the plugin.

Adding command-line options
...........................

If you want to add command-line options to the executable (particularly 
relevant e.g. to tune the parallelization level), you can pass each option 
as a string in a list, as follows::

  settings_dict = {  
      'cmdline': ['-option1', '-option2'],
  }

Note that very few user-level comand-line options (besides those
already inserted by AiiDA for MPI operation) are currently implemented.

Retrieving more files
.....................

If you know that your calculation is producing additional files that you want to
retrieve (and preserve in the AiiDA repository in the long term), you can add
those files as a list as follows::


  settings_dict = {  
    'additional_retrieve_list': ['aiida.EIG', 'aiida.ORB_INDX'],
  }


