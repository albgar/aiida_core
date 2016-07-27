siesta
++++++

Description
-----------

A plugin for Siesta. This is work in progress. Basic functionality is
implemented, but there remain some important points to address.

(These docs for version: `0.4` of the plugin)

Supported Siesta versions
-------------------------

In principle, all modern Siesta versions (>= 3.X) should work, but
some extra functionality might only be available for recent versions.
This will be documented more fully when v1.0 of the plugin is finalized.

Inputs
------

* **structure**, class :py:class:`StructureData <aiida.orm.data.structure.StructureData>`
    A structure. Siesta employs "species labels" to implement special
    conditions (such as basis set characteristics) to specific atoms
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

    
* **parameters**, class :py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>`
    A dictionary with (so far) scalar fdf variables. Units are
    specified for now as part of the value string. For example::

        {
          "mesh-cutoff": "200 Ry",
          "dm-tolerance": "0.0001",
        }


* **pseudo**, class :py:class:`PsfData <aiida.orm.data.psf.PsfData>`

  The PsfData class has been implemented along the lines of the Upf class for QE.

  One pseudopotential file per atomic element. Several species in the
  Siesta sense can share the same pseudopotential. For the example
  above::

    pseudo_file_to_species_map = [ ("C.psf", ['C', 'Cred']),
                              ("H.psf", 'H')
			    ]
  
  
  Alternatively, a pseudo for every atomic species can be set with the **use_pseudos_from_family**
  method, if a family of pseudopotentials has been installed.

* **basis**, class :py:class:`ParameterData  <aiida.orm.data.parameter.ParameterData>`
  
    A dictionary specifically intended for basis set
    information. Currently there is only a stub implementation
    targeting the `pao-basis-sizes` block, but soon there will be
    more functionality (most likely implemented as raw fdf-block
    content in a first round).

* **kpoints**, class :py:class:`KpointsData <aiida.orm.data.array.kpoints.KpointsData>`
  Reciprocal space points on which to build the wavefunctions. Can either be 
  a mesh or a list of points with/without weights

Outputs
-------
The output parser takes advantage of the structured output available
in Siesta as a CML file. The CML-writer functionality should be
compiled in and active in the run!

* **output_parameters** :py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>` 
  (accessed by ``calculation.res``)

    A dictionary with metadata and energy values. Units are specified
    by means of the second component of a tuple. Suggestions welcome::

        {
          "siesta:Version": "siesta-4.0-540",
          "siesta:E_fermi": ("-3.24", "eV"),
          "siesta:Free_EK": ("-6656.2343", "eV")
	}

* **output_array** :py:class:`ArrayData <aiida.orm.data.array.ArrayData>`

  Contains forces (eV/Angstrom) and stresses (GPa).
  

* **output_structure** :py:class:`StructureData <aiida.orm.data.structure.StructureData>`

  Present only if the calculation is moving the ions.
  Cell and ionic positions refer to the last configuration.

