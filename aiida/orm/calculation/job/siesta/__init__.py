# -*- coding: utf-8 -*-
import os

from aiida.orm.calculation.job import JobCalculation
from aiida.common.exceptions import InputValidationError
from aiida.orm import DataFactory
from aiida.common.datastructures import CalcInfo
from aiida.orm.data.psf import get_pseudos_from_structure
from aiida.common.utils import classproperty
from aiida.common.datastructures import CodeInfo

from aiida.orm.data.structure import StructureData
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.psf import PsfData
from aiida.orm.data.remote import RemoteData 

from aiida.common.constants import elements

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.7.0"
__contributors__ = "Victor M. Garcia-Suarez, Alberto Garcia, Emanuele Bosoni"

class SiestaCalculation(JobCalculation):
    """
    Add docs
    """
    _siesta_plugin_version = 'aiida-0.7--plugin-0.6.0'
    
    def _init_internal_params(self):
        super(SiestaCalculation, self)._init_internal_params()

        # Default Siesta output parser provided by AiiDA
        self._default_parser = "siesta"

        # Keywords that cannot be set
        # We need to canonicalize this!
        self._aiida_blocked_keywords = ['systemname','systemlabel']
        self._aiida_blocked_keywords.append('system-name')
        self._aiida_blocked_keywords.append('system-label')
        self._atom_blocked_keywords = ['numberofspecies','numberofatoms']
        self._atom_blocked_keywords.append('number-of-species')
        self._atom_blocked_keywords.append('number-of-atoms')
        #
        #
        self._atom_blocked_keywords.append('latticeconstant')
        self._atom_blocked_keywords.append('lattice-constant')
        self._atom_blocked_keywords.append('atomic-coordinates-format')
        self._atom_blocked_keywords.append('atomiccoordinatesformat')

        # Default input and output files                                        
        self._DEFAULT_INPUT_FILE = 'aiida.in'
        self._DEFAULT_OUTPUT_FILE = 'aiida.out'
        self._DEFAULT_XML_FILE = 'aiida.xml'
        self._DEFAULT_MESSAGES_FILE = 'MESSAGES'
	self._DEFAULT_BANDS_FILE = 'aiida.bands'

        self._PSEUDO_SUBFOLDER = './'
        self._OUTPUT_SUBFOLDER = './'
        self._PREFIX = 'aiida'
        self._INPUT_FILE_NAME = 'aiida.fdf'
        self._OUTPUT_FILE_NAME = 'aiida.out'
        self._XML_FILE_NAME = 'aiida.xml'
        self._MESSAGES_FILE_NAME = 'MESSAGES'
	self._BANDS_FILE_NAME = 'aiida.bands'

        # in restarts, it will copy from the parent the following
        # (fow now, just the density matrix file)
        self._restart_copy_from = os.path.join(self._OUTPUT_SUBFOLDER, '*.DM')

        # in restarts, it will copy the previous folder in the following one
        self._restart_copy_to = self._OUTPUT_SUBFOLDER


    @classproperty
    def _use_methods(cls):
        """
        Extend the parent _use_methods with further keys.
        """
        retdict = JobCalculation._use_methods

        retdict['kpoints'] = {
               'valid_types': KpointsData,
               'additional_parameter': None,
               'linkname': 'kpoints',
               'docstring': "Use the node defining the kpoint sampling to use",
               }
        retdict["structure"] = {
            'valid_types': StructureData,
            'additional_parameter': None,
            'linkname': 'structure',
            'docstring': "Choose the input structure to use",
            }
        retdict["basis"] = {
            'valid_types': ParameterData,
            'additional_parameter': None,
            'linkname': 'basis',
            'docstring': "Choose the input basis to use",
            }
        retdict["settings"] = {
            'valid_types': ParameterData,
            'additional_parameter': None,
            'linkname': 'settings',
            'docstring': "Use an additional node for special settings",
            }
        retdict["parameters"] = {
            'valid_types': ParameterData,
            'additional_parameter': None,
            'linkname': 'parameters',
            'docstring': ("Use a node that specifies the input parameters "
                          "for the namelists"),
            }
        retdict["parent_folder"] = {
            'valid_types': RemoteData,
            'additional_parameter': None,
            'linkname': 'parent_calc_folder',
            'docstring': ("Use a remote folder as parent folder (for "
                          "restarts and similar"),
            }
        retdict["pseudo"] = {
            'valid_types': PsfData,
            'additional_parameter': "kind",
            'linkname': cls._get_linkname_pseudo,
            'docstring': ("Use a node for the PSF pseudopotential of one of "
                          "the elements in the structure. You have to pass "
                          "an additional parameter ('kind') specifying the "
                          "name of the structure kind (i.e., the name of "
                          "the species) for which you want to use this "
                          "pseudo. You can pass either a string, or a "
                          "list of strings if more than one kind uses the "
                          "same pseudo"),
            }
	retdict['bandskpoints'] = {
               'valid_types': KpointsData,
               'additional_parameter': None,
               'linkname': 'bandskpoints',
               'docstring': "Use the node defining the kpoint sampling to use for bands calculation",
               }
        return retdict

    def _prepare_for_submission(self,tempfolder,
                                    inputdict):        
        """
        This is the routine to be called when you want to create
        the input files and related stuff with a plugin.
        
        :param tempfolder: a aiida.common.folders.Folder subclass where
                           the plugin should put all its files.
        :param inputdict: a dictionary with the input nodes, as they would
                be returned by get_inputdata_dict (without the Code!)
        """
        from aiida.common.utils import get_unique_filename, get_suggestion
        import re
        
        local_copy_list = []
        remote_copy_list = []

        # Process the settings dictionary first
        # Settings can be undefined, and defaults to an empty dictionary
        settings = inputdict.pop(self.get_linkname('settings'),None)
        if settings is None:
            settings_dict = {}
        else:
            if not isinstance(settings,  ParameterData):
                raise InputValidationError("settings, if specified, must be of "
                                           "type ParameterData")
            
            # Settings converted to UPPERCASE
            # Presumably to standardize the usage and avoid
            # ambiguities
            settings_dict = _uppercase_dict(settings.get_dict(),
                                            dict_name='settings')

        try:
            parameters = inputdict.pop(self.get_linkname('parameters'))
        except KeyError:
            raise InputValidationError("No parameters specified for this "
                "calculation")
        if not isinstance(parameters, ParameterData):
            raise InputValidationError("parameters is not of type "
                "ParameterData")

        # Basis can be undefined, and defaults to an empty dictionary,
        # Siesta will use default parameters
        basis = inputdict.pop(self.get_linkname('basis'),None)
        if basis is None:
            input_basis = {}
        else:
            if not isinstance(basis, ParameterData):
                raise InputValidationError("basis not of type ParameterData")
            input_basis=basis.get_dict()

        try:
            structure = inputdict.pop(self.get_linkname('structure'))
        except KeyError:
            raise InputValidationError("No structure specified for this "
                "calculation")
        if not isinstance(structure,  StructureData):
            raise InputValidationError("structure is not of type StructureData")

        # k-points
        # It is now possible to elide the kpoints node.
        #
        # Note also that a *different* set of k-points is needed if a band
        # calculation is carried out. This should be specified somehow in
        # the 'settings' dictionary (see QE example...)
        
        kpoints = inputdict.pop(self.get_linkname('kpoints'),None)
        if kpoints is None:
            # Do nothing. Assume it is a gamma-point calculation
            pass
        else:
            if not isinstance(kpoints,  KpointsData):
                raise InputValidationError("kpoints, if specified, must be of "
                                           "type KpointsData")

        bandskpoints = inputdict.pop(self.get_linkname('bandskpoints'),None)
        if bandskpoints is None:
            flagbands=False
        else:
	    flagbands=True
	    if not isinstance(bandskpoints,  KpointsData):
                raise InputValidationError("kpoints for bands is not of type KpointsData")



        
        pseudos = {}
        # I create here a dictionary that associates each kind name to a pseudo
        for link in inputdict.keys():
            if link.startswith(self._get_linkname_pseudo_prefix()):
                kindstring = link[len(self._get_linkname_pseudo_prefix()):]
                kinds = kindstring.split('_')
                the_pseudo = inputdict.pop(link)
                if not isinstance(the_pseudo, PsfData):
                    raise InputValidationError("Pseudo for kind(s) {} is not "
                    " of type PsfData".format(",".join(kinds)))
                #
                # Note that we can associate the same pseudo object to different
                # atom kinds
                #
                for kind in kinds:
                    pseudos[kind] = the_pseudo

        parent_calc_folder = inputdict.pop(self.get_linkname('parent_folder'),
            None)
        if parent_calc_folder is not None:
            if not isinstance(parent_calc_folder,  RemoteData):
                raise InputValidationError("parent_calc_folder, if specified,"
                    "must be of type RemoteData")

        try:
            code = inputdict.pop(self.get_linkname('code'))
        except KeyError:
            raise InputValidationError("No code specified for this calculation")
                                
        # Here, there should be no more parameters...
        if inputdict:
            raise InputValidationError("The following input data nodes are "
                "unrecognized: {}".format(inputdict.keys()))

        # Check structure, get species, check peudos
        kindnames = [k.name for k in structure.kinds]
        if set(kindnames) != set(pseudos.keys()):
            err_msg = ("Mismatch between the defined pseudos and the list of "
                       "kinds of the structure. Pseudos: {}; kinds: {}".format(
                        ",".join(pseudos.keys()), ",".join(list(kindnames))))
            raise InputValidationError(err_msg)

        ##############################
        # END OF INITIAL INPUT CHECK #
        ##############################


        # First-level keys as uppercase (i.e., namelist and card names)
        # and the second-level keys as lowercase
        # (deeper levels are unchanged)
        orig_input_params = _lowercase_dict(parameters.get_dict(),
            dict_name='parameters')

        # This is to be removed. Sanitization should be done
        # in the script, and the semicolon was an unfortunate choice.
        # substitute semicolon with .  ( AG: '-' is more legible)
        input_params = { k.replace(':','-') :v for k,v in
            orig_input_params.iteritems() }

        # look for blocked keywords here if present raise
        # add the keyword to the dictionary

        for blockedaiida in self._aiida_blocked_keywords:
            if blockedaiida in input_params:
                raise InputValidationError(
                    "You cannot specify explicitly the '{}' flag in the "
                        "input parameters".format(blockedaiida))
            input_params.update({blockedaiida: self._PREFIX})

        for blockedatm in self._atom_blocked_keywords:
            if blockedatm in input_params:
                raise InputValidationError(
                    "You cannot specify explicitly the '{}' flag in the "
                        "input parameters".format(blockedatm))
            
        input_params.update({'number-of-species': len(structure.kinds)})
        input_params.update({'number-of-atoms': len(structure.sites)})
        #
        # Regarding the lattice-constant parameter:
        # -- The variable "alat" is not typically kept anywhere, and
        # has already been used to define the vectors.
        # We need to specify that the units of these vectors are Ang...
        
        input_params.update({'lattice-constant': '1.0 Ang'})
        
        # Note that this  will break havoc with the band-k-points "pi/a"
        # option. The use of this option should be banned.

        # Note that the implicit coordinate convention of the Structure
        # class corresponds to the "Ang" convention in Siesta.
        # That is why the "atomic-coordinates-format" keyword is blocked
        # and reset.
        input_params.update({'atomic-coordinates-format': 'Ang'})

        # ============== Preparation of input data ===============
        #

        # ------------ CELL_PARAMETERS -----------
        cell_parameters_card = "%block lattice-vectors\n"
        for vector in structure.cell:
            cell_parameters_card += ("{0:18.10f} {1:18.10f} {2:18.10f}"
                                     "\n".format(*vector))
        cell_parameters_card += "%endblock lattice-vectors\n"

        # ------------- ATOMIC_SPECIES ------------
        # I create the subfolder that will contain the pseudopotentials
        tempfolder.get_subfolder(self._PSEUDO_SUBFOLDER, create=True)
        # I create the subfolder with the output data 
        tempfolder.get_subfolder(self._OUTPUT_SUBFOLDER, create=True)
                
        atomic_species_card_list = []

        # Dictionary to get the atomic number of a given element
        datmn=dict([(v['symbol'],k) for k,v in elements.iteritems()])

        spind={}
        spcount=0
        for kind in structure.kinds:

            ps = pseudos[kind.name]

            # I add this pseudo file to the list of files to copy,
            # with the appropiate name
            local_copy_list.append((ps.get_file_abs_path(),
                                    os.path.join(self._PSEUDO_SUBFOLDER,
                                                 kind.name+".psf")))
            spcount+=1
            spind[kind.name]=spcount
            atomic_species_card_list.append("{0:5} {1:5} {2:5}\n".format(
                spind[kind.name], datmn[kind.symbol], kind.name.rjust(6)))

        atomic_species_card_list = (["%block chemicalspecieslabel\n"] + 
            list(atomic_species_card_list))
        atomic_species_card = "".join(atomic_species_card_list)
        atomic_species_card += "%endblock chemicalspecieslabel\n"
        # Free memory
        del atomic_species_card_list

        # ------------ ATOMIC_POSITIONS -----------
        atomic_positions_card_list = [
            "%block atomiccoordinatesandatomicspecies\n"]
        countatm=0
        for site in structure.sites:
            countatm+=1
            atomic_positions_card_list.append(
            "{0:18.10f} {1:18.10f} {2:18.10f} {3:4} {4:6} {5:6}\n".format(
                site.position[0], site.position[1],
                site.position[2], spind[site.kind_name], 
                site.kind_name.rjust(6), countatm))
        atomic_positions_card = "".join(atomic_positions_card_list)
        del atomic_positions_card_list # Free memory
        atomic_positions_card +="%endblock atomiccoordinatesandatomicspecies\n"


        # --------------- K-POINTS ----------------
        if kpoints is not None:
            #
            # Get a mesh for sampling
            # NOTE that there is not yet support for the 'kgrid-cutoff'
            # option in Siesta.
            #
                try:
                    mesh,offset = kpoints.get_kpoints_mesh()
                    has_mesh = True
                except AttributeError:
                    raise InputValidationError("K-point sampling for scf "
                        "must be given in mesh form")
        
                kpoints_card_list = ["%block kgrid_monkhorst_pack\n"]
                #
                # This will fail if has_mesh is False (for the case of a list),
                # since in that case 'offset' is undefined.
                #
                if any( [ (i!=0. and i !=0.5) for i in offset] ):
                    raise InputValidationError("offset list must only be made "
                                               "of 0 or 0.5 floats")
                the_offset = [ 0 if i==0. else 1 for i in offset ]
                the_6_integers = list(mesh) + the_offset
                kpoints_card_list.append(
                "{0:6} {1:6} {2:6} {3:18.10f}\n".format(
                    mesh[0], 0, 0, the_offset[0]))
                kpoints_card_list.append(
                "{0:6} {1:6} {2:6} {3:18.10f}\n".format(
                    0, mesh[1], 0, the_offset[1]))
                kpoints_card_list.append(
                "{0:6} {1:6} {2:6} {3:18.10f}\n".format(
                    0, 0, mesh[2], the_offset[2]))
                
                kpoints_card = "".join(kpoints_card_list)
                kpoints_card += "%endblock kgrid_monkhorst_pack\n"
                del kpoints_card_list
            
 
        # --------------- K-POINTS-FOR-BANDS ----------------!
        #This part is computed only if flagbands=True
        #Two possibility are supported in Siesta: BandLines ad BandPoints
	#At the moment the user can't choose directly one of the two options
        #BandsLine is set automatically if bandskpoints has labels,
        #BandsPoints if bandskpoints has no labels
        #BandLinesScale =pi/a is not supported at the moment because currently 
        #a=1 always. BandLinesScale ReciprocalLatticeVectors is always set
        if flagbands:
            bandskpoints_card_list = ["BandLinesScale ReciprocalLatticeVectors\n"]
	    if bandskpoints.labels==None:
                bandskpoints_card_list.append("%block BandPoints\n")
		for s in bandskpoints.get_kpoints():
                    bandskpoints_card_list.append("{0:8.3f} {1:8.3f} {2:8.3f} \n".format(s[0], s[1], s[2]))
                fbkpoints_card = "".join(bandskpoints_card_list)
                fbkpoints_card += "%endblock BandPoints\n"
            else:
	        bandskpoints_card_list.append("%block BandLines\n")
	        savs=[]
                listforbands = bandskpoints.get_kpoints()
                for s,m in bandskpoints.labels:
    		    savs.append(s)                
    	        rawindex=0
	        for s,m in bandskpoints.labels:
    		    rawindex=rawindex+1
    		    x, y, z = listforbands[s]
    		    if rawindex==1:
		       	bandskpoints_card_list.append("{0:3} {1:8.3f} {2:8.3f} {3:8.3f} {4:1} \n".format(1, x, y, z, m))
    		    else:
                        bandskpoints_card_list.append("{0:3} {1:8.3f} {2:8.3f} {3:8.3f} {4:1} \n".format(s-savs[rawindex-2], x, y, z, m))
                fbkpoints_card = "".join(bandskpoints_card_list)
                fbkpoints_card += "%endblock BandLines\n"
            del bandskpoints_card_list

        # ================ Namelists and cards ===================
        
        input_filename = tempfolder.get_abs_path(self._INPUT_FILE_NAME)

        with open(input_filename,'w') as infile:
            # here print keys and values tp file
            
            for k, v in sorted(input_params.iteritems()):
                infile.write(get_input_data_text(k,v))
                # ,mapping=mapping_species))

            # Basis set info is processed just like the general
            # parameters section. Some discipline is needed to
            # put any basis-related parameters (including blocks)
            # in the basis dictionary in the input script.
            #
            if basis is not None:
              infile.write("#\n# -- Basis Set Info follows\n#\n")
              for k, v in input_basis.iteritems():
                  infile.write(get_input_data_text(k,v))

            # Write previously generated cards now
            infile.write("#\n# -- Structural Info follows\n#\n")
            infile.write(atomic_species_card)
            infile.write(cell_parameters_card)
            infile.write(atomic_positions_card)
            if kpoints is not None:
                infile.write("#\n# -- K-points Info follows\n#\n")
                infile.write(kpoints_card)
	    if flagbands:
		infile.write("#\n# -- Bandlines/Bandpoints Info follows\n#\n")
	        infile.write(fbkpoints_card)

            # Write max wall-clock time
            infile.write("#\n# -- Max wall-clock time block\n#\n")
            infile.write(
              "max.walltime {}".format(self.get_max_wallclock_seconds()))

        # ------------------------------------- END of fdf file creation
        
        # operations for restart
        # copy remote output dir, if specified
        
        if parent_calc_folder is not None:
            remote_copy_list.append(
                    (parent_calc_folder.get_computer().uuid,
                     os.path.join(parent_calc_folder.get_remote_path(),
                                  self._restart_copy_from),
                     self._restart_copy_to
                     ))

        calcinfo = CalcInfo()

        calcinfo.uuid = self.uuid
        #
        # Empty command line by default
        # Why use 'pop' ?
        cmdline_params = settings_dict.pop('CMDLINE', [])

        # Comment this paragraph better, if applicable to Siesta
        #
        #we commented calcinfo.stin_name and added it here in cmdline_params
        #in this way the mpirun ... pw.x ... < aiida.in 
        #is replaced by mpirun ... pw.x ... -in aiida.in
        # in the scheduler, _get_run_line, if cmdline_params is empty, it 
        # simply uses < calcinfo.stin_name
        
        if cmdline_params: 
            calcinfo.cmdline_params = list(cmdline_params)
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        
        calcinfo.stdin_name = self._INPUT_FILE_NAME
        calcinfo.stdout_name = self._OUTPUT_FILE_NAME
        calcinfo.xml_name = self._XML_FILE_NAME
        calcinfo.messages_name = self._MESSAGES_FILE_NAME

        #
        # Code information object
        #
        codeinfo = CodeInfo()
        codeinfo.cmdline_params = list(cmdline_params)
        codeinfo.stdin_name = self._INPUT_FILE_NAME
        codeinfo.stdout_name = self._OUTPUT_FILE_NAME
        codeinfo.xml_name = self._XML_FILE_NAME
        codeinfo.messages_name = self._MESSAGES_FILE_NAME
        codeinfo.code_uuid = code.uuid
        calcinfo.codes_info = [codeinfo]

        # Retrieve by default: the output file, the xml file, and the
        # messages file.
	# If flagbands=True we also add the bands file to the retrieve list!
	# This is extremely important because the parser parses the bands
	# only if aiida.bands is in the retrieve list!!
        
        calcinfo.retrieve_list = []         
        calcinfo.retrieve_list.append(self._OUTPUT_FILE_NAME)
        calcinfo.retrieve_list.append(self._XML_FILE_NAME)
        calcinfo.retrieve_list.append(self._MESSAGES_FILE_NAME)
        if flagbands:
	    calcinfo.retrieve_list.append(self._BANDS_FILE_NAME)

        # Any other files specified in the settings dictionary
        settings_retrieve_list = settings_dict.pop('ADDITIONAL_RETRIEVE_LIST',[])
        calcinfo.retrieve_list += settings_retrieve_list

        return calcinfo

    @classmethod
    def _get_linkname_pseudo_prefix(cls):
        """
        The prefix for the name of the link used for the pseudo for kind 'kind'
        
        :param kind: a string for the atomic kind for which we want
          to get the link name
        """
        return "pseudo_"

    @classmethod
    def _get_linkname_pseudo(cls, kind):
        """
        The name of the link used for the pseudo for kind 'kind'. 
        It appends the pseudo name to the pseudo_prefix, as returned by the
        _get_linkname_pseudo_prefix() method.
        
        :note: if a list of strings is given, the elements are appended
          in the same order, separated by underscores
        
        :param kind: a string (or list of strings) for the atomic kind(s) for 
            which we want to get the link name
        """
        # If it is a list of strings, and not a single string: join them
        # by underscore
        #
        # It might be better to use another character instead of '_'. As it
        # is now, it conflicts with species names of the form Symbol_extra.

        if isinstance(kind, (tuple, list)):
            suffix_string = "_".join(kind) 
        elif isinstance(kind, basestring):
            suffix_string = kind
        else:
            raise TypeError("The parameter 'kind' of _get_linkname_pseudo can "
                            "only be a string or a list of strings")
        return "{}{}".format(cls._get_linkname_pseudo_prefix(),suffix_string)

    def use_pseudos_from_family(self, family_name):
        """
        Set the pseudo to use for all atomic kinds, picking pseudos from the
        family with name family_name.
        
        :note: The structure must already be set.
        
        :param family_name: the name of the group containing the pseudos
        """
        from collections import defaultdict

        try:
          ##  structure = inputdict.pop(self.get_linkname('structure'))
          structure=self.get_inputs_dict()[self.get_linkname('structure')]
        except AttributeError:
            raise ValueError("Structure is not set yet! Therefore, the method "
                             "use_pseudos_from_family cannot automatically set "
                             "the pseudos")

        # A dict {kind_name: pseudo_object}
        kind_pseudo_dict = get_pseudos_from_structure(structure, family_name)
        
        # We have to group the species by pseudo, I use the pseudo PK
        # pseudo_dict will just map PK->pseudo_object
        pseudo_dict = {} 
        # Will contain a list of all species of the pseudo with given PK
        pseudo_species = defaultdict(list) 
        
        for kindname, pseudo in kind_pseudo_dict.iteritems():
            pseudo_dict[pseudo.pk] = pseudo
            pseudo_species[pseudo.pk].append(kindname)
                
        for pseudo_pk in pseudo_dict:
            pseudo = pseudo_dict[pseudo_pk]
            kinds = pseudo_species[pseudo_pk]
            # I set the pseudo for all species, sorting alphabetically
            self.use_pseudo(pseudo, sorted(kinds))

    def _set_parent_remotedata(self,remotedata):
        """
        Used to set a parent remotefolder in the restart of ph.
        """
        from aiida.common.exceptions import ValidationError
        
        if not isinstance(remotedata,RemoteData):
            raise ValueError('remotedata must be a RemoteData')
        
        # complain if another remotedata is already found
        input_remote = self.get_inputs(node_type=RemoteData)
        if input_remote:
            raise ValidationError("Cannot set several parent calculation to a "
                "{} calculation".format(self.__class__.__name__))

        self.use_parent_folder(remotedata)

    def create_restart(self,force_restart=False): 
        """
        Simple Function to restart a calculation that was not completed
        (for example, due to max walltime reached, or lack of convergence)
 
        This version requests that the density-matrix file be copied
        from the old calculation's output folder, and sets an fdf option to
        read it upon start. Other possibilites might be given by extra arguments
        in the future.

        No support for updated structures in variable-geometry runs yet.

        Returns a calculation c2, with all links prepared but not stored in DB.
        To submit it simply:
        c2.store_all()
        c2.submit()
        
        :param bool force_restart: restart also if parent is not in FINISHED 
        state (e.g. FAILED, IMPORTED, etc.). Default=False.
        """
        from aiida.common.datastructures import calc_states
        
        # Check the calculation's state using ``from_attribute=True`` to
        # correctly handle IMPORTED calculations.
        if self.get_state(from_attribute=True) != calc_states.FINISHED:
            if force_restart:
                pass
            else:
                raise InputValidationError("Calculation to be restarted must "
                    "be in the {} state. Otherwise, use the force_restart "
                    "flag".format(calc_states.FINISHED) )
        
        calc_inp = self.get_inputs_dict()

        old_inp_dict = calc_inp['parameters'].get_dict()
        
        # add the restart flag
        # In Siesta, this could be the option to read the DM, if
        # the restart is due to lack of convergence, but in general
        # one would need to pick up the latest structure if doing
        # geometry optimizations ...

        # Warning: we need a way to canonicalize the fdf options...
        # or in this case an "order-preserving" dictionary to put
        # this option at the beginning...
        
        old_inp_dict['dm-use-save-dm'] = True
        inp_dict = ParameterData(dict=old_inp_dict) 
        
        remote_folders = self.get_outputs(type=RemoteData)
        if len(remote_folders)!=1:
            raise InputValidationError("More than one output RemoteData found "
                                       "in calculation {}".format(self.pk))
        remote_folder = remote_folders[0]
        
        c2 = self.copy()
        
        # set the new links

        c2.use_code(calc_inp['code'])
        c2._set_parent_remotedata( remote_folder )

        # New fdf options
        c2.use_parameters(inp_dict)

        # This is too simple to deal with various kinds sharing
        # the same pseudo! As the pseudo is stored in the database,
        # we need to record the associated kinds in the calculation, somehow.
        
        for pseudo in self.get_inputs(node_type=PsfData):
            c2.use_pseudo(pseudo, kind=pseudo.element)

        # To meaningfully use the latest structure from
        # an aborted calculation, we need to implement trajectories
        # For now, test restarts of non-converged single-point calculations
        
        c2.use_structure(calc_inp['structure'])

        # These are optional...
        
        try:
            old_basis = calc_inp['basis']
        except KeyError:
            old_basis = None
        if old_basis is not None:
            c2.use_basis(old_basis)

        try:
            old_kpoints = calc_inp['kpoints']
        except KeyError:
            old_kpoints = None
        if old_kpoints is not None:
            c2.use_kpoints(old_kpoints)

        try:
            old_bandskpoints = calc_inp['bandskpoints']
        except KeyError:
            old_bandskpoints = None
        if old_bandskpoints is not None:
            c2.use_bandskpoints(old_bandskpoints)
            
        try:
            old_settings_dict = calc_inp['settings'].get_dict()
        except KeyError:
            old_settings_dict = {}
            
        if old_settings_dict: # if not empty dictionary
            settings = ParameterData(dict=old_settings_dict)
            c2.use_settings(settings)
        
        return c2

def get_input_data_text(key,val, mapping=None):
    """
    Given a key and a value, return a string (possibly multiline for arrays)
    with the text to be added to the input file.
    
    :param key: the flag name
    :param val: the flag value. If it is an array, a line for each element
            is produced, with variable indexing starting from 1.
            Each value is formatted using the conv_to_fortran function.
    :param mapping: Optional parameter, must be provided if val is a dictionary.
            It maps each key of the 'val' dictionary to the corresponding 
            list index. For instance, if ``key='magn'``, 
            ``val = {'Fe': 0.1, 'O': 0.2}`` and ``mapping = {'Fe': 2, 'O': 1}``,
            this function will return the two lines ``magn(1) = 0.2`` and
            ``magn(2) = 0.1``. This parameter is ignored if 'val' 
            is not a dictionary. 
    """
    from aiida.common.utils import conv_to_fortran
    # I check first the dictionary, because it would also match
    # hasattr(__iter__)
    if isinstance(val, dict):
        if mapping is None:
            raise ValueError("If 'val' is a dictionary, you must provide also "
                             "the 'mapping' parameter")

        list_of_strings = []
        for elemk, itemval in val.iteritems():
            try:
                idx = mapping[elemk]
            except KeyError:
                raise ValueError("Unable to find the key '{}' in the mapping "
                                 "dictionary".format(elemk))
            
            list_of_strings.append((idx,"  {0}({2}) = {1}\n".format(
                key, conv_to_fortran(itemval), idx)))
        
        # I first have to resort, then to remove the index from the first
        # column, finally to join the strings
        list_of_strings = zip(*sorted(list_of_strings))[1]
        return "".join(list_of_strings)                          
    elif hasattr(val,'__iter__'):
        # a list/array/tuple of values
        list_of_strings = [
            "{0}({2})  {1}\n".format(key, conv_to_fortran(itemval), idx+1)
            for idx, itemval in enumerate(val)]
        return "".join(list_of_strings)
    else:
        # single value
        if key[:6] == '%block':
            bname = key.split()[1]
            b1 = "{0}  {1}".format(key, my_conv_to_fortran(val))
            return b1 + "\n%endblock " + bname + "\n"
        else:
            return "{0}  {1}\n".format(key, my_conv_to_fortran(val))

def _lowercase_dict(d, dict_name):
    from collections import Counter
    
    if isinstance(d,dict):
        new_dict = dict((str(k).lower(), v) for k, v in d.iteritems())
        if len(new_dict) != len(d):
            num_items = Counter(str(k).lower() for k in d.keys())
            double_keys = ",".join([k for k, v in num_items if v > 1])
            raise InputValidationError(
                "Inside the dictionary '{}' there are the following keys that "
                "are repeated more than once when compared case-insensitively: "
                "{}."
                "This is not allowed.".format(dict_name, double_keys))
        return new_dict
    else:
        raise TypeError("_lowercase_dict accepts only dictionaries as argument")
    
def _uppercase_dict(d, dict_name):
    from collections import Counter
    
    if isinstance(d,dict):
        new_dict = dict((str(k).upper(), v) for k, v in d.iteritems())
        if len(new_dict) != len(d):
            
            num_items = Counter(str(k).upper() for k in d.keys())
            double_keys = ",".join([k for k, v in num_items if v > 1])
            raise InputValidationError(
                "Inside the dictionary '{}' there are the following keys that "
                "are repeated more than once when compared case-insensitively: "
                "{}."
                "This is not allowed.".format(dict_name, double_keys))
        return new_dict
    else:
        raise TypeError("_lowercase_dict accepts only dictionaries as argument")

def my_conv_to_fortran(val):
    """
    Special version to avoid surrounding strings with extra ' '. Otherwise the
    fdf tokenizer will not split values and units, for example.

    :param val: the value to be read and converted to a Fortran-friendly string.
    """
    # Note that bool should come before integer, because a boolean matches also
    # isinstance(...,int)
    if (isinstance(val, bool)):
        if val:
            val_str = '.true.'
        else:
            val_str = '.false.'
    elif (isinstance(val, (int, long))):
        val_str = "{:d}".format(val)
    elif (isinstance(val, float)):
        val_str = ("{:18.10e}".format(val)).replace('e', 'd')
    elif (isinstance(val, basestring)):
        val_str = "{!s}".format(val)
    else:
        raise ValueError("Invalid value passed, accepts only bools, ints, "
                         "floats and strings")

    return val_str

    
