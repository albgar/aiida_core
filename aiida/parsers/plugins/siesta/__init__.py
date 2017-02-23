# -*- coding: utf-8 -*-
from aiida.parsers.parser import Parser
from aiida.orm.calculation.job.siesta import SiestaCalculation
from aiida.orm.data.parameter import ParameterData

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.7.0"
__contributors__ = "Andrius Merkys, Giovanni Pizzi, Victor Garcia-Suarez, Alberto Garcia, Emanuele Bosoni"

# -*- coding: utf-8 -*-


# These auxiliary functions should be put in another module...
#
# List of scalar values from CML to be transferred to AiiDA
#
standard_output_list = [ 'siesta:FreeEK', 'siesta:Etot',
                         'siesta:Ebs', 'siesta:E_Fermi']

def get_parsed_xml_doc(xml_path):

     from xml.dom import minidom

     try:
          xmldoc = minidom.parse(xml_path)
     except:
          xmldoc = None

     return xmldoc

def get_dict_from_xml_doc(xmldoc):

     # Scalar items
     
     scalar_dict = {}
     
     # Metadata items
     itemlist = xmldoc.getElementsByTagName('metadata')
     for s in itemlist:
         name = s.attributes['name'].value
         value = s.attributes['content'].value
         scalar_dict[name] = value

     # Scalar output items
     # From the last module ("Finalization")
     # wrapped in <property> elements with a <scalar> child
     itemlist = xmldoc.getElementsByTagName('module')
     # 
     finalmodule = itemlist[-1]
     props = finalmodule.getElementsByTagName('property')

     for s in props:
       if s.attributes.has_key('dictRef'):
         name = s.attributes['dictRef'].value
         if name in standard_output_list:
             data = s.getElementsByTagName('scalar')[0]
             value = data.childNodes[0].nodeValue
             units = data.attributes['units'].value
             loc_colon = units.find(':')
             unit_name = units[loc_colon+1:]
             loc_colon = name.find(':')
             reduced_name = name[loc_colon+1:]
             # Put units in separate entries, as in QE
             scalar_dict[reduced_name] = value
             scalar_dict[reduced_name+"_units"] = unit_name

     scalar_dict['variable_geometry'] = is_variable_geometry(xmldoc)
     return scalar_dict


def is_variable_geometry(xmldoc):
     """
     Tries to guess whether the calculation involves changes in
     geometry. It needs at least one 'SCF finalization' item
     """
     
     itemlist = xmldoc.getElementsByTagName('module')
     geom_steps = 0
     for m in itemlist:
          if m.attributes.has_key('title'):
               if m.attributes['title'].value == "SCF Finalization":
                    geom_steps += 1

     if (geom_steps > 1):
          return True
     else:
          return False

def get_last_structure(xmldoc, input_structure):

    from aiida.orm import DataFactory

    # Final structure from the last module ("Finalization")

    itemlist = xmldoc.getElementsByTagName('module')
    # 
    finalmodule = itemlist[-1]
    atoms = finalmodule.getElementsByTagName('atom')
    cellvectors = finalmodule.getElementsByTagName('latticeVector')

    atomlist = []

    for a in atoms:
         kind = a.attributes['elementType'].value
         x = a.attributes['x3'].value
         y = a.attributes['y3'].value
         z = a.attributes['z3'].value
         atomlist.append([kind,[float(x),float(y),float(z)]])
    
    cell = []
    for l in cellvectors:
         data = l.childNodes[0].data.split()
         cell.append([float(s) for s in data])

    # Generally it is better to pass the input structure
    # and reset the data, since site 'names' are not handled by
    # the CML file (at least not in Siesta versions <= 4.0)
    #

    s = input_structure.copy()
    s.reset_cell(cell)
    new_pos = [atom[1] for atom in atomlist]
    s.reset_sites_positions(new_pos)
          
    return s

                        
def get_final_forces_and_stress(xmldoc):
 #
 # Extracts forces and stresses as lists of lists...
 # Where do I put them??
 #
 itemlist = xmldoc.getElementsByTagName('module')

 scf_final = None
 final = None
 for m in itemlist:
     if m.attributes.has_key('title'):
          # Get last scf finalization module
          if m.attributes['title'].value == "SCF Finalization":
               scf_final = m
          # Get (the only) finalization module
          if m.attributes['title'].value == "Finalization":
               final = m

 forces = None
 stress = None

 if scf_final is not None:
      props = scf_final.getElementsByTagName('property')

      for p in props:
        if p.attributes.has_key('dictRef'):
           if p.attributes['dictRef'].value=='siesta:forces':
                mat = p.getElementsByTagName('matrix')[0]
                # Get flat list and reshape as list of lists
                # using info on rows and columns in CML file
                rows = int(mat.attributes['rows'].value)
                cols = int(mat.attributes['columns'].value)
                f = mat.childNodes[0].data.split()
                f = [float(x) for x in f]
                forces = [ f[rows*i : rows*(i+1)] for i in range(cols)]

 if final is not None:
      props = final.getElementsByTagName('property')
      for p in props:
        if p.attributes.has_key('dictRef'):
           if p.attributes['dictRef'].value=='siesta:stress':
                mat = p.getElementsByTagName('matrix')[0]
                # Get flat list and reshape as list of lists
                # using info on rows and columns in CML file
                rows = int(mat.attributes['rows'].value)
                cols = int(mat.attributes['columns'].value)
                s = mat.childNodes[0].data.split()
                s = [float(x) for x in s]
                stress = [ s[rows*i : rows*(i+1)] for i in range(cols)]
          
 return forces, stress


#The parsing is different whether I have Bands or Points.
#I recognise this two situations looking at bandskpoints.label
#(like I did in the plugin)
def get_bands(self, bands_path):
    import numpy as np
    from aiida.common.exceptions import InputValidationError
    from aiida.common.exceptions import ValidationError
    tottx = []
    f=open(bands_path)
    tottx=f.read().split()
    ef = float(tottx[0])
    if self._calc.inp.bandskpoints.labels==None:
        minfreq,maxfreq = float(tottx[1]), float(tottx[2])
        nbands, nspins, nkpoints = int(tottx[3]), int(tottx[4]), int(tottx[5])
	spinup=np.zeros((nkpoints,nbands))
	spindown=np.zeros((nkpoints,nbands))
        for i in range(nkpoints):
            for j in range(nbands):
                spinup[i,j]=(float(tottx[i*(nbands*2+3)+6+j+3]))
                if (nspins==2):
                    spindown[i,j]=(float(tottx[i*(nbands*2+3)+6+j+3+nbands]))   
    else:
        mink,maxk = float(tottx[1]), float(tottx[2])
        minfreq,maxfreq = float(tottx[3]), float(tottx[4])
        nbands, nspins, nkpoints = int(tottx[5]), int(tottx[6]), int(tottx[7])
        spinup=np.zeros((nkpoints,nbands))
        spindown=np.zeros((nkpoints,nbands))
        for i in range(nkpoints):
            for j in range(nbands):
                spinup[i,j]=(float(tottx[i*(nbands*2+1)+8+j+1]))
       	        if (nspins==2):
            	    spindown[i,j]=(float(tottx[i*(nbands*2+1)+8+j+1+nbands]))
    if (nspins==2):
	bands = (spinup, spindown)
    elif (nspins==1):
	bands = spinup
    else:
	raise NotImplementedError('if nspins=4 could be a non collinear calculation: not implemented yet')
    return bands                     

def get_warnings_from_file(messages_path):
     """
     Generates a list of warnings from the 'MESSAGES' file, which
     contains a line per message, prefixed with 'INFO',
     'WARNING' or 'FATAL'.

     :param messages_path: 

     Returns a boolean indicating success (True) or failure (False)
     and a list of strings.
     """
     f=open(messages_path)
     lines=f.read().split('\n')   # There will be a final '' element

     import re
     
     # Search for 'FATAL:' messages and return immediately
     for line in lines:
          if re.match('^FATAL:.*$',line):
               return False, lines[:-1]  # Remove last (empty) element

     # Make sure that the job did finish (and was not interrupted
     # externally)

     if lines[-2] != 'INFO: Job completed':
          lines[-1] = 'FATAL: Job did not finish'
          return False, lines

     # (Insert any other "non-success" conditions before next section)
     # (e.g.: be very picky about (some) 'WARNING:' messages)
     
     # Return with success flag
     
     return True, lines[:-1]  # Remove last (empty) element

#----------------------------------------------------------------------
#----------------------------------------------------------------------

from aiida.parsers.exceptions import OutputParsingError

class SiestaOutputParsingError(OutputParsingError):
     pass
#---------------------------

class SiestaParser(Parser):
    """
    Parser for the output of siesta.
    """
    def __init__(self,calc):
        """
        Initialize the instance of SiestaParser
        """
        # check for valid input
        self._check_calc_compatibility(calc)
        super(SiestaParser, self).__init__(calc)

    def _check_calc_compatibility(self,calc):
        if not isinstance(calc,SiestaCalculation):
            raise SiestaOutputParsingError("Input calc must be a SiestaCalculation")

    def _get_output_nodes(self, output_path, messages_path, xml_path, bands_path):
        """
        Extracts output nodes from the standard output and standard error
        files. (And XML file)
        """
        from aiida.orm.data.array.trajectory import TrajectoryData
        import re

        parser_version = '0.6.0-warnings'
        parser_info = {}
        parser_info['parser_info'] = 'AiiDA Siesta Parser v{}'.format(parser_version)
        parser_info['parser_items'] = ['Metadata','Scalars','End Structure']
        parser_info['parser_warnings'] = []


        result_list = []

        if xml_path is None:
            self.logger.error("Could not find a CML file to parse")
            raise SiestaOutputParsingError("Could not find a CML file to parse")
        
        # We get everything from the CML file

        xmldoc = get_parsed_xml_doc(xml_path)
        if xmldoc is None:
            self.logger.error("Malformed CML file: cannot parse")
            raise SiestaOutputParsingError("Malformed CML file: cannot parse")
        
        # These are examples of how we can access input items
        #
        # Structure (mandatory)
        #
        in_struc = self._calc.get_inputs_dict()['structure']
        #
        # Settings (optional)
        #
        try:
             in_settings = self._calc.get_inputs_dict()['settings']
        except KeyError:
             in_settings = None

        result_dict = get_dict_from_xml_doc(xmldoc)

        # Add warnings
        successful = True
        if messages_path is None:
             # Perhaps using an old version of Siesta
             warnings_list = ['WARNING: No MESSAGES file...']
        else:
             successful, warnings_list = get_warnings_from_file(messages_path)

        result_dict["warnings"] = warnings_list
        
        # Add parser info dictionary
        parsed_dict = dict(result_dict.items() + parser_info.items())

        output_data = ParameterData(dict=parsed_dict)
        
        link_name = self.get_linkname_outparams()
        result_list.append((link_name,output_data))

        # If the structure has changed, save it
        if is_variable_geometry(xmldoc):
             # Get the input structure to copy its site names,
             # as the CML file traditionally contained only the
             # atomic symbols.
             #
             struc = get_last_structure(xmldoc,in_struc)
             result_list.append((self.get_linkname_outstructure(),struc))

        # Save forces and stress in an ArrayData object
        forces, stress = get_final_forces_and_stress(xmldoc)

        if forces is not None and stress is not None:
             from aiida.orm.data.array import ArrayData
             import numpy
             arraydata = ArrayData()
             arraydata.set_array('forces', numpy.array(forces))
             arraydata.set_array('stress', numpy.array(stress))
             result_list.append((self.get_linkname_outarray(),arraydata))

        # Parse band-structure information if available
        if bands_path is not None:
	    bands = get_bands(self,bands_path)
	    from aiida.orm.data.array.bands import BandsData
	    arraybands = BandsData()
            arraybands.set_kpoints(self._calc.inp.bandskpoints.get_kpoints(cartesian=True))
	    arraybands.set_bands(bands,units="eV")
	    result_list.append((self.get_linkname_bandsarray(),arraybands))
	
        return successful, result_list

    def parse_with_retrieved(self,retrieved):
        """
        Receives in input a dictionary of retrieved nodes.
        Does all the logic here.
        """
        
        from aiida.common.exceptions import InvalidOperation
        import os

        output_path = None
        messages_path  = None
        xml_path  = None
        bands_path = None
        try:
            output_path, messages_path, xml_path, bands_path = self._fetch_output_files(retrieved)
        except InvalidOperation:
            raise
        except IOError as e:
            self.logger.error(e.message)
            return False, ()

        if output_path is None and messages_path is None and xml_path is None:
            self.logger.error("No output files found")
            return False, ()

        successful, out_nodes = self._get_output_nodes(output_path,
                                           messages_path,
                                           xml_path,
                                           bands_path)
        
        return successful, out_nodes

    def _fetch_output_files(self, retrieved):
        """
        Checks the output folder for standard output and standard error
        files, returns their absolute paths on success.

        :param retrieved: A dictionary of retrieved nodes, as obtained from the
          parser.
        """
        from aiida.common.datastructures import calc_states
        from aiida.common.exceptions import InvalidOperation
        import os

        # check in order not to overwrite anything
#         state = self._calc.get_state()
#         if state != calc_states.PARSING:
#             raise InvalidOperation("Calculation not in {} state"
#                                    .format(calc_states.PARSING) )

        # Check that the retrieved folder is there
        try:
            out_folder = retrieved[self._calc._get_linkname_retrieved()]
        except KeyError:
            raise IOError("No retrieved folder found")

        list_of_files = out_folder.get_folder_list()

        output_path = None
        messages_path  = None
        xml_path  = None
	bands_path = None

        if self._calc._DEFAULT_OUTPUT_FILE in list_of_files:
            output_path = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_OUTPUT_FILE )
        if self._calc._DEFAULT_XML_FILE in list_of_files:
            xml_path = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_XML_FILE )
        if self._calc._DEFAULT_MESSAGES_FILE in list_of_files:
            messages_path  = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_MESSAGES_FILE )
        if self._calc._DEFAULT_BANDS_FILE in list_of_files:
            bands_path  = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_BANDS_FILE )

        return output_path, messages_path, xml_path, bands_path

    def get_linkname_outstructure(self):
        """
        Returns the name of the link to the output_structure
        Node exists if positions or cell changed.
        """
        return 'output_structure'

    def get_linkname_outarray(self):
        """                                                                     
        Returns the name of the link to the output_array                        
        In QE, Node may exist in case of calculation='scf'                             
        In Siesta, Node exists to hold the final forces and stress,
        pending the implementation of trajectory data.
        """
        return 'output_array'

    def get_linkname_bandsarray(self):
        """                                                                     
        Returns the name of the link to the bands_array                        
        In Siesta, Node exists to hold the bands,
        pending the implementation of trajectory data.
        """
        return 'bands_array'
                       
   
