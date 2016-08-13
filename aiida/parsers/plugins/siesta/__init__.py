# -*- coding: utf-8 -*-
from aiida.parsers.parser import Parser
from aiida.orm.calculation.job.siesta import SiestaCalculation
from aiida.orm.data.parameter import ParameterData

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.5.0"
__contributors__ = "Andrius Merkys, Giovanni Pizzi, Victor Garcia-Suarez, Alberto Garcia"

# -*- coding: utf-8 -*-


# These auxiliary functions should be put in another module...
#
# List of scalar values from CML to be transferred to AiiDA
#
standard_output_list = [ 'siesta:FreeEK', 'siesta:Etot',
                         'siesta:Ebs', 'siesta:E_Fermi']

def get_parsed_xml_doc(xml_path):

     from xml.dom import minidom

     xmldoc = minidom.parse(xml_path)
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

 for m in itemlist:
     if m.attributes.has_key('title'):
          # Get last scf finalization module
          if m.attributes['title'].value == "SCF Finalization":
               scf_final = m
          # Get (the only) finalization module
          if m.attributes['title'].value == "Finalization":
               final = m

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

                        
#----------------------------------------------------------------------
#----------------------------------------------------------------------
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
        from aiida.common.exceptions import ParsingError
        if not isinstance(calc,SiestaCalculation):
            raise ParsingError("Input calc must be a SiestaCalculation")

    def _get_output_nodes(self, output_path, error_path, xml_path):
        """
        Extracts output nodes from the standard output and standard error
        files. (And XML file)
        """
        from aiida.orm.data.array.trajectory import TrajectoryData
        import re

        parser_version = '0.5'
        parser_info = {}
        parser_info['parser_info'] = 'AiiDA Siesta Parser v{}'.format(parser_version)
        parser_info['parser_items'] = ['Metadata','Scalars','End Structure']
        parser_info['parser_warnings'] = []


        result_list = []

        if xml_path is None:
            self.logger.error("Could not find a CML file to parse")
            raise ParsingError("Could not find a CML file to parse")
        
        # We get everything from the CML file

        xmldoc = get_parsed_xml_doc(xml_path)
        
        # This is an example of how we can access input items
        in_struc = self._calc.get_inputs_dict()['structure']
        
        result_dict = get_dict_from_xml_doc(xmldoc)

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
        
        from aiida.orm.data.array import ArrayData
        import numpy
        arraydata = ArrayData()
        arraydata.set_array('forces', numpy.array(forces))
        arraydata.set_array('stress', numpy.array(stress))
        result_list.append((self.get_linkname_outarray(),arraydata))
        
        return result_list

    def parse_with_retrieved(self,retrieved):
        """
        Receives in input a dictionary of retrieved nodes.
        Does all the logic here.
        """
        
        from aiida.common.exceptions import InvalidOperation
        import os

        output_path = None
        error_path  = None
        xml_path  = None
        try:
            output_path, error_path, xml_path = self._fetch_output_files(retrieved)
        except InvalidOperation:
            raise
        except IOError as e:
            self.logger.error(e.message)
            return False, ()

        if output_path is None and error_path is None and xml_path is None:
            self.logger.error("No output files found")
            return False, ()

        return True, self._get_output_nodes(output_path, error_path, xml_path)

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
        error_path  = None
        xml_path  = None

        if self._calc._DEFAULT_OUTPUT_FILE in list_of_files:
            output_path = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_OUTPUT_FILE )
        if self._calc._DEFAULT_XML_FILE in list_of_files:
            xml_path = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_XML_FILE )
        if self._calc._DEFAULT_ERROR_FILE in list_of_files:
            error_path  = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_ERROR_FILE )

        return output_path, error_path, xml_path

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
                       
   
