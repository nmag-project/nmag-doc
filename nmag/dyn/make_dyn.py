#The next two items allow us to document less
#then the full functionality while containing
#the full documenation when creating the code
#documentation for developers using epydoc.

#anything folling this token in docstrings is not
#included in documentation produced with this tool.
developerinfo_token = "#dev-only#"

#if an object has an attribute __argdocstring__ then
#this (rather than introspection) will be used to
#distplay the function's (or constructor's or the method's)
#argument signature.
argumentdocstring_name = "__argdoc__"

argumentdocstringlong_name = "__argdoclong__"


import logging
logging.getLogger('').setLevel(logging.DEBUG)

from inspecthelp import pdef

def remove_left_space(list_of_strings):
    """
    Given a list of strings, this will remove any leading space on the
    left.

    Algorithm:

    - if there is only one line, remove leading spaces and return
      line as is

    - if the first line is empty, remove it (it should be empty
      for docstring following the epydoc convention)

    - find the lowest number of spaces for all lines starting
      from the second line and remove this number of spaces from
      all lines (with the exception of the first).

    :Parameters:
      ``list_of_strings : list of strings
         The input data to be processed

    :Returns:
      list of strings
        The processed list of strings
    """

    lines = list_of_strings #just shorthand variable name
    if len(lines) == 0:
        return lines

    elif len(lines) == 1:
        return [lines[0].lstrip()]

    # Otherwise, we have 2 or more lines. Find shortest space for all these:
    min_space = min(len(line) - len(line.lstrip(" "))
                    for line in lines[1:]
                    if len(line.strip()) > 0)

    # And remove from all lines
    newlines = map(lambda line: line[min_space:], lines)
    if len(lines[0].strip()) > 0:
        newlines.insert(0, lines[0])

    # Check there are no tabs
    for line in newlines:
        if line.startswith("\t"):
            raise ValueError("There is a tab in this line: '%s'. " % line +
                             "I cowardly refuse to delete this and suggest "
                             "that you replace it with spaces in the editor.")
    return newlines
        
        


def dyn_doc(outputfilename, things_to_doc):
    """

    :Parameters:
      `outputfilename` : string
         name of file to write to (existing file will be deleted)

      `modulename` : string
         name of module which contains the objects to be documented

      `moduleobjectstylelist` : list of 2-tuples = (string,char)
         list of 2-tuples, one for each object to document.
         The string carries the name of the object to document.
         The char carries the character to be used to tag the
         ((sub)sub)section for the documentation of this object.
         Typical char values are '-' or '~' or '=' etc.

    """
    logging.debug("About to Open '%s'" % outputfilename)

    f = open(outputfilename, 'w')
    for module_name, module_objsty_list in things_to_doc:
        dyn_doc_module(f, module_name, module_objsty_list)
    f.close()

def dyn_doc_module(f, modulename, moduleobjectstylelist):
    #import the module
    logging.debug("About to import module '%s'" % modulename)
    exec('import '+modulename)
    
    for objectname,headerstyle in moduleobjectstylelist:
        print "About to document '%s'" % objectname
        docstring = eval(modulename+'.'+objectname+'.__doc__')

        logging.debug('test')

        if docstring:
            logging.debug("Doc string has %d characters, starts with \n  %s" %
                      (len(docstring),docstring[:70]))
        else:
            logging.debug("-> no docstring found for '%s' " % objectname)
            docstring = "<undocumented XXX>\n\n"


        lines = docstring.split('\n')

        #remove everything after the developer token
        for i,line in enumerate(lines):
            if developerinfo_token in line:
                docstring = '\n'.join(lines[:i])

        objectrootname = objectname.split('.')[-1]

        headerline = len(objectrootname)*headerstyle

        f.write('.. _%s:\n\n' % objectrootname)
        f.write('%s\n' % objectrootname)
        f.write('%s\n' % headerline)
        f.write('\n')
        f.write('Module:\n\t``%s``\n\n' % modulename)
        f.write('Object:\n\t``%s``\n' % objectname)
        header,arguments = pdef( eval(modulename+'.'+objectname) )

        #Somtimes we need to introduce linebreaking into long
        #argument lists. Not nice but we need to be able to
        #read the documentation. (fangohr 10/10/2007 00:29)
        long_argument = False
        
        #see if __argdocstring__ is provided
        try:
            arguments = eval(modulename+'.'+objectname+'.'+argumentdocstring_name )
            print "Have found %s for '%s'" % (argumentdocstring_name,objectname)
        except AttributeError:
            pass

        try:
            arguments = eval(modulename+'.'+objectname+'.'+argumentdocstringlong_name )
            print "Have found %s for '%s'" % (argumentdocstringlong_name,objectname)
            long_argument = True
        except AttributeError:
            pass
        
        if arguments == None:
            #f.write("\n\n<Couldn't obtain signature>\n\n")
            print "Couldn't obtain signature for '%s' (okay if property)" % \
                  (modulename+'.'+objectname )
            f.write('%s\n\t``%s``\n' % (header,arguments))
        else:
            if long_argument == False:
                f.write('%s\n\t``%s``\n' % (header,arguments))
            else:
                f.write('%s\n\t::\n\n' % header)
                for line in arguments.split('\n'):
                    f.write('\t  %s\n' % line)
                f.write('\n')
        f.write('\n')


        #format dostring nicely (remove leading spaces and
        #leading empty lines):
        docstring = "\n".join(remove_left_space( docstring.split('\n') ))
        
        f.write('%s' % docstring)
        f.write('\n\n')

    exec("del "+modulename)


l1 = '-'
l2 = '~'

import nmag
from nmag.hlib import HMatrixSetup
nmag.HMatrixSetup = HMatrixSetup

if True:
    dyn_doc('nmag.txt',[('nmag',[('MagMaterial',l1),
                                 ('uniaxial_anisotropy',l2),
                                 ('cubic_anisotropy',l2),
                                 ('Simulation',l1),
                                 ('Simulation.advance_time',l2),
                                 ('Simulation.get_subfield',l2),
                                 ('Simulation.get_subfield_positions',l2),
                                 ('Simulation.get_subfield_sites',l2),
                                 ('Simulation.get_subfield_average',l2),
                                 ('Simulation.get_subfield_average_siv',l2),
                                 ('Simulation.probe_subfield',l2),
                                 ('Simulation.probe_subfield_siv',l2),
                                 ('Simulation.probe_H_demag_siv',l2),
                                 ('Simulation.hysteresis',l2),
                                 ('Simulation.load_mesh',l2),
                                 ('Simulation.load_m_from_h5file',l2),
                                 ('Simulation.save_restart_file',l2),
                                 ('Simulation.relax',l2),
                                 ('Simulation.save_data',l2),
                                 ('Simulation.set_m',l2),
                                 ('Simulation.set_H_ext',l2),
                                 ('Simulation.set_pinning',l2),
                                 ('Simulation.set_params',l2),
                                 ('get_subfield_from_h5file',l1),
                                 ('get_subfield_positions_from_h5file',l1),
                                 ('get_subfield_sites_from_h5file',l1),
                                 ('HMatrixSetup',l1),
                                 ('SI',l1),
                                 ('SI.value',l2),
                                 ('SI.units',l2),
                                 ('SI.in_units_of',l2),
                                 ('SI.is_compatible_with',l2),
                                 ('ipython',l1)
                                 ]),
                       ])

