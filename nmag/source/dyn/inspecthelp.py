
#-- stolen the next two functios from Ipython to quickly get
#the signature of a function. Surely this can be done neater
#(fangohr 20/04/2007)

import inspect,types

def __getargspec(obj):
    """Get the names and default values of a function's arguments.
    
    A tuple of four things is returned: (args, varargs, varkw, defaults).
    'args' is a list of the argument names (it may contain nested lists).
    'varargs' and 'varkw' are the names of the * and ** arguments or None.
    'defaults' is an n-tuple of the default values of the last n arguments.
    
    Modified version of inspect.getargspec from the Python Standard
    Library."""
    
    if inspect.isfunction(obj):
        func_obj = obj
    elif inspect.ismethod(obj):
        func_obj = obj.im_func
    else:
        raise TypeError, 'arg is not a Python function'
    args, varargs, varkw = inspect.getargs(func_obj.func_code)
    return args, varargs, varkw, func_obj.func_defaults

def __getdef(obj,oname=''):
    """Return the definition header for any callable object.
    
    If any exception is generated, None is returned instead and the
    exception is suppressed."""
    
    try:
        return oname + inspect.formatargspec(*__getargspec(obj))
    except:
        return None


def pdef(obj,oname=''):
    """Print the definition header for any callable object.
    
    If the object is a class, print the constructor information."""
    
    
    header = 'Arguments:'
    if type(obj) is types.ClassType:
        header = 'Class constructor information:'
        obj = obj.__init__
    elif type(obj) is types.InstanceType:
        obj = obj.__call__
    elif type(obj) is types.TypeType: #new class object
        header = 'Class constructor information:'
        obj = obj.__init__
    elif type(obj) is type(property()): #couldn't find types.Property type
        header = 'Property information'
        return header,None

    if not callable(obj):
        print 'Object is not callable.'
        return None,None

        
    output = __getdef(obj,oname)
    if output is None:
        return None,None
        #self.noinfo('definition header',oname)
        #print 'definition header',oname
    else:
        #print "header = ",header
        #print output
        return header,output
        #print >>Term.cout, header,self.format(output),

#    def pdoc(self,obj,oname='',formatter = None):
#        """Print the docstring for any object.
#
#        Optional:
#        -formatter: a function to run the docstring through for specially
#        formatted docstrings."""
#        
#        head = self.__head  # so that itpl can find it even if private
#        ds = getdoc(obj)
#        if formatter:
#            ds = formatter(ds)
#        if type(obj) is types.ClassType:
#            init_ds = getdoc(obj.__init__)
#            output = itpl('$head("Class Docstring:")\n'
#                          '$indent(ds)\n'
#                          '$head("Constructor Docstring"):\n'
#                          '$indent(init_ds)')
#        elif type(obj) is types.InstanceType and hasattr(obj,'__call__'):
#            call_ds = getdoc(obj.__call__)
#            if call_ds:
#                output = itpl('$head("Class Docstring:")\n$indent(ds)\n'
#                              '$head("Calling Docstring:")\n$indent(call_ds)')
#            else:
#                output = ds
#        else:
#            output = ds
#        if output is None:
#            self.noinfo('documentation',oname)
#            return
#        page(output)
