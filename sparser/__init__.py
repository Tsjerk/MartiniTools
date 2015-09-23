
"""

    Simple PARSER for command line options.
    ^      ^^^^^^

    (c)2014 T.A. Wassenaar

"""

import __main__ as main


class Usage(Exception):
    def __init__(self,msg,program=getattr(main,"__file__",None)):
        if program:
            self.msg = "%s: %s\nTry '%s --help' for more information."%(program,msg,program)
        else:
            self.msg = "Failed to parse options: %s"%msg
    def __str__(self):
        return self.msg


class Options:
    """
    A simple option parser. 
    All options are registered in the list __options as tuples consisting of:

        OPTION      - Option name.
        ATTRIBUTE   - Attribute on this class through which the option is available.
        TYPE        - Type of the value(s).
        NUMBER      - Number of arguments for option.
        DEFAULT     - Default value. 'None' if no default.
        MULTI       - True/False. If True, allow multiple instances (calls) of the option.
        DESCRIPTION - Description of the option.

    """

    def __init__(self, options, args=None):
        """Set up options object."""

        self.options = options

        # Make a dictionary from the option list
        self._optiondict = dict([(i[0],i[1:]) for i in options if not type(i) == str])

        # Set the options as attributes of this object
        for opt in self._optiondict.values():
            setattr(self,opt[0],([] if not opt[3] else [opt[3]]) if opt[4] else opt[3])

        # Parse the arguments, if given
        if args:
            self.parse(args)


    def __str__(self):
        """Make a string from the option list"""
        out = [main.__file__+"\n"]

        for thing in self.options:
            if type(thing) == str:
                out.append("     "+thing)
            else:
                out.append("     %10s   %s ( %s )"%(thing[0],thing[-1],str(getattr(self,thing[1]))))
            
        return "\n".join(out)+"\n"


    def parse(self,args):
        """Parse the (command-line) arguments."""

        self.help = False

        while args:
            opt = args.pop(0)
                        
            if opt in ("--help","-h"):
                self.help = True
                continue
            
            if not opt in self._optiondict:
                raise Usage("Unrecognized option '%s'"%opt)

            attr, typ, num, default, multi, description = self._optiondict[opt]

            if num > len(args):
                raise Usage("Option '%s' requires %d arguments"%(opt,num))

            if num:
                a = args.pop(0)
                try:
                    val = [typ(a) for i in range(num)]
                except ValueError:
                    raise Usage("Invalid argument to option '%s': %s"%(opt,repr(a)))
            else:
                # Boolean flag
                val = [True]

            if typ == bool:
                # A boolean option is simply set to True if given
                setattr(self,attr,True)
            elif multi:
                # A multi-option adds an item or a tuple to the list
                getattr(self,attr).append(val[0] if num == 1 else tuple(val))
            else:
                # Other options just set item or tuple
                setattr(self,attr,val[0] if num == 1 else tuple(val))

        return 
