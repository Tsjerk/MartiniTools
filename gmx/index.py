import collections

class Index(dict):
    """Class for reading a Gromacs index file and allow referencing by group name or number"""

    def __init__(self,*args):
        for item in args:
            if item.lower().endswith('.ndx') or item.lower().endswith(".ndx.gz"):
                self.from_file(item)
            elif item.lower().endswith('.select'):
                with open(item) as f:
                    for i in f:
                        self.fromExpression(i)
            else:
                self.fromExpression(i,ref)


    def __str__(self):
        return "\n".join(["%s (%d)"%(tag,len(atoms)) for tag,atoms in self.items()])


    def __getattr__(self,attr):
        attr = attr.lower()
        if attr in self.keys():
            return self[attr]
        else:
            for k in self.keys():
                if k.startswith(attr):
                    return self[k]
        raise IndexError("Index group not found: {}\nGroups:\n{}".format(attr, self))


    def from_file(self,filename):
        """Read index groups from GROMACS style index file"""

        # The GROMACS index files contains records with a [ header ]
        # followed by whitespace separated numbers.

        stuff = open(filename).read()

        #          name           numbers                                                           records
        data = [ (a.strip(),tuple(int(u)-1 for u in b.split())) for a,b in [i.split(']') for i in stuff.split('[')[1:]]]

        # Remove disallowed characters. Make everything lowercase
        for name, numbers in data:
            name = name.translate("{0:_>58}_{1:_^38}{1:_<159}".format("0123456789","abcdefghijklmnopqrstuvwxyz"))
            self[name] = numbers

    
    def bind(atoms):
        """Bind a structure to an index to allow selections based on atoms/residues/chains"""
        pass
