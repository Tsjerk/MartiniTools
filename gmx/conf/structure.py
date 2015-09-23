
import numpy

class Structure(numpy.ndarray):

    def __new__(cls, title=None, atoms=None, coord=None, q=None, b=None, m=1, box=None, center=False):

        if coord == None:
            stuff = zip(*atoms)
            names = stuff[1]
            atoms = numpy.asarray(zip(*stuff[:4]))
            coord = zip(*stuff[4:7])
            if not q and len(stuff) > 7:
                q = numpy.array(stuff[7])
            if not b and len(stuff) > 8:
                b = numpy.array(stuff[8])
            
        obj = numpy.asarray(coord).view(cls)

        obj.title = title
        obj.names = names
        obj.atoms = atoms
        obj.q     = q
        obj.b     = b
        obj.m     = m
        obj._box  = box
        
        if center:
            obj -= obj.mean()
            obj.centered = True
        else:
            obj.centered = False

        return obj


    def __array_finalize__(self,obj):
        if obj == None:
            return

        self.title    = getattr(obj,"title", None)
        self.atoms    = getattr(obj,"atoms", None)
        self.q        = getattr(obj,"q",     None)
        self.b        = getattr(obj,"b",     None)
        self.m        = getattr(obj,"m",     None)
        self._box     = getattr(obj,"_box",  None)
        self.centered = False


    def __array_wrap__(self,out,context):
        return numpy.ndarray.__array_wrap__(self,out,context)


    def __getslice__(self,i,j):
        return self.__getitem__(slice(i,j,None))


    def __getitem__(self,items):
        # Dealing with 'internal' calls.
        if numpy.rank(self) == 1:
            return numpy.ndarray.__getitem__(self,items)

        if type(items) in (int,slice):
            atoms, items = items, (items,slice(None,None,None))
        else:
            atoms = list(items[0])
            
        obj = numpy.ndarray.__getitem__(self,items)

        obj.title = self.title
        obj._box  = self._box

        obj.atoms = self.atoms[atoms,:]
        obj.q     = self.q[atoms] if self.q != None else None
        obj.b     = self.b[atoms] if self.b != None else None
        obj.m     = self.m[atoms] if self.m != 1    else 1

        obj.centered = False

        return obj

    
    def __iadd__(self,other):
        self.centered = False
        return numpy.ndarray.__iadd__(self,other)


    def x(self):
        # Return x similar to method for TRRFrames
        return self


    def box(self):
        # Return box similar to method for TRRFrames
        return self._box


    def mean(self,axis=0):
        return numpy.ndarray.mean(self,axis=axis)


    def radiusOfGyrationSquared(self):
        """Fetch the radius of gyration or determine and store the value and return it"""

        g = getattr(self,"_radiusOfGyrationSquared",None)

        if not g:
            if self.centered:
                g = numpy.sum(self**2)
            else:
                g = numpy.sum((self-self.mean())**2)

            self._radiusOfGyrationSquared = g

        return g
        

    def cov(self):
        """3x3 covariance matrix"""

        if self.centered:
            return numpy.dot(self.T,self).view(numpy.ndarray)
        else:
            s = self - self.mean()
            return numpy.dot(s.T,s).view(numpy.ndarray)


    def principal(self):
        """Principal axes and radii"""
        return numpy.linalg.eig(self.cov())


    def msd(self,target,align=True):
        "Mean square displacement and fitting"

        XX = self.radiusOfGyrationSquared
        YY = target.radiusOfGyrationSquared
        
        if self.centered or target.centered:
            XY = numpy.dot(self.T,target)
        else:
            XY = numpy.dot(self.T,(target-target.mean()))

        U,l,V = numpy.linalg.svd(XY)
        
