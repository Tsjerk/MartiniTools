import numpy

class CoordSet(numpy.ndarray):
    def __new__(cls, coord, time, box):
        obj      = numpy.asarray(coord).view(cls)
        obj.box  = numpy.asarray(box)
        obj.time = time
        return obj


class TRRFrame:

    def __init__(self,trr,nr,offset,time,lmb,box,x,v,f):
        self.nr     = 0
        self.offset = offset
        self.time   = time
        self.lmb    = lmb
        self._idx   = (box,x,v,f)
        self._trr   = trr
        self._box   = None
        self._x     = None
        self._v     = None
        self._f     = None


    def box(self):
        """Return box record from frame. Data is read on first request."""

        if not self._idx[0][0]:
            return None

        if self._box == None:
            if not self._trr.stream.tell() == self._idx[0][1]:
                self._trr.stream.seek(self._idx[0][1])
            self._box = numpy.fromfile(self._trr.stream,dtype=self._trr.dtypeE,count=self._idx[0][0]).reshape((3,3))

        return self._box


    def x(self,index=None):
        """Return coordinates from frame. Data is read on first request."""

        if not self._idx[1][0]:
            return None

        if self._x == None:
            if not self._trr.stream.tell() == self._idx[1][1]:
                self._trr.stream.seek(self._idx[1][1])
            self._x = numpy.fromfile(self._trr.stream,dtype=self._trr.dtypeE,count=self._idx[1][0])

        return CoordSet(self._x.reshape((-1,self._trr.dim)),time=self.time,box=self.box())


    def v(self,index=None):
        """Return velocities from frame. Data is read on first request."""

        if not self._idx[2][0]:
            return None

        if self._v == None:
            if not self._trr.stream.tell() == self._idx[2][1]:
                self._trr.stream.seek(self._idx[2][1])
            self._v = numpy.fromfile(self._trr.stream,dtype=self._trr.dtypeE,count=self._idx[2][0])

        return self._v.reshape((-1,self._trr.dim))


    def f(self,index=None):
        """Return forces from frame. Data is read on first request."""

        if not self._idx[3][0]:
            return None

        if self._f == None:
            if not self._trr.stream.tell() == self._idx[3][1]:
                self._trr.stream.seek(self._idx[3][1])
            self._f = numpy.fromfile(self._trr.stream,dtype=self._trr.dtypeE,count=self._idx[3][0])

        return self._f.reshape((-1,self._trr.dim))


    def clear(self):
        """Clear the box/coordinate/velocity/force data."""
        self._box = None
        self._x   = None
        self._v   = None
        self._f   = None


