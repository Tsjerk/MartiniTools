#!/usr/bin/env python

"""
Python implementation for reading Gromacs' TRR file format.
"""

__version__ = 0.1
__author__  = "Tsjerk A. Wassenaar"
__date__    = "02.12.2014"

import struct, sys, numpy


class TrrReadError(Exception):
    def __init__(self,msg):
        self.msg = msg
    def __str__(self):
        return self.msg


class TRRFrame:

    def __init__(self,trr,offset,time,lmb,box,x,v,f):
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

        return self._x.reshape((-1,self._trr.dim))


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




class TRR:
    """
    Class for reading Gromacs' TRR trajectory files
    """

    # First 24 bytes of each frame should be
    #              GMXMAGIC       LINEFEED       TITLELEN      TITLE(12)
    #          |-----1993-----||-----13-----||------12------||----------|
    # _tag = b'\x00\x00\x07\xc9\x00\x00\x00\r\x00\x00\x00\x0cGMX_trn_file'
    # _tagLen = len(_tag)

    def __init__(self,stream,offset=0,dim=3):
        if type(stream) == str:
            self.stream = open(stream,'rb',buffering=0)
        else:
            self.stream = stream


        # Check if the stream is opened correctly
        if not 'b' in self.stream.mode or hasattr(self.stream,"read1"):
            raise TrrReadError("TRR file stream should be opened with unbuffered binary mode.")


        # Check format and read some info about trajectory
        # Read in a chunk to make sure we have a complete title
        header = self.stream.read(1024)
        if header[:8] != b'\x00\x00\x07\xc9\x00\x00\x00\r': # 1993\r
            raise TrrReadError("Invalid magic number. Probably not a TRR file, or corrupted.")


        self.taglen  = 12 + struct.unpack('>l',header[8:12])[0]
        self.tag     = header[:self.taglen]
        stuff        = struct.unpack('>lllllllllllll',header[self.taglen:self.taglen+52])
        box          = stuff[2]
        self.atoms   = stuff[10]
        self.float   = box//(dim*dim) or (self.x or self.v or self.f)//self.block 
        self.dtype   = (self.float == 4 and "f") or (self.float == 8 and "d")
        self.dtypeE  = '>' + self.dtype # With Endianness
        self.hstr    = ">"+13*"l"+2*self.dtype # >lllllllllllllff
        self.hsize   = 52 + 2*self.float
        self.index   = []
        self.nframes = 0
        self.dim     = dim


        # Find the end
        self.stream.seek(0,2)
        self.size = self.stream.tell()

        
        # Wind the trajectory to the offset        
        self.stream.seek(offset,0)
        self.pos = offset
    

    def __del__(self):
        """Close file if iterator is stopped. Allows: ref = TRR().next()"""
        self.stream.close()


    def __iter__(self):
        return self


    def next(self):
        if self.pos >= self.size:
            raise StopIteration

	offset = self.pos

        if self.stream.tell() != self.pos:
            self.stream.seek(self.pos)

        hsize  = self.taglen+self.hsize
        header = self.stream.read(hsize)

        if header[:8] != b'\x00\x00\x07\xc9\x00\x00\x00\r': # 1993\r
            # Not a proper tag. Broken frame?
            raise StopIteration

        stuff      = struct.unpack(self.hstr,header[-self.hsize:])
        bytesize   = self.taglen + 52 + sum(stuff[:10]) + 2*self.float # Size of complete frame

        time,lmb   = stuff[-2:]
        
        # Lengths and positions of arrays
        box        = (self.dim**2, self.pos+hsize+sum(stuff[:2]))
        x          = (stuff[7]//self.float, self.pos+hsize+sum(stuff[:7]))
        v          = (stuff[8]//self.float, self.pos+hsize+sum(stuff[:8]))
        f          = (stuff[9]//self.float, self.pos+hsize+sum(stuff[:9]))

        self.pos += bytesize
        self.index.append(TRRFrame(self,offset=offset,time=time,lmb=lmb,box=box,x=x,v=v,f=f))

        # Go to next frame
        self.stream.seek(self.pos)

        return self.index[-1]



def processTrajectory(trj, offset=0):
    """Process a TRR file"""
    trr      = TRR(trj,offset=offset)
    frames   = 0
    for frame in trr:
        frames += 1
        t, x, v, f = frame.time, frame.x(), frame.v(), frame.f()
        print t, (x != None) and x.shape, (v != None) and v.shape, (f != None) and f.shape
    return frames


def main(argv=None):
    if argv == None:
        argv = sys.argv
    nframes = processTrajectory(argv[1])
    print 
    return 0


if __name__ == "__main__":
    sys.exit(main())
