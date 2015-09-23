
"""
Handle Gromacs' TRR file format - read in stuff as numpy arrays

(c)2014 Tsjerk A. Wassenaar
"""


import struct, sys, numpy


from frame import TRRFrame


class TrrReadError(Exception):
    def __init__(self,msg):
        self.msg = msg
    def __str__(self):
        return self.msg


class TRR:
    """Class for reading Gromacs' TRR trajectory files"""

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
    

    def __len__(self):
        if not self.nframes:
            start = self.pos
            for i in self:
                self.nframes += 1
            self.stream.seek(start,0)
            self.pos = start
            self.index = []
        return self.nframes


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
        self.index.append(TRRFrame(self,nr=len(self.index),offset=offset,time=time,lmb=lmb,box=box,x=x,v=v,f=f))

        # Go to next frame
        self.stream.seek(self.pos)

        return self.index[-1]


    def close(self):
        self.stream.close()
