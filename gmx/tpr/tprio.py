
import versions
import io
import struct
import numpy


class Integer(int):
    """Integer from a binary stream, with index"""
    def __new__(cls, stream):
        pos = stream.tell() 
        obj = int.__new__(cls, struct.unpack(">l", stream.read(4))[0])
        obj.position = pos
        return obj

    def pack(self,value):
        return struct.pack(">l", int(value))


class Unsigned(int):
    """Unsigned integer from a binary stream, with index"""
    def __new__(cls, stream):
        pos = stream.tell() 
        obj = int.__new__(cls, struct.unpack(">I", stream.read(4))[0])
        obj.position = pos
        return obj

    def pack(self,value):
        return struct.pack(">I", int(value))


class Long(int):
    """Long integer from a binary stream, with index"""
    def __new__(cls, stream):
        pos = stream.tell() 
        obj = int.__new__(cls, struct.unpack(">Q", stream.read(8))[0])
        obj.position = pos
        return obj

    def pack(self,value):
        return struct.pack(">Q", int(value))


class Char(str):
    """Character from a binary stream, with index"""
    def __new__(cls, stream):
        pos = stream.tell()
        obj = str.__new__(cls, chr(struct.unpack(">l", stream.read(4))[0]))
        obj.position = pos
        return obj

    def pack(self,value):
        return struct.pack(">l", ord(value))


class Float(float):
    """Float from a binary stream, with index"""
    def __new__(cls, stream, precision=4):
        pos = stream.tell()         
        if precision == 4:
            obj = float.__new__(cls, struct.unpack(">f", stream.read(4))[0])
        else:
            obj = float.__new__(cls, struct.unpack(">d", stream.read(8))[0])
        obj.position  = pos
        obj.precision = precision
        return obj

    def pack(self,value):
        if self.precision == 4:
            return struct.pack(">f", float(value))
        else:
            return struct.pack(">d", float(value))


class Double(Float):
    """Float from a binary stream, with index"""
    def __new__(cls, stream):
        obj = Float.__new__(cls, stream, precision="double")
        return obj


class Real(Float):
    """
    Float or double from a binary stream, with index
    
    This requires the stream to have an attribute 'precision', which
    determines whether floating points are encoded in single (4) or
    in double (!4) precision.
    """
    def __new__(cls, stream):
        return Float.__new__(cls, stream, precision=stream.precision)


class Vector(numpy.ndarray):
    """Vector from a binary stream, with index"""
    def __new__(cls, stream, precision=4):
        pos = stream.tell() 
        if precision == 4:
            buf = numpy.array(struct.unpack(">fff", stream.read(12)))
        else:
            buf = numpy.array(struct.unpack(">ddd", stream.read(24)))
        obj = numpy.ndarray.__new__(cls, shape=(3,), dtype="float", buffer=buf)
        obj.position  = pos
        obj.precision = precision
        return obj

    def pack(self, value):
        if issubclass(value, numpy.ndarray):
            value = value.tolist()
        if self.precision == 4:
            return struct.pack(">fff", value)
        else:
            return struct.pack(">ddd", value)


class RealVector(Vector):
    """
    Vector with floats/doubles, depending on precision of stream. With index.

    This requires the stream to have an attribute 'precision', which
    determines whether floating points are encoded in single (4) or
    in double (!4) precision.
    """
    def __new__(cls, stream):
        return Vector.__new__(cls, stream, precision=stream.precision)


class Matrix(numpy.ndarray):
    """Matrix from a binary stream, with index"""
    def __new__(cls, stream, precision=4):
        pos = stream.tell()
        if precision == 4:
            buf = numpy.array(struct.unpack(">fffffffff", stream.read(36)))
        else:
            buf = numpy.array(struct.unpack(">ddddddddd", stream.read(72)))
        obj = numpy.ndarray.__new__(cls, shape=(3,3), dtype="float", buffer=buf)
        obj.position  = pos
        obj.precision = precision
        return obj

    def pack(self, value):
        if isinstance(value, numpy.ndarray):
            value = value.flatten()
        if self.precision == 4:
            return struct.pack(">fffffffff", *value)
        else:
            return struct.pack(">ddddddddd", *value)


class RealMatrix(Matrix):
    """
    Matrix with floats/doubles, depending on precision of stream. With index.

    This requires the stream to have an attribute 'precision', which
    determines whether floating points are encoded in single (4) or
    in double (!4) precision.
    """
    def __new__(cls, stream):
        return Matrix.__new__(cls, stream, precision=stream.precision)


class Coordinates(numpy.ndarray):
    """Coordinate array from binary stream, with index"""
    def __new__(cls, stream, n, precision=4):
        pos = stream.tell()
        if precision == 4:
            buf = numpy.fromstring(stream.read(3*n*4),">f").astype("float")
        else:
            buf = numpy.fromstring(stream.read(3*n*8),">d").astype("float")
        obj = numpy.ndarray.__new__(cls, shape=(n,3), dtype="float", buffer=buf)
        obj.position  = pos
        obj.double    = precision != 4
        return obj
        
    def pack(self, arr):
        if not isinstance(arr, numpy.ndarray):
            arr = numpy.array(arr)
        if self.double:
            return arr.astype(">d").tostring()
        else:
            return arr.astype(">f").tostring()


class String(str):
    """String from a binary stream, with index"""
    def __new__(cls, stream):
        pos = stream.tell()
        n   = Integer(stream)+3 
        n  -= (n%4)             # Length of the string with padding
        ret = stream.read(4)    # Return '\x00\x00\x00\x07' (BELL) or '\x00\x00\x00\x0d' (CR) 
        obj = str.__new__(cls,stream.read(n).rstrip('\x00'))
        obj.position = pos
        return obj

    def pack(self):
        return self.ljust(len(self)+(4-len(self)%4),'\x00')


class Strings(tuple):
    """Set of strings from a binary stream, with index"""
    def __new__(cls, stream):
        pos = stream.tell()
        n   = Integer(stream)
        out = []
        for i in range(n):
            x  = Integer(stream)
            # Padding to full bytes
            m  = Integer(stream) + 3        
            m -= (m%4)
            out.append(stream.read(m))
        obj = tuple.__new__(cls, out)
        obj.position = pos
        return obj
  

class Tuple(tuple):
    """
    Tuple of somethings from binary stream, with index

    If n is an integer the corresponding 'typ' will be read n times, 
    giving a tuple of items.

    If n is a list/tuple, then for each number in it a tuple is read
    with the given 'typ', giving a nested tuple.

    If typ is a list/tuple, then the sequence will be read n times,
    giving a nested tuple. If n is a tuple/list, the result will be
    a double nested tuple. Further nesting is possible by providing
    a nested typ.
    """
    def __new__(cls, stream, n, typ):
        pos = stream.tell()
        if isinstance(n, int):
            if type(typ) in (list, tuple):
                if n == -1:
                    stuff = [ f(stream) for f in typ ]
                else:
                    stuff = [ Tuple(stream, -1, typ) for i in range(n) ]
            else:
                stuff = [ typ(stream) for i in range(n) ]
        else:
            stuff = [ Tuple(stream, num, typ) for num in n ]
        obj = tuple.__new__(cls, stuff)
        obj.position = pos
        return obj

    def pack(self, value):
        out = []
        if type(value) not in (tuple, list):
            for typ in self:
                out.append(typ.pack(value))
        else:
            for typ, val in zip(self, value):
                out.append(typ.pack(val))
        return "".join(out)


class Group(Tuple):
    """Like a Tuple, but reads the number of items from the stream."""
    def __new__(cls, stream, typ):
        pos = stream.tell()
        num = Integer(stream)
        obj = Tuple.__new__(cls, stream, num, typ)
        obj.position = pos
        obj.num = num
        return obj

    def pack(self, value, tp=None):
        if tp is None and len(self):
            tp = self[0]
        return "".join([ self.num.pack(len(value)) ] + [ tp.pack(val) for val in value ])


class ListWithNames(list):
    def __getattr__(self, attr):
        for item, value in self:
            if item == attr:
                return value
        raise IndexError("Item not found:",attr)

    def __getitem__(self, item):
        for i in self:
            if i[0] == item and i[1] != None:
                return i[1]
        raise IndexError        

    def __contains__(self, item):
        for i in self:
            if i[0] == item and i[1] != None:
                return True
        return False
    
    def search(self, item):
        if item in self:  #ListWithNames.__contains__(self,param):
            return self[item]
        else:
            for i in self:
                found = isinstance(i[1], ListWithNames) and i[1].search(item)
                if found is not None and found is not False:
                    return found
        return None


class TPRIO(io.BytesIO,versions.Test):

    def init(self):
        self.gmx        = String(self)
        self.precision  = Integer(self)
        self.version    = Integer(self)
        self.tag        = String(self)  if 77 <= self.version <= 80 else None
        self.generation = Integer(self) if       self.version >= 26 else 0
        self.tag        = String(self)  if       self.version >= 81 else self.tag
        self.seek(0)

    def __len__(self):
        """Get remaining length"""
        pos = self.tell()
        self.seek(0,2)
        end = self.tell()
        self.seek(pos)
        return end - pos

    def parse(self, pattern):

        parseDict = dict(
            S=self.readString, 
            I=self.readInteger, 
            F=self.readFloat, 
            D=self.readDouble, 
            M=self.readFloatMatrix, 
            A=self.readDoubleMatrix, 
            V=self.readFloatVector, 
            L=self.readDoubleVector,
            Q=self.readLongLong,
            U=self.readUnsigned,
            )

        if self.precision == 8:
            pattern = pattern.replace("R","D")
            pattern = pattern.replace("M","A")
            pattern = pattern.replace("V","L")
        else:
            pattern = pattern.replace("R","F")

        if len(pattern) == 1:
            return parseDict[pattern]()

        out = []
        for i in pattern:
            o = parseDict[i]()
            out.append(o)
        return out

    def do(self,stuff,target,version=None):
        for attr, tp, test in stuff:
            setattr(target,attr,target.test(test,version) and self.parse(tp))                            

    def readInteger(self):
        return struct.unpack(">l",self.read(4))[0]

    def readLongLong(self):
        return struct.unpack(">Q",self.read(8))[0]

    def readUnsigned(self):
        return struct.unpack(">I",self.read(4))[0]

    def readCharacter(self):
        return chr(self.readInteger())

    def readFloat(self):
        return struct.unpack(">f",self.read(4))[0]

    def readDouble(self):
        return struct.unpack(">d",self.read(8))[0]

    def readReal(self):
        if self.precision == 8:
            return struct.unpack(">d",self.read(8))[0]
        return struct.unpack(">f",self.read(4))[0]        

    def readFloatMatrix(self):
        #print self.tell(), len(self)
        return numpy.array(struct.unpack(">fffffffff",self.read(36))).reshape((3,3))
    
    def readDoubleMatrix(self):
        return numpy.array(struct.unpack(">ddddddddd",self.read(72))).reshape((3,3))

    def readFloatVector(self):
        return struct.unpack(">fff",self.read(12))

    def readDoubleVector(self):
        return struct.unpack(">ddd",self.read(24))

    def readString(self):
        n = self.readInteger()+3
        n -= (n%4)
        ret = self.read(4)
        return self.read(n)

    def readStrings(self):
        # Number of strings
        n = self.readInteger()
        o    = []
        for i in range(n):
            x = self.readInteger()
            # Padding to full bytes
            m = self.readInteger() + 3        
            m -= (m%4)         
            o.append((x,self.read(m)))
        return o
  

