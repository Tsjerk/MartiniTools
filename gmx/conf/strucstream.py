
import gzip, pdb, gro, itertools


def strucStream(stream):

    # First check whether we have have an open stream or a file
    # If it's a file, check whether it's zipped and open it
    if type(stream) == str:
        if stream.endswith("gz"):
            stream = gzip.open(stream)
        else:
            stream = open(stream)

    stored = [stream.next(),stream.next()]
    
    if stored[-1].strip().isdigit():
        # Must be a GRO file
        return gro.groFrameIterator(itertools.chain(stored,stream))
    else:
        # Then must be a PDB file
        return pdb.pdbFrameIterator(itertools.chain(stored,stream))


