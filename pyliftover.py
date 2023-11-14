# pyliftover.py

# Pyliftover was writtern by Konstantin Tretyakov
# https://github.com/konstantint/pyliftover

# Here three major scripts (intervaltree.py, chainfile.py, and liftover.py) 
# of Pyliftover were combined into one piece named pyliftover.py.
# For simplicity, many of the author's comments were removed and doc strings reduced.
# The program now only reads a provided chain file without downloading functionality.
# Finally, the convert_coordinate function in LiftOver.py was modified to output
# needed data in BED format.

# Please see the package for its original form and more imformation.

import sys
import gzip

# intervertree.py

class IntervalTree:
    '''
    Interval Tree data structure for indexing a set of 
    integer intervals of the form [start, end). 
    '''
    __slots__ = ['min', 'max', 'center', 'single_interval', 'left_subtree', 'right_subtree', 
                'mid_sorted_by_start', 'mid_sorted_by_end']

    def __init__(self, min, max):
        '''
        Creates a tree node for keeping intervals somewhere in the range [min...max).
        '''
        self.min = int(min)
        self.max = int(max)
        assert self.min < self.max
        self.center = (min + max)/2
        self.single_interval = None
        self.left_subtree = None
        self.right_subtree = None 
        self.mid_sorted_by_start = [] 
        self.mid_sorted_by_end = []
    
    def add_interval(self, start, end, data=None):
        '''
        Inserts an interval to the tree. 
        '''
        if (end - start) <= 0:
            return
        if self.single_interval is None:   # Empty tree
            self.single_interval = (start, end, data)
        elif self.single_interval == 0:   # Usual tree
            self._add_interval(start, end, data)
        else:   # Single interval. Convert to a usual tree.
            self._add_interval(*self.single_interval)
            self.single_interval = 0
            self._add_interval(start, end, data)
            
    def _add_interval(self, start, end, data=None):
        if end <= self.center:
            if self.left_subtree is None:
                self.left_subtree = IntervalTree(self.min, self.center)
            self.left_subtree.add_interval(start, end, data)
        elif start > self.center:
            if self.right_subtree is None:
                self.right_subtree = IntervalTree(self.center, self.max)
            self.right_subtree.add_interval(start, end, data)
        else:
            self.mid_sorted_by_start.append((start, end, data))
            self.mid_sorted_by_end.append((start, end, data))
    
    def sort(self):
        '''
        Must be invoked after all intevals have been added to sort mid_** arrays.
        '''
        if self.single_interval is None or self.single_interval != 0:
            return # Nothing to do for empty and leaf trees.
        self.mid_sorted_by_start.sort(key = lambda x: x[0])
        self.mid_sorted_by_end.sort(key = lambda x: x[1], reverse=True)
        if self.left_subtree is not None:
            self.left_subtree.sort()
        if self.right_subtree is not None:
            self.right_subtree.sort()
    
    def query(self, x):
        '''
        Returns all intervals in the tree, which overlap given point, i.e. all (start, end, data) 
        records, for which (start <= x < end).
        '''
        result = []
        self._query(x, result)
        return result
    
    def _query(self, x, result):
        '''
        Same as self.query, but uses a provided list to accumulate results into.
        '''
        if self.single_interval is None: # Empty
            return
        elif self.single_interval != 0:  # Single interval, just check whether x is in it
            if self.single_interval[0] <= x < self.single_interval[1]:
                result.append(self.single_interval)
        elif x < self.center:            # Normal tree, query point to the left of center
            if self.left_subtree is not None:
                self.left_subtree._query(x, result)
            for int in self.mid_sorted_by_start:
                if int[0] <= x:
                    result.append(int)
                else:
                    break
        else:  # Normal tree, query point to the right of center
            for int in self.mid_sorted_by_end:
                if int[1] > x:
                    result.append(int)
                else:
                    break
            if self.right_subtree is not None:
                self.right_subtree._query(x, result)

    def __len__(self):
        '''
        The number of intervals maintained in the tree.
        '''    
        if self.single_interval is None:
            return 0
        elif self.single_interval != 0:
            return 1
        else:
            size = len(self.mid_sorted_by_start)
            if self.left_subtree is not None:
                size += len(self.left_subtree)
            if self.right_subtree is not None:
                size += len(self.right_subtree)
            return size
            
    def __iter__(self):
        if self.single_interval is None:
            return
        elif self.single_interval != 0:
            yield self.single_interval
        else:
            if self.left_subtree is not None:
                for s in self.left_subtree:
                    yield s
            for s in self.mid_sorted_by_start:
                yield s
            if self.right_subtree is not None:
                for s in self.right_subtree:
                    yield s
                    
# chainfile.py

def open_liftover_chain_file(chain_file):
    '''
    Opens chain file for reading.
    '''
    if chain_file.endswith(".gz"):
        return gzip.open(chain_file, 'rb')
    else:
        return open(chain_file, 'rb')

class LiftOverChainFile:
    '''
    Loads and indexes USCS's .over.chain files.
    ''' 
    def __init__(self, f):
        '''
        Reads chain data from the file and initializes an interval index.
        '''
        self.chains = self._load_chains(f)
        self.chain_index = self._index_chains(self.chains)
        
    @staticmethod
    def _load_chains(f):
        '''
        Loads all LiftOverChain objects from a file into an array. Returns the result.
        '''
        chains = []
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith(b'#') or line.startswith(b'\n') or line.startswith(b'\r'):
                continue
            if line.startswith(b'chain'):  # header
                chains.append(LiftOverChain(line, f))
                continue
        return chains

    @staticmethod
    def _index_chains(chains, show_progress=False):
        '''
        Given a list of LiftOverChain objects, creates a
         dict: source_name --> 
            IntervalTree: <source_from, source_to> -->
                (target_from, target_to, chain)
        Returns the resulting dict.
        '''
        chain_index = {}
        source_size = {}
        target_size = {}
        for c in chains:
            # Verify that sizes of chromosomes are consistent over all chains
            source_size.setdefault(c.source_name, c.source_size)
            if source_size[c.source_name] != c.source_size:
                raise Exception("Chains have inconsistent specification of source chromosome size for %s (%d vs %d)" % (c.source_name, source_size[c.source_name], c.source_size))
            target_size.setdefault(c.target_name, c.target_size)
            if target_size[c.target_name] != c.target_size:
                raise Exception("Chains have inconsistent specification of target chromosome size for %s (%d vs %d)" % (c.target_name, target_size[c.target_name], c.target_size))
            chain_index.setdefault(c.source_name, IntervalTree(0, c.source_size))
            # Register all blocks from the chain in the corresponding interval tree
            tree = chain_index[c.source_name]
            for (sfrom, sto, tfrom) in c.blocks:
                tree.add_interval(sfrom, sto, (tfrom, c))
        # Sort all interval trees
        for k in chain_index:
            chain_index[k].sort()
        return chain_index

    def query(self, chromosome, position):
        '''
        Given a chromosome and position, returns all matching records from the chain index.
        Each record is an interval (source_from, source_to, data)
        where data = (target_from, target_to, chain). Note that depending on chain.target_strand, 
        the target values may need to be reversed (e.g. pos --> chain.target_size - pos).
        '''
        if type(chromosome).__name__ == 'bytes':
            chromosome = chromosome.decode('utf-8')
        if chromosome not in self.chain_index:
            return None
        else:
            return self.chain_index[chromosome].query(position)

class LiftOverChain:
    '''
    Represents a single chain from an .over.chain file.
    A chain basically maps a set of intervals from "source" coordinates to 
    corresponding coordinates in "target" coordinates. 
    '''
    __slots__ = ['score', 'source_name', 'source_size', 'source_start', 'source_end',
	             'target_name', 'target_size', 'target_strand', 'target_start', 'target_end', 'id', 'blocks']

    def __init__(self, header, f):
        '''
        Reads the chain from a stream given the first line and a file opened at all remaining lines.
        On error throws an exception.
        '''
        if sys.version_info >= (3, 0):
            header = header.decode('ascii') # In Python 2, work with usual strings.
        fields = header.split()
        if fields[0] != 'chain' and len(fields) not in [12, 13]:
            raise Exception("Invalid chain format. (%s)" % header)
        # chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1
        self.score = int(fields[1])        # Alignment score
        self.source_name = fields[2]       # E.g. chrY
        self.source_size = int(fields[3])  # Full length of the chromosome
        source_strand = fields[4]          # Must be +
        if source_strand != '+':
            raise Exception("Source strand in an .over.chain file must be +. (%s)" % header)
        self.source_start = int(fields[5]) # Start of source region
        self.source_end = int(fields[6])   # End of source region
        self.target_name = fields[7]       # E.g. chr5
        self.target_size = int(fields[8])  # Full length of the chromosome
        self.target_strand = fields[9]     # + or -
        if self.target_strand not in ['+', '-']:
            raise Exception("Target strand must be - or +. (%s)" % header)
        self.target_start = int(fields[10])
        self.target_end = int(fields[11])
        self.id = None if len(fields) == 12 else fields[12].strip()
        
        # Now read the alignment chain from the file and store it as a list (source_from, source_to) 
        # -> (target_from, target_to)
        sfrom, tfrom = self.source_start, self.target_start
        self.blocks = []

        fields = f.readline().decode('ascii').split()
        while len(fields) == 3:
            size, sgap, tgap = int(fields[0]), int(fields[1]), int(fields[2])
            self.blocks.append((sfrom, sfrom+size, tfrom))
            sfrom += size + sgap
            tfrom += size + tgap
            fields = f.readline().split()
        if len(fields) != 1:
            raise Exception("Expecting one number on the last line of alignments block. (%s)" % header)
        size = int(fields[0])
        self.blocks.append((sfrom, sfrom+size, tfrom))
        if (sfrom + size) != self.source_end  or (tfrom + size) != self.target_end:
            raise Exception("Alignment blocks do not match specified block sizes. (%s)" % header)

class LiftOver:
    def __init__(self, chain_file):
        '''
        '''       
        if chain_file.endswith('.gz'):
            f = gzip.open(chain_file)
        else:
            f = open(chain_file)
        f = open_liftover_chain_file(chain_file)
        self.chain_file = LiftOverChainFile(f)
        f.close()
       
    def convert_coordinate(self, chromosome, position, strand='+'):
        '''
        Returns a *list* of possible conversions for a given chromosome position.
        The list may be empty (no conversion), have a single element (unique conversion), or several 
        '''
        query_results = self.chain_file.query(chromosome, position)
        if query_results is None:
            return None
        else:
            # query_results contain the query point.
            results = []
            for (source_start, source_end, data) in query_results:
                target_start, chain = data
                result_position = target_start + (position - source_start)
                if chain.target_strand == '-':
                    result_position = chain.target_size - 1 - result_position
                result_strand = chain.target_strand if strand == '+' else ('+' if chain.target_strand == '-' else '-')
                results.append((chromosome, position, position+1, '+', chain.target_name, 
                                result_position, result_position+1, result_strand, source_end-source_start))
            # results.sort(key=lambda x: x[3], reverse=True)
            return results
