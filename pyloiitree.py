# pyloiitree.py

# Implementation of Heng Li's implicit interval tree (iitree) to Pyliftover
# for interval conversion.

# Pyliftover was writtern by Konstantin Tretyakov
# https://github.com/konstantint/pyliftover

# iitree.py was adapted from Heng Li's bedcov-py-iitree.py
# https://github.com/lh3/cgranges

# Please see the packages for their original form and more imformation.

import sys
import gzip

# iitree.py
# Adapted from bedcov-py-iitree.py by Heng Li

# Each node contains a interval (start, end)
# max_val is the largest interval end in the subtree, including this node
class node:
    def __init__(self, start, end, d=None):
        self.start = start
        self.end = end
        self.max_val = end
        self.data = d 

#  StackCell is used for checking overlapping intervals
class StackCell:
    def __init__(self, k, x, w):
        self.k = k # level of tree
        self.x = x # node index
        self.w = w # status, indicating whether left child has been processed

# This function index a given array a.
# Input:
#    a:  an array containing all the nodes
# Output:
#   k-1: max level of array 'a'
def index_core(a):
    # check for
    i,last_i,k = 0,0,1
    if (len(a) == 0): 
        return -1 # max level = -1 for empty array
    while(i < len(a)):
        last_i = i
        a[i].max_val = a[i].end
        last = a[i].max_val
        i +=2
    while(2**k < len(a)): # process internal nodes in the bottom-up order
        x = 2**(k-1) # x: index of the node
        i0 = (x*2) - 1 # i0 is the first node
        step = x*4  # x: index of the node
        i=i0
        # traverse all nodes at level k
        while(i <len(a)): 
            end_left = a[i - x].max_val; # max value of the left child
            end_right = a[i + x].max_val if i + x < len(a) else last # max value of the right child
            end = a[i].end
            a[i].max_val = max(end,end_left,end_right) # set the max value for node i to 
                                                       # max value of currect sub-tree
            i+=step
        # last_i now points to the index of parent of the original last_i
        last_i = last_i - x if last_i/2**k == 1 else last_i + x
        if last_i < len(a): # update 'last' accordingly
            if a[last_i].max_val > last:
                last = a[last_i].max_val # update max value for the whole tree
        k+=1
    return k - 1 # retruen total level of the array a

# This function checks all overlaping intervals in a given array
# Input:
#   a: an array containing interval nodes
#   max_level: max tree level of array 'a'
#   start, end : from input interval
# Output:
#   overlap_index: a list containing all the overlapping node index
def overlap(a, max_level, start, end):
    t = 0
    overlap_index = []
    # push the root; this is a top down traversal
    stack = [None]*64 # initialize an object list
    stack[t] = StackCell(max_level,2**max_level-1,0) # root, top-down traversal
    t+=1
    while t: # the following guarantees that numbers in "overlap_index" are always sorted
        z = stack[t-1] 
        t-=1
        # 1. if we are in a small subtree; traverse every node in this subtree
        if z.k <= 3:
            i0 = int(z.x /2**z.k) * 2**z.k # i0, start node index in the subtree
            i1 = i0 + 2**(z.k+1) - 1 # i1, maximum node index in subtree (next node at level k:i0+2^(k+1))
            if i1 >= len(a): 
                i1 = len(a)
            i = i0
            while (i < i1): 
                if (a[i].start < end) & (start < a[i].end): # if overlap, append to overlap_index[]
                    overlap_index.append(i)
                i+=1
        # 2. for a large tree, if left child not processed
        elif z.w == 0: 
            y = z.x - 2**(z.k-1) # the index of left child of z.x; 
                                 # NB: y may be out of range (i.e. y>=len(a))
            stack[t] = StackCell(z.k, z.x, 1); # re-add node z.x, but mark the left child having been processed
            t+=1
            if y >= len(a): # push the left child if y is out of range 
                stack[t] = StackCell(z.k - 1, y, 0)
                t+=1
            elif a[y].max_val > start: # push the left child if y  may overlap with the query
                stack[t] = StackCell(z.k - 1, y, 0)
                t+=1
        # 3. need to push the right child
        elif z.x < len(a):
            if ((a[z.x].start < end) & (start < a[z.x].end)): # test if z.x overlaps the query; if yes, append to overlap_index[]
                overlap_index.append(z.x)
            stack[t] = StackCell(z.k - 1, z.x + 2**(z.k-1), 0) # push the right child
            t+=1
    return overlap_index
                    
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
        self.chains = self._load_chains(f) # chain block headers
        self.chain_index, self.maxlevel_dict = self._index_chains(self.chains) # chain data in tree
        
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
            for (sfrom, sto, tfrom, tto) in c.blocks:
                if not c.source_name in chain_index: # check for chrom 
                    chain_index[c.source_name] = []
                # Add new node to the same chrom list
                chain_index[c.source_name].append(node(sfrom, sto, (tfrom, tto, c)))
        maxlevel_dict = {}
        for chrom in chain_index:
            chain_index[chrom]= sorted(chain_index[chrom], key=lambda l:l.start) # sort
            maxlevel_dict[chrom]=index_core(chain_index[chrom]) # append max level to another dictionary
        return chain_index, maxlevel_dict

    def query(self, chromosome, start, end):
        '''
        Given an interval, returns all matching records from the chain index.
        Each record is an interval (source_from, source_to, data)
        where data = (target_from, target_to, chain). Note that depending on chain.target_strand, 
        the target values may need to be reversed (e.g. pos --> chain.target_size - position).
        '''
        if type(chromosome).__name__ == 'bytes':
            chromosome = chromosome.decode('utf-8')
        if chromosome not in self.chain_index:
            return None
        else:
            # node index
            overlap_index = overlap(self.chain_index[chromosome], 
                                    self.maxlevel_dict[chromosome], 
                                    start, end)
            return overlap_index

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
            self.blocks.append((sfrom, sfrom+size, tfrom, tfrom+size))
            sfrom += size + sgap
            tfrom += size + tgap
            fields = f.readline().split()
        if len(fields) != 1:
            raise Exception("Expecting one number on the last line of alignments block. (%s)" % header)
        size = int(fields[0])
        self.blocks.append((sfrom, sfrom+size, tfrom, tfrom+size))
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
       
    def convert_coordinate(self, chromosome, start, end, strand='+'):
        '''
        Returns a *list* of possible conversions for a given chromosome position.
        The list may be empty (no conversion), have a single element (unique conversion), or several 
        '''
        query_results = self.chain_file.query(chromosome, start, end)
        if query_results is None:
            return None
        else:
            results = []
            for i in query_results:
                item = self.chain_file.chain_index[chromosome][i]
                query_start, query_end, data = item.start, item.end, item.data # sfrom, sto, data
                target_start, target_end, chain = data    # (tfrom, tto, c)
                segment_size = target_end - target_start
                if end - start == 1:  # point
                    target_start = target_start + (start - query_start)  # target positions
                    target_end = target_start + 1
                    query_start = query_start + (start - query_start)
                    query_end = query_start + 1                
                if chain.target_strand == '-':
                    target_start, target_end = chain.target_size-1-target_end, chain.target_size-1-target_start
                result_strand = chain.target_strand if strand == '+' else ('+' if chain.target_strand == '-' else '-')
                results.append((chromosome, query_start, query_end, '+', 
                                chain.target_name, target_start, target_end, 
                                result_strand, segment_size))
            return results
