# pyloAnnotate
A tool for conversion and annotation of genome coordinates

This command line tool consists of a modified version of [_Pyliftover_](https://github.com/konstantint/pyliftover) and a script to retrieve annotation for both query and target coordinates. It can convert a point coordinate like the original _pyliftover_, but it can also take a file with point coordinates or an interval. The result can be filtered with a threshold of aligned chain segment size. Annotation is enabled when annotation files (GFF3 or GTF) are provided. For a point file or interval input, the output is a BED if no annotation files provided or a tab-delimited file otherwise ready for downstream processing. 

Note here we refer UCSC's liftover (chain) terms "target" and "query" as "query" and "target" respectively for being more intuitive. _Pyliftover_ uses "source" for UCSC's "target" and  "target" for "query". Also note that coordinates in a chain file are 0-based, which _Pyliftover_ also uses. Single position (chromosome, position) corresponds to the _start_ position in a BED format (0-based).

(Update - as the result of a fun exercise, a new version (pyloAnnotateii.py) using implicit interval tree (pyloiitree.py, with tree script adapted from [Heng Li](https://github.com/lh3/cgranges)) was added for faster interval conversion.)

## **Dependencies**

The program was tested on Python version 3.11. The two required Python packages are:
1. [_Pybedtools_](https://daler.github.io/pybedtools/) (depending on [_bedtools_](https://bedtools.readthedocs.io/en/latest/index.html))
2. _Pandas_

## **Getting the program**

Download or copy the scripts into one place.

## **Usage**

```
python3 pyloAnnotate.py [options] [files] [-s segement_size] [-u] [-d distance]
```  
OR (assuming that the executables are in the PATH)
```  
pyloAnnotate.py [options] [files] [-s segement_size] [-u] [-d distance]
```

### **Required arguments**   
_-p_&emsp;_--point_&emsp;&emsp;&emsp;&emsp;point position (chromosome,position)  
_-f_&emsp;_--point_file  
&emsp;&emsp;--interval_file_&emsp;point positions or intervals in a file  
_-i_&emsp;_--interval_&emsp;&emsp;&emsp;interval positions (chromosome,start,end)  
_-c_&emsp;_--chain_&emsp;&emsp;&emsp;&emsp;chain file

( _-p_,_-f_, or _-i_ in one job) 

### **Optional arguments**

_-q_&emsp;_--query_gff_&emsp;&emsp;&emsp;query annotation file  
_-t_&emsp; _--target_gff_&emsp;&emsp;&emsp;target annotation file   
_-s_&emsp;_--segment_size_&emsp;length of the segment a point belongs to (an integer)  
_-u_&emsp;_--unlifted_&emsp;&emsp;&emsp;&emsp;annotate unlifted points (default: False)  
_-m_&emsp;_--merge_&emsp;&emsp;&emsp;&emsp;merge interval blocks (default: False)  
_-d_&emsp;_--distance_&emsp;&emsp;&emsp;&emsp;maximum distance between gaps of query segments allowed when merging (default: 0)

( _-q_ AND _-t_ in one job)  

Annotation files (GFF3 or GTF) can be from ESEMBL, GENCODE, and NCBI. Make sure the chromosome names in the annotation files are in UCSC style (chr1, chr2, etc.).

The _-s_ argument can be supplied with or without _-q_ and _-t_ arguments for a point file. It will reduce the size of output.

When _-u_ argument is given with _-q_ and _-t_ arguments, the program will annotate unlifted points besides lifted ones with slightly increased processing time.

The _-m_ argument is only available for the iitree version.

The _-d_ argument controls the distance between two features allowed to be merged in interval conversion. A larger distance produces a smaller number of entries in the final result. The default is 0. Please refer to [_bedtools_](https://bedtools.readthedocs.io/en/latest/index.html) merge command for more information. The distance can be seen as a gap between query segments in terms of chain, thus _-d_ setting will be the largest gap allowed. It can be supplied with or without annotation files. However, it may be better to merge annoation data in other programs for flexibilty and other considerations. Note that unlifted segment coordinates after a merge only indicate that there are unlifted segments in that interval.

### **Point coordinate format**
A point position on command line should be comma-separated without space:

chr1,60478

Multiple point positions in a file should be tab-delimited:

chr1&emsp;60478  
chr1&emsp;72291  
chr3&emsp;61745213  
chr3&emsp;82091822  

and can be in BED format, in which case the start positions will be used to search the target coordinates:  

chr1&emsp;60478&emsp;60479  
chr1&emsp;72291&emsp;72292  
chr3&emsp;61745213&emsp;61745214  
chr3&emsp;82091822&emsp;82091823

### **interval coordinate format**

Interval positions on command line should be comma-separated without space:  

chr1,60478,72291

### Examples  
(Note: option values for illustration only)  

With out _-q_ and _-t_ options:  
Single point:  
```  
pyloAnnotate.py -p chr21,33896581 -c hg38chr21ToMm39.over.chain.gz
```
Point file (with optional _-s_):  
``` 
pyloAnnotate.py -f 50hg38SNPs.bed -c hg38chr21ToMm39.over.chain.gz -s 20 >& err
```  
Interval:
```
pyloAnnotateii.py -i chr21,14143206,14143300 -c hg38chr21ToMm39.over.chain.gz
# Write out a BED file with all the segments in the interval
```
With annotation files (GENCODE chromosome GTF files with names changed):  
```
pyloAnnotate.py -p chr21,33896581 -c hg38chr21ToMm39.over.chain.gz -q hg38v44.gtf -t mm39vM33.gtf
``` 
```
pyloAnnotate.py -f 50hg38SNPs.bed -c hg38chr21ToMm39.over.chain.gz -q hg38v44.gtf -t mm39vM33.gtf -s 20 >& err
```
Annotate unlifted points as well (_-u_):
```
pyloAnnotate.py -p chr21,10428748 -c hg38chr21ToMm39.over.chain.gz -q hg38v44.gtf -t mm39vM33.gtf -u
``` 
```
pyloAnnotate.py -f 50hg38SNPs.bed -c hg38chr21ToMm39.over.chain.gz -q hg38v44.gtf -t mm39vM33.gtf -s 20 -u >& err
```  
Interval with annotation and merge:
```
pyloAnnotateii.py -i chr21,14143206,14143300 -c hg38chr21ToMm39.over.chain.gz -q hg38v44.gtf -t mm39vM33.gtf -u -m -d 20 >& err
```

## **_pyloAnnotate.py_ or _pyloAnnotateii.py_**
_pyloAnnotate.py_ for points file;  
_pyloAnnotateii.py_ for interval and interval file.

Changes in _pyloAnnotateii.py_:  
1. No _-p_ option. A point needs to be specified as an interval:  
 _-i_ chromosome,start,end.  
2. Optional _-m_ for merge operation. No by default.

## **Output**

If no GFF files are provided with _-q_ AND _-t_ options, the program will print out target coordinate(s) for a point like _Pyliftover_, or save the result as a BED file (.bed) for a group of points
in a file.

A query point without a target (unlifted point) will be written to _stdout_ and numbered in multiple point conversion. The points are saved internally for annotation when requested (_-u_ option).

Full output will be a tab-delimited text file (.txt) consists of 21 fields (columns), 10 for a query point, 10 for its target, and the last field (column) for the segment sizes where the points are in:  

1.	query chromosome
2.	query start  
3.  query end
4.	query strand
5.	query feature chromosome
6.	query feature start
7.	qurey feature end
8.	query feature strand
9.	query feature
10.	query gene name
11.	target chromosome
12.	target start  
13. target end  
14.	target strand
15.	target feature chromosome
16.	target feature start
17.	target feature end
18.	target feature strand
19.	target feature
20.	target gene name
21.	segment size

Multiple features and gene names will be respectively aggregated into comma-separated strings. Cells without annotation data will be _NaN_ or _NA_ after a tab-delimited file is read into a program (e.g., pandas with Jupyter Notebook).

The output of unlifted point (query) annotation will have ten columns for the query plus one for segement length.

The entries in a result of interval conversion depends on the distance (length of gaps) setting. At default value (0), the program will write out every query segment with corresponding target segments in that interval. The segment sizes are combined into one string.

## **Test data**  
A small point file (human SNP sites; for test purpose only) and human chromosome 21 to mouse chain (extracted from UCSC hg38ToMm39 chain) are provided. 

## **Acknowledgments**
The annotation functions were initially practice scripts for developing [_GAinSAW_](https://github.com/BrendelGroup/GAinSAW) package, which was an implementation of Professor [_Volker Brendel_](https://github.com/vpbrendel)'s approach to chain data.
