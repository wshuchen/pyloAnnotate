#!/usr/bin/env python

# pyloAnnotate.py by Wenshu Chen, 2023

# Based on Pyliftover (https://github.com/konstantint/pyliftover) and 
# Pybedtools (https://daler.github.io/pybedtools/),
# this Python program is a command line tool for
# liftover coordinate conversiont and annotation retrieval.
# It works for a single point, an interval as well as a file of points.

# Pyliftover uses 0-based coordinate as in chain file. A single point coordinate 
# (chromosome, position) will be the start position in a BED format.

# Here "query" and "target" refer to "target" and "query" respectively in UCSC chain-speak.
# Pyliftover refers to "query" as "source".

import os
import sys
import argparse
import subprocess
import pandas as pd
from pybedtools import BedTool
from pyliftover import LiftOver

def get_liftover_point(chain_file, query_point):
    '''
    Get the target search result for single point coordinate.
    '''
    unlifted = []
    lo = LiftOver(chain_file)
    query = query_point.split(",")
    chrom, pos = query[0], int(query[1])  
    result = lo.convert_coordinate(chrom, pos)
    if result is None:
        print(f"No chain data for query: {chrom}, {pos}")
        sys.exit()
    if result == []:
        print(f"No target for query: {chrom}, {pos}")
        unlifted.append((f"{chrom}\t{pos}\t{pos+1}\t+"))
    # In case of multiple targets, get the coordinates out of list.
    if len(result) > 1:
        result = [i for p in result for i in p]
    return result, unlifted

def get_liftover_file(chain_file, query_file):
    '''
    Get the target search result for multiple point coordinates
    in a tab-delimited file.
    '''
    lo = LiftOver(chain_file)
    points = []
    result = []
    unlifted = []
    count = 1
    with open(query_file) as f: 
        for line in f:
            if not line.startswith("chr"):
                print(f"Not valid coordiante in query file: {line}.")
                continue
            else:
                line = line.split("\t")
                if len(line) == 1:
                    print(f"Not tab-delimited: {line}")
                    continue
                points.append((line[0], int(line[1])))
        points = sorted(points)
    for point in points:
        chrom, pos = point[0], point[1]
        target = lo.convert_coordinate(chrom, pos)
        if target is None:
            print(f"No chain data for query: {chrom}, {pos}")
            continue
        if len(target) == 0:
            print(f"No target for query: {chrom}, {pos} ({count})")
            unlifted.append((f"{chrom}\t{pos}\t{pos+1}\t+"))
            count += 1
            continue
        result.append(target)
    # Get the tuples out of lists so the output contains a list of tuples.
    if len(result) > 1:
        result = [i for p in result for i in p]
    return result, unlifted

def write_interval_point(interval):
    """
    Turn an interval specified on command line into a list of every position.
    Could expand this to a file later.
    """
    ipoints = []
    pos = interval.split(",")
    chr, start, end = pos[0], int(pos[1]), int(pos[2])
    for i in range(start, end, 1):
        ipoints.append((chr, i))
    return ipoints

def get_liftover_interval(chain_file, query_interval):
    """
    Interval search by point search.
    """
    lo = LiftOver(chain_file)
    points = []
    result = []
    unlifted = []
    count = 1
    ipoints = write_interval_point(query_interval)
    for point in ipoints:
        chrom, pos = point[0], point[1]
        target = lo.convert_coordinate(chrom, pos)
        if target is None:
            print(f"No chain data for query: {chrom}, {pos}")
            continue
        if len(target) == 0:
            print(f"No target for query: {chrom}, {pos} ({count})")
            unlifted.append((f"{chrom}\t{pos}\t{pos+1}\t+"))
            count += 1
            continue
        result.append(target)
    # Get the tuples out of lists so the output contains a list of tuples.
    if len(result) > 1:
        result = [i for p in result for i in p]
    return result, unlifted

def clean_gff(gff_file):
    '''
    Write six columns of a genome annotation file to BED format.
    chr, start, end, strand, feature, gene_name (or associated parent attribute)
    chr5	143636995	143696005	+	mRNA	Cyth3
    Works for NCBI, ESEMBL, or GENCODE GFF3, or GENCODE GTF in current renditions.

    An inconsistences in the attributes of a annotation file may break this function,
    whick will in turn break the merge_feature function processing data frame
    (unmatched columns). A remedy is to add a print statement to identify the offensive line
    and edit the code accordingly.
    '''
    clean_annot = []
    name = ""
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.split("\t")
            if line[2] in ["region", "match", "cDNA_match", "chromosome", 
                            "supercontig", "biological_region", "."]:
                continue
            if line[8].startswith("ID"):
                name = line[8].split("gene=")  # ncbi
                if len(name) == 1:  # cannot be split
                    name = line[8].split("Name=")  # esembl
                    if len(name) == 1:
                        name = line[8].split("gene_name=")  # gencode
                        if len(name) == 1:
                            name = line[8].split(":") # esembl CDS
                        else:
                            print("Annotation style not recognized.")
                            name = line[8][0].strip()
                            break
                name = name[1].split(";")[0]               
            elif line[8].startswith("Parent"):  # canon-gff3 ncbi or esembl
                name = line[8].split("=")  # rna-* or transcript:*
                if len(name) > 2:  # esemble gff3
                    name = line[8].split("Name=")[1].split(";")[0]
                else:
                    name = name[1].strip()
            elif line[8].startswith("gene_id"):  # GTF
                name = line[8].split(" gene_name ")
                if len(name) > 1:  # Not necessay but just in case
                    name = name[1].split(";")[0].replace('"', '')
                else:
                    name = line[8][0].strip()
            clean_annot.append((line[0], int(line[3])-1, line[4], line[6], line[2], name))
    return sorted(list(set(clean_annot)))

def get_feature(x):
    '''
    Collect the features into one list, remove duplicates, set priority,
    and join items in the last tier into one comma separated string.
    '''
    # Collect the features into a list.
    features = ','.join(x).split(',')
    # Remove duplicates.
    features = list(sorted(set(features)))
    # Set priority.
    # features = [i for i in features if i not in ["exon", "mRNA", "gene"]]   
    if "pseudogene" in features:
        feature = "pseudogene"
    elif "CDS" in features:
        feature = "CDS"
    elif "five_prime_UTR" in features:
        feature = "five_prime_UTR"
    elif "five_prime_utr" in features:
        feature = "five_prime_utr"
    elif "three_prime_UTR" in features:
        feature = "three_prime_UTR"
    elif "three_prime_utr" in features:
        feature = "three_prime_utr"
    elif "UTR" in features:
        feature = "UTR"
    else:
        feature = ','.join(features)
    return feature

def join_gene(x):
    '''
    Collect the gene names into one string.
    '''
    # Collect gene and RNA names into a list and remove duplicates.
    genes = ",".join(x).split(",")
    genes = list((sorted(set(genes))))
    # Keep only gene names. Names will be in alphabetical order,
    # which may not be in the order of features in final data frame. 
    genes = ",".join([i for i in genes if 'rna-' not in i])
    return genes

def merge_feature(annotations):
    '''
    Merge overlapping coordinates and collect features into one string.
    The new start is the minimum of original starts, 
    The new end is the maximum of original ends.
    Potential gaps in the orginal coordiantes will not be considered.
    '''
    annotation = annotations.to_dataframe()  # pybedtools method
    # "a" (point) and "b" (GFF) as -a and -b respectively in bedtools intersect.
    annotation.columns = ["a_chr", "a_start", "a_end", "a_strand", 
                "b_chr", "b_start", "b_end", "b_strand", "feature", "gene"]
    agg_how = {'a_chr':'first', 'a_start':'first', 'a_end':'first', 'a_strand':'first', 
                'b_chr': 'first', "b_start":'min', 'b_end':'max', 'b_strand':'first',
                'feature':get_feature, 'gene':join_gene}
    merged_annotation = annotation.groupby(["a_chr", "a_start", "a_end"], 
                as_index=False, sort=False).aggregate(agg_how).reindex()
    merged_annotation = merged_annotation.sort_values(by=["a_chr", "a_start"])  
    return merged_annotation

def get_annotation(liftover_result, q_gff_file, t_gff_file):
    '''
    Intersect for features of query and target coordinates.
    '''
    q_gff = BedTool(clean_gff(q_gff_file))
    t_gff = BedTool(clean_gff(t_gff_file))
    # Seperate the points from the search result using command line cut.
    # Add a PID tag to prevent overwriting and removal in parallel runs.
    temp = f"lo{os.getpid()}.bed"
    BedTool(liftover_result).saveas(temp)
    q_points = subprocess.run(["cut", "-f1-4", temp], 
                            capture_output=True, text=True)
    t_points = subprocess.run(["cut", "-f5-8", temp], 
                            capture_output=True, text=True)
    q_points = BedTool(q_points.stdout, from_string=True)
    t_points = BedTool(t_points.stdout, from_string=True).sort()
    subprocess.run(["rm", temp])
    # Intersect for features.
    q_annot = q_points.intersect(q_gff, wa=True, wb=True)     
    t_annot = t_points.intersect(t_gff, wa=True, wb=True)
    # What may happen with single point search
    # No annotation at all.
    if len(q_annot) == 0 and len(t_annot) == 0:
        for i in liftover_result:
            print("\t".join(map(str, i)))
        print("No annotation available.\n")
        sys.exit()
    # Only query or target gets annotation.
    if len(q_annot) == 0 and len(t_annot) != 0:
        q_point = "\t".join(map(str, liftover_result[0][0:4]))
        q_annot = f"{q_point}\tNA\tNA\tNA\tNA\tNA\tNA"
        t_annot = merge_feature(t_annot).astype(str).to_string(index=False, header=False)
        # Simply put them together.  
        print(f"{q_annot}\t{t_annot}\t{liftover_result[0][8]}\n")
        sys.exit()
    if len(q_annot) != 0 and len(t_annot) == 0:
        q_annot = merge_feature(q_annot).astype(str).to_string(index=False, header=False)
        t_point = "\t".join(map(str, liftover_result[0][4:8]))
        t_annot = f"{t_point}\tNA\tNA\tNA\tNA\tNA\tNA"
        print(f"{q_annot}\t{t_annot}\t{liftover_result[0][8]}\n")
        sys.exit()
    # All well; and, typically with a file of multiple points.
    else:   
        q_annot = merge_feature(q_annot).astype(str)
        t_annot = merge_feature(t_annot).astype(str)
    q_annot.columns = ["q_chr", "q_start", "q_end", "q_strand", 
                    "qf_chr", "qf_start", "qf_end", "qf_strand","q_feature", "q_gene"]
    t_annot.columns = ["t_chr", "t_start", "t_end", "t_strand", 
                    "tf_chr", "tf_start", "tf_end", "tf_strand", "t_feature", "t_gene"]
    # Merge the data frames guided by the search result.
    lo_df = BedTool(liftover_result).to_dataframe().astype(str)
    lo_df.columns = ["q_chr", "q_start", "q_end", "q_strand",
                     "t_chr", "t_start", "t_end", "t_strand", "segment_size"]
    df = pd.merge(lo_df, q_annot, how="outer", 
                    on=["q_chr", "q_start", "q_end", "q_strand"], sort=False)
    df = pd.merge(df, t_annot, how="outer", 
                    on=["t_chr", "t_start", "t_end", "t_strand"], sort=False)
    annotation = pd.concat([df.iloc[:, 0:4], df.iloc[:, 9:15], 
                    df.iloc[:, 4:8], df.iloc[:, 15:21], df.iloc[:, 8]], 
                    axis=1).sort_values(by=["q_chr", "q_start"])
    return annotation

def annotate_unlifted(unlifted_points, q_gff_file):
    '''
    Intersect for features of unlifted (query without target) points.
    '''
    q_gff = BedTool(clean_gff(q_gff_file))
    unlifted = BedTool(unlifted_points)
    unlifted_annot = unlifted.intersect(q_gff, wa=True, wb=True)
    # Single point, no annotation      
    if len(unlifted_annot) == 0:
        print("No annotation available.\n")
        sys.exit()
    # Single point with annotation, or multiple points.
    else:   
        unlifted_annot = merge_feature(unlifted_annot).astype(str)
    unlifted_annot.columns = ["u_chr", "u_start", "u_end", "u_strand", 
                "uf_chr", "uf_start", "uf_end", "uf_strand","u_feature", "u_gene"]
    # Merge with point file to preserve the points without annotation.
    unlifted_point_df = BedTool(unlifted_points).to_dataframe().astype(str)
    unlifted_point_df.columns = ["u_chr", "u_start", "u_end", "u_strand"]
    unlifted_annot = pd.merge(unlifted_point_df, unlifted_annot, how="outer", 
                on=["u_chr", "u_start", "u_end", "u_strand"], 
                sort=False).sort_values(by=["u_chr", "u_start"]) 
    return unlifted_annot  

def main():
    parser = argparse.ArgumentParser(description=
            "pyloAnnotate converts genome coordinates with annotation")
    parser.add_argument("-c", "--chain_file", type=str, 
                        help="chain file, uncompressed or gzipped")
    parser.add_argument("-p", "--query_point", type=str, 
                        help="a query point in the form of 'chrom,position'")
    parser.add_argument("-f", "--query_file", type=str, 
                        help="query points in a file, delimited with tab or in BED format")
    parser.add_argument("-i", "--query_interval", type=str, 
                        help="query interval in the form of 'chrom,start,end'")
    parser.add_argument("-q", "--query_gff", type=str, help="query annotation file")
    parser.add_argument("-t", "--target_gff", type=str, help="target annotation file")
    parser.add_argument("-s", "--segment_size", type=int, 
                        help="chain segment (aligned block) size)")
    parser.add_argument("-u", "--unlifted", action="store_true", 
                        help="annotate unlifted points")
    parser.add_argument("-d", "--distance", type=int, default=0, 
                        help="maximum distance allowed for merging interval blocks")               

    args = parser.parse_args()
    
    chain_file     = args.chain_file
    query_point    = args.query_point
    query_file     = args.query_file
    query_interval = args.query_interval
    q_gff_file     = args.query_gff
    t_gff_file     = args.target_gff
    segment_size   = args.segment_size

    if len(sys.argv)==1:   
        parser.print_help()
        print("\nPlease provide minimum a query position and a chain file.\n")
        sys.exit(1)

    # If both annotation files are provided, annotation will be done;
    # Otherwise, the program is similar to orginal Pyliftover 
    # but also takes an interval or point file.
    # Result for a single point input will be printed to stdout; 
    # result for an interval or a file will be saved as a BED (no annotation)
    # or a tab-delimited file (with annotation).
    # Annotation of unlifted point is optional.  
    
    if args.query_point:
        print("\nMapping query coordinates to target")
        print(f"Query point: {args.query_point}")
        print(f"Chain file: {args.chain_file}\n")
        result, unlifted = get_liftover_point(chain_file, query_point)
        if not args.query_gff and not args.target_gff:
            for i in result:
                print("\t".join(map(str, i)))
            print("")
        if args.query_gff and args.target_gff:
            if len(result) > 0:
                print("\nRetrieving features for lifted point")
                print(f"Query annotation files: {args.query_gff}")
                print(f"Target annotation file: {args.target_gff}\n")
                annotation = get_annotation(result, q_gff_file, t_gff_file)       
                print(annotation.to_string(header=False, index=False))
                print("")
            else:
                if args.unlifted:
                    if len(unlifted) == 0:
                        print("\nNo unlifted point.")
                        sys.exit()
                    print("\nRetrieving features for unlifted point")
                    print(f"Query annotation files: {args.query_gff}")
                    print(f"Target annotation file: {args.target_gff}\n")
                    unlifted_annot = annotate_unlifted(unlifted, args.query_gff)
                    print(unlifted_annot.to_string(header=False, index=False))
                    print("")

    if args.query_file:
        print("\nMapping query coordinates to target")
        print(f"Query file: {args.query_file}")
        print(f"Chain file: {args.chain_file}\n")
        query_name = os.path.basename(args.query_file).split(".")[0]
        chain_name = os.path.basename(args.chain_file).split(".")[0]
        out_name = f"{query_name}_{chain_name}"
        result, unlifted = get_liftover_file(chain_file, query_file)
        if args.segment_size:  # Filter the results by segment_size if given
            result = [i for i in result if i[-1] >= int(segment_size)]
        if not args.query_gff and not args.target_gff:  
            out_file = f"{out_name}.bed"
            print(f"\nDone. The result was saved as {out_file}\n")
            BedTool(result).saveas(out_file)
        if args.query_gff and args.target_gff:
            print("\nRetrieving features for lifted points")
            print(f"Query annotation files: {args.query_gff}")
            print(f"Target annotation file: {args.target_gff}\n")
            out_file = f"{out_name}.txt"    
            annotation = get_annotation(result, q_gff_file, t_gff_file)
            annotation.to_csv(out_file, header=False, index=False, sep="\t")
            print(f"Done. The result was saved as {out_file}\n")
            if args.unlifted:
                if len(unlifted) == 0:
                    print("\nNo unlifted points.")
                    sys.exit()
                print("Retrieving features for unlifted points")
                out_file = f"{out_name}_unlifted.txt"
                unlifted_annot = annotate_unlifted(unlifted, args.query_gff)
                unlifted_annot.to_csv(out_file, header=False, index=False, sep="\t")
                print(f"\nDone. The result was saved as {out_file}\n")

    if args.query_interval:
        print("\nMapping query coordinates to target")
        print(f"Query interval: {args.query_interval}")
        print(f"Chain file: {args.chain_file}\n")
        query_name = "_".join((args.query_interval).split(","))
        chain_name = os.path.basename(args.chain_file).split(".")[0]
        out_name = f"{query_name}_{chain_name}"
        result, unlifted = get_liftover_interval(chain_file, query_interval)
        if not args.query_gff and not args.target_gff:  
            out_file = f"{out_name}.bed"
            print(f"\nDone. The result was saved as {out_file}\n")
            BedTool(result).merge(
                d=args.distance, 
                c="4,5,6,7,8,9", 
                o="distinct,distinct,min,max,distinct,distinct").saveas(out_file)
        if args.query_gff and args.target_gff:
            print("\nRetrieving features for interval")
            print(f"Query annotation files: {args.query_gff}")
            print(f"Target annotation file: {args.target_gff}\n")
            out_file = f"{out_name}.txt"    
            annotation = get_annotation(result, q_gff_file, t_gff_file)
            BedTool.from_dataframe(annotation).merge(
                d=args.distance, 
                c="4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21", 
                o="distinct,distinct,min,max,distinct,distinct,distinct,distinct,min,max,distinct,distinct,min,max,distinct,distinct,distinct,distinct"
                ).saveas(out_file)
            print(f"Done. The result was saved as {out_file}\n")
            if args.unlifted:
                print("Retrieving features for unlifted positions")
                out_file = f"{out_name}_unlifted.txt"
                if len(unlifted) == 0:
                    print("\nNo unlifted points.")
                    sys.exit()
                unlifted_annot = annotate_unlifted(unlifted, args.query_gff)
                BedTool.from_dataframe(unlifted_annot).merge(
                    d=args.distance, 
                    c="4,5,6,7,8,9,10", 
                    o="distinct,distinct,min,max,distinct,distinct,distinct",
                    ).saveas(out_file)
                print(f"\nDone. The result was saved as {out_file}\n")
                       
if __name__ == "__main__":
    main()
