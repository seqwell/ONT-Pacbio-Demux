#!/usr/bin/env python3
"""
extract_bam_tags.py  <input.bam>  <output.tsv>

Writes TSV:  <read_name>\t <space><SAM-tags space-separated>

The leading space after the tab ensures reheader_reads.py produces:
  @uuid <tag1> <tag2> ...
with a proper space separating UUID from tags in the FASTQ header line.
"""
import pysam, sys

bam_path = sys.argv[1]
out_path = sys.argv[2]

def tag_to_sam(key, val, type_char):
    """Serialize a pysam tag to SAM-style string using the actual type char."""
    if type_char in ('i', 'c', 'C', 's', 'S', 'I'):
        return f"{key}:i:{int(val)}"
    elif type_char == 'f':
        return f"{key}:f:{val}"
    elif type_char in ('Z', 'A'):
        return f"{key}:{type_char}:{val}"
    elif type_char == 'H':
        return f"{key}:H:{val.hex() if isinstance(val, bytes) else val}"
    elif type_char == 'B':
        if isinstance(val, tuple) and len(val) == 2:
            subtype, arr = val
            vals_str = ",".join(str(v) for v in arr)
            return f"{key}:B:{subtype},{vals_str}"
        else:
            vals_str = ",".join(str(v) for v in val)
            return f"{key}:B:i,{vals_str}"
    else:
        return f"{key}:Z:{val}"

written = 0
with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam, \
     open(out_path, "w") as out:
    for read in bam.fetch(until_eof=True):
        if read.query_name is None:
            continue
        tags_str = " ".join(
            tag_to_sam(k, v, t) for k, v, t in read.get_tags(with_value_type=True)
        )
        # NOTE: tab then SPACE before tags — reheader_reads.py stores the value
        # after the tab as `extra` and writes @uuid{extra}, so the space here
        # ensures the FASTQ header reads "@uuid qs:f:... " not "@uuidqs:f:..."
        out.write(f"{read.query_name}\t {tags_str}\n")
        written += 1

print(f"[extract_bam_tags] Wrote {written} entries to {out_path}", file=sys.stderr)
