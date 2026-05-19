#!/usr/bin/env python3
import gzip
import sys
import argparse

def parse_args():
    p = argparse.ArgumentParser(
        description="Restore ONT header tags to demuxed FASTQ using UUID lookup TSV")
    p.add_argument("fastq",          help="Input demuxed FASTQ (gzipped)")
    p.add_argument("output",         help="Output reheadered FASTQ (gzipped)")
    p.add_argument("tsv", nargs="+", help="One or more uuid_tags.tsv files")
    return p.parse_args()

def build_lookup(tsv_files):
    lookup = {}
    for f in tsv_files:
        with open(f) as fh:
            for line in fh:
                line = line.rstrip()
                if not line:
                    continue
                parts = line.split("\t", 1)
                if len(parts) == 2:
                    lookup[parts[0]] = parts[1]
    print(f"[reheader_reads] Loaded {len(lookup)} UUID entries from {len(tsv_files)} TSV(s)",
          file=sys.stderr)
    return lookup

def reheader(fastq_in, fastq_out, lookup):
    missing = 0
    total   = 0
    with gzip.open(fastq_in, "rt") as inp, \
         gzip.open(fastq_out, "wt") as out:
        while True:
            # Read exactly 4 lines per FASTQ record — avoids misidentifying
            # quality score lines that start with '@' as header lines
            header   = inp.readline()
            if not header:
                break           # EOF
            seq      = inp.readline()
            plus     = inp.readline()
            qual     = inp.readline()

            if not header.startswith("@"):
                print(f"[reheader_reads] WARNING: unexpected header line: {header.rstrip()}",
                      file=sys.stderr)
                continue

            header_body = header[1:].rstrip("\n")
            parts = header_body.split(None, 1)   # split on first whitespace only
            uuid  = parts[0] if parts else ""

            if not uuid:
                # Empty header — write as-is and continue
                out.write(header)
                out.write(seq)
                out.write(plus)
                out.write(qual)
                continue

            total += 1
            extra  = lookup.get(uuid)

            if extra:
                out.write(f"@{uuid}{extra}\n")
            else:
                missing += 1
                out.write(header)

            out.write(seq)
            out.write(plus)
            out.write(qual)

    if missing:
        print(f"[reheader_reads] WARNING: {missing}/{total} reads had no tag match",
              file=sys.stderr)
    else:
        print(f"[reheader_reads] {total} reads reheadered -> {fastq_out}",
              file=sys.stderr)

def main():
    args   = parse_args()
    lookup = build_lookup(args.tsv)
    reheader(args.fastq, args.output, lookup)

if __name__ == "__main__":
    main()
