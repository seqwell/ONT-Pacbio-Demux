#!/usr/bin/env python3
import pysam, gzip, sys, re, os, shutil, subprocess

fastq_path = sys.argv[1]
out_path   = sys.argv[2]
tsv_files  = sys.argv[3:]

sample_id     = os.path.basename(out_path).replace(".bam", "")
barcode_label = sample_id.split(".")[0]

TAG_RE = re.compile(r'([A-Za-z][A-Za-z0-9]):([AifZHBc]):(.+)')

def open_fastq(path):
    if path.endswith(".gz") and shutil.which("pigz"):
        proc = subprocess.Popen(
            ["pigz", "-dc", "-p", "2", path],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        return proc.stdout
    return gzip.open(path, "rb")

# --- Pass 1: collect needed UUIDs from FASTQ (~instant) ---
needed = set()
fh = open_fastq(fastq_path)
readline = fh.readline
while True:
    h = readline()
    if not h:
        break
    sp = h.find(b' ')
    uuid = (h[1:sp] if sp != -1 else h[1:]).rstrip().decode()
    needed.add(uuid)
    readline(); readline(); readline()  # skip seq, +, qual
fh.close()

print(f"[fastq_to_bam] {len(needed)} reads in FASTQ", file=sys.stderr)

def parse_tag(tag_str):
    m = TAG_RE.match(tag_str)
    if not m:
        return None
    key, typ, val = m.group(1), m.group(2), m.group(3)
    if typ == 'i':
        return (key, int(val))
    elif typ == 'f':
        return (key, float(val))
    elif typ in ('Z', 'A'):
        return (key, val)
    elif typ == 'H':
        return (key, bytes.fromhex(val))
    elif typ in ('B', 'c'):
        parts = val.split(",")
        subtype = parts[0]
        try:
            if subtype in ('c', 'C', 's', 'S', 'i', 'I'):
                return (key, (subtype, [int(x) for x in parts[1:]]))
            else:
                return (key, (subtype, [float(x) for x in parts[1:]]))
        except Exception:
            return None
    return (key, val)

# --- Pass 2: scan TSVs but only parse rows we actually need ---
lookup = {}
for tsv in tsv_files:
    if len(lookup) == len(needed):
        break                           # found everything, stop early
    with open(tsv, buffering=1 << 20) as fh:
        for line in fh:
            tab = line.find('\t')
            if tab == -1:
                continue
            uuid = line[:tab]
            if uuid not in needed:
                continue                # skip parse entirely — just a dict lookup
            parsed = [t for t in (parse_tag(t) for t in line[tab+1:].split()) if t is not None]
            parsed = [(k, v) for k, v in parsed if k != 'BC']
            parsed.append(('BC', barcode_label))
            lookup[uuid] = parsed
            if len(lookup) == len(needed):
                break                   # found all — stop scanning this file

print(f"[fastq_to_bam] Loaded {len(lookup)}/{len(needed)} UUIDs from TSVs", file=sys.stderr)

EMPTY_TAGS = [('BC', barcode_label)]

sam_header = {
    "HD": {"VN": "1.6", "SO": "unknown"},
    "PG": [{"ID": "fastq_to_bam", "PN": "fastq_to_bam"}]
}

# --- Pass 3: write BAM ---
written = missing = 0
with pysam.AlignmentFile(out_path, "wb", header=sam_header) as bam:
    bam_header = bam.header
    fh = open_fastq(fastq_path)
    readline = fh.readline
    try:
        while True:
            h_raw = readline()
            if not h_raw:
                break
            seq_raw  = readline()
            readline()
            qual_raw = readline()

            h_body = h_raw[1:].rstrip()
            sp     = h_body.find(b' ')
            uuid   = (h_body[:sp] if sp != -1 else h_body).decode()
            seq    = seq_raw.rstrip().decode()
            qual   = qual_raw.rstrip().decode()

            a = pysam.AlignedSegment(bam_header)
            a.query_name      = uuid
            a.query_sequence  = seq
            a.flag            = 4
            a.query_qualities = pysam.qualitystring_to_array(qual)

            if uuid in lookup:
                tags = lookup[uuid]
                written += 1
            else:
                tags = EMPTY_TAGS
                missing += 1

            try:
                a.set_tags(tags)
            except Exception as e:
                safe = [(k, v) for k, v in tags if not isinstance(v, tuple)]
                a.set_tags(safe)
                print(f"[fastq_to_bam] WARNING: tag error on {uuid}: {e}", file=sys.stderr)

            bam.write(a)
    finally:
        fh.close()

print(f"[fastq_to_bam] {written} reads with tags, {missing} without", file=sys.stderr)
