
import argparse
import subprocess
# from pyfasta import Fasta
from pyfaidx import Fasta
import tempfile


#create argument parser
parser = argparse.ArgumentParser(
    description="Automatically identify primers around a SNP position."
)
parser.add_argument("fasta_file", type=str, help="Path to FASTA file")
parser.add_argument(
    "chromosome",
    type=str,
    help="Chromosome that must match FASTA header."
)
parser.add_argument("position", type=int, help="Position to amplify")
args = parser.parse_args()


#open fasta file and read sequence
genome = Fasta(args.fasta_file, as_raw=True, sequence_always_upper=True)
if args.chromosome not in genome:
    available = ", ".join(list(genome.keys())[:20])
    raise KeyError(f"Chromosome '{args.chromosome}' not found. Available: {available}")

chrom_seq = str(genome[args.chromosome])
chrom_len = len(chrom_seq)


# Identify sequence 500 bp upstream and downstream of requested position and store as string
pos0 = args.position - 1 if args.position >= 1 else args.position
if pos0 < 0 or pos0 >= chrom_len:
    raise ValueError(f"Position {args.position} is out of bounds.")

flank = 500
start = max(0, pos0 - flank)
end = min(chrom_len, pos0 + flank)
template_seq = chrom_seq[start:end]
snp_offset = pos0 - start


# Add braces around SNP position for visualization
vis_left = max(0, snp_offset - 100)
vis_right = min(len(template_seq), snp_offset + 100)
visualized_seq = (
    template_seq[:vis_left]
    + "{"
    + template_seq[vis_left:vis_right]
    + "}"
    + template_seq[vis_right:]
)


# Create input file containing the DNA sequence and specify a product length of 600-800 bp
left_start = 0
left_len = max(0, snp_offset - 100)
right_start = min(len(template_seq), snp_offset + 100)
right_len = max(0, len(template_seq) - right_start)

primer3_input = "\n".join(
    [
        f"SEQUENCE_ID={args.chromosome}:{args.position}",
        f"SEQUENCE_TEMPLATE={template_seq}",
        (
            "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST="
            f"{left_start},{left_len},{right_start},{right_len}"
        ),
        "PRIMER_TASK=generic",
        "PRIMER_NUM_RETURN=1",
        "PRIMER_OPT_SIZE=20",
        "PRIMER_MIN_SIZE=18",
        "PRIMER_MAX_SIZE=25",
        "PRIMER_OPT_TM=60.0",
        "PRIMER_MIN_TM=57.0",
        "PRIMER_MAX_TM=63.0",
        "PRIMER_MIN_GC=20.0",
        "PRIMER_MAX_GC=80.0",
        "PRIMER_PRODUCT_SIZE_RANGE=600-800",
        "=",
        "",
    ]
)

with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".p3") as tf:
    tf.write(primer3_input)
    primer3_input_path = tf.name


# Run primer3_core
with open(primer3_input_path, "r") as fh:
    result = subprocess.run(
        ["primer3_core"],
        stdin=fh,
        capture_output=True,
        text=True,
    )

if result.returncode != 0:
    raise RuntimeError(
        "primer3_core failed.\n"
        f"STDERR:\n{result.stderr}\n"
        f"STDOUT:\n{result.stdout}\n"
    )

primer3_stdout = result.stdout


# Parse primer3 output to extract primer sequences
out = {}
for line in primer3_stdout.splitlines():
    line = line.strip()
    if not line or line == "=":
        continue
    if "=" not in line:
        continue
    k, v = line.split("=", 1)
    out[k] = v

left_primer = out.get("PRIMER_LEFT_0_SEQUENCE")
right_primer = out.get("PRIMER_RIGHT_0_SEQUENCE")
print(left_primer.lower())
print(right_primer.lower())