#!/usr/bin/env python3
"""
breakdancer2vcf.py
Convert BreakDancer.ctx output to VCF format.
Usage: breakdancer2vcf.py --input <ctx> --output <vcf> --sample <id> --reference <hg38|hg19>
"""

import argparse
import sys
from datetime import date

SVTYPE_MAP = {
    "DEL": "DEL", "INS": "INS", "INV": "INV",
    "ITX": "DUP", "CTX": "BND", "Unknown": "BND"
}

def parse_args():
    p = argparse.ArgumentParser(description="Convert BreakDancer ctx to VCF")
    p.add_argument("--input",     required=True,  help="BreakDancer .ctx file")
    p.add_argument("--output",    required=True,  help="Output VCF file")
    p.add_argument("--sample",    required=True,  help="Sample ID")
    p.add_argument("--reference", default="hg38", help="Genome build (hg38 or hg19)")
    return p.parse_args()

def write_header(fh, sample, reference):
    fh.write("##fileformat=VCFv4.2\n")
    fh.write(f"##fileDate={date.today().strftime('%Y%m%d')}\n")
    fh.write(f"##source=BreakDancer\n")
    fh.write(f"##reference={reference}\n")
    fh.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
    fh.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
    fh.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
    fh.write('##INFO=<ID=SCORE,Number=1,Type=Integer,Description="BreakDancer score">\n')
    fh.write('##INFO=<ID=NUM_READS,Number=1,Type=Integer,Description="Supporting reads">\n')
    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fh.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")

def convert(args):
    sv_count = 0
    with open(args.input) as fin, open(args.output, 'w') as fout:
        write_header(fout, args.sample, args.reference)
        for line in fin:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chr1, pos1, ori1, chr2, pos2, ori2, svtype, size, score = parts[:9]
            num_reads = parts[9] if len(parts) > 9 else "."
            svtype_vcf = SVTYPE_MAP.get(svtype, "BND")
            try:
                sv_len = int(size)
                end    = int(pos1) + abs(sv_len)
                qual   = int(score)
            except ValueError:
                continue
            alt = f"<{svtype_vcf}>"
            info = f"SVTYPE={svtype_vcf};SVLEN={sv_len};END={end};SCORE={score};NUM_READS={num_reads}"
            flt = "PASS" if qual >= 30 else "LowQual"
            fout.write(f"{chr1}\t{pos1}\tBreakDancer_{sv_count}\tN\t{alt}\t{qual}\t{flt}\t{info}\tGT\t./.\n")
            sv_count += 1
    print(f"Converted {sv_count} BreakDancer calls to VCF", file=sys.stderr)

if __name__ == "__main__":
    convert(parse_args())
