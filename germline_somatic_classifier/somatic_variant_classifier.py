#!/usr/bin/env python3

import argparse
import math
import os
import sys
import csv
from multiprocessing import Pool, cpu_count
from cyvcf2 import VCF
import pysam
import re
import warnings
import time

warnings.filterwarnings("ignore", message="no intervals found")

GNOMAD = None
COSMIC = None
CONFIG = None


####################################
# Utility Logging (stderr only)
####################################

def log(msg):
    print(msg, file=sys.stderr, flush=True)


####################################
# Chromosome Normalization
####################################

def strip_chr(chrom: str) -> str:
    return chrom.replace("chr", "")

def add_chr(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else "chr" + chrom


####################################
# Scoring
####################################

def compute_germline_score(af, an, an_threshold):
    if af is None or an is None:
        return 0.0
    confidence = min(an / an_threshold, 1.0)
    return af * confidence

def compute_somatic_score(cosmic_count):
    if cosmic_count is None:
        return 0.0
    return math.log10(cosmic_count + 1)


####################################
# Reference Accessors
####################################

class GnomadAccessor:
    def __init__(self, path):
        self.vcf = VCF(path)
        first_contig = next(iter(self.vcf.seqnames))
        self.has_chr = first_contig.startswith("chr")

    def normalize(self, chrom):
        chrom = strip_chr(chrom)
        return add_chr(chrom) if self.has_chr else chrom

    def query(self, chrom, pos, ref, alt):
        for var in self.vcf(f"{self.normalize(chrom)}:{pos}-{pos}"):

            if strip_chr(var.CHROM) != strip_chr(chrom):
                continue
            if var.REF != ref:
                continue

            for i, allele in enumerate(var.ALT):
                if allele == alt:
                    af = var.INFO.get("AF")
                    an = var.INFO.get("AN")
                    if isinstance(af, (list, tuple)):
                        af = af[i]
                    return af, an

        return None, None


class CosmicAccessor:
    def __init__(self, path):
        self.tabix = pysam.TabixFile(path)

    def query(self, chrom, pos, ref, alt):
        try:
            records = self.tabix.fetch(add_chr(strip_chr(chrom)), pos-1, pos)
        except Exception:
            try:
                records = self.tabix.fetch(strip_chr(chrom), pos-1, pos)
            except Exception:
                return None, ""

        total = 0
        tissues = set()

        for line in records:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue

            if (strip_chr(parts[0]) == strip_chr(chrom) and
                int(parts[1]) == pos and
                parts[3] == ref and
                parts[4] == alt):

                matches = re.findall(r"(\d+)\(([^)]+)\)", parts[5])
                for count, tissue in matches:
                    total += int(count)
                    tissues.add(tissue)

        if total == 0:
            return None, ""

        return total, ",".join(sorted(tissues))


####################################
# Worker Init
####################################

def init_worker(gnomad_path, cosmic_path, config):
    global GNOMAD, COSMIC, CONFIG
    GNOMAD = GnomadAccessor(gnomad_path)
    COSMIC = CosmicAccessor(cosmic_path)
    CONFIG = config


####################################
# Classification
####################################

def classify_variant(g_score, s_score, g_thr, s_thr):
    if g_score >= g_thr and s_score < s_thr:
        return "LIKELY_GERMLINE"
    if g_score < g_thr and s_score >= s_thr:
        return "LIKELY_SOMATIC"
    if g_score >= g_thr and s_score >= s_thr:
        return "CONFLICTING"
    return "UNKNOWN"


def process_variant(var_dict):

    chrom = var_dict["chrom"]
    pos = var_dict["pos"]
    ref = var_dict["ref"]
    alt = var_dict["alt"]

    g_af, g_an = GNOMAD.query(chrom, pos, ref, alt)
    c_count, c_tissues = COSMIC.query(chrom, pos, ref, alt)

    g_score = compute_germline_score(g_af, g_an, CONFIG["an_threshold"])
    s_score = compute_somatic_score(c_count)

    classification = classify_variant(
        g_score, s_score,
        CONFIG["g_thr"],
        CONFIG["s_thr"]
    )

    var_dict.update({
        "gnomad_af": g_af if g_af is not None else "",
        "gnomad_an": g_an if g_an is not None else "",
        "cosmic_count": c_count if c_count is not None else "",
        "cosmic_tissues": c_tissues,
        "germline_score": round(g_score, 4),
        "somatic_score": round(s_score, 4),
        "classification": classification
    })

    return var_dict


####################################
# Main
####################################

def main():
    try:

        parser = argparse.ArgumentParser()
        parser.add_argument("--sample-vcf", required=True)
        parser.add_argument("--gnomad-vcf", required=True)
        parser.add_argument("--cosmic-tsv", required=True)
        parser.add_argument("--output", required=True)
        parser.add_argument("--threads", type=int, default=cpu_count())
        parser.add_argument("--germline-threshold", type=float, default=0.1)
        parser.add_argument("--somatic-threshold", type=float, default=1.0)
        parser.add_argument("--an-confidence-threshold", type=int, default=100000)

        args = parser.parse_args()

        log("Starting Somatic vs Germline Variant Classifier")

        config = {
            "g_thr": args.germline_threshold,
            "s_thr": args.somatic_threshold,
            "an_threshold": args.an_confidence_threshold
        }

        os.makedirs(os.path.dirname(args.output), exist_ok=True)

        sample_vcf = VCF(args.sample_vcf)

        with open(args.output, "w", newline="") as out_f:

            fieldnames = [
                "chrom","pos","ref","alt",
                "filter","sample_dp","sample_af",
                "gnomad_af","gnomad_an",
                "cosmic_count","cosmic_tissues",
                "germline_score","somatic_score",
                "classification"
            ]

            writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()

            jobs = []

            for var in sample_vcf:

                filter_status = var.FILTER if var.FILTER not in (None, "PASS", ".") else ""

                dp = var.format("DP")
                af = var.format("AF")
                ad = var.format("AD")

                sample_dp = dp[0][0] if dp is not None else ""

                sample_af = ""
                if af is not None:
                    sample_af = af[0][0]
                elif ad is not None and len(ad[0]) > 1:
                    ref_count = ad[0][0]
                    alt_count = ad[0][1]
                    total = ref_count + alt_count
                    if total > 0:
                        sample_af = alt_count / total

                for alt in var.ALT:
                    jobs.append({
                        "chrom": strip_chr(var.CHROM),
                        "pos": var.POS,
                        "ref": var.REF,
                        "alt": alt,
                        "filter": filter_status,
                        "sample_dp": sample_dp,
                        "sample_af": sample_af
                    })

            start = time.time()

            with Pool(args.threads, initializer=init_worker,
                      initargs=(args.gnomad_vcf, args.cosmic_tsv, config)) as pool:

                for result in pool.imap(process_variant, jobs, chunksize=200):
                    writer.writerow(result)

            elapsed = round(time.time() - start, 2)

        log(f"Completed successfully in {elapsed} seconds. Exit Code: 0")
        sys.exit(0)

    except Exception as e:
        log(f"ERROR: {str(e)}. Exit Code: 1")
        sys.exit(1)


if __name__ == "__main__":
    main()
