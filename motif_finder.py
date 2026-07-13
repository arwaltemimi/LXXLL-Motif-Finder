"""
LXXLL-like Motif Finder
------------------------
Scans protein sequences for LXXLL-like nuclear-receptor-interaction motifs
(the "NR box" and known variants), based on the motif logic described in:

  "Conserved LXXLL-like Motifs in Viral Proteins as Potential Interactors
   of Host Nuclear Receptors: Protein-Virus-Cancer Connections"

Usage examples:
    python motif_finder.py --sequence MEEPQSDPSVEPPLSQETFSDLWKLL
    python motif_finder.py --fasta my_proteins.fasta --plot
    python motif_finder.py --uniprot P03126 --out results.csv --plot
"""

import re
import csv
import argparse
import sys
from dataclasses import dataclass
from typing import List, Optional, Tuple
from Bio import SeqIO

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False


# Known LXXLL-like motif patterns. "." matches any single amino acid (regex).
# These cover the canonical NR box and the most commonly reported variants
# in the nuclear receptor / coactivator literature.
MOTIF_PATTERNS = {
    "LXXLL (canonical NR box)": r"L..LL",
    "FXXLF": r"F..LF",
    "LXXLL variant (LXXIL)": r"L..[IL]L",
    "LXXML": r"L..ML",
    "LXXHL": r"L..HL",
}

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class MotifMatch:
    protein_id: str
    motif_type: str
    matched_sequence: str
    start: int  # 1-indexed, inclusive
    end: int    # 1-indexed, inclusive


def clean_sequence(raw: str) -> str:
    """Strip whitespace/newlines and uppercase a raw sequence string."""
    seq = re.sub(r"\s+", "", raw).upper()
    return seq


def find_motifs(sequence: str, protein_id: str = "") -> List[MotifMatch]:
    """Find all LXXLL-like motif matches in a protein sequence.

    Overlapping matches of the SAME pattern are all reported (uses lookahead
    so e.g. LLLLL doesn't hide a second overlapping motif).
    """
    sequence = clean_sequence(sequence)
    matches: List[MotifMatch] = []

    for motif_name, pattern in MOTIF_PATTERNS.items():
        # lookahead capture group lets us find overlapping matches
        lookahead_pattern = f"(?=({pattern}))"
        for m in re.finditer(lookahead_pattern, sequence):
            matched_seq = m.group(1)
            start = m.start() + 1  # 1-indexed
            end = start + len(matched_seq) - 1
            matches.append(MotifMatch(
                protein_id=protein_id,
                motif_type=motif_name,
                matched_sequence=matched_seq,
                start=start,
                end=end,
            ))

    matches.sort(key=lambda x: (x.start, x.motif_type))
    return matches


def parse_fasta(filepath: str) -> List[Tuple[str, str]]:
    """Parse FASTA using Biopython SeqIO."""
    records = []

    for record in SeqIO.parse(filepath, "fasta"):
        records.append((record.id, str(record.seq)))

    return records

def fetch_uniprot_sequence(accession: str) -> Optional[str]:
    """Fetch a protein sequence from UniProt REST API by accession number."""
    if not HAS_REQUESTS:
        print("The 'requests' library is not installed. Run: pip install requests")
        return None

    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        resp = requests.get(url, timeout=15)
        resp.raise_for_status()
        lines = resp.text.strip().split("\n")
        if not lines or not lines[0].startswith(">"):
            print(f"Unexpected response for accession {accession}.")
            return None
        seq = "".join(lines[1:])
        return seq
    except Exception as e:
        print(f"Error fetching UniProt accession '{accession}': {e}")
        return None


def write_csv_report(all_matches: List[MotifMatch], output_path: str) -> None:
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Protein ID", "Motif Type", "Matched Sequence", "Start", "End"])
        for m in all_matches:
            writer.writerow([m.protein_id, m.motif_type, m.matched_sequence, m.start, m.end])


def plot_motifs(matches: List[MotifMatch], seq_length: int, protein_id: str, output_path: str) -> None:
    if not HAS_MATPLOTLIB:
        print("matplotlib is not installed; skipping plot. Run: pip install matplotlib")
        return

    fig, ax = plt.subplots(figsize=(10, 2.8))

    # backbone
    ax.plot([1, seq_length], [0, 0], color="#d9d9d9", linewidth=8,
            solid_capstyle="round", zorder=1)

    colors = {
        "LXXLL (canonical NR box)": "#d62728",
        "FXXLF": "#1f77b4",
        "LXXLL variant (LXXIL)": "#2ca02c",
        "LXXML": "#9467bd",
        "LXXHL": "#ff7f0e",
    }

    seen_labels = set()
    for m in matches:
        color = colors.get(m.motif_type, "black")
        label = m.motif_type if m.motif_type not in seen_labels else None
        ax.plot([m.start, m.end], [0, 0], color=color, linewidth=12,
                solid_capstyle="butt", zorder=2, label=label)
        seen_labels.add(m.motif_type)

    ax.set_xlim(1, max(seq_length, 2))
    ax.set_ylim(-1, 1)
    ax.set_yticks([])
    ax.set_xlabel("Residue position")
    ax.set_title(f"LXXLL-like motif map: {protein_id}  ({seq_length} aa, {len(matches)} motif hits)")

    if seen_labels:
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.35),
                   ncol=3, fontsize=8, frameon=False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Plot saved to: {output_path}")


def run(records: List[Tuple[str, str]], out_csv: str, make_plot: bool) -> None:
    all_matches: List[MotifMatch] = []

    for protein_id, seq in records:
        seq_clean = clean_sequence(seq)
        invalid_chars = set(seq_clean) - VALID_AA
        if invalid_chars:
            print(f"Warning: '{protein_id}' contains non-standard characters: "
                  f"{sorted(invalid_chars)} (still scanned as-is)")

        matches = find_motifs(seq_clean, protein_id)
        all_matches.extend(matches)

        print(f"\n>{protein_id}  ({len(seq_clean)} aa)  -  {len(matches)} motif hit(s)")
        for m in matches:
            print(f"   [{m.start:>4}-{m.end:<4}] {m.motif_type:<28} {m.matched_sequence}")

    write_csv_report(all_matches, out_csv)
    print(f"\nFull report written to: {out_csv}")

    if make_plot and records:
        first_id, first_seq = records[0]
        first_seq_clean = clean_sequence(first_seq)
        first_matches = [m for m in all_matches if m.protein_id == first_id]
        plot_path = out_csv.rsplit(".", 1)[0] + "_plot.png"
        plot_motifs(first_matches, len(first_seq_clean), first_id, plot_path)


def main():
    parser = argparse.ArgumentParser(
        description="Find LXXLL-like nuclear-receptor-interaction motifs in protein sequences."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--fasta", help="Path to a FASTA file (can contain multiple sequences)")
    group.add_argument("--uniprot", help="UniProt accession number to fetch and scan (e.g. P03126)")
    group.add_argument("--sequence", help="Raw protein sequence given directly on the command line")

    parser.add_argument("--out", default="motif_report.csv", help="Output CSV path (default: motif_report.csv)")
    parser.add_argument("--plot", action="store_true",
                         help="Generate a PNG map of motif positions (first sequence only)")
    args = parser.parse_args()

    if args.fasta:
        records = parse_fasta(args.fasta)
    elif args.uniprot:
        seq = fetch_uniprot_sequence(args.uniprot)
        records = [(args.uniprot, seq)] if seq else []
    else:
        records = [("input_sequence", args.sequence)]

    if not records:
        print("No sequences to process. Exiting.")
        sys.exit(1)

    run(records, args.out, args.plot)


if __name__ == "__main__":
    main()
