#!/usr/bin/env python3
import os
import subprocess
import tempfile
import statistics
from collections import Counter
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment

def IUPAC_ambiguity(base_counts):
    IUPAC_CODES = {
        frozenset(["A"]): "A", frozenset(["C"]): "C", frozenset(["G"]): "G", frozenset(["T"]): "T",
        frozenset(["A", "G"]): "R", frozenset(["C", "T"]): "Y",
        frozenset(["G", "C"]): "S", frozenset(["A", "T"]): "W",
        frozenset(["G", "T"]): "K", frozenset(["A", "C"]): "M",
        frozenset(["A", "C", "G"]): "V", frozenset(["A", "C", "T"]): "H",
        frozenset(["A", "G", "T"]): "D", frozenset(["C", "G", "T"]): "B",
        frozenset(["A", "C", "G", "T"]): "N",
    }
    bases = frozenset(base_counts.keys())
    return IUPAC_CODES.get(bases, "N")

def create_consensus_with_ambiguity(aligned_records, threshold=0.7, gap_fraction_cutoff=0.6):
    alignment = MultipleSeqAlignment(aligned_records)
    consensus = []
    ambiguity_count = 0
    length = alignment.get_alignment_length()

    for i in range(length):
        column = alignment[:, i]
        base_counts = Counter(column)
        total = sum(base_counts.values())
        gap_count = base_counts.get("-", 0)

        # Check gap fraction
        if (gap_count / total) > gap_fraction_cutoff:
            # Skip this column from the consensus
            continue

        non_gap_total = total - gap_count

        if non_gap_total == 0:
            # No non-gap bases
            consensus.append("N")
            ambiguity_count += 1
            continue

        # Identify bases that meet the threshold
        passed = [(base, count) for base, count in base_counts.items()
                  if base in "ACGT" and (count / non_gap_total) >= threshold]

        if len(passed) == 1:
            # One base meets threshold
            consensus.append(passed[0][0])
        elif len(passed) > 1:
            # Multiple bases meet threshold, use ambiguity over these bases
            filtered_counts = {b: c for b, c in base_counts.items() if b in {x[0] for x in passed} and b in "ACGT"}
            
            # Apply 5% filter: only keep bases >5% of total reads
            filtered_counts = {b: c for b, c in filtered_counts.items() if (c / total) > 0.05}

            if not filtered_counts:
                # If no bases remain after 5% filter, use N
                consensus.append("N")
                ambiguity_count += 1
            else:
                code = IUPAC_ambiguity(filtered_counts)
                if code not in "ACGT":
                    ambiguity_count += 1
                consensus.append(code)
        else:
            # No single base meets threshold
            # Use ambiguity over all ACGT present
            filtered_counts = {b: c for b, c in base_counts.items() if b in "ACGT"}
            
            # Apply 5% filter here as well
            filtered_counts = {b: c for b, c in filtered_counts.items() if (c / total) > 0.05}

            if not filtered_counts:
                # No bases remain after filtering, use N
                consensus.append("N")
                ambiguity_count += 1
            else:
                code = IUPAC_ambiguity(filtered_counts)
                if code not in "ACGT":
                    ambiguity_count += 1
                consensus.append(code)

    return "".join(consensus), ambiguity_count, [i for i in range(length) if i < len(consensus)]

def find_fastq_file():
    for f in os.listdir('.'):
        if f.endswith('.fastq'):
            return f
    raise FileNotFoundError("No FASTQ file found in the current directory.")

def run_muscle(input_fasta, output_fasta):
    cmd = ["muscle", "-in", input_fasta, "-out", output_fasta]
    subprocess.run(cmd, check=True)

def align_qualities_to_alignment(aligned_records, original_records):
    quality_map = {rec.id: rec.letter_annotations["phred_quality"] for rec in original_records.values()}
    aligned_qualities = []
    for rec in aligned_records:
        orig_qual = quality_map[rec.id]
        aligned_seq = str(rec.seq)
        qual_index = 0
        new_qual_list = []
        for char in aligned_seq:
            if char == '-':
                new_qual_list.append(None)
            else:
                new_qual_list.append(orig_qual[qual_index])
                qual_index += 1
        aligned_qualities.append((rec.id, aligned_seq, new_qual_list))
    return aligned_qualities

def compute_column_averages(aligned_qualities):
    alignment_length = len(aligned_qualities[0][1]) if aligned_qualities else 0
    column_quals = []
    for col in range(alignment_length):
        column_values = [seq_qual[col] for (_, _, seq_qual) in aligned_qualities if seq_qual[col] is not None]
        avg = statistics.mean(column_values) if column_values else None
        column_quals.append(avg)
    return column_quals

def quality_to_bin(q, binsize=5, max_qual=40):
    if q is None:
        return "gap"
    if q < 0: q = 0
    if q > max_qual: q = max_qual
    bin_index = q // binsize
    return f"q{bin_index}"

def compute_column_base_counts(aligned_qualities):
    """
    Compute base counts per column excluding gaps.
    Returns a list of dictionaries: column_base_counts[i] = {base: count, ...}
    """
    alignment_length = len(aligned_qualities[0][1]) if aligned_qualities else 0
    column_base_counts = []
    for col in range(alignment_length):
        counts = Counter()
        for (_, seq, _) in aligned_qualities:
            base = seq[col]
            if base in "ACGT":
                counts[base] += 1
        column_base_counts.append(counts)
    return column_base_counts

def write_html(aligned_qualities, column_quals, consensus_full_line, consensus_full_qual, 
               output_file="alignment.html", binsize=5, max_qual=40, threshold=0.3):
    alignment_length = len(aligned_qualities[0][1]) if aligned_qualities else 0

    # Precompute column base counts for highlighting mismatches
    column_base_counts = compute_column_base_counts(aligned_qualities)

    def should_highlight(col_index, base):
        if base not in "ACGT":
            return False
        counts = column_base_counts[col_index]
        non_gap_total = sum(counts.values())
        if non_gap_total == 0:
            return False
        count_base = counts[base]
        fraction_of_base = count_base / non_gap_total
        return fraction_of_base < threshold  # less than 30%, means >=70% differ

    # Compute highlight counts per column (for all sequences except the consensus row)
    highlight_counts = [0]*alignment_length
    for (seq_id, aligned_seq, q_list) in aligned_qualities:
        if seq_id == "Consensus":
            # Skip consensus for counting highlights
            continue
        for col_i, base in enumerate(aligned_seq):
            if base in "ACGT" and should_highlight(col_i, base):
                highlight_counts[col_i] += 1

    css_rules = []
    num_bins = (max_qual // 5) + 1
    # Correction: binsize is defined above as a parameter. Use that instead of hard coding 5.
    num_bins = (max_qual // binsize) + 1
    for bin_i in range(num_bins):
        q_val = bin_i * binsize
        fraction = q_val / max_qual if max_qual > 0 else 0
        fraction = max(min(fraction, 1), 0)
        r = int(255 - fraction * 255)
        g = int(fraction * 255)
        b = 0
        color = f"#{r:02X}{g:02X}{b:02X}"
        css_rules.append(f".q{bin_i} {{ background-color: {color}; }}")

    css_rules.append(".gap { background-color: #CCCCCC; }")
    css_rules.append(".highlight { border: 2px solid yellow; }")

    html = []
    html.append("<html><head>")
    html.append("<meta charset='UTF-8'>")
    html.append("<style>")
    html.append("body { font-family: monospace; }")
    html.append("table { border-collapse: collapse; }")
    html.append("th, td { border: 1px solid #ccc; padding: 2px; text-align: center; }")
    html.append("\n".join(css_rules))
    html.append(".avg-row { font-weight: bold; }")
    html.append("</style></head><body>")

    html.append("<h1>Aligned Sequences with Quality Coloring</h1>")
    html.append("<p>Top row: Average Q per column</p>")
    html.append("<p>Second row: Highlight counts (number of cells in that column that differ from >=70% of the column)</p>")
    html.append("<p>Following rows: Sequences labeled Seq1, Seq2, ...</p>")
    html.append("<p>Last row: Consensus. '-' means column not included in final consensus FASTA.</p>")
    html.append("<p>Yellow outline indicates a cell's nucleotide differs from >=70% of that column.</p>")

    html.append("<table>")

    # Header row for average Q
    html.append("<tr class='avg-row'>")
    html.append("<th></th>")
    for avg_q in column_quals:
        if avg_q is None:
            html.append("<th>-</th>")
        else:
            html.append(f"<th>{int(round(avg_q))}</th>")
    html.append("</tr>")

    # Row for highlight counts
    html.append("<tr class='avg-row'>")
    html.append("<th>Highlights</th>")
    for hc in highlight_counts:
        html.append(f"<th>{hc}</th>")
    html.append("</tr>")

    # Rows for sequences
    # We'll write all non-consensus first
    non_consensus = [r for r in aligned_qualities if r[0] != "Consensus"]
    for (seq_id, aligned_seq, q_list) in non_consensus:
        html.append("<tr>")
        html.append(f"<td>{seq_id}</td>")
        for col_i, (base, q) in enumerate(zip(aligned_seq, q_list)):
            if base == '-':
                qclass = "gap"
                highlight = False
            else:
                qclass = quality_to_bin(q, binsize=binsize, max_qual=max_qual)
                highlight = should_highlight(col_i, base)

            cell_class = qclass + (" highlight" if highlight else "")
            html.append(f"<td class='{cell_class}'>{base}</td>")
        html.append("</tr>")

    # Consensus row
    consensus_row = [r for r in aligned_qualities if r[0] == "Consensus"]
    if consensus_row:
        # We know there's only one consensus row if it exists
        html.append("<tr>")
        html.append("<td>Consensus</td>")
        for col_i, (base, q) in enumerate(zip(consensus_full_line, consensus_full_qual)):
            if base == '-':
                qclass = "gap"
                highlight = False  # no highlight for consensus
            else:
                qclass = quality_to_bin(q, binsize=binsize, max_qual=max_qual)
                highlight = False
            cell_class = qclass + (" highlight" if highlight else "")
            html.append(f"<td class='{cell_class}'>{base}</td>")
        html.append("</tr>")

    html.append("</table>")
    html.append("</body></html>")

    # Write using UTF-8 encoding
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("\n".join(html))
    print(f"HTML alignment visualization written to {output_file}")

if __name__ == "__main__":
    fastq_file = find_fastq_file()
    print(f"Found FASTQ file: {fastq_file}")

    original_records = {}
    for record in SeqIO.parse(fastq_file, "fastq"):
        original_records[record.id] = record

    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp_in:
        input_fasta = tmp_in.name
        with open(input_fasta, "w") as fh:
            for rec_id, rec in original_records.items():
                fh.write(f">{rec_id}\n{str(rec.seq)}\n")

    output_fasta = input_fasta + "_aligned.fasta"

    print("Running MUSCLE alignment...")
    run_muscle(input_fasta, output_fasta)

    aligned_records = list(SeqIO.parse(output_fasta, "fasta"))
    aligned_qualities = align_qualities_to_alignment(aligned_records, original_records)
    column_quals = compute_column_averages(aligned_qualities)

    consensus_seq, ambiguity_count, kept_columns = create_consensus_with_ambiguity(
        aligned_records, threshold=0.7, gap_fraction_cutoff=0.4
    )

    alignment_length = len(aligned_qualities[0][1]) if aligned_qualities else 0
    kept_col_set = set(kept_columns)
    col_to_consensus_char = {}
    for idx, c in enumerate(kept_columns):
        col_to_consensus_char[c] = consensus_seq[idx]

    consensus_full_line = []
    consensus_full_qual = []
    for i in range(alignment_length):
        if i in kept_col_set:
            char = col_to_consensus_char[i]
            q = column_quals[i]
            q_val = int(round(q)) if q is not None else 0
            consensus_full_line.append(char)
            consensus_full_qual.append(q_val)
        else:
            consensus_full_line.append('-')
            consensus_full_qual.append(None)

    # Rename sequences as Seq1, Seq2, ...
    for i in range(len(aligned_qualities)):
        seq_id, aligned_seq, q_list = aligned_qualities[i]
        aligned_qualities[i] = (f"Seq{i+1}", aligned_seq, q_list)

    if not kept_columns:
        # No kept columns, add empty consensus row
        aligned_qualities.append(("Consensus", "", []))
        write_html(aligned_qualities, column_quals, consensus_full_line, consensus_full_qual, "alignment.html")
        print("No columns passed the gap threshold. No consensus fasta generated.")
        print("Done.")
        exit(0)

    # Add the consensus row at the end
    # Actually, we already handle consensus row in write_html by checking if it exists in aligned_qualities.
    # Let's append the consensus row now.
    aligned_qualities.append(("Consensus", ''.join(consensus_full_line), consensus_full_qual))

    write_html(aligned_qualities, column_quals, consensus_full_line, consensus_full_qual, "alignment.html")

    consensus_filename = f"consensus-{ambiguity_count}ambig.fasta"
    with open(consensus_filename, "w", encoding="utf-8") as cf:
        cf.write(f">Consensus\n{consensus_seq}\n")
    print(f"Consensus FASTA file written to {consensus_filename}")
    print("Done.")
