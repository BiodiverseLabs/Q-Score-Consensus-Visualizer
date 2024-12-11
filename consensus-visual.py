#!/usr/bin/env python3
import os
import subprocess
import tempfile
import statistics
from collections import Counter
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

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

def create_consensus_with_ambiguity(aligned_records, threshold=0.7, gap_fraction_cutoff=0.4):
    alignment = MultipleSeqAlignment(aligned_records)
    consensus = []
    ambiguity_count = 0
    length = alignment.get_alignment_length()
    kept_columns = []

    for i in range(length):
        column = alignment[:, i]
        base_counts = Counter(column)
        total = sum(base_counts.values())
        gap_count = base_counts.get("-", 0)
        n_count = base_counts.get("N", 0)

        # Check gap fraction
        if (gap_count / total) > gap_fraction_cutoff:
            continue

        kept_columns.append(i)
        non_gap_total = total - gap_count - n_count

        if non_gap_total == 0:
            # No ACGT bases, use N
            consensus.append("N")
            ambiguity_count += 1
            continue

        # Identify bases that meet threshold
        passed = [(base, count) for base, count in base_counts.items()
                  if base in "ACGT" and (count / non_gap_total) >= threshold]

        def apply_5_percent_filter(bases_dict):
            return {b:c for b,c in bases_dict.items() if (c/total) > 0.05}

        if len(passed) == 1:
            b = passed[0][0]
            # 5% filter
            if (base_counts[b]/total) <= 0.05:
                consensus.append("N")
                ambiguity_count += 1
            else:
                consensus.append(b)
        elif len(passed) > 1:
            filtered_counts = {b: c for b,c in base_counts.items() if b in {x[0] for x in passed} and b in "ACGT"}
            filtered_counts = apply_5_percent_filter(filtered_counts)
            if not filtered_counts:
                consensus.append("N")
                ambiguity_count += 1
            else:
                code = IUPAC_ambiguity(filtered_counts)
                if code not in "ACGT":
                    ambiguity_count += 1
                consensus.append(code)
        else:
            # No base meets threshold
            filtered_counts = {b:c for b,c in base_counts.items() if b in "ACGT"}
            filtered_counts = apply_5_percent_filter(filtered_counts)
            if not filtered_counts:
                consensus.append("N")
                ambiguity_count += 1
            else:
                code = IUPAC_ambiguity(filtered_counts)
                if code not in "ACGT":
                    ambiguity_count += 1
                consensus.append(code)

    return "".join(consensus), ambiguity_count, kept_columns

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
    if not aligned_qualities:
        return []
    alignment_length = len(aligned_qualities[0][1])
    column_quals = []
    for col in range(alignment_length):
        column_values = [q_list[col] for (_, _, q_list) in aligned_qualities if q_list[col] is not None]
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
    if not aligned_qualities:
        return []
    alignment_length = len(aligned_qualities[0][1])
    column_base_counts = []
    for col in range(alignment_length):
        counts = Counter()
        for (_, seq, _) in aligned_qualities:
            base = seq[col]
            if base in "ACGT":
                counts[base] += 1
        column_base_counts.append(counts)
    return column_base_counts

def records_from_aligned_qualities(aligned_qualities):
    records = []
    for (sid, s, ql) in aligned_qualities:
        records.append(SeqRecord(Seq(s), id=sid, description=""))
    return records

def write_html(aligned_qualities, column_quals, consensus_full_line, consensus_full_qual, highlight_counts,
               top_consensus_full_line, top_consensus_full_qual,
               highlight_map,  # highlight_map is a set of (row,col) that should be highlighted
               output_file="alignment.html", binsize=5, max_qual=40):
    alignment_length = len(aligned_qualities[0][1]) if aligned_qualities else 0

    css_rules = []
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
    html.append("<p>Second row: Highlight counts (cells differing ≥70%)</p>")
    html.append("<p>Third row: Consensus(top) after N conversions</p>")
    html.append("<p>Then sequences: Seq1, Seq2, ...</p>")
    html.append("<p>Last row: Consensus(bottom). '-' means column excluded. 'N' means forced changes or no suitable base.</p>")
    html.append("<p>Yellow outline: cell differs from ≥70% or was changed to N from highlight+Q<12 condition.</p>")

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

    # Top consensus row
    html.append("<tr>")
    html.append("<td>Consensus(top)</td>")
    for col_i, (base, q) in enumerate(zip(top_consensus_full_line, top_consensus_full_qual)):
        qclass = quality_to_bin(q, binsize=binsize, max_qual=max_qual) if q is not None else "gap"
        highlight = False
        cell_class = qclass + (" highlight" if highlight else "")
        html.append(f"<td class='{cell_class}'>{base}</td>")
    html.append("</tr>")

    # Non-consensus sequences
    non_consensus = [r for r in aligned_qualities if r[0] != "Consensus"]
    for row_i, (seq_id, aligned_seq, q_list) in enumerate(non_consensus):
        html.append("<tr>")
        html.append(f"<td>{seq_id}</td>")
        for col_i, (base, q) in enumerate(zip(aligned_seq, q_list)):
            if q is None:
                qclass = "gap"
            else:
                qclass = quality_to_bin(q, binsize=binsize, max_qual=max_qual)
            highlight = (row_i, col_i) in highlight_map
            cell_class = qclass + (" highlight" if highlight else "")
            html.append(f"<td class='{cell_class}'>{base}</td>")
        html.append("</tr>")

    # Bottom consensus row
    consensus_row = [r for r in aligned_qualities if r[0] == "Consensus"]
    if consensus_row:
        html.append("<tr>")
        html.append("<td>Consensus</td>")
        for col_i, (base, q) in enumerate(zip(consensus_full_line, consensus_full_qual)):
            qclass = quality_to_bin(q, binsize=binsize, max_qual=max_qual) if q is not None else "gap"
            # bottom consensus highlight = no special highlight needed for consensus row itself
            highlight = False
            cell_class = qclass + (" highlight" if highlight else "")
            html.append(f"<td class='{cell_class}'>{base}</td>")
        html.append("</tr>")

    html.append("</table>")
    html.append("</body></html>")

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("\n".join(html))
    print(f"HTML alignment visualization written to {output_file}")

if __name__ == "__main__":
    fastq_file = find_fastq_file()
    print(f"Found FASTQ file: {fastq_file}")

    original_records = {r.id:r for r in SeqIO.parse(fastq_file, "fastq")}

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

    highlight_threshold = 0.3
    column_base_counts = compute_column_base_counts(aligned_qualities)

    # Detect highlight cells pre-transform
    pre_highlight_cells = []
    n_converted_cells = set()  # to store cells changed to N
    for row_i, (seq_id, seq, q_list) in enumerate(aligned_qualities):
        for col_i, (base, q) in enumerate(zip(seq, q_list)):
            if base in "ACGT":
                counts = column_base_counts[col_i]
                non_gap_total = sum(counts.values())
                if non_gap_total > 0:
                    fraction_of_base = counts[base]/non_gap_total
                    if fraction_of_base < highlight_threshold:
                        pre_highlight_cells.append((row_i, col_i))

    # Convert highlight+Q<12 → N
    for (r, c) in pre_highlight_cells:
        seq_id, seq, q_list = aligned_qualities[r]
        q = q_list[c]
        if q is not None and q < 12:
            seq = list(seq)
            seq[c] = 'N'
            seq = ''.join(seq)
            q_list[c] = None
            aligned_qualities[r] = (seq_id, seq, q_list)
            n_converted_cells.add((r, c))  # mark this cell as N-converted

    # Recompute after transformations
    column_base_counts = compute_column_base_counts(aligned_qualities)

    # Recompute highlight counts after N conversion
    highlight_counts = [0]*len(aligned_qualities[0][1])
    # Remember we don't count N cells in highlight counts (they act like gaps for highlight frequency)
    for row_i, (seq_id, seq, q_list) in enumerate(aligned_qualities):
        if seq_id == "Consensus":
            continue
        for col_i, (base, q) in enumerate(zip(seq, q_list)):
            if base in "ACGT":
                counts = column_base_counts[col_i]
                non_gap_total = sum(counts.values())
                if non_gap_total > 0:
                    fraction_of_base = counts[base]/non_gap_total
                    if fraction_of_base < highlight_threshold:
                        highlight_counts[col_i] += 1
            # N or '-' do not count towards highlight counts now

    # Recompute column averages ignoring N (q=None)
    column_quals = compute_column_averages(aligned_qualities)

    # Re-run consensus after changes
    final_records = records_from_aligned_qualities(aligned_qualities)
    consensus_seq, ambiguity_count, kept_columns = create_consensus_with_ambiguity(
        final_records, threshold=0.7, gap_fraction_cutoff=0.4
    )

    alignment_length = len(aligned_qualities[0][1]) if aligned_qualities else 0
    kept_col_set = set(kept_columns)

    consensus_full_line = []
    consensus_full_qual = []
    for i in range(alignment_length):
        if i in kept_col_set:
            idx = kept_columns.index(i)
            char = consensus_seq[idx]
            q = column_quals[i]
            q_val = int(round(q)) if q is not None else None
            consensus_full_line.append(char)
            consensus_full_qual.append(q_val)
        else:
            consensus_full_line.append('-')
            consensus_full_qual.append(None)

    # Rename sequences as Seq1, Seq2, ...
    for i in range(len(aligned_qualities)):
        seq_id, aligned_seq, q_list = aligned_qualities[i]
        aligned_qualities[i] = (f"Seq{i+1}", aligned_seq, q_list)

    # Create top consensus same as final:
    top_consensus_full_line = consensus_full_line[:]
    top_consensus_full_qual = consensus_full_qual[:]

    # Add final consensus row at bottom
    aligned_qualities.append(("Consensus", ''.join(consensus_full_line), consensus_full_qual))

    # Now determine final highlighting:
    # - N-converted cells remain highlighted.
    # - Check all cells again. If ACGT <30%, highlight. If cell is N due to conversion, highlight.
    highlight_map = set()
    final_column_base_counts = compute_column_base_counts(aligned_qualities)
    # We know "Consensus" row does not get highlighted
    non_consensus = [r for r in aligned_qualities if r[0] != "Consensus"]
    for row_i, (seq_id, seq, q_list) in enumerate(non_consensus):
        for col_i, (base, q) in enumerate(zip(seq, q_list)):
            if (row_i, col_i) in n_converted_cells:
                # remain highlighted
                highlight_map.add((row_i, col_i))
            else:
                # Check highlight conditions for ACGT bases
                if base in "ACGT":
                    counts = final_column_base_counts[col_i]
                    non_gap_total = sum(counts.values())
                    if non_gap_total > 0:
                        fraction_of_base = counts[base]/non_gap_total
                        if fraction_of_base < highlight_threshold:
                            highlight_map.add((row_i, col_i))
                # If base is N or '-' but not n_converted_cells, no highlight unless it was from Q<12

    # Write HTML
    write_html(
        aligned_qualities, column_quals, consensus_full_line, consensus_full_qual, highlight_counts,
        top_consensus_full_line, top_consensus_full_qual,
        highlight_map,
        "alignment.html"
    )

    # Write consensus fasta
    if kept_columns:
        consensus_filename = f"consensus-{ambiguity_count}ambig.fasta"
        with open(consensus_filename, "w", encoding="utf-8") as cf:
            cf.write(f">Consensus\n{consensus_seq}\n")
        print(f"Consensus FASTA file written to {consensus_filename}")
    else:
        print("No columns passed the gap threshold. No consensus fasta generated.")

    print("Done.")
