# Q-Score Consensus Visualizer

# UNDER DEVELOPMENT 

This repository contains a Python script that:
1. Reads a FASTQ file from the current directory.
2. Runs a MUSCLE alignment on all sequences.
3. Computes a consensus sequence, excluding columns with >60% gaps from the final consensus and applying IUPAC ambiguity codes where applicable.
4. Uses a 5% occurrence filter to exclude rare bases (≤5% frequency) from ambiguity calculations.
5. Creates an HTML visualization of the alignment, showing:
   - Average Q scores per column.
   - A highlight count row indicating how many cells in each column are considered outliers.
   - Each nucleotide cell colored by its Q score.
   - Outlier cells (those that differ from ≥70% of the column) outlined in yellow.
   - Converts bases with quality scores ≤ 12 and are highlighted to 'N' for robust consensus generation.
   - A consensus row at the top and bottom with `'-'` indicating excluded columns.
6. Outputs a consensus FASTA file named `consensus-{X}ambig.fasta`, where `{X}` is the count of ambiguous nucleotides in the consensus.

**Key Features:**
- **Column Filtering:** Columns with more than 60% gaps are not included in the final consensus FASTA (they appear as `'-'` in the HTML consensus row).
- **Ambiguity Codes:** Ambiguity codes are computed using only bases that appear in >5% of the reads in that column.
- **Quality Score Coloring:** Each base cell is colored based on its phred quality score, ranging from red (low) to green (high).
- **Outlier Highlighting:** Cells that differ from ≥70% of the column are outlined in yellow.
- **Renumbering Sequences:** Original sequence IDs are replaced by `Seq1`, `Seq2`, etc. for readability.

## Requirements

- Python 3.7+  
- [Biopython](https://biopython.org/) (`pip install biopython`)
- [MUSCLE](https://www.drive5.com/muscle/) installed and available on your PATH
- `matplotlib` is not required since we removed image generation in the final version.

## Usage

1. Place the script and your FASTQ file in the same directory.
2. Ensure MUSCLE is installed and accessible (e.g., `muscle` command works in your terminal).
3. Run:
   ```bash
   python consensus-visual.py
The script will:
Find the FASTQ file in the current directory.
Run MUSCLE to align the sequences.
Compute column averages, consensus sequence, and filters columns.
Produce alignment.html which shows the alignment with color-coded Q scores and highlight outlines.
Produce a consensus-{X}ambig.fasta file if columns are kept. If no columns are kept, no consensus FASTA is produced.
Output Files
alignment.html:
Displays a table with:

Top row: Average Q score per column.
Second row: Highlight counts per column.
Subsequent rows: Sequences labeled Seq1, Seq2, etc., with colored and possibly highlighted nucleotides.
Last row: Consensus sequence line showing the final consensus, with '-' for excluded columns.
consensus-{X}ambig.fasta:
A FASTA file containing the final consensus sequence, where {X} is the number of ambiguous nucleotides in the consensus.

Interpretation
Colored Cells: Green indicates high-quality bases (around Q40), and red indicates low quality (close to Q0).
Yellow Outlines: If a cell is outlined in yellow, that base is uncommon in that column, differing from ≥70% of the column.
'-' in Consensus: Columns that were gap-heavy (>40% gaps) are not included in the final consensus FASTA. These appear as '-' in the HTML visualization's consensus row.
Troubleshooting
Unicode Errors:
The script writes using UTF-8 encoding. If you encounter any encoding errors, ensure your environment supports UTF-8.
Index Errors / Mismatches:
Make sure MUSCLE runs correctly and produces a proper alignment. All sequences in the alignment should have the same length.
No Consensus Columns:
If all columns are gap-heavy, no columns are kept. In that case, the consensus row is empty and no consensus.fasta is created.

