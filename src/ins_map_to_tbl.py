#!/usr/bin/env python3
""""
Summarize INS map hits using the RepBase classification hierarchy
Prerequsite: have src/repbase_repeat_classes.tsv, 2026-10-12 version has the epbase_repeat_classes.tsv from RepBase v31.09, if RepBase change their classification hierarchy, you cann use extract_repbase_classes.py to generate new tsv file.
How this works?
It applied RepeatMasker ProcessRepeats.txt to deal with overlapping region.



Execute in terminal:
(
   PROJECT_DIR=/home/luan/t2t_data/snkm_x_carp
   SAMPLE_DIR="${PROJECT_DIR}/output/(after_changed)_test_input"
   OUTPUT_DIR="${PROJECT_DIR}/output/annotate_genome_ins_test"
   LOG_DIR="${OUTPUT_DIR}/logs"

   SCRIPT_FILE="${PROJECT_DIR}/src/ins_map_to_tbl.py"
   MAP_FILE="${OUTPUT_DIR}/test_input.ins.map"
   GENOME_FILE="${SAMPLE_DIR}/03_seqer/test_input.mfa"
   REPBASE_FILE="${PROJECT_DIR}/lib/Vertebrates_RepBase31.09.fa"
   TE_LIBRARY_FILE="${SAMPLE_DIR}/11_generate_annotated_repeat_lib/test_input_Denovo_TE_Library.fasta"
   CLASSES_FILE="${PROJECT_DIR}/src/repbase_repeat_classes.tsv"
   OUTPUT_FILE="${OUTPUT_DIR}/test_input.mfa.tbl"
   LOG_FILE="${LOG_DIR}/test_input.tbl.log"

   mkdir -p "${LOG_DIR}" || exit 1
   exec > >(tee "${LOG_FILE}") 2>&1

   echo Running on `hostname`
   echo Start time is `date`
   trap 'echo End time is `date`' EXIT

   python3 "${SCRIPT_FILE}" \
      --map "${MAP_FILE}" \
      --genome "${GENOME_FILE}" \
      --repbase "${REPBASE_FILE}" \
      --te_library "${TE_LIBRARY_FILE}" \
      --classes "${CLASSES_FILE}" \
      --output "${OUTPUT_FILE}"
)
"""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from Bio import SeqIO


def read_classes(path):
   """Read the TSV at path and return rows and label-to-ancestor paths"""
   rows = []
   paths = {}

   with path.open(encoding="utf-8", newline="") as handle:
      reader = csv.DictReader(handle, delimiter="\t")
      for row in reader:
         label = row["Repeat_label"]
         rows.append(row)
         paths[label] = row["Full_path"].split(" > ")

   if not rows:
      raise ValueError("[ERROR] Classification TSV contains no data rows")

   return rows, paths


def read_repbase(path):
   """Read reference IDs and classification labels from the RepBase FASTA at path"""
   labels = {}

   with path.open(encoding="utf-8") as handle:
      for line in handle:
         if not line.startswith(">"):
            continue

         # Read the actual tab separators used in the current RepBase FASTA
         header = line[1:].rstrip()
         fields = header.split("\t")
         sequence_id = fields[0].split()[0]
         if len(fields) > 1:
            label = fields[1].strip() or None
         else:
            label = None

         if sequence_id in labels and labels[sequence_id] != label:
            raise ValueError(f"[ERROR] Conflicting RepBase labels: {sequence_id}")

         labels[sequence_id] = label

   return labels


def get_labels_from_TE_library_header(path, labels):
   """Read the TE FASTA at path and add inherited classes to the labels dictionary"""
   with path.open(encoding="utf-8") as handle:
      for line in handle:
         if not line.startswith(">"):
            continue

         # Keep the ID before whitespace, excluding the annotation description
         sequence_id = line[1:].split()[0]
         label = None

         if sequence_id.endswith("#Unclassified"):
            label = "Unclassified"
         elif "#" not in sequence_id and ":" in sequence_id:
            target = sequence_id.split(":", 1)[1]
            label = labels.get(target)

         if sequence_id in labels and labels[sequence_id] != label:
            raise ValueError(f"[ERROR] Conflicting reference labels: {sequence_id}")

         labels[sequence_id] = label


def get_genome_lengths(path):
   """Return ID-to-length mappings from the genome FASTA at path, including N bases"""
   lengths = {}

   with path.open(encoding="utf-8") as handle:
      for record in SeqIO.parse(handle, "fasta"):
         if record.id in lengths:
            raise ValueError(f"[ERROR] Duplicate genome ID: {record.id}")
         lengths[record.id] = len(record.seq)

   if not sum(lengths.values()):
      raise ValueError("[ERROR] Genome FASTA has no sequence bases")

   return lengths


def union_bp(intervals_by_contig):
   """Count the union of ascending, 1-based inclusive intervals per sequence"""
   total = 0

   for intervals in intervals_by_contig.values():
      if not intervals:
         continue

      ordered = sorted(intervals)
      left, right = ordered[0]

      for start, end in ordered[1:]:
         if start <= right + 1:
            right = max(right, end)
         else:
            total += right - left + 1
            left, right = start, end

      total += right - left + 1

   return total


def main():
   parser = argparse.ArgumentParser(description="Summuary INS .map file to .tbl")
   parser.add_argument("--map", required=True, type=Path, help="The .map file from INS")
   parser.add_argument("--genome", required=True, type=Path, help="The INS query genome FASTA")
   parser.add_argument("--repbase", required=True, type=Path, help="The RepBase FASTA used by INS")
   parser.add_argument("--te_library", required=True, type=Path, help="The CARP de novo TE FASTA")
   parser.add_argument("--classes", required=True, type=Path, help="The RepBase hierarchy .tsv file")
   parser.add_argument("--output", required=True, type=Path, help="The output .tbl file")
   args = parser.parse_args()

   for input_path in (
      args.map, args.genome, args.repbase, args.te_library, args.classes,
   ):
      if args.output.resolve() == input_path.resolve():
         raise ValueError("[ERROR] Output must differ from every input")

   class_rows, class_paths = read_classes(args.classes)
   reference_labels = read_repbase(args.repbase)
   get_labels_from_TE_library_header(args.te_library, reference_labels)
   genome_lengths = get_genome_lengths(args.genome)
   genome_bp = sum(genome_lengths.values())

   all_intervals = defaultdict(list)

   # A new category automatically receives a sequence-to-interval-list dictionary
   category_intervals = defaultdict(lambda: defaultdict(list))
   hit_counts = Counter()
   unmatched_labels = Counter()
   total_hits = 0
   map_hits = []

   with args.map.open(encoding="utf-8") as handle:
      for line_number, line in enumerate(handle, start=1):
         fields = line.split()
         if not fields:
            continue
         if len(fields) != 10:
            raise ValueError(f"[ERROR] Expected 10 map columns at line {line_number}")
         # Example of .ins.map file first row
         # test_contig_1	134	158	ConCri-5.22	137	161	d	0.92	NA	37.603
         contig = fields[0]
         start = int(fields[1])
         end = int(fields[2])
         reference_id = fields[3]
         score = float(fields[9])

         if contig not in genome_lengths:
            raise ValueError(f"[ERROR] Genome ID not found: {contig}")
         if not 1 <= start <= end <= genome_lengths[contig]:
            raise ValueError(f"[ERROR] Invalid genome interval at line {line_number}")
         if reference_id not in reference_labels:
            raise ValueError(f"[ERROR] Reference ID not found: {reference_id}")

         interval = (start, end)
         all_intervals[contig].append(interval)
         total_hits += 1
         label = reference_labels[reference_id]

         # Only explicit Unknown or Unclassified labels enter Unclassified
         if isinstance(label, str) and label.split("/", 1)[0].casefold() in (
            "unknown", "unclassified",
         ):
            nodes = ["Unclassified"]
         elif label in class_paths:
            nodes = class_paths[label]
         else:
            raw_label = label if label is not None else "Missing_header_classification"
            unmatched_labels[raw_label] += 1
            nodes = [
               "Unmatched_RepBase_class",
               f"Unmatched_RepBase_class > {raw_label}",
            ]

         map_hits.append((contig, start, end, score, nodes))
         for node in nodes:
            hit_counts[node] += 1
   
   ##############################
   # Handle overlapping region ##
   ##############################
   
   
   # Follow the approach used in ProcessRepeats.txt to handle overlapping regions
   # Sorting column Sequence and start ascend; end and BLAST bit score descend
   """
   | Item | ProcessRepeats code you provided |
   |---|---|
   | Overlapping hits within the same class | Subtracts coverage already counted |
   | Overlapping hits from different classes | Assigns newly covered bases to the class of the current annotation |
   | Fully contained hits | May contribute 0 bp but still contribute to the element count |
   | Element counting | Checks `getLeftLinkedHit()` to avoid counting linked fragments as separate elements |
   | Classification rules | Groups annotations using regular expressions applied to class names |
   | `Unclassified` | Includes `Unknown`, `Composite`, and `Other`; classes matching no summary rule also enter `OTHER` |
   | Percentage denominator | Supports `-excln` |
   | Preprocessing before summarization | Includes multiple rounds of trimming, joining, and reclassification |
   """
   # ====================
   # Compare sequence IDs
   # ====================
   # If equal, compare starts in ascending order
   # If starts are also equal, compare ends in descending order
   # Compare scores only when all preceding values are equal
   map_hits.sort(key=lambda hit: (hit[0], hit[1], -hit[2], -hit[3]))

   past_contig = None
   past_hit_start = 0
   past_hit_end = 0
   past_nodes = None

   for contig, current_hit_start, current_hit_end, score, nodes in map_hits:
      if contig != past_contig:
         past_contig = contig
         past_hit_start = 0
         past_hit_end = 0
         past_nodes = None
      # =================================================================
      # Condition 1: overlapping hits within the same class
      # =================================================================
      # The current hit starts inside the previous hit and has the same class
      """
         Condition 1-1: overlapping hits within the same class, and past hit contain the current hit
               |----------------------------------------------------------------|
         past_hit_start                                                   past_hit_end
                                 |--------------------------------|
                           current_hit_start               current_hit_end
         Coverage assignment: Take the larger region which is the past hit in here, and skip the rest of this loop iteration and process the next hit
      """
      if (
         past_hit_start <= current_hit_start <= past_hit_end
         and nodes == past_nodes
      ):
         if current_hit_end <= past_hit_end:
            continue
      
         """
         Condition 1-2: overlapping hits within the same class, and past hit and current hit is partial overlap
               |--------------------------------------|
               10                                    100
         past_hit_start                          past_hit_end
                                 |----------------------------------------------|
                                 80                                             150
                           current_hit_start               current_hit_end
         Coverage assignment: assigned_start will be 101 which is past_hit_end + 1, and assigned_interval = (101, 150)
         """ 
         assigned_start = past_hit_end + 1

      # Condition 2: overlapping hits from different classes
      
         """
         Condition 2-1: overlapping hits from different classes, and past hit contain the current hit
         Previous hit: L1
                  |---------------------------------------------|
                  10                                            100
            past_hit_start                                 past_hit_end

         Current hit: SINE
                              |-------------------------|
                              30                         80
                     current_hit_start           current_hit_end

         Coverage assignment:
                  |---------------------------------------------|
                  10                                            100
                  <---------------- L1: 91 bp ------------------>

                              |-------------------------|
                              30                         80
                              Fully contained in previous coverage
                              SINE adds 0 bp
                              Hit record count is retained
         Condition 2-2: overlapping hits from different classes, and past hit and current hit is partial overlap
         Previous hit: L1
                  |-----------------------------------|
                  10                                  100
            past_hit_start                       past_hit_end

         Current hit: SINE
                                    |---------------------------|
                                    80                          150
                           current_hit_start             current_hit_end

         Coverage assignment:
                     |-----------------------------------| |---------------|
                     10                                  100 101           150
                     <--------- L1: 91 bp ---------------> <─ SINE: 50 bp ─>
         """
      elif (
               past_hit_start <= current_hit_start <= past_hit_end
               and nodes != past_nodes
            ):
               if current_hit_end <= past_hit_end:
                  continue

               assigned_start = past_hit_end + 1

      # Condition 3: no overlap with previously counted coverage
      else:
         assigned_start = current_hit_start

      assigned_interval = (assigned_start, current_hit_end)

      # Add the assigned segment to the current class and all its ancestors
      for node in nodes:
         category_intervals[node][contig].append(assigned_interval)

      # A contained hit skips these updates and preserves the earlier boundary
      past_hit_start = current_hit_start
      past_hit_end = current_hit_end
      past_nodes = nodes
   
   ########################
   # Coverage calculation #
   ########################
   covered_bp = union_bp(all_intervals)
   covered_percent = 100.0 * covered_bp / genome_bp
   table_rows = []

   for row in class_rows:
      label = row["Repeat_label"]
      depth = len(class_paths[label]) - 1
      display_name = "   " * depth + label
      bases = union_bp(category_intervals[label])
      table_rows.append((display_name, hit_counts[label], bases, depth))

   unclassified_bp = union_bp(category_intervals["Unclassified"])
   table_rows.append((
      "Unclassified", hit_counts["Unclassified"], unclassified_bp, 0,
   ))

   group = "Unmatched_RepBase_class"
   unmatched_bp = union_bp(category_intervals[group])
   table_rows.append((
      "Unmatched RepBase class", hit_counts[group], unmatched_bp, 0,
   ))
   # Iterate over unmatched label names in alphabetical order
   # Example: BDDF1, LTR, Repeat region
   for label in sorted(unmatched_labels):
      key = f"{group} > {label}"
      bases = union_bp(category_intervals[key])
      table_rows.append(("   " + label, unmatched_labels[label], bases, 1))

   label_width = 28
   for name, count, bases, depth in table_rows:
      label_width = max(label_width, len(name))
   table_width = label_width + 12 + 14 + 12 + 3

   with args.output.open("w", encoding="utf-8") as output:
      output.write(
         "INS map classification summary\n"
         f"Sequence records: {len(genome_lengths)}\n"
         f"Input genome length: {genome_bp} bp\n\n"
      )
      output.write("=" * table_width + "\n")
      output.write(
         f"{'':<{label_width}} {'number of':>12} "
         f"{'length':>14} {'percentage':>12}\n"
         f"{'Repeat category':<{label_width}} {'hit records*':>12} "
         f"{'occupied':>14} {'of sequence':>12}\n"
      )
      output.write("-" * table_width + "\n")

      for index, (name, count, bases, depth) in enumerate(table_rows):
         if index > 0 and (
            depth == 0 or (depth == 1 and name.strip() in class_paths)
         ):
            output.write("\n")

         percentage = 100.0 * bases / genome_bp
         output.write(
            f"{name:<{label_width}} {count:>12d} "
            f"{bases:>11d} bp {percentage:>10.2f} %\n"
         )

      output.write("\n" + "-" * table_width + "\n")
      output.write(
         f"{'Total annotated (union)':<{label_width}} {total_hits:>12d} "
         f"{covered_bp:>11d} bp {covered_percent:>10.2f} %\n"
      )
      output.write("=" * table_width + "\n")
      output.write(
         "* Hit records are alignment records, not repeat copies\n"
         "Parent rows include descendants\n"
         "Genomic overlap is removed before category totals are calculated\n"
         "Do not add parent and descendant rows together\n"
         "Coverage is calculated from map intervals, not a masked FASTA\n"
      )

   print(f"Genome bases: {genome_bp}")
   print(f"Hit records: {total_hits}")
   print(f"Annotated union bases: {covered_bp}")
   print(f"Annotated union percent: {covered_percent:.4f}")
   print(f"Output: {args.output}")


if __name__ == "__main__":
   main()