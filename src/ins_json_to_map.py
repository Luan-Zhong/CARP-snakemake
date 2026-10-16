#!/usr/bin/env python3
"""Convert ins JSON output to the map fields consumed by CARP.

This is a CARP compatibility format, not authentic Censor output.
Column 9 is unused by the current Java readers.
Column 10 contains the BLAST bit score, not a Censor score.
"""

import argparse
import json
from pathlib import Path


def read_hits(path):
   # ins emits consecutive JSON objects, not a JSON array or JSONL.
   text = path.read_text(encoding="utf-8")
   decoder = json.JSONDecoder()
   position = 0

   while position < len(text):
      if text[position].isspace():
         position += 1
         continue

      hit, position = decoder.raw_decode(text, position)
      if not isinstance(hit, dict):
         raise ValueError("Expected a JSON object for each hit")
      yield hit


def interval(start, end):
   # Expect the integer coordinates stored by the current INS implementation
   if type(start) is not int or type(end) is not int:
      raise ValueError("[ERROR] Coordinates must be integers")

   if start <= end:
      # Forward example: start=0, end=84 represents bases 1 through 84
      left = start + 1
      right = end
   else:
      # INS subtracts one from the first BLAST endpoint before storing it
      # Reverse example: start=83, end=1 represents bases 1 through 84
      left = end
      right = start + 1

   if left < 1 or right < left:
      raise ValueError(f"[ERROR] Invalid interval: {start}, {end}")

   return left, right


def identifier(value):
   if not isinstance(value, str) or not value:
      raise ValueError("Missing sequence identifier")
   if any(character.isspace() for character in value):
      raise ValueError(f"Whitespace in sequence identifier: {value!r}")
   return value


def convert_hit(hit):
   # In ins JSON, Subject is the input consensus sequence;
   # Query is the reference-library sequence.
   consensus = identifier(hit["SubjectAccVer"])
   reference = identifier(hit["QueryAccVer"])

   start, end = interval(hit["SubjectStart"], hit["SubjectEnd"])
   ref_start, ref_end = interval(hit["QueryStart"], hit["QueryEnd"])

   strand = hit["Strand"]
   if strand not in (-1, 1):
      raise ValueError(f"Unsupported strand: {strand}")

   return [
      consensus,
      start,
      end,
      reference,
      ref_start,
      ref_end,
      "c" if strand == -1 else "d",
      float(hit["PctIdentity"]) / 100.0,
      "NA",
      float(hit["BitScore"]),
   ]


def main():
   parser = argparse.ArgumentParser(description=__doc__)
   parser.add_argument("input_json", type=Path)
   parser.add_argument("output_map", type=Path)
   args = parser.parse_args()

   if args.input_json.resolve() == args.output_map.resolve():
      parser.error("Input and output paths must be different")

   # Validate all records before opening the output file.
   rows = [convert_hit(hit) for hit in read_hits(args.input_json)]

   # Match the deterministic ordering used in the compatibility test.
   # This does not select a biologically preferred hit.
   rows.sort(key=lambda row: (row[0], row[1], row[2], row[3]))

   with args.output_map.open("w", encoding="utf-8") as output:
      for row in rows:
         output.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
   main()