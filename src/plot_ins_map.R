# 2026-10-11
# Visualize CARP consensus annotation results
# This script is called by the Snakemake pipeline
# Execute in terminal:
library(Gviz)
library(grid)

############################################
# PART 1: Define shared annotation settings
############################################

# Keep database sections in the same order in every PNG
DATABASE_ORDER <- c(
   "RepBase / our-known",
   "Swiss-Prot",
   "GB_TE",
   "Retrovirus"
)

# Define left-inclusive BLAST bit-score intervals
# Example: 50 belongs to the 50-<80 interval
SCORE_BREAKS <- c(
   0,
   40,
   50,
   80,
   200,
   Inf
)

# Display one label for each bit-score interval
SCORE_LABELS <- c(
   "<40",
   "40-<50",
   "50-<80",
   "80-<200",
   ">=200"
)

# Match each bit-score interval to one plotting color
SCORE_COLORS <- c(
   "black",
   "blue",
   "green3",
   "magenta",
   "red"
)


# Return the shared column structure used by every annotation reader
#
# Every reader must return these columns in this exact order
# This allows data from different databases to be combined with rbind()
empty_annotation_data <- function() {
   data.frame(
      Consensus = character(),
      query_from = numeric(),
      query_to = numeric(),
      reference = character(),
      reference_from = numeric(),
      reference_to = numeric(),
      direction = character(),
      identity = numeric(),
      bit_score = numeric(),
      evalue = numeric(),
      database = character(),
      stringsAsFactors = FALSE
   )
}


############################################
# PART 2: Read consensus lengths from FASTA
############################################

# file_path:
#   Path to the consensus FASTA associated with the annotation results
#   Expected input is one or more DNA FASTA records
#   Sequence IDs are taken from the first word of each FASTA header
#   Example: ">family000000_consensus (2 members)"
#            becomes "family000000_consensus"
read_consensus_lengths <- function(file_path) {
   if (!file.exists(file_path)) {
      stop("[ERROR] FASTA file does not exist: ", file_path)
   }

   # Read one sequence length for each FASTA record
   # seqtype = "DNA" validates DNA sequence symbols
   # use.names = TRUE keeps each complete FASTA header as a vector name
   sequence_lengths <- tryCatch(
      Biostrings::fasta.seqlengths(
         filepath = file_path,
         seqtype = "DNA",
         use.names = TRUE
      ),
      error = function(error) {
         stop(
            "[ERROR] Could not read FASTA: ",
            conditionMessage(error),
            call. = FALSE
         )
      }
   )

   if (length(sequence_lengths) == 0L) {
      stop("[ERROR] FASTA contains no sequences")
   }

   # Remove the description after the first whitespace
   # Example: "family000000_consensus (2 members)"
   #          becomes "family000000_consensus"
   consensus_ids <- sub(
      "[[:space:]].*$",
      "",
      names(sequence_lengths)
   )

   if (any(!nzchar(consensus_ids))) {
      stop("[ERROR] FASTA contains an empty sequence ID")
   }

   duplicate_ids <- unique(
      consensus_ids[duplicated(consensus_ids)]
   )

   if (length(duplicate_ids) > 0L) {
      stop(
         "[ERROR] Duplicate FASTA ID: ",
         paste(duplicate_ids, collapse = ", ")
      )
   }

   if (any(sequence_lengths == 0L)) {
      stop("[ERROR] FASTA contains a sequence with zero length")
   }

   data.frame(
      Consensus = consensus_ids,
      length = unname(sequence_lengths),
      stringsAsFactors = FALSE
   )
}


############################################
# PART 3: Read the repeat-library map
############################################

# file_path:
#   Path to the 10-column repeat-library map
#
# Expected columns:
#   1  consensus ID
#   2  consensus start
#   3  consensus end
#   4  reference ID
#   5  reference start
#   6  reference end
#   7  direction using d or c
#   8  identity as a proportion from 0 to 1
#   9  unused value
#   10 BLAST bit score
read_map_file <- function(file_path) {
   if (!file.exists(file_path)) {
      stop("[ERROR] Map file does not exist: ", file_path)
   }

   map_lines <- readLines(
      file_path,
      warn = FALSE
   )

   map_lines <- map_lines[
      nzchar(trimws(map_lines))
   ]

   # An empty map is valid because a run may have no repeat-library hits
   if (length(map_lines) == 0L) {
      return(empty_annotation_data())
   }

   fields <- strsplit(
      trimws(map_lines),
      "[[:space:]]+"
   )

   bad_rows <- which(
      lengths(fields) != 10L
   )

   if (length(bad_rows) > 0L) {
      stop(
         "[ERROR] Expected 10 map columns at non-empty row: ",
         bad_rows[1]
      )
   }

   map_data <- read.table(
      text = paste(map_lines, collapse = "\n"),
      header = FALSE,
      colClasses = "character",
      comment.char = "",
      quote = "",
      stringsAsFactors = FALSE,
      na.strings = "NA"
   )

   names(map_data) <- c(
      "Consensus",
      "query_from",
      "query_to",
      "reference",
      "reference_from",
      "reference_to",
      "direction",
      "identity",
      "unused",
      "bit_score"
   )

   numeric_columns <- c(
      "query_from",
      "query_to",
      "reference_from",
      "reference_to",
      "identity",
      "bit_score"
   )

   for (column in numeric_columns) {
      map_data[[column]] <- suppressWarnings(
         as.numeric(map_data[[column]])
      )

      if (any(!is.finite(map_data[[column]]))) {
         stop(
            "[ERROR] Invalid numeric value in map column: ",
            column
         )
      }
   }

   coordinate_columns <- c(
      "query_from",
      "query_to",
      "reference_from",
      "reference_to"
   )

   for (column in coordinate_columns) {
      values <- map_data[[column]]

      if (any(values < 1 | values != floor(values))) {
         stop(
            "[ERROR] Map coordinates must be positive integers: ",
            column
         )
      }
   }

   if (any(map_data$query_from > map_data$query_to)) {
      stop(
         "[ERROR] A map consensus interval starts after it ends"
      )
   }

   if (any(map_data$reference_from > map_data$reference_to)) {
      stop(
         "[ERROR] A map reference interval starts after it ends"
      )
   }

   if (any(
      map_data$identity < 0 |
      map_data$identity > 1
   )) {
      stop(
         "[ERROR] Map identity must be in [0, 1]"
      )
   }

   if (any(map_data$bit_score < 0)) {
      stop(
         "[ERROR] Map bit score must be non-negative"
      )
   }

   if (any(!map_data$direction %in% c("d", "c"))) {
      stop(
         "[ERROR] Map direction must be d or c"
      )
   }

   # Convert map direction symbols into plotting strand symbols
   direction <- ifelse(
      map_data$direction == "d",
      "+",
      "-"
   )

   data.frame(
      Consensus = map_data$Consensus,
      query_from = map_data$query_from,
      query_to = map_data$query_to,
      reference = map_data$reference,
      reference_from = map_data$reference_from,
      reference_to = map_data$reference_to,
      direction = direction,
      identity = map_data$identity,
      bit_score = map_data$bit_score,

      # The repeat-library map does not contain BLAST E-values
      evalue = rep(
         NA_real_,
         nrow(map_data)
      ),

      # The map combines both repeat libraries without source provenance
      database = rep(
         "RepBase / our-known",
         nrow(map_data)
      ),
      stringsAsFactors = FALSE
   )
}


############################################
# PART 4: Read BLAST outfmt 6 annotations
############################################

# file_path:
#   Path to a standard 12-column BLAST outfmt 6 file
#   An empty file is accepted and represents no database hits
#
# database:
#   Database name displayed beside each hit
#   Expected value is Swiss-Prot, GB_TE, or Retrovirus
#   Example: database = "Swiss-Prot"
read_blast_hits <- function(
    file_path,
    database) {
   if (!file.exists(file_path)) {
      stop(
         "[ERROR] BLAST result file does not exist: ",
         file_path
      )
   }

   if (!is.character(database) ||
      length(database) != 1L ||
      is.na(database) ||
      !nzchar(database)) {
      stop(
         "[ERROR] Database name must be one non-empty string"
      )
   }

   if (!database %in% DATABASE_ORDER) {
      stop(
         "[ERROR] Unsupported database name: ",
         database
      )
   }

   blast_lines <- readLines(
      file_path,
      warn = FALSE
   )

   blast_lines <- blast_lines[
      nzchar(trimws(blast_lines))
   ]

   # An empty file means that this database produced no hits
   if (length(blast_lines) == 0L) {
      return(empty_annotation_data())
   }

   fields <- strsplit(
      trimws(blast_lines),
      "[[:space:]]+"
   )

   bad_rows <- which(
      lengths(fields) != 12L
   )

   if (length(bad_rows) > 0L) {
      stop(
         "[ERROR] Expected 12 BLAST columns in ",
         database,
         " at non-empty row: ",
         bad_rows[1]
      )
   }

   blast_data <- read.table(
      text = paste(blast_lines, collapse = "\n"),
      header = FALSE,
      colClasses = "character",
      comment.char = "",
      quote = "",
      stringsAsFactors = FALSE
   )

   names(blast_data) <- c(
      "query_id",
      "subject_id",
      "percent_identity",
      "alignment_length",
      "mismatches",
      "gap_opens",
      "query_start",
      "query_end",
      "subject_start",
      "subject_end",
      "evalue",
      "bit_score"
   )

   numeric_columns <- c(
      "percent_identity",
      "alignment_length",
      "mismatches",
      "gap_opens",
      "query_start",
      "query_end",
      "subject_start",
      "subject_end",
      "evalue",
      "bit_score"
   )

   for (column in numeric_columns) {
      blast_data[[column]] <- suppressWarnings(
         as.numeric(blast_data[[column]])
      )

      if (any(!is.finite(blast_data[[column]]))) {
         stop(
            "[ERROR] Invalid numeric value in ",
            database,
            " column: ",
            column
         )
      }
   }

   coordinate_columns <- c(
      "query_start",
      "query_end",
      "subject_start",
      "subject_end"
   )

   for (column in coordinate_columns) {
      values <- blast_data[[column]]

      if (any(values < 1 | values != floor(values))) {
         stop(
            "[ERROR] Coordinates must be positive integers in ",
            database,
            " column: ",
            column
         )
      }
   }

   if (any(
      blast_data$percent_identity < 0 |
      blast_data$percent_identity > 100
   )) {
      stop(
         "[ERROR] BLAST percent identity must be in [0, 100]"
      )
   }

   if (any(
      blast_data$alignment_length < 1 |
      blast_data$alignment_length !=
         floor(blast_data$alignment_length)
   )) {
      stop(
         "[ERROR] BLAST alignment length must be a positive integer"
      )
   }

   if (any(
      blast_data$mismatches < 0 |
      blast_data$mismatches != floor(blast_data$mismatches)
   )) {
      stop(
         "[ERROR] BLAST mismatch count must be a non-negative integer"
      )
   }

   if (any(
      blast_data$gap_opens < 0 |
      blast_data$gap_opens != floor(blast_data$gap_opens)
   )) {
      stop(
         "[ERROR] BLAST gap-open count must be a non-negative integer"
      )
   }

   if (any(blast_data$evalue < 0)) {
      stop(
         "[ERROR] BLAST E-value must be non-negative"
      )
   }

   if (any(blast_data$bit_score < 0)) {
      stop(
         "[ERROR] BLAST bit score must be non-negative"
      )
   }

   if (any(!nzchar(blast_data$query_id))) {
      stop(
         "[ERROR] BLAST contains an empty query ID"
      )
   }

   if (any(!nzchar(blast_data$subject_id))) {
      stop(
         "[ERROR] BLAST contains an empty subject ID"
      )
   }

   # Descending query coordinates indicate a reverse-strand alignment
   direction <- ifelse(
      blast_data$query_start <= blast_data$query_end,
      "+",
      "-"
   )

   data.frame(
      Consensus = blast_data$query_id,
      query_from = pmin(
         blast_data$query_start,
         blast_data$query_end
      ),
      query_to = pmax(
         blast_data$query_start,
         blast_data$query_end
      ),
      reference = blast_data$subject_id,
      reference_from = pmin(
         blast_data$subject_start,
         blast_data$subject_end
      ),
      reference_to = pmax(
         blast_data$subject_start,
         blast_data$subject_end
      ),
      direction = direction,

      # Convert a percentage into the map identity scale
      # Example: 88.333 becomes 0.88333
      identity = blast_data$percent_identity / 100,
      bit_score = blast_data$bit_score,
      evalue = blast_data$evalue,
      database = rep(
         database,
         nrow(blast_data)
      ),
      stringsAsFactors = FALSE
   )
}


############################################
# PART 5: Build one Gviz track for one hit
############################################

# hit:
#   One row from the combined annotation data frame
#   The row describes one alignment on one consensus
#
# hit_color:
#   Color selected from the BLAST bit-score interval
#   Example: a bit score of 250 uses red
create_hit_track <- function(
    hit,
    hit_color) {
   database <- hit$database[[1]]
   reference <- hit$reference[[1]]
   identity <- hit$identity[[1]]
   bit_score <- hit$bit_score[[1]]
   evalue <- hit$evalue[[1]]

   if (is.na(evalue)) {
      statistics_text <- sprintf(
         "Bit score: %.1f | Identity: %.1f%%",
         bit_score,
         100 * identity
      )
   } else {
      statistics_text <- sprintf(
         "Bit score: %.1f | Identity: %.1f%% | E-value: %.2g",
         bit_score,
         100 * identity,
         evalue
      )
   }

   # Example:
   # Swiss-Prot | sp|P11369|LORF2_MOUSE
   # Bit score: 119.0 | Identity: 61.6% | E-value: 9.48e-30
   track_name <- sprintf(
      "%s | %s\n%s",
      database,
      reference,
      statistics_text
   )

   # Use dark text on bright green and magenta bars
   text_color <- if (
      hit_color %in% c("green3", "magenta")
   ) {
      "black"
   } else {
      "white"
   }

   AnnotationTrack(
      name = "",
      start = hit$query_from[[1]],
      end = hit$query_to[[1]],
      strand = hit$direction[[1]],
      chromosome = "chrNA",
      id = track_name,
      group = track_name,
      showFeatureId = FALSE,
      showId = TRUE,
      just.group = "below",
      fontsize.group = 12,
      fontcolor.group = "black",
      rotation.group = 0,
      showTitle = FALSE,
      stacking = "squish",
      collapse = FALSE,
      shape = "fixedArrow",
      arrowHeadWidth = 10,
      stackHeight = 0.20,
      fill = hit_color,
      col = hit_color,
      fontsize.item = 9,
      fontcolor.item = text_color,
      rotation = 0
   )
}


############################################
# PART 6: Generate one PNG per consensus
############################################

# data:
#   Combined annotation data returned by read_map_file()
#   and read_blast_hits()
#
# consensus_lengths:
#   Data frame returned by read_consensus_lengths()
#   Expected columns are Consensus and length
#
# output_dir:
#   Directory for PNG files and plot_manifest.tsv
#   The directory is created if it does not already exist
#
# width_in:
#   PNG width in inches
#   Default 14 gives long annotation labels more horizontal space
#
# dpi:
#   PNG resolution in pixels per inch
#   Default 300 produces report-ready images
generate_plots <- function(
    data,
    consensus_lengths,
    output_dir,
    width_in = 14,
    dpi = 300) {
   if (!is.numeric(width_in) ||
      length(width_in) != 1L ||
      !is.finite(width_in) ||
      width_in <= 0) {
      stop(
         "[ERROR] width_in must be one positive number"
      )
   }

   if (!is.numeric(dpi) ||
      length(dpi) != 1L ||
      !is.finite(dpi) ||
      dpi <= 0) {
      stop(
         "[ERROR] dpi must be one positive number"
      )
   }

   required_columns <- names(
      empty_annotation_data()
   )

   missing_columns <- setdiff(
      required_columns,
      names(data)
   )

   if (length(missing_columns) > 0L) {
      stop(
         "[ERROR] Annotation data is missing columns: ",
         paste(missing_columns, collapse = ", ")
      )
   }

   unknown_databases <- setdiff(
      unique(data$database),
      DATABASE_ORDER
   )

   if (length(unknown_databases) > 0L) {
      stop(
         "[ERROR] Unsupported databases: ",
         paste(unknown_databases, collapse = ", ")
      )
   }

   unknown_ids <- setdiff(
      unique(data$Consensus),
      consensus_lengths$Consensus
   )

   if (length(unknown_ids) > 0L) {
      stop(
         "[ERROR] Annotation ID is absent from FASTA: ",
         paste(unknown_ids, collapse = ", ")
      )
   }

   # Match sequence lengths by consensus ID rather than row order
   data$length <- consensus_lengths$length[
      match(
         data$Consensus,
         consensus_lengths$Consensus
      )
   ]

   bad_rows <- which(
      data$query_from < 1 |
      data$query_from > data$query_to |
      data$query_to > data$length
   )

   if (length(bad_rows) > 0L) {
      bad_row <- bad_rows[1]

      stop(
         "[ERROR] Annotation interval is outside the consensus: ",
         data$Consensus[bad_row],
         "; interval=",
         data$query_from[bad_row],
         "-",
         data$query_to[bad_row],
         "; length=",
         data$length[bad_row]
      )
   }

   # Restrict IDs to characters that are safe in PNG filenames
   if (any(!grepl(
      "^[A-Za-z0-9][A-Za-z0-9_.-]*$",
      consensus_lengths$Consensus
   ))) {
      stop(
         "[ERROR] A consensus ID contains unsafe filename characters"
      )
   }

   dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
   )

   if (!dir.exists(output_dir)) {
      stop(
         "[ERROR] Could not create output directory: ",
         output_dir
      )
   }

   all_plots <- vector(
      "list",
      nrow(consensus_lengths)
   )

   for (
      consensus_index in seq_len(nrow(consensus_lengths))
   ) {
      consensus <- consensus_lengths$Consensus[
         consensus_index
      ]

      consensus_length <- consensus_lengths$length[
         consensus_index
      ]

      subset_data <- data[
         data$Consensus == consensus,
         ,
         drop = FALSE
      ]

      # Group hits by database and then place higher scores first
      subset_data <- subset_data[
         order(
            match(
               subset_data$database,
               DATABASE_ORDER
            ),
            -subset_data$bit_score,
            subset_data$query_from,
            subset_data$query_to,
            subset_data$reference
         ),
         ,
         drop = FALSE
      ]

      n_hits <- nrow(subset_data)

      # right = FALSE creates left-inclusive intervals
      # Example: 50 belongs to [50, 80)
      score_bins <- cut(
         subset_data$bit_score,
         breaks = SCORE_BREAKS,
         labels = SCORE_LABELS,
         right = FALSE
      )

      hit_colors <- SCORE_COLORS[
         as.integer(score_bins)
      ]

      track_list <- vector(
         "list",
         n_hits
      )

      for (hit_index in seq_len(n_hits)) {
         hit <- subset_data[
            hit_index,
            ,
            drop = FALSE
         ]

         track_list[[hit_index]] <- create_hit_track(
            hit = hit,
            hit_color = hit_colors[hit_index]
         )
      }

      # Increase image height as the number of hits increases
      height_in <- max(
         3.5,
         2.8 + 0.5 * n_hits
      )

      title_text <- sprintf(
         "%s | %d bp | %d matches",
         consensus,
         consensus_length,
         n_hits
      )

      png_name <- paste0(
         consensus,
         ".png"
      )

      png_path <- file.path(
         output_dir,
         png_name
      )

      png(
         filename = png_path,
         width = width_in,
         height = height_in,
         units = "in",
         res = dpi,
         bg = "white"
      )

      # Always close the PNG device if drawing raises an error
      tryCatch({
         grid.newpage()

         # Row 1 contains the bit-score legend
         # Row 2 contains the consensus and all annotation tracks
         pushViewport(viewport(
            layout = grid.layout(
               nrow = 2,
               ncol = 1,
               heights = unit(
                  c(
                     0.55,
                     height_in - 0.55
                  ),
                  "inches"
               )
            )
         ))

         ####################################
         # Draw the bit-score legend
         ####################################

         pushViewport(viewport(
            layout.pos.row = 1
         ))

         grid.text(
            "BLAST bit score",
            x = 0.02,
            y = 0.82,
            just = "left",
            gp = gpar(
               fontsize = 10
            )
         )

         legend_x <- seq(
            0.08,
            0.84,
            length.out = length(SCORE_LABELS)
         )

         for (
            legend_index in seq_along(SCORE_LABELS)
         ) {
            grid.rect(
               x = legend_x[legend_index],
               y = 0.29,
               width = 0.035,
               height = 0.25,
               gp = gpar(
                  fill = SCORE_COLORS[legend_index],
                  col = NA
               )
            )

            grid.text(
               SCORE_LABELS[legend_index],
               x = legend_x[legend_index] + 0.025,
               y = 0.35,
               just = "left",
               gp = gpar(
                  fontsize = 10
               )
            )
         }

         popViewport()

         ####################################
         # Draw the consensus and hit tracks
         ####################################

         pushViewport(viewport(
            layout.pos.row = 2,
            clip = "off"
         ))

         # The consensus and hit panels use the same horizontal width
         pushViewport(viewport(
            width = 0.94,
            height = 0.90,
            clip = "off",
            layout = grid.layout(
               nrow = 2,
               ncol = 1,
               heights = unit.c(
                  unit(
                     0.75,
                     "inches"
                  ),
                  unit(
                     1,
                     "null"
                  )
               )
            )
         ))

         # Add coordinate padding around both ends of the consensus
         coordinate_padding <- max(
            10,
            ceiling(
               consensus_length * 0.03
            )
         )

         view_from <- 1 - coordinate_padding
         view_to <- consensus_length +
            coordinate_padding

         ####################################
         # Draw the consensus bar
         ####################################

         pushViewport(viewport(
            layout.pos.row = 1
         ))

         pushViewport(viewport(
            xscale = c(
               view_from,
               view_to
            ),
            yscale = c(
               0,
               1
            ),
            clip = "off"
         ))

         grid.rect(
            x = unit(
               (consensus_length + 1) / 2,
               "native"
            ),
            y = 0.72,
            width = unit(
               consensus_length,
               "native"
            ),
            height = 0.32,
            gp = gpar(
               fill = "#52C3C5",
               col = NA
            )
         )

         grid.text(
            title_text,
            x = unit(
               (consensus_length + 1) / 2,
               "native"
            ),
            y = 0.72,
            gp = gpar(
               fontsize = 12,
               fontface = "bold"
            )
         )

         # Select readable coordinate tick positions
         tick_positions <- pretty(
            c(
               1,
               consensus_length
            ),
            n = 5
         )

         # Remove interior ticks that are too close to either endpoint
         tick_positions <- tick_positions[
            tick_positions >
               1 + 0.04 * consensus_length &
            tick_positions <
               consensus_length -
               0.04 * consensus_length
         ]

         # Always include the first and final consensus positions
         tick_positions <- sort(unique(c(
            1,
            tick_positions,
            consensus_length
         )))

         grid.segments(
            x0 = unit(
               tick_positions,
               "native"
            ),
            x1 = unit(
               tick_positions,
               "native"
            ),
            y0 = 0.56,
            y1 = 0.36,
            gp = gpar(
               col = "black",
               lwd = 1.5
            )
         )

         grid.text(
            format(
               tick_positions,
               trim = TRUE,
               scientific = FALSE
            ),
            x = unit(
               tick_positions,
               "native"
            ),
            y = 0.29,
            just = c(
               "center",
               "top"
            ),
            gp = gpar(
               fontsize = 12
            )
         )

         # Leave the native-coordinate and consensus-row viewports
         popViewport(2)

         ####################################
         # Draw every annotation track
         ####################################

         pushViewport(viewport(
            layout.pos.row = 2,
            clip = "off"
         ))

         if (n_hits > 0L) {
            plotTracks(
               track_list,
               sizes = rep(
                  1,
                  length(track_list)
               ),
               chromosome = "chrNA",
               from = view_from,
               to = view_to,
               panel.only = TRUE,
               showTitle = FALSE,
               background.panel = "white"
            )
         } else {
            grid.text(
               "No matches",
               gp = gpar(
                  fontsize = 12
               )
            )
         }

         # Leave the annotation-row viewport
         popViewport()

         # Leave the shared plotting and body-row viewports
         popViewport(2)

         # Leave the outer layout viewport
         popViewport()
      }, finally = {
         dev.off()
      })

      all_plots[[consensus_index]] <- data.frame(
         Consensus = consensus,
         length = consensus_length,
         n_matches = n_hits,
         png_file = png_name,
         stringsAsFactors = FALSE
      )
   }

   plot_manifest <- do.call(
      rbind,
      all_plots
   )

   # Store PNG basenames for later use in the HTML report
   write.table(
      plot_manifest,
      file = file.path(
         output_dir,
         "plot_manifest.tsv"
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
   )

   return(plot_manifest)
}


############################################
# PART 7: Read pipeline arguments and run
############################################

main <- function() {
   # Expected argument order:
   # 1 consensus FASTA
   # 2 repeat-library map
   # 3 Swiss-Prot BLAST outfmt 6
   # 4 GB_TE BLAST outfmt 6
   # 5 Retrovirus BLAST outfmt 6
   # 6 output directory
   # 7 comma-separated consensus IDs
   arguments <- commandArgs(
      trailingOnly = TRUE
   )

   if (length(arguments) != 7L) {
      stop(
         paste(
            "[ERROR] Usage:",
            "Rscript src/plot_ins_map.R",
            "<consensus_fasta>",
            "<map_file>",
            "<swissprot_blast>",
            "<gb_te_blast>",
            "<retrovirus_blast>",
            "<output_dir>",
            "<consensus_ids>"
         ),
         call. = FALSE
      )
   }

   consensus_fasta <- arguments[[1]]
   map_file <- arguments[[2]]
   swissprot_file <- arguments[[3]]
   gb_te_file <- arguments[[4]]
   retrovirus_file <- arguments[[5]]
   output_dir <- arguments[[6]]
   consensus_ids_text <- arguments[[7]]

   # Split the seventh argument at literal commas
   # Example: "family1,family2" becomes c("family1", "family2")
   selected_consensus_ids <- strsplit(
      consensus_ids_text,
      split = ",",
      fixed = TRUE
   )[[1]]

   # Remove accidental whitespace surrounding an ID
   selected_consensus_ids <- trimws(
      selected_consensus_ids
   )

   if (any(!nzchar(selected_consensus_ids))) {
      stop(
         "[ERROR] Selected consensus IDs must not be empty",
         call. = FALSE
      )
   }

   if (anyDuplicated(selected_consensus_ids)) {
      stop(
         "[ERROR] Selected consensus IDs must not contain duplicates",
         call. = FALSE
      )
   }

   consensus_lengths <- read_consensus_lengths(
      consensus_fasta
   )

   # Replace the all keyword with every consensus ID from the FASTA file
   # Example: a FASTA with 31 records produces 31 selected IDs
   if (identical(selected_consensus_ids, "all")) {
      selected_consensus_ids <- consensus_lengths$Consensus
   }

   # Find requested IDs that do not exist in the consensus FASTA
   missing_consensus_ids <- setdiff(
      selected_consensus_ids,
      consensus_lengths$Consensus
   )

   if (length(missing_consensus_ids) > 0L) {
      stop(
         "[ERROR] Selected consensus ID is absent from FASTA: ",
         paste(missing_consensus_ids, collapse = ", "),
         call. = FALSE
      )
   }

   # Keep only selected records and preserve their YAML order
   consensus_lengths <- consensus_lengths[
      match(
         selected_consensus_ids,
         consensus_lengths$Consensus
      ),
      ,
      drop = FALSE
   ]

   map_hits <- read_map_file(
      map_file
   )

   swissprot_hits <- read_blast_hits(
      file_path = swissprot_file,
      database = "Swiss-Prot"
   )

   gb_te_hits <- read_blast_hits(
      file_path = gb_te_file,
      database = "GB_TE"
   )

   retrovirus_hits <- read_blast_hits(
      file_path = retrovirus_file,
      database = "Retrovirus"
   )

   # Every reader returns the same columns in the same order
   annotation_data <- rbind(
      map_hits,
      swissprot_hits,
      gb_te_hits,
      retrovirus_hits
   )

   # Keep annotation rows belonging to selected consensus IDs
   annotation_data <- annotation_data[
      annotation_data$Consensus %in% selected_consensus_ids,
      ,
      drop = FALSE
   ]

   plot_manifest <- generate_plots(
      data = annotation_data,
      consensus_lengths = consensus_lengths,
      output_dir = output_dir,
      width_in = 14,
      dpi = 300
   )

   print(plot_manifest)
}


main()