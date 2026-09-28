#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Требуется 2 аргумента: primer3_output, primer_results_output", call. = FALSE)
}

primer3_output <- args[1]
primer_results_output <- args[2]

# Читаем результаты primer3
lines <- readLines(primer3_output, warn = FALSE)

# Обрабатываем результаты
results <- list()
current_seq <- NULL
current_primers <- list()
current_pair <- 0
in_primer_section <- FALSE

for (line in lines) {
  line <- trimws(line)
  if (line == "" || line == "=") {
    # Сохраняем результаты для текущей последовательности
    if (!is.null(current_seq)) {
      if (length(current_primers) > 0 && !is.null(current_primers$left_seq) && !is.null(current_primers$right_seq)) {
        results[[length(results) + 1]] <- data.frame(
          sequence_id = current_seq,
          left_primer = current_primers$left_seq,
          right_primer = current_primers$right_seq,
          product_size = current_primers$product_size,
          left_tm = current_primers$left_tm,
          right_tm = current_primers$right_tm,
          left_gc = current_primers$left_gc,
          right_gc = current_primers$right_gc,
          penalty = current_primers$penalty,
          offtarget_score = current_primers$offtarget_score,
          stringsAsFactors = FALSE
        )
      }
    }
    # Сбрасываем переменные
    current_seq <- NULL
    current_primers <- list()
    current_pair <- 0
    in_primer_section <- FALSE
    next
  }

  if (startsWith(line, "SEQUENCE_ID=")) {
    current_seq <- str_replace(line, "SEQUENCE_ID=", "")
    next
  }

  if (startsWith(line, "PRIMER_PAIR_NUM_RETURNED=0") || startsWith(line, "PRIMER_ERROR=")) {
    # Нет праймеров для этой последовательности
    if (!is.null(current_seq)) {
      results[[length(results) + 1]] <- data.frame(
        sequence_id = current_seq,
        left_primer = NA,
        right_primer = NA,
        product_size = NA,
        left_tm = NA,
        right_tm = NA,
        left_gc = NA,
        right_gc = NA,
        penalty = NA,
        offtarget_score = NA,
        stringsAsFactors = FALSE
      )
    }
    next
  }

  # Для первого набора праймеров (пара 0)
  if (startsWith(line, "PRIMER_LEFT_0_SEQUENCE=")) {
    current_primers$left_seq <- str_replace(line, "PRIMER_LEFT_0_SEQUENCE=", "")
  } else if (startsWith(line, "PRIMER_RIGHT_0_SEQUENCE=")) {
    current_primers$right_seq <- str_replace(line, "PRIMER_RIGHT_0_SEQUENCE=", "")
  } else if (startsWith(line, "PRIMER_LEFT_0_TM=")) {
    current_primers$left_tm <- as.numeric(str_replace(line, "PRIMER_LEFT_0_TM=", ""))
  } else if (startsWith(line, "PRIMER_RIGHT_0_TM=")) {
    current_primers$right_tm <- as.numeric(str_replace(line, "PRIMER_RIGHT_0_TM=", ""))
  } else if (startsWith(line, "PRIMER_LEFT_0_GC_PERCENT=")) {
    current_primers$left_gc <- as.numeric(str_replace(line, "PRIMER_LEFT_0_GC_PERCENT=", ""))
  } else if (startsWith(line, "PRIMER_RIGHT_0_GC_PERCENT=")) {
    current_primers$right_gc <- as.numeric(str_replace(line, "PRIMER_RIGHT_0_GC_PERCENT=", ""))
  } else if (startsWith(line, "PRIMER_PAIR_0_PRODUCT_SIZE=")) {
    current_primers$product_size <- as.numeric(str_replace(line, "PRIMER_PAIR_0_PRODUCT_SIZE=", ""))
  } else if (startsWith(line, "PRIMER_PAIR_0_PENALTY=")) {
    current_primers$penalty <- as.numeric(str_replace(line, "PRIMER_PAIR_0_PENALTY=", ""))
  } else if (startsWith(line, "PRIMER_PAIR_0_LIBRARY_MISPRIMING=")) {
    score_str <- str_replace(line, "PRIMER_PAIR_0_LIBRARY_MISPRIMING=", "")
    # Может содержать дополнительную информацию через запятую
    score_val <- str_split(score_str, ",")[[1]][1]
    current_primers$offtarget_score <- as.numeric(score_val)
  }
}

# Преобразуем в data frame
if (length(results) > 0) {
  results_df <- bind_rows(results)

  # Добавляем информацию о специфичности
  results_df <- results_df %>%
    mutate(
      specificity = case_when(
        is.na(offtarget_score) ~ "Unknown",
        offtarget_score < 8 ~ "High",
        offtarget_score < 12 ~ "Medium",
        TRUE ~ "Low"
      ),
      left_tm_diff = ifelse(!is.na(left_tm), abs(left_tm - 60), NA),
      right_tm_diff = ifelse(!is.na(right_tm), abs(right_tm - 60), NA),
      tm_diff = ifelse(!is.na(left_tm) & !is.na(right_tm), abs(left_tm - right_tm), NA)
    ) %>%
    arrange(desc(!is.na(penalty)), penalty)

  # Сохраняем результаты
  write_tsv(results_df, primer_results_output)

  cat(sprintf("Найдено %d наборов праймеров\n", nrow(results_df)))
  print(results_df %>% select(sequence_id, left_primer, right_primer, product_size, specificity))
} else {
  warning("Не найдено праймеров для последовательностей")
  # Создаем пустой файл с заголовками
  empty_df <- data.frame(
    sequence_id = character(0),
    left_primer = character(0),
    right_primer = character(0),
    product_size = numeric(0),
    left_tm = numeric(0),
    right_tm = numeric(0),
    left_gc = numeric(0),
    right_gc = numeric(0),
    penalty = numeric(0),
    offtarget_score = numeric(0),
    specificity = character(0),
    left_tm_diff = numeric(0),
    right_tm_diff = numeric(0),
    tm_diff = numeric(0)
  )
  write_tsv(empty_df, primer_results_output)
}

cat("Результаты primer3 успешно обработаны и сохранены\n")
