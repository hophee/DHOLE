#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# --genome
# --genome_annotation
# --gene_name
# --target_plasmid
# --output_dir
# --cas_plasmid (optional; reserved for downstream use)
# --annotation_format (optional: bakta, gff; default: bakta)

args <- commandArgs(trailingOnly = TRUE)
script_path <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
legacy_dir <- dirname(normalizePath(script_path))
source(file.path(legacy_dir, "pipeline_utils.R"))
input <- parse_named_args(
  args,
  required = c("genome", "genome_annotation", "gene_name", "target_plasmid", "output_dir"),
  optional = c("cas_plasmid", "annotation_format")
)

result_dir <- input$output_dir
genome_path <- input$genome
pTarget_path <- input$target_plasmid
cas_plasmid_path <- input$cas_plasmid
annotation_format <- if (is.null(input$annotation_format)) "bakta" else input$annotation_format

genome_table <- read_genome_annotation(input$genome_annotation, annotation_format)
g_idx <- find_target_feature(genome_table, input$gene_name)

# Загружаем функцию callPrimer3 из GitHub
# suppressMessages(devtools::source_url("https://gist.githubusercontent.com/IdoBar/5e78ae7a5cc7277a04b126ce6f595d6e/raw/45c60662f3479f41765bce839835c4988a7e5b36/callPrimer3.R"))
source(file.path(dirname(legacy_dir), "callPrimer3.R"))

#
is_range_contained <- function(range1, range2) {
  start1 <- min(range1)
  end1 <- max(range1)
  start2 <- min(range2)
  end2 <- max(range2)

  return(start2 <= start1 && end1 <= end2)
}
# check_primer_table <- function(primer_table) {
#   if (is.null(primer_table) || !is.data.frame(primer_table) || nrow(primer_table) == 0) {
#     warning("Пропускаем текущие границы вырезания — нет праймеров для правого плеча")
#     cut_region_range <- cut_region_range + c(1, -1)
#
#     if (!is_range_contained(gRNA_ranges, cut_region_range) & exists("ha_primers")) {
#       cat("Не удалось подобрать оптимальную пару праймеров без выхода за границы гена!\n")
#       break
#     }
#     if ((((diff(cut_region_range) + 1) / gene_length) < 0.3) & exists("ha_primers")) {
#       cat("Не удалось подобрать оптимальную пару праймеров без выхода за границы гена!\n")
#       break
#     }
#     next
#   }
# }

# Загружаем геном
genome <- readDNAStringSet(genome_path, nrec = 1)[[1]]
genome_length <- length(genome)

# Подготовка данных из аннотации
on_plus_strand <- genome_table$strand[g_idx] == "+"
gene_start <- genome_table$start[g_idx]
gene_end <- genome_table$stop[g_idx]
gene_center <- (gene_end + gene_start) %/% 2
gene_name <- genome_table$gene[g_idx]
gene_length <- (gene_end - gene_start + 1)
cat(paste0("Название гена: ", gene_name, "\n"))
cat(paste0("Длина гена: ", gene_length, "\n"))
cat(paste0("Границы гена: ", gene_start, "-", gene_end, "\n"))

nearest_genes <- genome_table %>%
  slice(g_idx - 1, g_idx + 1) %>%
  select(start, stop) %>%
  as.matrix() %>%
  sort()
nearest_genes <- nearest_genes[2:3]

# Читаем таблицу с результатами ChopChop
gRNA_table_path <- file.path(result_dir, "n20_table.tsv")
gRNA_table <- read_tsv(gRNA_table_path, show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  mutate(genomic_start = as.numeric(gsub("contig_\\d:", "", genomic_location))) %>%
  mutate(genomic_end = genomic_start + 23 - 1) %>%
  mutate(n20_mid = (genomic_start + genomic_end) %/% 2) %>%
  mutate(mid_closeness = (abs(gene_center - n20_mid)) / gene_length) %>%
  filter(self_complementarity == 0) %>%
  filter(mm0 == 0) # %>%
# filter(mid_closeness <= 0.18)

# Проверяем, что таблица не пустая
if (nrow(gRNA_table) == 0) {
  stop("Таблица с результатами ChopChop пустая", call. = FALSE)
}

# Выбираем лучшие gRNA, чтобы они были не на одной стрэнде
plus_strand <- gRNA_table %>%
  filter(strand == "+") %>%
  slice_head(n = 2)
minus_strand <- gRNA_table %>%
  filter(strand == "-") %>%
  slice_head(n = 2)
best_gRNA <- bind_rows(plus_strand, minus_strand) %>%
  arrange(mid_closeness) %>%
  slice_head(n = 3)
best_gRNA %>%
  select(1:11) %>%
  write_tsv(file.path(result_dir, "selected_n20_table.tsv"))
selected_gRNA_set <- best_gRNA$target_sequence %>%
  substring(1, 20) %>%
  DNAStringSet()
names(selected_gRNA_set) <- paste0(gene_name, "_n20_", 1:3)
writeXStringSet(selected_gRNA_set,
  filepath = file.path(result_dir, paste0(gene_name, "_selected_n20.fasta"))
)
writeXStringSet(selected_gRNA_set,
  filepath = file.path(result_dir, paste0(gene_name, "_selected_n20.txt"))
)

if (nrow(best_gRNA) < 3) {
  stop("Недостаточно gRNA для определения границ вырезания (менее 2)", call. = FALSE)
}

gRNA_ranges <- best_gRNA %>%
  select(genomic_start, genomic_end) %>%
  range() + c(-40, 40)

# Определяем границы вырезания и извлекаем последовательности плеч гомологии
lha_length <- 300
rha_length <- 400

#
# if ((gene_end - gene_start + 1) > 500) {
#   range_borders <- 30
#   cut_region_range <- range(best_gRNA$genomic_location) + c(-range_borders, range_borders)
#   cut_region_range[1] <- max(1, cut_region_range[1]) # genome_table$start[g_idx]
#   cut_region_range[2] <- min(genome_length, cut_region_range[2]) # genome_table$stop[g_idx]
# }else{
#   range_borders <- (gene_end - gene_start + 1)%/%4
#   cut_region_range <- c(gene_start, gene_end) + c(range_borders, -range_borders)
# }

if (gene_length <= 500) {
  gap <- round(gene_length / 10, -1) * 2
  cut_start <- gene_start - gap
  cut_end <- gene_end + gap

  cut_start <- max(cut_start, nearest_genes[1])
  cut_end <- min(cut_end, nearest_genes[2])
} else if (gene_length > 500 & gene_length < 1500) {
  side <- 250
  mid <- gene_start + gene_length %/% 2
  cut_start <- mid - side
  cut_end <- mid + side
} else {
  cut_start <- floor((gene_end - gene_start) / 3) + gene_start
  cut_end <- gene_end - floor((gene_end - gene_start) / 3)
}
cut_region_range <- c(cut_start, cut_end)
gap <- 50
if (on_plus_strand) {
  left_arm_start <- max(1, cut_region_range[1] - (lha_length + gap))
  left_arm_end <- cut_region_range[1] - 1
  right_arm_start <- cut_region_range[2] + 1
  right_arm_end <- min(genome_length, cut_region_range[2] + (rha_length + gap))
} else {
  right_arm_start <- max(1, cut_region_range[1] - (rha_length + gap))
  right_arm_end <- cut_region_range[1] - 1
  left_arm_start <- cut_region_range[2] + 1
  left_arm_end <- min(genome_length, cut_region_range[2] + (lha_length + gap))
}
ticker <- 0
repeat{
  ticker <- ticker + 1
  if (on_plus_strand) {
    # left_arm_end <- cut_region_range[1] - 1
    # right_arm_start <- cut_region_range[2] + 1

    left_arm_end <- max((cut_region_range[1] - 1), nearest_genes[1])
    right_arm_start <- min((cut_region_range[2] + 1), nearest_genes[2])

    left_arm <- genome[left_arm_start:left_arm_end]
    right_arm <- genome[right_arm_start:right_arm_end]
  } else {
    right_arm_end <- max((cut_region_range[1] - 1), nearest_genes[1])
    left_arm_start <- min((cut_region_range[2] + 1), nearest_genes[2])
    left_arm <- genome[left_arm_start:left_arm_end] %>%
      reverse() %>%
      complement()
    right_arm <- genome[right_arm_start:right_arm_end] %>%
      reverse() %>%
      complement()
  }


  # Проверяем, что последовательности не пустые
  if (length(left_arm) == 0) {
    stop(sprintf("Левое плечо пустое! Границы: %d-%d", left_arm_start, left_arm_end), call. = FALSE)
  }
  if (length(right_arm) == 0) {
    stop(sprintf("Правое плечо пустое! Границы: %d-%d", right_arm_start, right_arm_end), call. = FALSE)
  }

  # Создаем DNAStringSet для плеч гомологии
  homology_arms <- DNAStringSet(list(
    left_arm_seq = left_arm,
    right_arm_seq = right_arm
  ))
  names(homology_arms) <- c(
    paste0("left_arm_", left_arm_start, "_", left_arm_end),
    paste0("right_arm_", right_arm_start, "_", right_arm_end)
  )

  # Сохраняем плечи гомологии в FASTA
  homology_arms_path <- file.path(result_dir, "homology_arms_before_primer_search.fa")
  writeXStringSet(homology_arms, homology_arms_path)
  homology_arms_path <- normalizePath(homology_arms_path)
  # cat(sprintf("Последовательности плеч гомологии сохранены в: %s\n", homology_arms_path))


  # === ВЫЗОВ PRIMER3 ДЛЯ ГЕНЕРАЦИИ ПРАЙМЕРОВ ===
  # cat("\n=== ГЕНЕРАЦИЯ ПРАЙМЕРОВ С ПОМОЩЬЮ PRIMER3 ===\n")

  # Определяем пути к файлам и настройкам
  primer3_path <- "primer3/src/primer3_core"
  thermo_params_path <- "primer3/src/primer3_config/"

  min_prod_size <- min(length(right_arm), length(left_arm))
  max_prod_size <- max(length(right_arm), length(left_arm))
  product_range <- round(c(min_prod_size, max_prod_size * 1.5))

  # Создаем файл настроек для primer3
  settings_path <- file.path(result_dir, "primer3_settings.txt")
  writeLines(c(
    "Primer3 File - http://primer3.org",
    "P3_FILE_TYPE=settings",
    "",
    "PRIMER_TASK=generic",
    "PRIMER_PICK_LEFT_PRIMER=1",
    "PRIMER_PICK_RIGHT_PRIMER=1",
    "PRIMER_NUM_RETURN=10",
    "PRIMER_MIN_SIZE=18",
    "PRIMER_OPT_SIZE=21",
    "PRIMER_MAX_SIZE=27",
    "PRIMER_MIN_TM=59.0",
    "PRIMER_OPT_TM=60.0",
    "PRIMER_MAX_TM=61.0",
    "PRIMER_PAIR_MAX_DIFF_TM=8.0",
    "PRIMER_MIN_GC=40.0",
    "PRIMER_OPT_GC_PERCENT=50.0",
    "PRIMER_MAX_GC=60.0",
    "PRIMER_MAX_SELF_ANY=12.0",
    "PRIMER_MAX_SELF_END=8.0",
    "PRIMER_PAIR_MAX_COMPL_ANY=12.0",
    "PRIMER_PAIR_MAX_COMPL_END=8.0",
    "PRIMER_MAX_HAIRPIN_TH=47.0",
    "PRIMER_MAX_POLY_X=5",
    "PRIMER_MAX_NS_ACCEPTED=0",
    "PRIMER_SALT_CORRECTIONS=1",
    "PRIMER_DNA_CONC=50.0",
    "PRIMER_WT_SIZE_LT=0.5",
    "PRIMER_WT_SIZE_GT=0.5",
    "PRIMER_WT_TM_LT=0.5",
    "PRIMER_WT_TM_GT=0.5",
    "PRIMER_PAIR_WT_DIFF_TM=0.2",
    #  paste0("PRIMER_MISPRIMING_LIBRARY=", offtarget_path),
    #  "PRIMER_MAX_LIBRARY_MISPRIMING=15.0",
    #  "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=30.0",
    paste0("PRIMER_PRODUCT_SIZE_RANGE=", product_range[1], "-", product_range[2]),
    "PRIMER_EXPLAIN_FLAG=1",
    "PRIMER_FIRST_BASE_INDEX=1",
    "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0",
    "PRIMER_MAX_END_STABILITY=9.0",
    "PRIMER_PAIR_MAX_LIB_SIM=12.0",
    "PRIMER_LIBERAL_BASE=2",
    "="
  ), settings_path)


  # Проверяем существование primer3_core
  if (!file.exists(primer3_path)) {
    stop(sprintf("primer3_core не найден по пути: %s", primer3_path), call. = FALSE)
  }
  # cat("\nГенерирую праймеры для плеч гомологии...\n")

  # Вызываем функцию callPrimer3 для левого плеча
  l_l <- c(lha_length, length(left_arm))
  left_primers <- tryCatch(
    {
      callPrimer3(
        seq = as.character(left_arm),
        size_range = paste0(l_l[1], "-", l_l[2]),
        Tm = c(60.5, 61.0, 62.5),
        Tm_diff = 2.0,
        name = names(homology_arms)[1],
        primer_num = 10,
        primer3 = primer3_path,
        thermo.param = thermo_params_path,
        report = file.path(result_dir, "left_arm_report.txt"),
        settings = settings_path
      )
    },
    error = function(e) {
      warning(sprintf("Primer3 не нашёл праймеров для левого плеча: %s", e$message))
      NULL
    }
  )
  if (is.null(left_primers) || !is.data.frame(left_primers) || nrow(left_primers) == 0) {
    warning("Пропускаем текущие границы вырезания — нет праймеров для левого плеча")
    cut_region_range <- cut_region_range + c(1, -1)
    cut_region_range[1] <- min(gRNA_ranges[1], cut_region_range[1])
    cut_region_range[2] <- max(gRNA_ranges[2], cut_region_range[2])

    if (!is_range_contained(gRNA_ranges, cut_region_range) & exists("ha_primers")) {
      cat("Не удалось подобрать оптимальную пару праймеров без выхода за границы гена!\n")
      break
    }
    if ((((diff(cut_region_range) + 1) / gene_length) < 0.3) & exists("ha_primers")) {
      cat("Не удалось подобрать оптимальную пару праймеров без выхода за границы гена!\n")
      break
    }
    next
  }

  l_r <- c(rha_length, length(right_arm))
  # Вызываем функцию callPrimer3 для правого плеча
  right_primers <- tryCatch(
    {
      callPrimer3(
        seq = as.character(right_arm),
        size_range = paste0(l_r[1], "-", l_r[2]),
        Tm = c(60.5, 61.0, 62.5),
        Tm_diff = 2.0,
        name = names(homology_arms)[2],
        primer_num = 10,
        primer3 = primer3_path,
        thermo.param = thermo_params_path,
        report = file.path(result_dir, "right_arm_report.txt"),
        settings = settings_path
      )
    },
    error = function(e) {
      warning(sprintf("Primer3 не нашёл праймеров для правого плеча: %s", e$message))
      NULL
    }
  )
  if (is.null(right_primers) || !is.data.frame(right_primers) || nrow(right_primers) == 0) {
    warning("Пропускаем текущие границы вырезания — нет праймеров для правого плеча")
    cut_region_range <- cut_region_range + c(1, -1)
    cut_region_range[1] <- min(gRNA_ranges[1]-1, cut_region_range[1])
    cut_region_range[2] <- max(gRNA_ranges[2]+1, cut_region_range[2])

    if (!is_range_contained(gRNA_ranges, cut_region_range) & exists("ha_primers")) {
      cat("Не удалось подобрать оптимальную пару праймеров без выхода за границы гена!\n")
      break
    }
    if ((((diff(cut_region_range) + 1) / gene_length) < 0.3) & exists("ha_primers")) {
      cat("Не удалось подобрать оптимальную пару праймеров без выхода за границы гена!\n")
      break
    }
    next
  }

  if (on_plus_strand) {
    left_primers <- left_primers %>%
      mutate(
        genome_start = (PRIMER_LEFT_pos + left_arm_start - 1),
        genome_end = (PRIMER_RIGHT_pos + left_arm_start - 1)
      ) %>%
      filter(genome_end <= nearest_genes[1])
    right_primers <- right_primers %>%
      mutate(
        genome_start = (PRIMER_LEFT_pos + right_arm_start - 1),
        genome_end = (PRIMER_RIGHT_pos + right_arm_start - 1)
      ) %>%
      filter(genome_start >= nearest_genes[2])
  } else {
    left_primers <- left_primers %>%
      mutate(
        genome_start = left_arm_end - PRIMER_RIGHT_pos - 1,
        genome_end   = left_arm_end - PRIMER_LEFT_pos - 1
      ) %>%
      filter(genome_start >= nearest_genes[2])
    right_primers <- right_primers %>%
      mutate(
        genome_start = right_arm_end - PRIMER_RIGHT_pos + 1,
        genome_end   = right_arm_end - PRIMER_LEFT_pos + 1
      ) %>%
      filter(genome_end <= nearest_genes[2])
  }


  # cat("\nПодбираю лучшую пару праймеров...\n")
  ha_primers <- bind_rows(left_primers, right_primers)

  if (on_plus_strand) {
    lha_ends <- left_primers %>% pull(genome_end)
    rha_starts <- right_primers %>% pull(genome_start)

    ranges_matrix <- outer(rha_starts, lha_ends, "-") - 1
    inframe_indexes <- which(((ranges_matrix %% 3) == 0) & (ranges_matrix < gene_length))

    if (length(inframe_indexes) == 0) {
      # cat(paste0("Не найдены in-frame пары праймеров для гена ", gene_name))
      best_ha_primers <- bind_rows(
        left_primers[(left_primers$PRIMER_RIGHT_pos %>% which.max()), ],
        right_primers[(right_primers$PRIMER_LEFT_pos %>% which.min()), ]
      )
    } else {
      selected_index <- which.min(ranges_matrix[inframe_indexes])
      if (is.na(selected_index)) {
        stop("Ошибка при выборе минимального индекса.")
      }

      final_index <- inframe_indexes[selected_index]
      row_col_indices <- arrayInd(final_index, dim(ranges_matrix))

      best_ha_primers <- bind_rows(
        left_primers[row_col_indices[2], ],
        right_primers[row_col_indices[1], ]
      )
    }
  } else {
    lha_starts <- left_primers %>% pull(genome_start)
    rha_ends <- right_primers %>% pull(genome_end)

    ranges_matrix <- outer(lha_starts, rha_ends, "-") - 1
    inframe_indexes <- which(((ranges_matrix %% 3) == 0) & (ranges_matrix < gene_length))

    if (length(inframe_indexes) == 0) {
      # cat(paste0("Не найдены in-frame пары праймеров для гена ", gene_name))
      best_ha_primers <- bind_rows(
        left_primers[(left_primers$PRIMER_LEFT_pos %>% which.min()), ],
        right_primers[(right_primers$PRIMER_RIGHT_pos %>% which.max()), ]
      )
    } else {
      selected_index <- which.min(ranges_matrix[inframe_indexes])
      if (is.na(selected_index)) {
        stop("Ошибка при выборе минимального индекса.")
      }

      final_index <- inframe_indexes[selected_index]
      row_col_indices <- arrayInd(final_index, dim(ranges_matrix))

      best_ha_primers <- bind_rows(
        left_primers[row_col_indices[1], ],
        right_primers[row_col_indices[2], ]
      )
    }
  }

  if (nrow(best_ha_primers) != 2) {
    stop(paste0("Не удалось найти праймеры для гена:", genome_table$gene[g_idx]))
  }


  homology_arms_range <- best_ha_primers %>%
    select(genome_start, genome_end) %>%
    as.matrix()
  offtarget_range <- homology_arms_range %>% range()
  ha_ticks <- homology_arms_range %>% sort()
  if (is_range_contained(ha_ticks[2:3], c(gene_start, gene_end))) {
    break
  } else {
    cut_region_range <- cut_region_range + c(1, -1)
    if (!is_range_contained(gRNA_ranges, cut_region_range)) {
      cat("Не удалось подобрать оптимальную пару праймеров без выхода за границы гена!\n")
      break
    }
    if (((diff(cut_region_range) + 1) / gene_length) < 0.3) {
      cat("Не удалось подобрать оптимальную пару праймеров без выхода за границы гена!\n")
      break
    }
  }
  
  if (ticker>1000){
    cat("Не удалось подобрать оптимальную пару праймеров за 1000 попыток!\n")
    break
    }
}
if (((ha_ticks[3] - ha_ticks[2] - 1) %% 3) == 0) {
  in_frame <- "in-frame"
}

pontis <- "ATGACTGCCCGCAAG"
revers_pontis <- DNAString(pontis) %>%
  reverse() %>%
  complement() %>%
  as.character()
if (((ha_ticks[3] - ha_ticks[2] - 1) %% 3) == 1) {
  pontis <- substr(pontis, 1, 14)
  revers_pontis <- substr(revers_pontis, 2, 15)
  in_frame <- "изменена длина вставки - 14"
} else if (((ha_ticks[3] - ha_ticks[2] - 1) %% 3) == 2) {
  pontis <- substr(pontis, 1, 13)
  revers_pontis <- substr(revers_pontis, 3, 15)
  in_frame <- "изменена длина вставки - 13"
}

DNAStringSet(c(
  left_homology_arm = left_arm[best_ha_primers$PRIMER_LEFT_pos[1]:best_ha_primers$PRIMER_RIGHT_pos[1]],
  right_homology_arm = right_arm[best_ha_primers$PRIMER_LEFT_pos[2]:best_ha_primers$PRIMER_RIGHT_pos[2]]
)) %>%
  writeXStringSet(file.path(result_dir, "homology_arms.fasta"))

cat("\nСобираю оффтаргет...\n")

# Создаем последовательность для проверки офф-таргетов
offtarget_path <- file.path(result_dir, "offtarget.fa")
offtg1 <- genome[1:offtarget_range[1]]
offtg2 <- genome[offtarget_range[2]:length(genome)]

DNAStringSet(list(offtarget1 = offtg1, offtarget2 = offtg2)) %>% writeXStringSet(offtarget_path)

# Скрининг генома праймеры
screening_range <- offtarget_range + c(-200, 200)
screening_seq <- genome[screening_range[1]:screening_range[2]]
genome_primers <- tryCatch(
  {
    callPrimer3(
      seq = as.character(screening_seq),
      size_range = paste0(length(screening_seq) - 100, "-", length(screening_seq)),
      Tm = c(60.5, 61.0, 62.5),
      Tm_diff = 2.0,
      name = paste0("genome_screening:", screening_range[1], "-", screening_range[2]),
      primer_num = 5,
      primer3 = primer3_path,
      thermo.param = thermo_params_path,
      report = file.path(result_dir, "genome_screening_report.txt"),
      settings = settings_path
    )
  },
  error = function(e) {
    warning(sprintf("Ошибка при генерации праймеров для скрининга: %s", e$message))
    return(NA)
  }
) %>%
  mutate(genome_start = (PRIMER_LEFT_pos + screening_range[1]), genome_end = (PRIMER_RIGHT_pos + screening_range[1]))
genome_screening_f <- genome_primers[1, 3] %>% DNAString()
genome_screening_r <- genome_primers[1, 4] %>% DNAString()
genome_screening_primers <- DNAStringSet(list(
  genome_screening_f,
  genome_screening_r
))
names(genome_screening_primers) <- paste0(gene_name, "_", "scr", c("F", "R"))
genome_screening_tm_f <- genome_primers[1, 5] %>% round(1)
genome_screening_tm_r <- genome_primers[1, 6] %>% round(1)
genome_screening_product_size_before <- genome_primers[1, 26]
genome_screening_start <- genome_primers[1, 28]
genome_screening_old_end <- genome_primers[1, 29]

bind_rows(ha_primers, genome_primers) %>% write_tsv(file.path(result_dir, "primer3_table.tsv"))


cat("\nОбъединяю результаты...\n")

# Результаты
gRNAs_forward_primers <- DNAStringSet(paste0(
  "ACGACTAGT",
  substring(best_gRNA$target_sequence, 1, 20),
  "GTTTTAGAGCTAGAAATAGCAAGTTaaaataaggct"
))
gRNA_revers_primer <- DNAStringSet("AGTTGACGCTAAAAAAAGCACCGACTCGGTGCC")
names(gRNA_revers_primer) <- paste0(gene_name, "_sgR")
names(gRNAs_forward_primers) <- paste0(gene_name, "_sgF", 1:3)

left_homology_arm_primer <- DNAStringSet(c(
  paste0("AGCGTCAACT", best_ha_primers$PRIMER_LEFT_SEQUENCE[1]),
  paste0(revers_pontis, best_ha_primers$PRIMER_RIGHT_SEQUENCE[1])
))
names(left_homology_arm_primer) <- paste0(
  gene_name, "_L", c("F", "R"),
  round(c(best_ha_primers$PRIMER_LEFT_TM[1], best_ha_primers$PRIMER_RIGHT_TM[1]))
)
right_homology_arm_primer <- DNAStringSet(c(
  paste0(pontis, best_ha_primers$PRIMER_LEFT_SEQUENCE[2]),
  paste0("ACGCTGCAG", best_ha_primers$PRIMER_RIGHT_SEQUENCE[2])
))
names(right_homology_arm_primer) <- paste0(
  gene_name, "_R", c("F", "R"),
  round(c(best_ha_primers$PRIMER_LEFT_TM[2], best_ha_primers$PRIMER_RIGHT_TM[2]))
)

c(
  gRNAs_forward_primers,
  gRNA_revers_primer,
  left_homology_arm_primer,
  right_homology_arm_primer,
  genome_screening_primers
) %>%
  writeXStringSet(file.path(result_dir, "all_primers.fasta"))

primers_without_service_seq <- DNAStringSet(
  c(
    (best_gRNA$target_sequence %>% substring(1, 20)),
    best_ha_primers[1, 3],
    best_ha_primers[1, 4],
    best_ha_primers[2, 3],
    best_ha_primers[2, 4]
  )
)
names(primers_without_service_seq) <- c(
  paste0(gene_name, "_sgF", 1:3),
  names(left_homology_arm_primer),
  names(right_homology_arm_primer)
)
writeXStringSet(primers_without_service_seq,
  filepath = file.path(result_dir, "primers_without_service_sequences.fasta")
)
writeXStringSet(primers_without_service_seq,
  filepath = file.path(result_dir, "primers_without_service_sequences.txt")
)

virtual_pcr_config3 <- file.path(result_dir, "pcr_config3.conf")
writeLines(c(
  paste0("targets_path=", genome_path),
  paste0("output_path=", file.path(result_dir, "genome_offtarget_non_service_seq.txt")),
  paste0("primers_path=", file.path(result_dir, "primers_without_service_sequences.fasta")),
  "type=primer",
  #  "ShowOnlyAmplicons=true",
  "ShowPCRProducts=true",
  "ShowPrimerAlignment=true",
  "ShowPrimerAlignmentPCRproduct=false",
  "primerstatistic=true"
), virtual_pcr_config3)


virtual_pcr_config <- file.path(result_dir, "pcr_config.conf")
writeLines(c(
  paste0("targets_path=", genome_path),
  paste0("output_path=", file.path(result_dir, "genome_offtarget_check.txt")),
  paste0("primers_path=", file.path(result_dir, "all_primers.fasta")),
  "type=primer",
  #  "ShowOnlyAmplicons=true",
  "ShowPCRProducts=true",
  "ShowPrimerAlignment=true",
  "ShowPrimerAlignmentPCRproduct=false",
  "primerstatistic=true"
), virtual_pcr_config)


writeXStringSet(genome_screening_primers, file.path(result_dir, "screening_primers.fasta"))
virtual_pcr_config2 <- file.path(result_dir, "pcr_config2.conf")
writeLines(c(
  paste0("targets_path=", pTarget_path),
  paste0("output_path=", file.path(result_dir, "screening_offtarget_check.txt")),
  paste0("primers_path=", file.path(result_dir, "screening_primers.fasta")),
  "type=primer",
  #  "ShowOnlyAmplicons=true",
  "ShowPCRProducts=true",
  "ShowPrimerAlignment=true",
  "ShowPrimerAlignmentPCRproduct=false",
  "primerstatistic=true"
), virtual_pcr_config2)


# pTarget res
pTargert_seq <- readDNAStringSet(pTarget_path)[[1]]
ins_start <- regexpr("ACTAGT", pTargert_seq, ignore.case = T)[1]
ins_end <- regexpr("CTGCAG", pTargert_seq, ignore.case = T)[1]

pcr_products_gRNA <- sapply(gRNAs_forward_primers, function(x) {
  c(x, DNAString("agtccgttatcaacttgaaaaagtggcaccgagtcggtgctttttttAGCGTCAACT"))
})
names(pcr_products_gRNA) <- paste0("gRNA_F", 1:3, "_pcr_product")

starts <- best_ha_primers$PRIMER_LEFT_pos + best_ha_primers$PRIMER_LEFT_len
ends <- best_ha_primers$PRIMER_RIGHT_pos - best_ha_primers$PRIMER_RIGHT_len

pcr_product_lha <- c(
  left_homology_arm_primer[[1]],
  left_arm[starts[1]:ends[1]],
  (left_homology_arm_primer[[2]] %>% complement() %>% reverse())
)
pcr_product_rha <- c(
  right_homology_arm_primer[[1]],
  right_arm[starts[2]:ends[2]],
  (right_homology_arm_primer[[2]] %>% complement() %>% reverse())
)

concatinate_pcr_result <- sapply(pcr_products_gRNA, function(x) {
  c(x, pcr_product_lha[-(1:10)], pcr_product_rha[-(1:15)])
})
names(concatinate_pcr_result) <- paste0("concatinate_pcr_product_", 1:3)

PCR_product_sequnces <- DNAStringSet(c(pcr_products_gRNA,
  "pcr_product_lha" = pcr_product_lha,
  "pcr_product_lha" = pcr_product_rha,
  concatinate_pcr_result
))
writeXStringSet(PCR_product_sequnces, file.path(result_dir, "PCR_product_sequnces.fasta"))

pTarget_sequences <- sapply(PCR_product_sequnces[6:8], function(x) {
  c(
    pTargert_seq[1:ins_start],
    subseq(x, start = 4, end = (length(x) - 3)),
    pTargert_seq[ins_end:length(pTargert_seq)]
  )
}) %>% DNAStringSet()
names(pTarget_sequences) <- paste0("pTarget_", 1:3)

writeXStringSet(pTarget_sequences, file.path(result_dir, "pTargets.fasta"))

edited_genome <- DNAStringSet(c(edited_genome = c(
  offtg1,
  left_arm[starts[1]:ends[1]],
  right_homology_arm_primer$right_homology_arm_primer_f,
  right_arm[starts[2]:ends[2]],
  offtg2
)))
writeXStringSet(edited_genome, file.path(result_dir, "edited_genome.fasta"))

pTarget_start <- regexpr("GCATCTGTGCGGTATTTCAC", pTargert_seq)[1]
old_end <- (DNAString("TGCTTATGGAGCTGCACATG") %>%
  complement() %>%
  reverse() %>%
  regexpr(pTargert_seq))
new_end <- (DNAString("TGCTTATGGAGCTGCACATG") %>%
  complement() %>%
  reverse() %>%
  regexpr(pTarget_sequences[[1]]))
primer_screning_before <- (old_end[1] + attr(old_end, "match.length")) - pTarget_start
primer_screning_after <- (new_end[1] + attr(new_end, "match.length")) - pTarget_start

genome_screening_new_end <- genome_screening_r %>%
  complement() %>%
  reverse() %>%
  regexpr(edited_genome)
genome_screening_product_size_after <- (genome_screening_new_end[1] + attr(genome_screening_new_end, "match.length")) - genome_screening_start


report_path <- file.path(result_dir, "report.txt")
writeLines(c(
  paste0("Температура отжига прямого праймера для левого плеча гомологии=", round(best_ha_primers[1, 5], 1)),
  paste0("Температура отжига обратного праймера для левого плеча гомологии=", round(best_ha_primers[1, 6], 1)),
  paste0("Температура отжига прямого праймера для правого плеча гомологии=", round(best_ha_primers[2, 5], 1)),
  paste0("Температура отжига обратного праймера для правого плеча гомологии=", round(best_ha_primers[2, 6], 1)),
  "",
  paste0("Длина продукта скрининга плазмиды до модификации=", primer_screning_before),
  paste0("Длина продукта скрининга плазмиды после модификации=", primer_screning_after),
  paste0("Длина продукта скрининга генома до модификации=", genome_screening_product_size_before),
  paste0("Длина продукта скрининга генома после модификации=", genome_screening_product_size_after),
  "",
  paste0("Температура отжига прямого праймера для скрининга генома=", genome_screening_tm_f),
  paste0("Температура отжига обратного праймера для скрининга генома=", genome_screening_tm_r),
  paste0("in-frame stutus: ", in_frame)
), report_path)
