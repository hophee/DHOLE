#!/usr/bin/bash
#source: https://github.com/hophee/2PAC

eval "$(conda shell.bash hook)"
conda activate oligo_design
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Использование:
  ./chopa_ncRNA.sh --genome GENOME --genome_annotation ANNOTATION --gene_name GENE \
    --target_plasmid PTARGET --output_dir OUTPUT_DIR \
    [--cas_plasmid CAS] [--annotation_format bakta|gff]
EOF
}

genome=""
genome_annotation=""
gene_name=""
target_plasmid=""
output_dir=""
cas_plasmid=""
annotation_format="bakta"

while [ "$#" -gt 0 ]; do
  case "$1" in
    --genome) genome="$2" ;;
    --genome_annotation|--genome-annotation) genome_annotation="$2" ;;
    --gene_name|--gene-name) gene_name="$2" ;;
    --target_plasmid|--target-plasmid) target_plasmid="$2" ;;
    --output_dir|--output-dir) output_dir="$2" ;;
    --cas_plasmid|--cas-plasmid) cas_plasmid="$2" ;;
    --annotation_format|--annotation-format) annotation_format="$2" ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Неизвестный аргумент: $1" >&2; usage >&2; exit 2 ;;
  esac
  if [ "$#" -lt 2 ]; then
    echo "Для аргумента $1 не задано значение" >&2
    exit 2
  fi
  shift 2
done

for required_name in genome genome_annotation gene_name target_plasmid output_dir; do
  if [ -z "${!required_name}" ]; then
    echo "Не задан обязательный аргумент --${required_name}" >&2
    usage >&2
    exit 2
  fi
done

output_dir=${output_dir,,}
genome_name=$(basename "$genome" | cut -d. -f1)
genome_dir="$(cd -- "$(dirname -- "$genome")" && pwd)"

# 2bit
if [ ! -f "$genome_dir/${genome_name}.2bit" ]; then
  faToTwoBit "$genome" "$genome_dir/${genome_name}.2bit"
fi

# Проверка bowtie индексов
if [ ! -f "$genome_dir/${genome_name}.1.ebwt" ]; then
  bowtie-build "$genome" "$genome_dir/${genome_name}"
fi

annotation_args=(
  --genome "$genome"
  --genome_annotation "$genome_annotation"
  --gene_name "$gene_name"
  --annotation_format "$annotation_format"
)

TARGET_REGION=$(Rscript "$script_dir/get_cords2.R" "${annotation_args[@]}")

# Проверяем, что TARGET_REGION не пуст
if [ -z "$TARGET_REGION" ]; then
  echo -e "Ошибка: Rscript не вернул координаты." >&2
  exit 1
fi


echo -e "Координаты таргета: $TARGET_REGION"

mkdir -p "$output_dir"

# Запуск chopchop
echo 'Запуск chopchop'
python chopchop/chopchop.py \
  -Target "$TARGET_REGION" \
  -G "$genome_name" \
  -M NGG \
  -T 1 \
  -g 20 \
  -m 3 \
  --padSize 0 \
  -O 50 \
  --scoringMethod DOENCH_2016 \
  -o "$output_dir" > "$output_dir"/n20_table.tsv

mkdir -p "$output_dir"/offtargets_n20
mv "$output_dir"/*offtargets "$output_dir"/offtargets_n20/

echo -e "Извлекаю последовательности плеч гомологии..."
grna_args=(
  "${annotation_args[@]}"
  --target_plasmid "$target_plasmid"
  --output_dir "$output_dir"
)
if [ -n "$cas_plasmid" ]; then
  grna_args+=(--cas_plasmid "$cas_plasmid")
fi
Rscript "$script_dir/get_gRNA_place2.R" "${grna_args[@]}"

if [ $? -ne 0 ]; then
  echo -e "Ошибка при извлечении последовательностей плеч гомологии" >&2
  exit 1
fi


echo -e "\\n Проверка на оффтаргет для праймеров"
java -jar virtualPCR/dist/virtualPCR.jar "$output_dir"/pcr_config.conf &> /dev/null
java -jar virtualPCR/dist/virtualPCR.jar "$output_dir"/pcr_config2.conf &> /dev/null
java -jar virtualPCR/dist/virtualPCR.jar "$output_dir"/pcr_config3.conf &> /dev/null

rm "$output_dir"/sequence.fa "$output_dir"/gene_file.fa
mkdir -p "$output_dir"/additional_files
mv "$output_dir"/bowtie.err "$output_dir"/additional_files/
mv "$output_dir"/twoBitToFa.err "$output_dir"/additional_files/
mv "$output_dir"/right_arm_report.txt "$output_dir"/additional_files/
mv "$output_dir"/left_arm_report.txt "$output_dir"/additional_files/
mv "$output_dir"/pcr_config.conf "$output_dir"/additional_files/
mv "$output_dir"/pcr_config2.conf "$output_dir"/additional_files/
mv "$output_dir"/primer3_settings.txt "$output_dir"/additional_files/
mv "$output_dir"/genome_screening_report.txt "$output_dir"/additional_files/
mv "$output_dir"/primer3_table.tsv "$output_dir"/additional_files/

mkdir -p "$output_dir"/sequences
mv "$output_dir"/offtarget.fa "$output_dir"/sequences/
mv "$output_dir"/homology_arms.fasta "$output_dir"/sequences/
mv "$output_dir"/homology_arms_before_primer_search.fa "$output_dir"/sequences/

cp "$output_dir"/all_primers.fasta "$output_dir"/all_primers.txt

echo -e " \\n \\n!!Результаты сохранены в директории" "$(readlink -f $output_dir)!! \n"
