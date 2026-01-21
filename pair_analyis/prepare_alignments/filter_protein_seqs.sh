#!/bin/bash
# Создаем целевую папку

DIRS=(
GCA_014839755.1
GCF_015227805.2
GCF_036417665.1
GCF_047830755.1
GCF_009650955.1
GCA_014706295.1
GCF_023634155.1
GCF_015220075.1
# GCA_964417175.1
# GCA_963932325.2
GCF_947461875.1
GCF_964188355.1
)

mkdir -p longest_isoforms_proteomes_clean

# Обрабатываем каждый файл в исходной папке
for id in ${DIRS[@]}; do
    if [[ -f longest_isoforms_proteomes/${id}_proteins.fa ]]; then
        # Извлекаем имя файла
        #filename=$(basename "$file")
        
        # Фильтруем последовательности без точек и сохраняем в новую папку
        awk '
        /^>/ {
            if (seq && seq !~ /\./) {
                print header
                print seq
            }
            header = $0
            seq = ""
            next
        }
        {
            seq = seq $0
        }
        END {
            if (seq && seq !~ /\./) {
                print header
                print seq
            }
        }' longest_isoforms_proteomes/${id}_proteins.fa > "longest_isoforms_proteomes_clean/${id}_proteins.fa"
        
        echo "Обработан: $id"
    fi
done

echo "Готово! Отфильтрованные файлы в папке longest_isoforms_proteomes_clean"