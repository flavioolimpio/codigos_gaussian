#!/bin/bash

# Defina o diretório de entrada aqui
input_dir="."

# Verifique se o diretório de entrada existe
if [ ! -d "$input_dir" ]; then
    echo "Diretório de entrada não encontrado: $input_dir"
    exit 1
fi

# Loop em todos os arquivos log do Gaussian no diretório de entrada
for input_file in "$input_dir"/*.log; do
    # Verifique se o arquivo de entrada existe
    if [ ! -f "$input_file" ]; then
        continue
    fi

    # Encontre a última ocorrência da orientação padrão no arquivo de entrada
    start_line=$(grep -n 'Standard orientation:' "$input_file" | tail -1 | cut -d: -f1)

    # Verifique se a orientação padrão foi encontrada
    if [ -z "$start_line" ]; then
        echo "Orientação padrão não encontrada no arquivo de entrada: $input_file"
        continue
    fi

    # Encontre a linha de início dos dados da tabela
    start_line=$((start_line + 5))

    # Encontre a linha final dos dados da tabela
    end_line=$(sed -n "${start_line},\$p" "$input_file" | grep -n ' ----' | head -1 | cut -d: -f1)
    end_line=$((end_line + start_line - 2))

    # Extraia os valores das colunas Atomic Number e coordenadas X, Y e Z
    data=$(sed -n "${start_line},${end_line}p" "$input_file" | awk '{print $2, $4, $5, $6}')

    # Crie os arquivos de saída com base no nome do arquivo de entrada
    base_name="${input_file%.*}"
    output_file_1="${base_name}_N.gjf"
    output_file_2="${base_name}_N+1.gjf"
    output_file_3="${base_name}_N-1.gjf"

    # Escreva os dados nos arquivos de saída com as linhas adicionais especificadas
    input=`echo $input_file | sed 's/.log//g' | sed 's|^./||'`
	echo -e "%nprocs=4\n%mem=8GB\n%chk=${input}_N.chk\n# m062x/6-31+g(d) scrf=smd out=wfn\n\n${input}\n" > "$output_file_1"
	echo -e "%nprocs=4\n%mem=8GB\n%chk=${input}_N+1.chk\n# m062x/6-31+g(d) scrf=smd out=wfn\n\n${input}\n" > "$output_file_2"
	echo -e "%nprocs=4\n%mem=8GB\n%chk=${input}_N-1.chk\n# m062x/6-31+g(d) scrf=smd out=wfn\n\n${input}\n" > "$output_file_3"
	
	echo "0 1" >> "$output_file_1"
    echo "-1 2" >> "$output_file_2"
    echo "1 2" >> "$output_file_3"
	
	echo "$data" >> "$output_file_1"
	echo " " >> "$output_file_1"
    echo "$data" >> "$output_file_2"
	echo " " >> "$output_file_2"
    echo "$data" >> "$output_file_3"
	echo " " >> "$output_file_3"
	
	echo "/home/valter/Flavio/elisa_ueg/Fukui/${input}_N.wfn" >> "$output_file_1"
	echo "/home/valter/Flavio/elisa_ueg/Fukui/${input}_N+1.wfn" >> "$output_file_2"
	echo "/home/valter/Flavio/elisa_ueg/Fukui/${input}_N-1.wfn" >> "$output_file_3"
	
done
