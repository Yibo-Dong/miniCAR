#!/bin/bash

fail_list=()
safe_list=()
unsafe_list=()
timeout_list=()

res_path=$1
aig_path=$2

for res_file in ${res_path}/*.res; do
    name=$(basename "$res_file" .res)
    aig_file="${aig_path}/$name.aig"
    if [ -f "$aig_file" ]; then
        ./tools/aigsim -c "$aig_file" "$res_file"
        status=$?
        first_line=$(head -n 1 "$res_file")
        if [ $status -ne 0 ]; then
            fail_list+=("$name")
        elif [ "$first_line" = "0" ]; then
            safe_list+=("$name")
        elif [ "$first_line" = "1" ]; then
            unsafe_list+=("$name")
        else
            timeout_list+=("$name")
        fi
    else
        echo "Error: $aig_file not found."
    fi
done

echo ""
echo "Failed files: ${#fail_list[@]}"
echo "${fail_list[@]}"
echo ""
echo "Safe files: ${#safe_list[@]}"
echo "${safe_list[@]}"
echo ""
echo "Unsafe files: ${#unsafe_list[@]}"
echo "${unsafe_list[@]}"
echo ""
echo "Timeout files: ${#timeout_list[@]}"
echo "${timeout_list[@]}"

