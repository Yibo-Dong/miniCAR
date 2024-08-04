#!/bin/bash

# 统计当前目录下以.res结尾的非空文件个数

count=0

if [ -z "$1" ]
then
  dir="."
else
  dir=$1
fi


# 循环查找以.res结尾的文件
for file in "${dir}"/*.res; do
    # 检查文件是否非空
    if [ -s "$file" ]; then
      first_line=$(head -n 1 "$file")
      if [ "$first_line" = "1" ]; then
        count=$((count+1))
      fi
    fi
done

echo "当前目录下非空的.res文件个数为: $count"


