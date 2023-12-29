!#/bin/bash


for file in *.fq.gz; do
    new_name=$(echo "$file" | sed 's/Read/_Read/g')
    echo $new_name
    mv "$file" "$new_name"
done