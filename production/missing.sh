ls -l $1/*lhe  | awk '{print $9}' | xargs -n1 -P8 -i bash get_missing.sh {} $2 > missing_files.txt

sed -i 's/\/\//\//g' missing_files.txt

#. get_missing.sh $2
