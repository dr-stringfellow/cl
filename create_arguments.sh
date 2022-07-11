while read f;
do
    while read v;
    do
	echo $f,$v | sed 's/ /,/g'
    done < $2
done < $1
