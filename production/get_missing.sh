file=$1
augmentations=$2

#echo $file
    
flatfile="${file/lhes/flat}" 
#echo $flatfile
while read v;
do
    variation=`echo $v | sed 's/ .*//g'`
    setting=`echo $v | sed 's/.* //g'`
    
    #echo $variation
    flatflat="${flatfile/.lhe/_${variation}_flat.root}"
    #echo $flatflat
    
    # check if file exists:
    if [ -f "$flatflat" ]; then
	#echo YES FILE $flatflat exists!
	continue
    else
	echo $file
    fi
done < $2
