read -p 'Enter root filename:' rootname

for f in $rootname*
do   
    echo "deleted $f"
    rm $f
done


