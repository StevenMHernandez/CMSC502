make clean

make


for ((i=10;i<=25;i++)); do
    sh gen_cities.sh $i > tmp.txt
    echo "run $i: \c"
    ./sequential tmp.txt
done

echo "done"
