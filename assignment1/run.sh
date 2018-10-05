make clean

make sequential
make threaded
make mpi

# Time tests (sequential)
echo "\c" > tmp/seq_runs.txt
for ((i=10;i<=25;i++)); do
    sh gen_cities.sh $i > tmp.txt
    echo "run $i: \c" >> tmp/seq_runs.txt
    ./sequential tmp.txt >> tmp/seq_runs.txt
done

# Time tests (thread)
echo "\c" > tmp/thread_runs.txt
for ((i=10;i<=250;i++)); do
    sh gen_cities.sh $i > tmp.txt
    echo "run $i: \c" >> tmp/thread_runs.txt
    ./threaded tmp.txt >> tmp/thread_runs.txt
done

## Time tests (mpi)
#echo "\c" > tmp/mpi_runs.txt
#for ((i=10;i<=100;i++)); do
#    sh gen_cities.sh $i > tmp.txt
#    echo "run $i: \c" >> tmp/mpi_runs.txt
#    ./mpi tmp.txt >> tmp/mpi_runs.txt
#done

echo "done"
