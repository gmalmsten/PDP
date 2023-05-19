echo -e '\n'
echo New run >> strong.txt
for i in 2 4 8 16;
do
echo -n -e $i '\t' >> strong.txt
mpirun -np $i mc 1000000 >> strong.txt
done