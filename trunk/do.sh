#sh exe.sh images/lena.png 10 10 10 0 64
#sh exe.sh images/niagarafalls_aerial_photo.jpg 10 10 10 0 100
p=1;
p=$1;
for (( i=1; i<=p; i++ ))
do
    export OMP_NUM_THREADS=$i
    #sh exe.sh images/lena.png 10 10 10 0 64
    sh exe.sh images/niagarafalls_aerial_photo.jpg 10 10 10 0 100
done
