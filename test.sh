cp /home/travis/build/s-sajid-ali/pyFFTW/test/* .
for i in $(ls | grep test_*); do
  coverage run $i
done    
