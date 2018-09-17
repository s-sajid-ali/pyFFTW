for i in $(ls /home/travis/build/s-sajid-ali/pyFFTW/test/ | grep test_*); do
  coverage run i
done    
