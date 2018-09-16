for i in $(ls /test/ | grep test_*); do
  coverage run i
done    
