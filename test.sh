for i in $(ls $pwd/test/ | grep test_*); do
  coverage run i
done    
