for i in $(ls $HOME/test/ | grep test_*); do
  coverage run i
done    
