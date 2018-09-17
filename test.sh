for i in $(ls | grep test_*)
do
  echo $i
  coverage run $i
done    
converage setup.py test
coverage combine
