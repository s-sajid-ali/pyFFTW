for i in $(ls test/ | grep test_*)
do
	echo $i
	coverage run test/$i
done
coverage run setup.py test
coverage combine

