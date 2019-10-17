cd test/
for i in $(ls | grep test_*)
do
	echo $i
	coverage run $i
done
coverage setup.py test
coverage combine

