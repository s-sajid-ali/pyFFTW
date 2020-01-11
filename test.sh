#for i in $(ls test/ | grep test_*)
#do
#	echo $i
#	coverage run test/$i
#done

cd test
for i in $(ls | grep test_*)
do
	echo $i
	coverage run test/$i
done
cd ../
coverage run setup.py test
coverage combine
