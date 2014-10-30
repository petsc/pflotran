#!/bin/sh

clear

cur_dir=`echo $PWD`

(cd ../src/pflotran; make pflotran;)

for fname in `find ./ -name '*.in' | grep -v obsolete | grep -v scaling | cut -c 4- `
do

  fname2=`echo $fname| awk -F '.' '{print $1}'`

  cp $fname ${fname2}_test.in

  final_time=`cat ${fname2}_test.in | grep "FINAL_TIME" | grep -v \#`
  sed -i .bk "s/${final_time}/FINAL_TIME 1.d0 s/g" ${fname2}_test.in

  dir_path="${fname%/*}"

  echo "==========================================================="
  echo "Inputfile  : " $fname

  (cd $dir_path; $cur_dir/../src/pflotran/pflotran -pflotranin ${cur_dir}/${fname2}_test.in > logfile);
  exit_status=$?

  echo "Exit status: " $exit_status

  if [ "$exit_status" -ne 86 ]
  then
    echo "Tail of logfile: "
    tail -5 $dir_path/logfile
  fi

  rm -f $dir_path/logfile
  rm -rf ${fname2}_test.in ${fname2}_test.in.bk
done

# Do cleanup
for fname in `find ./ -name '*test.out' `
do
  rm -rf $fname
done

for fname in `find ./ -name '*test.h5' `
do
  rm -rf $fname
done

for fname in `find ./ -name '*test.tec' `
do
  rm -rf $fname
done
