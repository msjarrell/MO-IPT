# Clean the current directory
#set -x
#set -e 

rm -f *.dat HK_constructed HK_data  log_file out_file
#
#
# Get user input
echo -n "Which sample do you want to run [enter a number between 1 and 5] ? "
read ns
echo
echo -n "Path for mpirun - Press enter to use `which mpirun` "
read mpr
echo
if [ -z $mpr ]
then
	mpr=mpirun
fi
echo
echo "The execution is going to use "`which $mpr`
echo
if [ $ns -le 3 ]
then
	np=1
else
echo -n "How many processors do you want to use? "
read np
fi
#
#
# Copy relevant files for execution
#
rm -rf output/sample_"$ns"_output
mkdir -p output/sample_"$ns"_output
pushd output/sample_"$ns"_output;

cp -rf ../../sample_"$ns"/input/* .
#
# Now execute and store the screen output in log_file
#
if [[ ! -f ../../src/MO_IPT.run ]]; then
    echo "The executable doesn't exist! Compiling now ..."
    pushd ../../..
    make
    popd
fi

rm -rf MO_IPT.run
ln -s ../../../MO_IPT.run .

date > test_outfile  
$mpr -np "$np" ./MO_IPT.run | tee log_file
date  >> out_file
popd;
#
#
echo
echo
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "Now a comparison may be made between the output of this run"
echo "and the one saved in data/sample_"$ns"/output."
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
