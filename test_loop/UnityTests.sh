#!/bin/bash
#
#-----------------------------Operation and Case Options------------------------#
#
echo "========================================================================="
echo "This script offers several options, 1, 2 or 3"
echo " "
echo "1 - To install the CFD solver CHAPSim"
echo "2 - To run the test cases"
echo "3 - To compare the results to an existing reference"
echo " "
echo "Type 1, 2 or 3"
read scriptoption
echo " "
#
if [[ $scriptoption == 1 ]]; then
	echo "You have picked the option 1 and the code will be installed"
elif [[  $scriptoption == 2 ]]; then
	echo "You have picked the option 2 and the tests will be run."
	echo " "
	echo "There are several test cases, 1, 2, 3, 4"
	echo " "
	echo "1 - the channel flow with Cartesian Coordinates."
	echo "2 - the pipe flow with Cylindrical Coordinates."
	echo "3 - the annular flow with Cylindrical Coordinates."
	echo "4 - the Taylor Green Vortex flow in a box with Cartesian Coordinates."
	echo " "
	echo "Type 1, 2, 3 or 4"
	read caseoption
	echo " "
elif [[  $scriptoption == 3 ]];then
	echo "You have picked the option 3 and the comparisons will be carried out"
else
	echo "The only available options are 1, 2 or 3"
	exit
fi

#-----------------------------Compile the code----------------------------------#
#
cd ../
#
export ROOTDIR_CHAPSIM=`pwd`
export BINDIR_CHAPSIM=$ROOTDIR_CHAPSIM/bin
#
if [[ $scriptoption == 1 ]]; then
	make all
	echo "  "
	echo "The root directory is $ROOTDIR_CHAPSIM"
	echo "and CHAPSim will be installed under $BINDIR_CHAPSIM"
	echo " "
	echo "The installation is complete"
	echo " "
	echo "Please restart the script with the second option for the tests" 
	echo " "
elif [[ $scriptoption == 2 ]]; then
#
	if [[  $caseoption == 1 ]]; then
		allcases='
		TC11_Channel_TG_Ret180_flow_only_08
		TC12_Channel_TG_Ret180_flow_only_16
		TC21_Channel_IO_Ret180_flow_only_08
		TC22_Channel_IO_Ret180_flow_only_16
		TC31_Channel_TG_IO_Ret180_flow_only_08
		TC32_Channel_TG_IO_Ret180_flow_only_16
		TC41_Channel_IO_Ret180_thermal_08
		TC42_Channel_IO_Ret180_thermal_16
		TC51_Channel_TG_IO_Ret180_thermal_08
		TC52_Channel_TG_IO_Ret180_thermal_16'
	elif [[  $caseoption == 2 ]]; then 
		allcases='
		TP11_Pipe_TG_Ret180_flow_only_08
		TP12_Pipe_TG_Ret180_flow_only_16
		TP21_Pipe_IO_Ret180_flow_only_08
		TP22_Pipe_IO_Ret180_flow_only_16
		TP31_Pipe_TG_IO_Ret180_flow_only_08
		TP32_Pipe_TG_IO_Ret180_flow_only_16
		TP41_Pipe_IO_Ret180_thermal_08
		TP42_Pipe_IO_Ret180_thermal_16
		TP51_Pipe_TG_IO_Ret180_thermal_08
		TP52_Pipe_TG_IO_Ret180_thermal_16'
	elif [[  $caseoption == 3 ]]; then
		allcases='
		TA11_Annular_TG_Ret180_flow_only_08
		TA12_Annular_TG_Ret180_flow_only_16
		TA21_Annular_IO_Ret180_flow_only_08
		TA22_Annular_IO_Ret180_flow_only_16
		TA31_Annular_TG_IO_Ret180_flow_only_08
		TA32_Annular_TG_IO_Ret180_flow_only_16
		TA41_Annular_IO_Ret180_thermal_08
		TA42_Annular_IO_Ret180_thermal_16
		TA51_Annular_TG_IO_Ret180_thermal_08
		TA52_Annular_TG_IO_Ret180_thermal_16'
	elif [[  $caseoption == 4 ]]; then
		allcases='
		TT21_TGV_IO_Ret180_flow_only_08
		TT22_TGV_IO_Ret180_flow_only_16
		TT41_TGV_IO_Ret180_thermal_08
		TT42_TGV_IO_Ret180_thermal_16'
	else
		exit
	fi
	#
	for testcase in $allcases; do
		#
		nbr_procs=${testcase:(-2)}
		#
		echo "==> The case $testcase is run on $nbr_procs processors"
		#
		cd $ROOTDIR_CHAPSIM/test_cases/$testcase
		#
		./clean_CHAPSim_case.sh
		#
		mpirun -np $nbr_procs $BINDIR_CHAPSIM/CHAPSim
		#. 
		echo " "
		#
		tail -4 $ROOTDIR_CHAPSIM/test_cases/$testcase/0_log_monitors/CHA*log
		#
		echo " "
		echo "The case $testcase is complete"
		echo " "
		#
	done
	#
	echo " "
	echo "All the tests are done.  Please restart the script with the third option to check the validity of the results"
	echo " "
	#
elif [[ $scriptoption == 3 ]]; then
	for testcase in $allcases; do
	  #
	  cd $ROOTDIR_CHAPSIM/test_cases/$testcase
	  #
	  diff 0_log_monitors/monitor*log monitor_reference.log
	  #
	  echo " "
	  echo "The difference with the reference for testcase $testcase is complete"
	  echo " "
	  #
	done
else
  	exit
fi
#
