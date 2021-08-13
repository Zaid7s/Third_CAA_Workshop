cd TimeStepSol/
rm *.dat
cd ../
#read -p "Before running, Make sure you have the right path to write-out the unsteady solution. This is found in 'Main.f90' line 600ish (look for OPEN(57...). The path can be found by typing pwd in terminal. Rememeber it must end with '/......../TimeStepSol/' which is a directory where the unsteady solution is located at (it also includes a matlab script to view the   solutions. Type Any key to run, or control+C to terminate and fix path"     
#rm *.mod
#rm *.o
rm *.dat
./run.exe
