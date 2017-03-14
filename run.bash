#! /bin/bash

# Make the thing
make anne

# Setup result directory
dir=RESLT
if [ -e $dir ]; then
    echo " "
    echo " "
    echo "Directory " $dir "already exists -- please move it out of the way."
    echo " "
    echo " "
    exit
fi
mkdir $dir

# Copy driver code in for future reference
cp anne.cc $dir

# Run the bastard
echo "Running..."
./anne --y_bl 3 --x_mid 4.5 --use_oomph_gmres > $dir/OUTPUT 
echo "...done"

# Post-process
cd $dir
oomph-convert -z soln*.dat
makePvd soln soln.pvd
cd ..

echo "=============================================================="
echo " "
echo "Done! Go to directory $dir and animate the flow field with "
echo " "
echo " paraview --state=../flow.pvsm"
echo " "
echo "=============================================================="
