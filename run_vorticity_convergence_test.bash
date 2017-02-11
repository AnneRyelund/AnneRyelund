#! /bin/bash

# Make the thing
make anne

# Setup result directory
dir=RESLT
if [ -e $dir ]; then
    echo " "
    echo " "
    echo "Directory " $main_dir "already exists -- please move it out of the way."
    echo " "
    echo " "
    exit
fi
mkdir $dir

# Copy driver code in for future reference
cp anne.cc $dir

# Run the bastard
echo "Running..."
./anne --validate_projection > $dir/OUTPUT 
echo "...done"

# Post-process
cd $dir
oomph-convert -z analytical_vorticity_and_indicator*.dat 
makePvd analytical_vorticity_and_indicator analytical_vorticity_and_indicator.pvd
oomph-convert soln*.dat 
makePvd soln soln.pvd

cd ..

echo "=============================================================="
echo " "
echo "Done! Go to directory $dir and animate the validation with "
echo " "
echo "   paraview --state=../vorticity_convergence.pvsm"
echo " "
echo "and/or "
echo " "
echo "   tecplot ../vorticity_convergence.lay"
echo " "
echo "=============================================================="
