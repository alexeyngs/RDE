blockMesh
setFields
decomposePar -force
mpirun -np 4 RDE -parallel
reconstructPar
rm -r processor*