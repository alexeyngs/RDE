mv ./system/setFieldsDict ./system/setFieldsDict-temp
mv ./system/setFieldsDict0 ./system/setFieldsDict
setFields
blockMesh
mv ./system/setFieldsDict ./system/setFieldsDict0
mv ./system/setFieldsDict-temp ./system/setFieldsDict
setFields
decomposePar -force
mpirun -np 4 RDE -parallel
reconstructPar
rm -r processor*