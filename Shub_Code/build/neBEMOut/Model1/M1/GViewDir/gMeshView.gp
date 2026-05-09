set title "neBEM mesh in gnuplot VIEWER"
#set pm3d
#set style data pm3d
#set palette model CMY
set hidden3d
set nokey
set xlabel "X"
set ylabel "Y"
set zlabel "Z"
set view 70, 335, 1, 1

splot \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim1.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim1.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim2.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim2.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim3.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim3.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim4.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim4.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim5.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim5.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim6.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim6.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim7.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim7.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim8.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim8.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim9.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim9.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim10.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim10.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim11.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim11.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim12.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim12.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim13.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim13.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim14.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim14.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim15.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim15.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim16.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim16.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim17.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim17.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim18.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim18.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim19.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim19.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim20.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim20.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim21.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim21.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim22.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim22.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim23.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim23.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim24.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim24.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim25.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim25.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim26.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim26.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim27.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim27.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim28.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim28.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim29.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim29.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim30.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim30.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim31.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim31.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim32.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim32.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim33.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim33.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim34.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim34.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim35.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim35.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim36.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim36.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim37.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim37.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim38.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim38.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim39.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim39.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim40.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim40.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim41.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim41.out' w p ps 1, \
 'neBEMOut/Model1/M1/GViewDir/gpMeshOnPrim42.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim42.out' w p ps 1

pause-1