set title "neBEM elements in gnuplot VIEWER"
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
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim1.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim2.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim3.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim4.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim5.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim6.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim7.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim8.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim9.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim10.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim11.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim12.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim13.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim14.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim15.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim16.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim17.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim18.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim19.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim20.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim21.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim22.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim23.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim24.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim25.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim26.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim27.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim28.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim29.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim30.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim31.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim32.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim33.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim34.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim35.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim36.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim37.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim38.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim39.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim40.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim41.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim42.out' w p

pause-1