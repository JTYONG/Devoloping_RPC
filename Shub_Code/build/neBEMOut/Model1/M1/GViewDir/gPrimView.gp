set title "neBEM primitives in gnuplot VIEWER"
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
 'neBEMOut/Model1/M1/GViewDir/gpPrim1.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim2.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim3.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim4.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim5.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim6.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim7.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim8.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim9.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim10.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim11.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim12.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim13.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim14.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim15.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim16.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim17.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim18.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim19.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim20.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim21.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim22.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim23.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim24.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim25.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim26.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim27.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim28.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim29.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim30.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim31.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim32.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim33.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim34.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim35.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim36.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim37.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim38.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim39.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim40.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim41.out' w l, \
 'neBEMOut/Model1/M1/GViewDir/gpPrim42.out' w l

pause-1