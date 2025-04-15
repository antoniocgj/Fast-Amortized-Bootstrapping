echo "##### BINARY #####"
make -B PARAM=SET_2_3 && ./main
make -B PARAM=SET_4_5 && ./main
make -B PARAM=SET_6_7 && ./main
make -B PARAM=SET_8_9 && ./main
make -B PARAM=SET_8_9_b && ./main
make -B PARAM=SET_8_9_HIGH_FR && ./main
echo "##### TERNARY #####"
make -B KEY=TERNARY PARAM=SET_2_3 && ./main
make -B KEY=TERNARY PARAM=SET_4_5 && ./main
make -B KEY=TERNARY PARAM=SET_6_7 && ./main
make -B KEY=TERNARY PARAM=SET_8_9 && ./main
make -B KEY=TERNARY PARAM=SET_8_9_HIGH_FR && ./main
echo "##### ARBITRARY #####"
make -B KEY=ARBITRARY PARAM=SET_A2 && ./main
make -B KEY=ARBITRARY PARAM=SET_A4 && ./main
make -B KEY=ARBITRARY PARAM=SET_A6 && ./main
echo "##### LARGE #####"
make -B KEY=TERNARY PARAM=SET_L0 && ./main
make -B KEY=TERNARY PARAM=SET_L1 && ./main
make -B KEY=TERNARY PARAM=SET_L2 && ./main
make -B KEY=TERNARY PARAM=SET_L3 && ./main
make -B KEY=ARBITRARY PARAM=SET_LA1 && ./main
make -B KEY=ARBITRARY PARAM=SET_LA2 && ./main
make -B KEY=ARBITRARY PARAM=SET_LA3 && ./main
make -B KEY=ARBITRARY PARAM=SET_LA4 && ./main