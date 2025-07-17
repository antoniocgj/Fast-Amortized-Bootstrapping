echo "##### BINARY #####"
make -B PARAM=SET_2_3_2048 && ./main
make -B PARAM=SET_2_3_4096 && ./main
make -B PARAM=SET_2_3_8192 && ./main
make -B PARAM=SET_4_5_2048 && ./main
make -B PARAM=SET_4_5_4096 && ./main
make -B PARAM=SET_4_5_8192 && ./main
make -B PARAM=SET_6_7_4096 && ./main
make -B PARAM=SET_6_7_8192 && ./main
make -B PARAM=SET_8_9_4096 && ./main
make -B PARAM=SET_8_9_8192 && ./main
make -B PARAM=SET_8_9_HIGH_FR && ./main
echo "##### TERNARY #####"
make -B KEY=TERNARY PARAM=SET_2_3_2048 && ./main
make -B KEY=TERNARY PARAM=SET_2_3_4096 && ./main
make -B KEY=TERNARY PARAM=SET_2_3_8192 && ./main
make -B KEY=TERNARY PARAM=SET_4_5_2048 && ./main
make -B KEY=TERNARY PARAM=SET_4_5_4096 && ./main
make -B KEY=TERNARY PARAM=SET_4_5_8192 && ./main
make -B KEY=TERNARY PARAM=SET_6_7_4096 && ./main
make -B KEY=TERNARY PARAM=SET_6_7_8192 && ./main
make -B KEY=TERNARY PARAM=SET_8_9_4096 && ./main
make -B KEY=TERNARY PARAM=SET_8_9_8192 && ./main
echo "##### ARBITRARY #####"
make -B KEY=ARBITRARY PARAM=SET_A2 && ./main
make -B KEY=ARBITRARY PARAM=SET_A3 && ./main
make -B KEY=ARBITRARY PARAM=SET_A4 && ./main
make -B KEY=ARBITRARY PARAM=SET_A5 && ./main
