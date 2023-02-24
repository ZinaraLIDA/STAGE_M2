#!/bin/bash

for i in b1_1 b1_2 b1_3 b1_4 b1_5 b1_6 b1_7 b1_8 b1_9 b1_10 b1_11 b1_12 b2_1 b2_2 b2_3 b2_4 b2_5 b2_6 b2_7 b2_8 b2_9 b2_10 b2_11 b2_12 b2_13 b2_14 b2_15;
	do sbatch 6-scripts_purecn/purecn_$i.sh;
done
