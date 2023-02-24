#!/bin/bash

for ((i=1; i<=27; i++));
	do sbatch polysolver$i.sh;
done
