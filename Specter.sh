#!/bin/bash

spark-submit --driver-memory $1 Specter_Spark.py $2 $3 $4 $5 $6 $7 $8

Rscript SpecterQuant.R $2 $3