#!/usr/bin/bash


matrixFile=GSE131617-GPL5175_expr_TC_Baak_stage.txt
script=anova.py
python=~/Biosoft/miniconda2/envs/py2/bin/python3
prefix="TC"
$python $script $matrixFile $prefix
