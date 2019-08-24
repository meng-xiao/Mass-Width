#!/bin/bash
for ch in 2e2mu #4e 4mu
	do
	for cat_vbf in 0 #1
	do
	for onshell in 1 #1
	do
		./job.lsf $ch $cat_vbf $onshell
		#bsub -q cmscaf1nd job.lsf $ch $cat_vbf $onshell
	done
	done
	done
