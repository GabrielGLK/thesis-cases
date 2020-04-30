#!/bin/bash


awk 'NR == 1{p = $4;next}
	{print $1, ($4 - p)/0.01; p = $4}
	END {print $1,p}' data> data-1
