#!/bin/bash
#this script generates a detailed assignment of each read
featureCounts -a ../refgenome/Mus_musculus.GRCm38.91.gtf -s 2 -o fcount -Q 10 -R *bam
