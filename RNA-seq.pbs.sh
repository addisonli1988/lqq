#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -pe mpi 4
#$ -q main.q
#$ -N RNA-seq_13.yangyili_20170624
sh /gluster/home/yuanxiao/scripts/RNA-seq.sh /gluster/home/yuanxiao/RNA-seq/13.yangyili_20170624 hg19 PE
