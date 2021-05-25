cd /data/huebner2/ANALYSES/20171204_RatHeartTranslatome/R/RI_panel/2b_find_correct_SDP

qsub -l h_vmem=60G -b y /usr/bin/Rscript Lv_PolyA.R
qsub -l h_vmem=60G -b y /usr/bin/Rscript Liver_PolyA.R
qsub -l h_vmem=60G -b y /usr/bin/Rscript Lv_Ribo.R
qsub -l h_vmem=60G -b y /usr/bin/Rscript Liver_Ribo.R
qsub -l h_vmem=60G -b y /usr/bin/Rscript Lv_RiboRnaResid.R
qsub -l h_vmem=60G -b y /usr/bin/Rscript Liver_RiboRnaResid.R
