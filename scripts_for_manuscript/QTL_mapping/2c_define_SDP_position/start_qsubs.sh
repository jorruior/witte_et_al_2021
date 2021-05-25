cd RatHeartTranslatome/R/RI_panel/2c_define_SDP_position

###lv polyA (475020 unknowns)
for (( i=0; i < 23; i++ ))
do
cd RatHeartTranslatome/results/QTL_mapping/Result_tables/20180213_Lv_polyA
value=20000
start=$(expr $i \* $value)
start=$(expr $start + 1)
end=$(expr $i + 1)
end=$(expr $end \* $value)
id=$(expr $i + 1)
qsub -N lvpolya -l h_vmem=40G -b y /usr/bin/Rscript RatHeartTranslatome/R/RI_panel/2c_define_SDP_position/Lv_PolyA.R $id $start $end
done;

qsub -N lvpolya -l h_vmem=40G -b y /usr/bin/Rscript RatHeartTranslatome/R/RI_panel/2c_define_SDP_position/Lv_PolyA.R 24 460001 475020

###lv ribo

for (( i=0; i < 23; i++ ))
do
cd RatHeartTranslatome/results/QTL_mapping/Result_tables/20180213_Lv_Ribo
value=20000
start=$(expr $i \* $value)
start=$(expr $start + 1)
end=$(expr $i + 1)
end=$(expr $end \* $value)
id=$(expr $i + 1)
qsub -N lvribo -l h_vmem=40G -b y /usr/bin/Rscript RatHeartTranslatome/R/RI_panel/2c_define_SDP_position/Lv_Ribo.R $id $start $end
done;

qsub -N lvribo -l h_vmem=40G -b y /usr/bin/Rscript RatHeartTranslatome/R/RI_panel/2c_define_SDP_position/Lv_Ribo.R 24 460001 475020


###lv RiboRnaResid

for (( i=0; i < 23; i++ ))
do
cd RatHeartTranslatome/results/QTL_mapping/Result_tables/20180213_Lv_Ribo_RNA_resid
value=20000
start=$(expr $i \* $value)
start=$(expr $start + 1)
end=$(expr $i + 1)
end=$(expr $end \* $value)
id=$(expr $i + 1)
qsub -N lvresid -l h_vmem=40G -b y /usr/bin/Rscript RatHeartTranslatome/R/RI_panel/2c_define_SDP_position/Lv_RiboRnaResid.R $id $start $end
done;

qsub -N lvresid -l h_vmem=40G -b y /usr/bin/Rscript RatHeartTranslatome/R/RI_panel/2c_define_SDP_position/Lv_RiboRnaResid.R 24 460001 475020



#liver (420840)

###liver polyA 
for (( i=0; i < 21; i++ ))
do
cd RatHeartTranslatome/results/QTL_mapping/Result_tables/20180213_Liver_polyA
value=20000
start=$(expr $i \* $value)
start=$(expr $start + 1)
end=$(expr $i + 1)
end=$(expr $end \* $value)
id=$(expr $i + 1)
qsub -N lipolya -l h_vmem=40G -b y /usr/bin/Rscript RatHeartTranslatome/R/RI_panel/2c_define_SDP_position/Liver_PolyA.R $id $start $end
done;

qsub -N lipolya -l h_vmem=40G -b y /usr/bin/Rscript RatHeartTranslatome/R/RI_panel/2c_define_SDP_position/Liver_PolyA.R 22 420001 420840

###liver ribo

for (( i=0; i < 21; i++ ))
do
cd RatHeartTranslatome/results/QTL_mapping/Result_tables/20180213_Liver_Ribo
value=20000
start=$(expr $i \* $value)
start=$(expr $start + 1)
end=$(expr $i + 1)
end=$(expr $end \* $value)
id=$(expr $i + 1)
qsub -N liribo -l h_vmem=40G -b y /usr/bin/Rscript RatHeartTranslatome/R/RI_panel/2c_define_SDP_position/Liver_Ribo.R $id $start $end
done;

qsub -N liribo -l h_vmem=40G -b y /usr/bin/Rscript RatHeartTranslatome/R/RI_panel/2c_define_SDP_position/Liver_Ribo.R 22 420001 420840


###liver RiboRnaResid

for (( i=0; i < 21; i++ ))
do
cd RatHeartTranslatome/results/QTL_mapping/Result_tables/20180403_Liver_Ribo_RNA_resid
value=20000
start=$(expr $i \* $value)
start=$(expr $start + 1)
end=$(expr $i + 1)
end=$(expr $end \* $value)
id=$(expr $i + 1)
qsub -N liresid -l h_vmem=40G -b y /usr/bin/Rscript RatHeartTranslatome/R/RI_panel/2c_define_SDP_position/Liver_RiboRnaResid.R $id $start $end
done;

qsub -N liresid -l h_vmem=40G -b y /usr/bin/Rscript RatHeartTranslatome/R/RI_panel/2c_define_SDP_position/Liver_RiboRnaResid.R 22 420001 420840

