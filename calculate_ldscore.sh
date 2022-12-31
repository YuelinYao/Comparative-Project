t=DVG
m=toplast_10%

for chr in {1..22}
do
../ldscore/ldsc/ldsc.py \
--l2 \
--bfile ../ldscore/1kg_eur/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr \
--ld-wind-cm 1 \
--annot  "annotation/$t""_$m/$t""_$m.$chr.annot.gz" \
--out  "annotation/$t""_$m/$t""_$m.$chr" \
--print-snps ../ldscore/ldsc/hapmap3_snps/hm.$chr.snp
done
