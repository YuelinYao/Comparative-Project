t=DVG
while read p
do
for chr in {22..1}
do
../ldscore/ldsc/make_annot.py \
--gene-set-file GeneList/$t/Top_Last_10%/Top_10%/$p-top10%.txt \
--gene-coord-file Cattle_Human_17k_Orthologous_Genes_GRCh37 \
--windowsize 50000 \
--bimfile ../ldscore/1kg_eur/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr.bim \
--annot-file "annotation/$t""_toplast_10%/$p-top10%.$chr"

../ldscore/ldsc/make_annot.py \
--gene-set-file GeneList/$t/Top_Last_10%/Last_10%/$p-last10%.txt \
--gene-coord-file Cattle_Human_17k_Orthologous_Genes_GRCh37 \
--windowsize 50000 \
--bimfile ../ldscore/1kg_eur/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr.bim \
--annot-file "annotation/$t""_toplast_10%/$p-last10%.$chr"
done
done < tissue.txt
