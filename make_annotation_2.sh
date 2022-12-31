t=DVG
for chr in {22..1}
do
../ldscore/ldsc/make_annot.py \
--gene-set-file GeneList/$t/Overlap_Across_Tissues/conserved.txt \
--gene-coord-file Cattle_Human_17k_Orthologous_Genes_GRCh37 \
--windowsize 50000 \
--bimfile ../ldscore/1kg_eur/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr.bim \
--annot-file "annotation/$t""_overlap_across_tissues/conserved.$chr"

../ldscore/ldsc/make_annot.py \
--gene-set-file GeneList/$t/Overlap_Across_Tissues/diverged.txt \
--gene-coord-file Cattle_Human_17k_Orthologous_Genes_GRCh37 \
--windowsize 50000 \
--bimfile ../ldscore/1kg_eur/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr.bim \
--annot-file "annotation/$t""_overlap_across_tissues/diverged.$chr"

../ldscore/ldsc/make_annot.py \
--gene-set-file GeneList/$t/Overlap_Across_Tissues/overlap.txt \
--gene-coord-file Cattle_Human_17k_Orthologous_Genes_GRCh37 \
--windowsize 50000 \
--bimfile ../ldscore/1kg_eur/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr.bim \
--annot-file "annotation/$t""_overlap_across_tissues/overlap.$chr"
done
