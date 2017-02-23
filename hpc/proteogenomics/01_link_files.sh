ln -s ../trinity/transdecoder.pep transdecoder.fasta
ln -s ../trinity/Trinity.fasta.transdecoder_dir/longest_orfs.gff3 transdecoder.gff3
ln -s ../trinity/Trinity.fasta transcriptome_na.fasta

for f in ../msdata/Elite/MC038/*.mzML; do ln -s $f .;done
for f in ../msdata/QE/MC383/*.mzML; do ln -s $f .;done

# This is so we only use the one good technical rep here
# but keep names standardised
mv Alvaro_S1A_150810090743.mzML Alvaro_S1A.mzML