# Downstream analysis of novel proteins identified by PG

List all copies of novel peptides in both novel and known databases

```bash
	cat ../proteogenomics/novel.fasta ../proteogenomics/known.fasta > merged.fasta
	bioawk -c 'fastx' '{printf(">%s\n%s\n",$name,$seq)}' merged.fasta  > merged_greppable.fasta
	cat really_novel_pg.gff3 | ./peptides_from_gff.rb | sort -u | ./peptide_uniqueness.rb > novel_peptides_copies.txt
```

Use the peptide uniqueness file to annotate the gff

```bash
	./score_by_uniqueness.rb > really_novel_uniqueness_pg.gff3
```

Extract corresponding ORFs

samtools faidx novel.fasta $(cat really_novel_uniqueness_pg.gff3 |awk -F '\t' '{printf("lcl|%s\n", substr($9,7+index($9,"Parent=")))}' | sort -u | tr '\n' ' ') > really_novel_uniqueness_pg.fasta

