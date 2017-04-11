# for Aim 1
# remove the known peptides if it occurs on the same strand

#use bedtools subtract Remove intervals based on overlaps b/w two files
# 	-s	Require same strandedness.  That is, only report hits in B
		that overlap A on the _same_ strand.
		- By default, overlaps are reported without respect to strand.

    #-f	Minimum overlap required as a fraction of A.
		- Default is 1E-9 (i.e., 1bp).
		- FLOAT (e.g. 0.50)
#set overlap to 1.0 means the fraction of overlap is 100% (overlaps completely)
###

```bash
bedtools subtract -s -f 1.0 -a novel_pg.gff3 -b known_pg.gff3
```

# count the number of lines (966) (reduced from 2773 from novel_pg.gff3 as obtained by (wc -l novel_pg.gff3))

```bash
bedtools subtract -s -f 1.0 -a novel_pg.gff3 -b known_pg.gff3 | wc -l
```

#put this output in a different file called really.novel.gff3


```bash
bedtools subtract -s -f 1.0 -a novel_pg.gff3 -b known_pg.gff3 > really_novel.gff3
```

#Cut the first column containing file names

```bash
cat really_novel.gff3 | cut -f1
```

#sort them by newest first

```bash
cat really_novel.gff3 | cut -f1 | sort -u
```

#copy transcriptome fasta file to current dir

```bash
cp ../proteogenomics/
```

```bash
cp ../proteogenomics/transcriptome_na.fasta .
```

#Use samtools fastaindex (faidx) to match novel peptide names to their respective fasta transcript on separate lines

```bash
$ samtools faidx transcriptome_na.fasta $(cat really_novel.gff3 | cut -f1 | sort -u | tr '\n' ' ')
```

#put this output in a fasta index file called really_novel.fasta

```bash
samtools faidx transcriptome_na.fasta $(cat really_novel.gff3 | cut -f1 | sort -u | tr '\n' ' ') > really_novel.fasta
```

#Copy ruby files from proteogenomics folder to current dir to get the ruby program filter_gff.rb

```bash
cp ../proteogenomics/
```

```bash
cp ../proteogenomics/*.rb .
```

#cut the first column with id names in it from really_novel.gff3

```bash
cat really_novel.gff3 | cut -f1
```

#Sort them according to newest first

```bash
cat really_novel.gff3 | cut -f1 | sort -u
```

#Put results in .txt file

```bash
cat really_novel.gff3 | cut -f1 | sort -u > really_novel_ids.txt
```

#run filter_gff ruby program to mix known and novels

```bash
./filter_gff.rb
```

#place results into a gff file

```bash
./filter_gff.rb > known_with_novel.gff3
```
