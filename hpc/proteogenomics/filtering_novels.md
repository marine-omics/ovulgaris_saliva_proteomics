#for aim 1 - to make database
  #have known and novel gff files but some overlapped so wanted to remove the knowns from the novels
     # use bedtools to take out the known peptide if it occurs on the same strand (i.e it points the same direction "5-3" "3-5" are two strands)


###from bedtools help
         #bedtools   subtract      Remove intervals based on overlaps b/w two files.
    #Usage:   bedtools subtract [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>
# 	-s	Require same strandedness.  That is, only report hits in B
		that overlap A on the _same_ strand.
		- By default, overlaps are reported without respect to strand.

    #-f	Minimum overlap required as a fraction of A.
		- Default is 1E-9 (i.e., 1bp).
		- FLOAT (e.g. 0.50)
#set overlap to 1.0 means the fraction of overlap is 100% (overlaps completely)
###



```bash
bedtools subtract -f 1.0 -s -a novel_pg.gff3 -b known.gff3
```

#count the number of lines (resulted in 689)

```bash
bedtools subtract -f 1.0 -s -a novel_pg.gff3 -b known.gff3 | wc -l
```

#put this output in a different file called really.novel.gff3

```bash
bedtools subtract -f 1.0 -s -a novel_pg.gff3 -b known.gff3 > really_novel.gff3
```

### then used awk and samtools transform this file into a FASTA index on separate lines

```bash
samtools faidx transcriptome_na.fasta `cat really_novel.gff3 | awk '{print $1}' | sort -u  | tr '\n' ' '`
[fai_load] build FASTA index.
```

#put these results in a new file called really.novel.FASTA

```bash
samtools faidx transcriptome_na.fasta `cat really_novel.gff3 | awk '{print $1}' | sort -u  | tr '\n' ' '` > really_novel.fasta
```

#created ruby program ./filter_gff.rb to ... ???

### as in filter_gff.gb file

#! /usr/bin/env ruby

require 'set'

really_novel_ids = Set.new()

File.open('really_novel_ids.txt').each do |line|
  really_novel_ids.add(line.chomp)
end

#require 'byebug';byebug

File.open('known.gff3').each do |line|
  line.chomp!
  line_id = line.split(/\s/)[0]
#  puts line_id
  puts line if really_novel_ids.include?(line_id)
end
