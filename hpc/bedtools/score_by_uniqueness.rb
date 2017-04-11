#!/usr/bin/env ruby

require 'bio'

peptide_uniqueness = {}

peptide_is_covered_by_known = {}

File.foreach("novel_peptides_copies.txt") do |line|  
	line.chomp!
	lparts = line.split(" ")
	if peptide_uniqueness[lparts[0]]
		peptide_uniqueness[lparts[0]] += 1
	else
		peptide_uniqueness[lparts[0]] = 1
		peptide_is_covered_by_known[lparts[0]] = false
	end

	if lparts[1] !~ /frame/
		peptide_is_covered_by_known[lparts[0]] = true
	end

end



File.foreach("really_novel_pg.gff3") do |line|  
	line.chomp!
	lparts = line.split(/\s/)
	peptide = lparts[8].match(/[0-9]+\.([A-Z]+)\.[0-9]+\;/).captures[0]
	lparts[5] = 1/(peptide_uniqueness[peptide]*1.0)

	unless peptide_is_covered_by_known[peptide]
		puts lparts.join("\t")
	end
end