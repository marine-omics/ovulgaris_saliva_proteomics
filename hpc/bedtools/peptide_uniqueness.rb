#!/usr/bin/env ruby

require 'bio'

ARGF.each do |item|
	item.chomp!
	res = %x(grep -B1 #{item} merged_greppable.fasta)


	res_items = res.split("\n")
	(0...res_items.count).each { |x| 
		ri = res_items[x]
		if ri =~ /^>/
			puts "#{item} #{ri}"
		end
	}
end
