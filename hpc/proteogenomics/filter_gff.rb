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
