#!/usr/bin/env ruby

require 'bio'

ARGF.each do |item|
	puts item.match(/[0-9]+\.([A-Z]+)\.[0-9]+\;/).captures
end