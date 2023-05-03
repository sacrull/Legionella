''' Adds missing fields to list of taxonomy strings. Use: python fix_taxonomy.py my_taxa.txt > new_taxa.txt '''

import sys

with open(sys.argv[1], "r") as f: # read in file
		for line in f: # for each line in file
			x = line.rstrip("\n") # assign line to object 'x', strip new line character from object
			y = x.count(";") # count the number of occurances of ';', assign to object 'y'
			if y == 0: # if only one taxonomic level assigned (e.g., unknown or kingdom level)
				new = (x.split(";")[0] + ";") + ("Bacteria_unknown;") * 5 + "Bacteria_unknown"
				print(new) # print the new taxonomy string
			elif y < 6: # else, if there are fewer than eight fields
				add = 6 - y # add that many fields to the taxonomy string
				new = (";" + x.split(";")[-1] + ("_unknwon")) * add
				print(x + new) # print the old plus new taxonomy string
			else: # if the taxonomy string is full (e.g., all taxonomic levels)
				print(x) # don't do anything