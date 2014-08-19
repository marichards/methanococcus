# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 17:01:09 2014

@author: Matt
"""


#Bring in RegEx
import re
#Load file to parse
infile = open('all_met_formulas.txt','r')
#Load file to write to
outfile = open('formatted_formulas.txt','w')
#Write the first line
outfile.write('model.metFormulas={...\n')
#Cycle through the lines and 
for line in infile:
    #If it has a formula
    if re.match("\w",line):
        #For each line
        new = re.sub(r"(\w+)",r"'\1'",line)
        #Print the line
        outfile.write('%s;...\n' %new.strip())
    #If it doesn't, just print it    
    else:
        outfile.write('%s;...\n' %line.strip())
            
            