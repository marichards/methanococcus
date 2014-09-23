# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 18:31:00 2014

@author: Matt
"""

#Import the packages I'll need
import urllib2
from bs4 import BeautifulSoup
import re
import contextlib
#Should take in the list of reaction IDs
txt = open('rxns_for_subs.txt','r')
#Also open a new file to write to
new = open('pulled_subs.txt','w')

for item in txt:
    
    #Strip the e0 or c0
    rxn = re.sub('_e\n0|_c0\n','',item)
    
    #For each ID in reaction IDs, open the URL for that one on the SEED
    with contextlib.closing(urllib2.urlopen("http://seed-viewer.theseed.org/seedviewer.cgi?page=ReactionViewer&reaction=%s&model=" %rxn)) as page:

        #Make it into the soup
        soup = BeautifulSoup(page)

        #Search for the Subsystems tag and Kegg Maps tag
        #For each one of them, check the tag itself, then pull out the next string
        #Just do for subsystems now
        for tag in soup.find_all(name='th'):
            if tag.string=='Subsystems':
                stuff=tag.next_sibling.string
            
            #Now write it to the new file
        new.write('%s\t%s\n' %(item,stuff))