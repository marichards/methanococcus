# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 18:31:00 2014

@author: Matt
"""

#Import the packages I'll need
import url
from bs4 import BeautifulSoup
import sys

#Should take in the list of reaction IDs

#For each ID in reaction IDs, open the URL for that one on the SEED
#"http://seed-viewer.theseed.org/seedviewer.cgi?page=ReactionViewer&reaction=rxn00006&model="

#Used this before, I don't know that I need it now
f = file(sys.argv[1],'r')
txt = f.read()


soup = BeautifulSoup(txt)
#print soup.prettify()

tags = soup.findAll(name='p', attrs={'class' : 'title'})

for p in tags:

    # tag = '<p>...</p>'
    # tag.contents[0] = '<a href='...'> </a>'
    
    a = p.contents[0]
    print "a = %s" % a
    print "a.attrs = %s" % a.attrs
    print "link = %s" % a.attrs[0][1]
    print "title = %s" % a.contents[0]