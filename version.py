'''
Created on 21 Oct 2018

@author: thomasgumbricht
'''
__version__ = '0.3.1'
VERSION = tuple( int(x) for x in __version__.split('.') )
metadataD = { 'name':'export', 
             'author':'Thomas Gumbricht', 
             'author_email':'thomas.gumbricht@gmail.com',
             'title':'Export spatial data either for layout or backup.', 
             'label':'Processes for exporting and backing up spatial data, either for publication or backup storage.',
             'prerequisites':'The data to export must exist in the Framework; export to layout requires that translation to binary (0-255) range is defined.',
             'image':'avg-trmm-3b43v7-precip_3B43_trmm_2001-2016_A'}