# -*- coding: utf-8 -*-
#Oscar Gerardo Lazo Arjona

from commands import *
from os import system
def get_files(path):
	t=getoutput('ls -l '+path).split('\n')
	files=[]; folders=[]
	for i in t:

		f=i.split()[-1]
		if i[0]=='d':
			folders+=[path+'/'+f]
		elif i[0]=='-':
			files+=[path+'/'+f]
	return files,folders

def files_to_delete(files):
	ff=[]
	for f in files:
		#print f
		if (f[-3:] != '.py' and f[-5:] != '.data' and f[-4:] != '.mod' and f[-4:] != '.pbs' and f[-5:] != 'ipynb'):
			band=True
			for ex in file_exceptions:
				if ex in f:
					band=False
			if band:
				ff+=[f]
	return ff

def find_files_to_delete(folder):
	files,folders=get_files(folder)
	#print folders
	
	to_delete=files_to_delete(files)
	#print folder, to_delete

	
	for f in to_delete:
		com= 'rm '+f
		print com
		system(com)
	#print 222,folders
	for folder in folders:
		if folder not in folder_exceptions:
			find_files_to_delete(folder)
	
folder_exceptions=[]
file_exceptions=['CHANGELOG','README']
folder='.'

find_files_to_delete(folder)







def count_lines(f):
	f=file(f,'r')
	n_lines=len(f.readlines())
	f.close()
	return n_lines

#~ for f in files:
	#~ n=count_lines(f)
	#~ s0+=n
	#~ print n,f
#~ print s0,'lines for the model.'
#~ print
#~ s1=0
#~ files=get_python_files(path+'examples/')
#~ for f in files:
	#~ n=count_lines(f)
	#~ s1+=n
	#~ print n,f
#~ print s1,'lines for examples.'
#~ 
#~ print s0+s1,'total lines.'
