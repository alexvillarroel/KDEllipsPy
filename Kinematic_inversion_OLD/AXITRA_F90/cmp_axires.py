#!/usr/bin/python

from struct import *
file1="axi.res"
file2="axi.res_new"
file=open(file1,'rb')
file_ori=open(file2,'rb')
donnes=range(0,14)
donnes2=range(0,14)
nfreq=128
print "compare",file1,' with ',file2
for k in range(1,nfreq):
	a1=file.read(4)
	b1=file_ori.read(4)
	a=file.read(96)
	b=file_ori.read(96)
	if k<4 or k>124:
		(length,)=unpack('i',a1)
		(length_ori,)=unpack('i',b1)
		print "k,largo = ",k,length,length_ori
		(donnes)=unpack('12d',a)
		(donnes2)=unpack('12d',b)
		for i in range(0,12,2):
			print "( %12.5e,%12.5e)=?=( %12.5e,%12.5e)"%(donnes[i],donnes[i+1],donnes2[i],donnes2[i+1])
	a=file.read(4)
	b=file_ori.read(4)
