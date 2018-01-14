# -*- coding: utf-8 -*-
import sys

"""
headercreat
在新的csv中创建表头
参数一：新的文件的描述符
参数二：第一个文件的表头所在行
参数三：第二个文件的表头所在行
"""
def headercreate(filename,file1header,file2header):
	datafinal=file1header.split(',')
	datafinallen = len(datafinal)
	for j in range(0,datafinallen-1):
		filename.write(datafinal[j])
		filename.write(',')
	data2=datafinal[datafinallen-1]
	data2=data2[:-1]
	print data2
	filename.write(data2)
	filename.write(',')

	print file2header
	datanorm=file2header.split(',')
	datanormlen = len(datanorm)
	print datanormlen
	for j in range(1,datanormlen-1):
		filename.write(datanorm[j])
		filename.write(',')
	filename.write(datanorm[datanormlen-1])

"""
mergerdata
合并两张表中共有的参数
参数1：新的csv文件的文件描述符
参数2：第一张表中数据
参数3：第一张表中数据开始行
参数4：第一张表中数据结束行
参数5：第二张表中数据
参数6：第二张表中数据开始行
参数7：第二张表中数据结束行
"""
def mergerdata(filename,file1data,file1start,file1end,file2data,file2start,file2end):
	okflag = 0;	
	for i in range(file1start,file1end):
		datafinal=file1data[i].split(',')
		datafinallen = len(datafinal)
		
		for k in range(file2start,file2end):
			datanorm=file2data[k].split(',')
			datanormlen = len(datanorm)
			if datafinal[0] == datanorm[0]:
				okflag=1
				break

		if okflag == 1:
			for j in range(0,datafinallen-1):
				filename.write(datafinal[j])
				filename.write(',')
			data2=datafinal[datafinallen-1]
			data2=data2[:-1]
		#	print data2
			filename.write(data2)
			filename.write(',')

			for j in range(1,datanormlen-1):
				filename.write(datanorm[j])
				filename.write(',')
			filename.write(datanorm[datanormlen-1])

def main():
	ffros=open('final_results_of_MT_vs_PB.csv','rb')
	finaldata = ffros.read()
	ffros.close()
	finallist = finaldata.split('\n')
	finallistlen = len(finallist)

	print finallistlen

	fnorm=open('normcounts.csv','rb')
	normdata = fnorm.read()
	fnorm.close()
	normlist = normdata.split('\n')
	normlistlen = len(normlist)

	print normlistlen
	
	fnew=open('new.csv','wb')
	
	headercreate(fnew,finallist[0],normlist[0]);#header在csv表格的第几行，方括号中的数字就改为几
	mergerdata(fnew,finallist,1,finallistlen,normlist,1,normlistlen);#合并两张csv的数据
	
	fnew.close()
	print "ok" 
		
		

if __name__=="__main__":
	main()
