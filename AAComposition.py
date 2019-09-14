# -*- coding: utf-8 -*-
"""
###############################################################################

The module is used for computing the composition of amino acids, dipetide and 

3-mers (tri-peptide) for a given protein sequence. You can get 8420 descriptors 

for a given protein sequence. You can freely use and distribute it. If you hava 

any problem, you could contact with us timely!

References:

[1]: Reczko, M. and Bohr, H. (1994) The DEF data base of sequence based protein

fold class predictions. Nucleic Acids Res, 22, 3616-3619.

[2]: Hua, S. and Sun, Z. (2001) Support vector machine approach for protein

subcellular localization prediction. Bioinformatics, 17, 721-728.


[3]:Grassmann, J., Reczko, M., Suhai, S. and Edler, L. (1999) Protein fold class

prediction: new methods of statistical classification. Proc Int Conf Intell Syst Mol

Biol, 106-112.

Authors: Dongsheng Cao and Yizeng Liang.

Date: 2012.3.27

Email: oriental-cds@163.com

###############################################################################
"""

import re
import sys
import numpy as np
#from function_data import *
AALetter=["1","2","3","4","5","6","7"]
#############################################################################################
def CalculateAAComposition(ProteinSequence):

	"""
	########################################################################
	Calculate the composition of Amino acids 
	
	for a given protein sequence.
	
	Usage:
	
	result=CalculateAAComposition(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing the composition of 
	
	20 amino acids.
	########################################################################
	"""
	LengthSequence=len(ProteinSequence)
	Result={}
	for i in AALetter:
		Result[i]=round(float(ProteinSequence.count(i))/LengthSequence*100,3)
	return Result

#############################################################################################
def CalculateDipeptideComposition(ProteinSequence):
	"""
	########################################################################
	Calculate the composition of dipeptidefor a given protein sequence.
	
	Usage: 
	
	result=CalculateDipeptideComposition(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing the composition of 
	
	400 dipeptides.
	########################################################################
	"""

	LengthSequence=len(ProteinSequence)
	Result={}
	for i in AALetter:
		for j in AALetter:
			Dipeptide=i+j
			Result[Dipeptide]=round(float(ProteinSequence.count(Dipeptide))/(LengthSequence-1)*100,2)
	return Result



#############################################################################################

def Getkmers():
	"""
	########################################################################
	Get the amino acid list of 3-mers. 
	
	Usage: 
	
	result=Getkmers()
	
	Output: result is a list form containing 8000 tri-peptides.
	
	########################################################################
	"""
	kmers=list()
	for i in AALetter:
		for j in AALetter:
			for k in AALetter:
				kmers.append(i+j+k)
	return kmers

#############################################################################################
def GetSpectrumDict(proteinsequence):
	"""
	########################################################################
	Calcualte the spectrum descriptors of 3-mers for a given protein.
	
	Usage: 
	
	result=GetSpectrumDict(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing the composition values of 8000
	
	3-mers.
	########################################################################
	"""
	result={}
	kmers=Getkmers()
	for i in kmers:
		result[i]=len(re.findall(i,proteinsequence))
	return result

#############################################################################################
def CalculateAADipeptideComposition(ProteinSequence):

	"""
	########################################################################
	Calculate the composition of AADs, dipeptide and 3-mers for a 
	
	given protein sequence.
	
	Usage:
	
	result=CalculateAADipeptideComposition(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing all composition values of 
	
	AADs, dipeptide and 3-mers (8420).
	########################################################################
	"""

	result={}
	result.update(CalculateAAComposition(ProteinSequence))
	result.update(CalculateDipeptideComposition(ProteinSequence))
	result.update(GetSpectrumDict(ProteinSequence))

	return result

def rep(text):
	p1=re.compile(r'[A|G|V]')
	p2=re.compile(r'[I|L|P|F]')
	p3=re.compile(r'[Y|M|T|S]')
	p4=re.compile(r'[H|Q|N|W]')
	p5=re.compile(r'[R|K]')
	p6=re.compile(r'[D|E]')
	p7=re.compile(r'C')
	text=re.sub(p1,'1',text)
	text=re.sub(p2,'2',text)
	text=re.sub(p3,'3',text)
	text=re.sub(p4,'4',text)
	text=re.sub(p5,'5',text)
	text=re.sub(p6,'6',text)
	text=re.sub(p7,'7',text)
	return text

def norm_special(array):
	s1=array.shape[0]
	s2=array.shape[1]
	b=np.zeros((s1,s2))
	for i in range(s1):
		for j in range(s2):
			b[i,j]=float(array[i,j])/float(np.sum(array[i,:]))
	return b


#############################################################################################
###############################################################################################
def GetSequenceFromTxt(data):
	text1='\t'.join(data)
	text3=rep(text1)
	list1=text3.split('\t')
	a=calculate(list1)
	return a

def calculate(list1):
	ll1=[]
	for index,i in enumerate(list1):
		itrim=str.strip(i)
		if itrim == "":
			continue
		else:
			Des=GetSpectrumDict(itrim)
			a1=list(Des.values())
			ll1.append(a1)
	ll2=np.array(ll1)
	ll3=norm_special(ll2)
	return ll3


if __name__=="__main__":
	GetSequenceFromTxt("/home/lhlai_pkuhpc/lustre1/tlsun/drug_prediction/large_data/TA/",sys.argv[1],sys.argv[2])
