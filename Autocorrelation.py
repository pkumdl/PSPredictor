# -*- coding: utf-8 -*-
"""
##########################################################################################

This module is used for computing the Autocorrelation descriptors based different

 properties of AADs.You can also input your properties of AADs, then it can help you

to compute Autocorrelation descriptors based on the property of AADs. Currently, You 

can get 720 descriptors for a given protein sequence based on our provided physicochemical

properties of AADs. You can freely use and distribute it. If you hava  any problem, 

you could contact with us timely!

References:

[1]: http://www.genome.ad.jp/dbget/aaindex.html

[2]:Feng, Z.P. and Zhang, C.T. (2000) Prediction of membrane protein types based on

the hydrophobic index of amino acids. J Protein Chem, 19, 269-275.

[3]:Horne, D.S. (1988) Prediction of protein helix content from an autocorrelation

analysis of sequence hydrophobicities. Biopolymers, 27, 451-477.

[4]:Sokal, R.R. and Thomson, B.A. (2006) Population structure inferred by local

spatial autocorrelation: an Usage from an Amerindian tribal population. Am J

Phys Anthropol, 129, 121-131.

Authors: Dongsheng Cao and Yizeng Liang.

Date: 2010.11.22

Email: oriental-cds@163.com

##########################################################################################
"""

import math,string
import numpy as np
import sys
import re

AALetter=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

_Hydrophobicity={"A":0.62,"R":-2.53,"N":-0.78,"D":-0.9,"C":0.29,"Q":-0.85,"E":-0.74,"G":0.48,"H":-0.4,"I":1.38,"L":1.06,"K":-1.5,"M":0.64,"F":1.19,"P":0.12,"S":-0.18,"T":-0.05,"W":0.81,"Y":0.26,"V":1.08}

_Hydrophilicity={"A":-0.5,"R":3,"N":2,"D":3,"C":-1,"Q":0.2,"E":3,"G":0,"H":-0.5,"I":-1.8,"L":-1.8,"K":3,"M":-1.3,"F":-2.5,"P":0,"S":0.3,"T":-0.4,"W":-3.4,"Y":-2.3,"V":-1.5}

_Polarizability={"A":0.046,"R":0.291,"N":0.134,"D":0.105,"C":0.128,"Q":0.180,"E":0.151,"G":0.000,"H":0.230,"I":0.186,"L":0.186,"K":0.219,"M":0.221,"F":0.290,"P":0.131,"S":0.062,"T":0.108,"W":0.409,"Y":0.298,"V":0.140}

_SidechainV={"A":27.5,"R":105,"N":58.7,"D":40,"C":44.6,"Q":80.7,"E":62,"G":0,"H":79,"I":93.5,"L":93.5,"K":100,"M":94.1,"F":115.5,"P":41.9,"S":29.3,"T":51.3,"W":145.5,"Y":117.3,"V":71.5}

_ResidueASA={"A":1.181,"R":2.56,"N":1.655,"D":1.587,"C":1.461,"Q":1.932,"E":1.862,"G":0.881,"H":2.025,"I":1.81,"L":1.931,"K":2.258,"M":2.034,"F":2.228,"P":1.468,"S":1.298,"T":1.525,"W":2.663,"Y":2.368,"V":1.645}

_Polarity={"A":8.1,"R":10.5,"N":11.6,"D":13,"C":5.5,"Q":10.5,"E":12.3,"G":9,"H":10.4,"I":5.2,"L":4.9,"K":11.3,"M":5.7,"F":5.2,"P":8,"S":9.2,"T":8.6,"W":5.4,"Y":6.2,"V":5.9}

_NCI={"A":0.007187,"R":0.043587,"N":0.005392,"D":-0.02382,"C":-0.03661,"Q":0.049211,"E":0.006802,"G":0.179052,"H":-0.01069,"I":0.021631,"L":0.051672,"K":0.017708,"M":0.002683,"F":0.037552,"P":0.239531,"S":0.004627,"T":0.003352,"W":0.037977,"Y":0.023599,"V":0.057004}

#_Mutability={"A":100.0,"R":65.0,"N":134.0,"D":106.0,"C":20.0,"Q":93.0,"E":102.0,"G":49.0,"H":66.0,"I":96.0,"L":40.0,"K":-56.0,"M":94.0,"F":41.0,"P":56.0,"S":120.0,"T":97.0,"W":18.0,"Y":41.0,"V":74.0}


###You can continuely add other properties of AADs to compute the descriptors of protein sequence.


_AAProperty=(_Hydrophobicity,_Hydrophilicity,_Polarizability,_SidechainV,_ResidueASA,_Polarity,_NCI)

_AAPropertyName=('_Hydrophobicity','_Hydrophilicity','_Polarizability','_SidechainV','_ResidueASA','_Polarity','_NCI')			 

##################################################################################################
def _mean(listvalue):
	"""
	The mean value of the list data.
	"""
	return sum(listvalue)/len(listvalue)
##################################################################################################
def _std(listvalue,ddof=1):
	"""
	The standard deviation of the list data.
	"""
	mean=_mean(listvalue)
	temp=[math.pow(i-mean,2) for i in listvalue]
	res=math.sqrt(sum(temp)/(len(listvalue)-ddof))
	return res
##################################################################################################

def NormalizeEachAAP(AAP):
	"""
	####################################################################################
	All of the amino acid indices are centralized and 
	
	standardized before the calculation.
	
	Usage:
	
	result=NormalizeEachAAP(AAP)
	
	Input: AAP is a dict form containing the properties of 20 amino acids.
	
	Output: result is the a dict form containing the normalized properties 
	
	of 20 amino acids.
	####################################################################################
	"""
	if len(AAP.values())!=20:
		print('You can not input the correct number of properities of Amino acids!')
	else:
		Result={}
		for i,j in AAP.items():
			Result[i]=(j-_mean(AAP.values()))/_std(AAP.values(),ddof=0)

	return Result

def CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,AAP,AAPName):
		
	"""
	####################################################################################
	you can use the function to compute MoreauBrotoAuto
	
	descriptors for different properties based on AADs.
	
	Usage:
	
	result=CalculateEachNormalizedMoreauBrotoAuto(protein,AAP,AAPName)
	
	Input: protein is a pure protein sequence.
	
	AAP is a dict form containing the properties of 20 amino acids (e.g., _Hydrophilicity).
	
	AAPName is a string used for indicating the property (e.g., '_Hydrophilicity'). 
	
	Output: result is a dict form containing 30 Normalized Moreau-Broto autocorrelation
	
	descriptors based on the given property.
	####################################################################################
	"""
		
	AAPdic=NormalizeEachAAP(AAP)

	Result={}
	for i in range(1,31):
		temp=0
		for j in range(len(ProteinSequence)-i):
			temp=temp+AAPdic[ProteinSequence[j]]*AAPdic[ProteinSequence[j+1]]
		if len(ProteinSequence)-i==0:
			Result['MoreauBrotoAuto'+AAPName+str(i)]=round(temp/(len(ProteinSequence)),3)
		else:
			Result['MoreauBrotoAuto'+AAPName+str(i)]=round(temp/(len(ProteinSequence)-i),3)

	return Result


def CalculateEachMoranAuto(ProteinSequence,AAP,AAPName):

	"""
	####################################################################################
	you can use the function to compute MoranAuto
	
	descriptors for different properties based on AADs.
	
	Usage:
	
	result=CalculateEachMoranAuto(protein,AAP,AAPName)
	
	Input: protein is a pure protein sequence.
	
	AAP is a dict form containing the properties of 20 amino acids (e.g., _Hydrophilicity).
	
	AAPName is a string used for indicating the property (e.g., '_Hydrophilicity'). 
	
	Output: result is a dict form containing 30 Moran autocorrelation
	
	descriptors based on the given property.
	####################################################################################
	"""

	AAPdic=NormalizeEachAAP(AAP)

	cds=0
	for i in AALetter:
		cds=cds+(ProteinSequence.count(i))*(AAPdic[i])
	Pmean=cds/len(ProteinSequence)

	cc=[]
	for i in ProteinSequence:
		cc.append(AAPdic[i])

	K=(_std(cc,ddof=0))**2

	Result={}
	for i in range(1,31):
		temp=0
		for j in range(len(ProteinSequence)-i):
				
			temp=temp+(AAPdic[ProteinSequence[j]]-Pmean)*(AAPdic[ProteinSequence[j+i]]-Pmean)
		if len(ProteinSequence)-i==0:
			Result['MoranAuto'+AAPName+str(i)]=round(temp/(len(ProteinSequence))/K,3)
		else:
			Result['MoranAuto'+AAPName+str(i)]=round(temp/(len(ProteinSequence)-i)/K,3)

	return Result


def CalculateEachGearyAuto(ProteinSequence,AAP,AAPName):

	"""
	####################################################################################
	you can use the function to compute GearyAuto
	
	descriptors for different properties based on AADs.
	
	Usage:
	
	result=CalculateEachGearyAuto(protein,AAP,AAPName)
	
	Input: protein is a pure protein sequence.
	
	AAP is a dict form containing the properties of 20 amino acids (e.g., _Hydrophilicity).
	
	AAPName is a string used for indicating the property (e.g., '_Hydrophilicity'). 
	
	Output: result is a dict form containing 30 Geary autocorrelation
	
	descriptors based on the given property.
	####################################################################################
	"""

	AAPdic=NormalizeEachAAP(AAP)

	cc=[]
	for i in ProteinSequence:
		cc.append(AAPdic[i])

	K=((_std(cc))**2)*len(ProteinSequence)/(len(ProteinSequence)-1)
	Result={}
	for i in range(1,31):
		temp=0
		for j in range(len(ProteinSequence)-i):
				
			temp=temp+(AAPdic[ProteinSequence[j]]-AAPdic[ProteinSequence[j+i]])**2
		if len(ProteinSequence)-i==0:
			Result['GearyAuto'+AAPName+str(i)]=round(temp/(2*(len(ProteinSequence)))/K,3)
		else:
			Result['GearyAuto'+AAPName+str(i)]=round(temp/(2*(len(ProteinSequence)-i))/K,3)
	return Result


##################################################################################################

def CalculateNormalizedMoreauBrotoAuto(ProteinSequence,AAProperty,AAPropertyName): 

	"""
	####################################################################################
	A method used for computing MoreauBrotoAuto for all properties.
	
	Usage:
	
	result=CalculateNormalizedMoreauBrotoAuto(protein,AAP,AAPName)
	
	Input: protein is a pure protein sequence.
	
	AAProperty is a list or tuple form containing the properties of 20 amino acids (e.g., _AAProperty).
	
	AAPName is a list or tuple form used for indicating the property (e.g., '_AAPropertyName'). 
	
	Output: result is a dict form containing 30*p Normalized Moreau-Broto autocorrelation
	
	descriptors based on the given properties.
	####################################################################################
	
	"""
	Result={}
	for i in range(len(AAProperty)):
		Result[AAPropertyName[i]]=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,AAProperty[i],AAPropertyName[i])


	return Result


def CalculateMoranAuto(ProteinSequence,AAProperty,AAPropertyName):  
	"""
	####################################################################################
	A method used for computing MoranAuto for all properties
	
	Usage:
	
	result=CalculateMoranAuto(protein,AAP,AAPName)
	
	Input: protein is a pure protein sequence.
	
	AAProperty is a list or tuple form containing the properties of 20 amino acids (e.g., _AAProperty).
	
	AAPName is a list or tuple form used for indicating the property (e.g., '_AAPropertyName'). 
	
	Output: result is a dict form containing 30*p Moran autocorrelation
	
	descriptors based on the given properties.
	####################################################################################
	"""
	Result={}
	for i in range(len(AAProperty)):
		Result[AAPropertyName[i]]=CalculateEachMoranAuto(ProteinSequence,AAProperty[i],AAPropertyName[i])

	return Result



def CalculateGearyAuto(ProteinSequence,AAProperty,AAPropertyName):  
	"""
	####################################################################################
	A method used for computing GearyAuto for all properties
	
	Usage:
	
	result=CalculateGearyAuto(protein,AAP,AAPName)
	
	Input: protein is a pure protein sequence.
	
	AAProperty is a list or tuple form containing the properties of 20 amino acids (e.g., _AAProperty).
	
	AAPName is a list or tuple form used for indicating the property (e.g., '_AAPropertyName'). 
	
	Output: result is a dict form containing 30*p Geary autocorrelation
	
	descriptors based on the given properties.
	####################################################################################
	"""
	Result={}
	for i in range(len(AAProperty)):
		Result[AAPropertyName[i]]=CalculateEachGearyAuto(ProteinSequence,AAProperty[i],AAPropertyName[i])

	return Result


########################NormalizedMoreauBorto##################################
def CalculateNormalizedMoreauBrotoAutoHydrophobicity(ProteinSequence):

	"""
	####################################################################################
	Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on
	
	hydrophobicity.
	
	Usage:
	
	result=CalculateNormalizedMoreauBrotoAutoHydrophobicity(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Normalized Moreau-Broto Autocorrelation
	
	descriptors based on Hydrophobicity.
	####################################################################################
	"""
	
	result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
	return result


def CalculateNormalizedMoreauBrotoAutoHydrophilicity(ProteinSequence):

	"""
	####################################################################################
	Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on
	
	Hydrophilicity.
	
	Usage:
	
	result=CalculateNormalizedMoreauBrotoAutoHydrophilicity(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Normalized Moreau-Broto Autocorrelation
	
	descriptors based on Hydrophilicity.
	####################################################################################
	"""
	
	result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_Hydrophilicity,'_Hydrophilicity')
	return result


def CalculateNormalizedMoreauBrotoAutoPolarizability(ProteinSequence):

	"""
	####################################################################################
	Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on
	
	Polarizability.
	
	Usage:
	
	result=CalculateNormalizedMoreauBrotoAutoPolarizability(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Normalized Moreau-Broto Autocorrelation
	
	descriptors based on Polarizability.
	####################################################################################
	"""
	
	result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_Polarizability,'_Polarizability')
	return result


def CalculateNormalizedMoreauBrotoAutoSidechainV(ProteinSequence):

	"""
	####################################################################################
	Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on
	
	SidechainV.
	
	Usage:
	
	result=CalculateNormalizedMoreauBrotoAutoSidechainV(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Normalized Moreau-Broto Autocorrelation
	
	descriptors based on SidechainV.
	####################################################################################
	"""
	
	result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_SidechainV,'_SidechainV')
	return result



def CalculateNormalizedMoreauBrotoAutoResidueASA(ProteinSequence):

	"""
	####################################################################################
	Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on
	
	ResidueASA.
	
	Usage:
	
	result=CalculateNormalizedMoreauBrotoAutoResidueASA(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Normalized Moreau-Broto Autocorrelation
	
	descriptors based on ResidueASA.
	####################################################################################
	"""
	
	result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_ResidueASA,'_ResidueASA')
	return result


def CalculateNormalizedMoreauBrotoAutoPolarity(ProteinSequence):

	"""
	####################################################################################
	Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on
	
	Polarity.
	
	Usage:
	
	result=CalculateNormalizedMoreauBrotoAutoPolarity(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Normalized Moreau-Broto Autocorrelation
	
	descriptors based on Polarity.
	####################################################################################
	"""
	
	result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_Polarity,'_Polarity')
	return result
	
def CalculateNormalizedMoreauBrotoAutoNCI(ProteinSequence):

	"""
	####################################################################################
	Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on NCI.
	
	Usage:
	
	result=CalculateNormalizedMoreauBrotoAutoNCI(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Normalized Moreau-Broto Autocorrelation
	
	descriptors based on NCI.
	####################################################################################
	"""
	
	result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_NCI,'_NCI')
	return result


def CalculateNormalizedMoreauBrotoAutoMutability(ProteinSequence):

	"""
	####################################################################################
	Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on Mutability.
	
	Usage:
	
	result=CalculateNormalizedMoreauBrotoAutoMutability(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Normalized Moreau-Broto Autocorrelation
	
	descriptors based on Mutability.
	####################################################################################
	"""
	
	result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_Mutability,'_Mutability')
	return result
############################################################################

##############################MoranAuto######################################
def CalculateMoranAutoHydrophobicity(ProteinSequence):

	"""
	####################################################################################
	Calculte the MoranAuto Autocorrelation descriptors based on hydrophobicity.
	
	Usage:
	
	result=CalculateMoranAutoHydrophobicity(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Moran Autocorrelation
	
	descriptors based on hydrophobicity.
	####################################################################################
	"""
	
	result=CalculateEachMoranAuto(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
	return result
	

def CalculateMoranAutoHydrophilicity(ProteinSequence):

	"""
	####################################################################################
	Calculte the MoranAuto Autocorrelation descriptors based on
	
	Hydrophilicity.
	
	Usage:
	
	result=CalculateMoranAutoHydrophilicity(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Moran Autocorrelation
	
	descriptors based on Hydrophilicity.
	####################################################################################
	"""
	
	result=CalculateEachMoranAuto(ProteinSequence,_Hydrophilicity,'_Hydrophilicity')
	return result


def CalculateMoranAutoPolarizability(ProteinSequence):

	"""
	####################################################################################
	Calculte the MoranAuto Autocorrelation descriptors based on
	
	Polarizability.
	
	Usage:
	
	result=CalculateMoranAutoPolarizability(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Moran Autocorrelation
	
	descriptors based on Polarizability.
	####################################################################################
	"""
	
	result=CalculateEachMoranAuto(ProteinSequence,_Polarizability,'_Polarizability')
	return result


def CalculateMoranAutoSidechainV(ProteinSequence):

	"""
	####################################################################################
	Calculte the MoranAuto Autocorrelation descriptors based on
	
	SidechainV.
	
	Usage:
	
	result=CalculateMoranAutoSidechainV(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Moran Autocorrelation
	
	descriptors based on SidechainV.
	####################################################################################
	"""
	
	result=CalculateEachMoranAuto(ProteinSequence,_SidechainV,'_SidechainV')
	return result



def CalculateMoranAutoResidueASA(ProteinSequence):

	"""
	####################################################################################
	Calculte the MoranAuto Autocorrelation descriptors based on
	
	ResidueASA.
	
	Usage:
	
	result=CalculateMoranAutoResidueASA(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Moran Autocorrelation
	
	descriptors based on ResidueASA.
	####################################################################################
	"""
	
	result=CalculateEachMoranAuto(ProteinSequence,_ResidueASA,'_ResidueASA')
	return result


def CalculateMoranAutoPolarity(ProteinSequence):

	"""
	####################################################################################
	Calculte the MoranAuto Autocorrelation descriptors based on
	
	Polarity.
	
	Usage:
	
	result=CalculateMoranAutoPolarity(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Moran Autocorrelation
	
	descriptors based on Polarity.
	####################################################################################
	"""
	
	result=CalculateEachMoranAuto(ProteinSequence,_Polarity,'_Polarity')
	return result
	
def CalculateMoranAutoNCI(ProteinSequence):

	"""
	####################################################################################
	Calculte the MoranAuto Autocorrelation descriptors based on
	
	AutoNCI.
	
	Usage:
	
	result=CalculateMoranAutoNCI(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Moran Autocorrelation
	
	descriptors based on AutoNCI.
	####################################################################################
	"""
	
	result=CalculateEachMoranAuto(ProteinSequence,_NCI,'_NCI')
	return result


def CalculateMoranAutoMutability(ProteinSequence):

	"""
	####################################################################################
	Calculte the MoranAuto Autocorrelation descriptors based on
	
	Mutability.
	
	Usage:
	
	result=CalculateMoranAutoMutability(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Moran Autocorrelation
	
	descriptors based on Mutability.
	####################################################################################
	"""
	
	result=CalculateEachMoranAuto(ProteinSequence,_Mutability,'_Mutability')
	return result
############################################################################


################################GearyAuto#####################################
def CalculateGearyAutoHydrophobicity(ProteinSequence):

	"""
	####################################################################################
	Calculte the GearyAuto Autocorrelation descriptors based on
	
	hydrophobicity.
	
	Usage:
	
	result=CalculateGearyAutoHydrophobicity(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Geary Autocorrelation
	
	descriptors based on hydrophobicity.
	####################################################################################
	"""
	
	result=CalculateEachGearyAuto(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
	return result
	

def CalculateGearyAutoHydrophilicity(ProteinSequence):

	"""
	####################################################################################
	Calculte the GearyAuto Autocorrelation descriptors based on
	
	Hydrophilicity.
	
	Usage:
	result=CalculateGearyAutoHydrophilicity(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Geary Autocorrelation
	
	descriptors based on Hydrophilicity.
	####################################################################################
	"""
	
	result=CalculateEachGearyAuto(ProteinSequence,_Hydrophilicity,'_Hydrophilicity')
	return result


def CalculateGearyAutoPolarizability(ProteinSequence):

	"""
	####################################################################################
	Calculte the GearyAuto Autocorrelation descriptors based on
	
	Polarizability.
	
	Usage:
	
	result=CalculateGearyAutoPolarizability(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Geary Autocorrelation
	
	descriptors based on Polarizability.
	####################################################################################
	"""
	
	result=CalculateEachGearyAuto(ProteinSequence,_Polarizability,'_Polarizability')
	return result


def CalculateGearyAutoSidechainV(ProteinSequence):
	
	"""
	####################################################################################
	Calculte the GearyAuto Autocorrelation descriptors based on
	
	SidechainV.
	
	Usage:
	
	result=CalculateGearyAutoSidechainV(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Geary Autocorrelation
	
	descriptors based on SidechainV.
	####################################################################################
	"""
	
	result=CalculateEachGearyAuto(ProteinSequence,_SidechainV,'_SidechainV')
	return result



def CalculateGearyAutoResidueASA(ProteinSequence):
	
	"""
	####################################################################################
	Calculte the GearyAuto Autocorrelation descriptors based on
	
	ResidueASA.
	
	Usage:
	
	result=CalculateGearyAutoResidueASA(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Geary Autocorrelation
	
	descriptors based on ResidueASA.
	####################################################################################
	"""
	
	result=CalculateEachGearyAuto(ProteinSequence,_ResidueASA,'_ResidueASA')
	return result


def CalculateGearyAutoPolarity(ProteinSequence):
	
	"""
	####################################################################################
	Calculte the GearyAuto Autocorrelation descriptors based on
	
	Polarity.
	
	Usage:
	
	result=CalculateGearyAutoPolarity(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Geary Autocorrelation
	
	descriptors based on Polarity.
	####################################################################################
	"""
	
	result=CalculateEachGearyAuto(ProteinSequence,_Polarity,'_Polarity')
	return result
	
def CalculateGearyAutoNCI(ProteinSequence):
	
	"""
	####################################################################################
	Calculte the GearyAuto Autocorrelation descriptors based on
	
	NCI.
	
	Usage:
	
	result=CalculateGearyAutoNCI(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Geary Autocorrelation
	
	descriptors based on NCI.
	####################################################################################
	"""
	
	result=CalculateEachGearyAuto(ProteinSequence,_NCI,'_NCI')
	return result


def CalculateGearyAutoMutability(ProteinSequence):

	"""
	####################################################################################
	Calculte the GearyAuto Autocorrelation descriptors based on
	
	Mutability.
	
	Usage:
	
	result=CalculateGearyAutoMutability(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30 Geary Autocorrelation
	
	descriptors based on Mutability.
	####################################################################################
	"""
	
	result=CalculateEachGearyAuto(ProteinSequence,_Mutability,'_Mutability')
	return result
##################################################################################################

def CalculateNormalizedMoreauBrotoAutoTotal(ProteinSequence):
	"""
	####################################################################################
	A method used for computing normalized Moreau Broto autocorrelation descriptors based 
	
	on 8 proterties of AADs.
	
	Usage:
	
	result=CalculateNormalizedMoreauBrotoAutoTotal(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30*8=240 normalized Moreau Broto 
	
	autocorrelation descriptors based on the given properties(i.e., _AAPropert).
	#################################################################################### 
	"""
	result={}
	result.update(CalculateNormalizedMoreauBrotoAutoHydrophobicity(ProteinSequence))
	result.update(CalculateNormalizedMoreauBrotoAutoHydrophilicity(ProteinSequence))
	result.update(CalculateNormalizedMoreauBrotoAutoPolarizability(ProteinSequence))
	result.update(CalculateNormalizedMoreauBrotoAutoSidechainV(ProteinSequence))
	result.update(CalculateNormalizedMoreauBrotoAutoResidueASA(ProteinSequence))
	result.update(CalculateNormalizedMoreauBrotoAutoPolarity(ProteinSequence))
	result.update(CalculateNormalizedMoreauBrotoAutoNCI(ProteinSequence))
	#result.update(CalculateNormalizedMoreauBrotoAutoMutability(ProteinSequence))
	return result

def CalculateMoranAutoTotal(ProteinSequence):
	"""
	####################################################################################
	A method used for computing Moran autocorrelation descriptors based on 8 properties of AADs.
	
	Usage:
	
	result=CalculateMoranAutoTotal(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30*8=240 Moran
	
	autocorrelation descriptors based on the given properties(i.e., _AAPropert).
	####################################################################################
	"""
	result={}
	result.update(CalculateMoranAutoHydrophobicity(ProteinSequence))
	result.update(CalculateMoranAutoHydrophilicity(ProteinSequence))
	result.update(CalculateMoranAutoPolarizability(ProteinSequence))
	result.update(CalculateMoranAutoSidechainV(ProteinSequence))
	result.update(CalculateMoranAutoResidueASA(ProteinSequence))
	result.update(CalculateMoranAutoPolarity(ProteinSequence))
	result.update(CalculateMoranAutoNCI(ProteinSequence))
#	result.update(CalculateMoranAutoMutability(ProteinSequence))
	return result

def CalculateGearyAutoTotal(ProteinSequence):
	"""
	####################################################################################
	A method used for computing Geary autocorrelation descriptors based on 8 properties of AADs.
	
	Usage:
	
	result=CalculateGearyAutoTotal(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30*8=240 Geary
	
	autocorrelation descriptors based on the given properties(i.e., _AAPropert).
	####################################################################################
	"""
	result={}
	result.update(CalculateGearyAutoHydrophobicity(ProteinSequence))
	result.update(CalculateGearyAutoHydrophilicity(ProteinSequence))
	result.update(CalculateGearyAutoPolarizability(ProteinSequence))
	result.update(CalculateGearyAutoSidechainV(ProteinSequence))
	result.update(CalculateGearyAutoResidueASA(ProteinSequence))
	result.update(CalculateGearyAutoPolarity(ProteinSequence))
	result.update(CalculateGearyAutoNCI(ProteinSequence))
	result.update(CalculateGearyAutoMutability(ProteinSequence))
	return result

##################################################################################################
def CalculateAutoTotal(ProteinSequence):
	"""
	####################################################################################
	A method used for computing all autocorrelation descriptors based on 8 properties of AADs.
	
	Usage:
	
	result=CalculateGearyAutoTotal(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing 30*8*3=720 normalized Moreau Broto, Moran, and Geary
	
	autocorrelation descriptors based on the given properties(i.e., _AAPropert).
	####################################################################################
	"""
	result={}
	result.update(CalculateNormalizedMoreauBrotoAutoTotal(ProteinSequence))
	result.update(CalculateMoranAutoTotal(ProteinSequence))
	result.update(CalculateGearyAutoTotal(ProteinSequence))
	return result

##################################################################################################
#if __name__=="__main__":
#	GetSequenceFromTxt("/lustre1/lhlai_pkuhpc/tlsun/data/propy-1.0/src/propy/",sys.argv[1],sys.argv[2])	
	#protein="ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
	#temp1=CalculateNormalizedMoreauBrotoAuto(protein,AAProperty=_AAProperty,AAPropertyName=_AAPropertyName)
	#print temp1
	#temp2=CalculateMoranAutoMutability(protein)
	#print temp2
	#print len(CalculateAutoTotal(protein))

##################################################################################################
def GetSequenceFromTxt1(data):
	#aaa=f2.read()
	ll1=[]
	#flags=re.search(r'[^ACDEFGHIKLMNPQRSTVWY\n]',aaa,re.I)
	#if flags:
	#	print("please check")
	#else:
	for index,i in enumerate(data):
		itrim=str.strip(i)
	   # flags=re.search(r'[^ACDEFGHIKLMNPQRSTVWY]',itrim,re.I)
	   # if flags:
	   #		 print("please check the sequences for characters other than 20 amino acids")
	   #		 sys.exit(1)
		#else:
		try:
			D=CalculateMoranAutoTotal(itrim)
		except:
			print(i)
			continue
		print("finish %d" %(index+1))
		tt=list(D.values())
		ll1.append(tt)
	ll2=np.array(ll1)
	return ll2
		#f1=file(path+savefile,'wb')
		#a=np.loadtxt(path+savefile)
		#c,d=a.shape
		#print(c,d)
		#b=np.empty([c,d])
		#b[:][:]=(a[:][:]-a[:][:].min())/(a[:][:].max()-a[:][:].min())
		#print(b)
		#np.savetxt('result',b)
if __name__=="__main__":
	GetSequenceFromTxt("/lustre1/lhlai_pkuhpc/tlsun/drug_prediction/large_data/neg3")












				

