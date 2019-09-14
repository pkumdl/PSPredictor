import sys
from CTD_zong_all import calculate
import hyper as hp
import numpy as np
from sklearn import metrics,preprocessing
from AAComposition1 import GetSequenceFromTxt
from Autocorrelation3 import GetSequenceFromTxt1
from random import shuffle
import pandas as pd
from embed1 import ol_gram,embedding,hla
import random
hla_weight=hla('my_wordvec_model_trembl_size_200')

def norm(x):
	min_max_scaler = preprocessing.MinMaxScaler()
	x_n= min_max_scaler.fit_transform(x)
	#print(min_max_scaler.get_params())

	return x_n

def process(method,r,name):
	f_p=open(hp.train_p)
	f_n=open(hp.train_n)
	
	p=[i.strip() for i in f_p.readlines()]
	n=[i.strip() for i in f_n.readlines()]

	s1=len(p)
	s2=len(n)
	
	if s1 > s2:
		if r==0:
			p1=p
			n1=n
		else:
			shuffle(p)
			p1=p[:int(len(n)*r)]
			n1=n
	else:
		shuffle(n)
		n1=n[:len(p)*r]
		p1=p

	p2=pd.DataFrame(p1)
	n2=pd.DataFrame(n1)

	p_label=pd.DataFrame(np.ones([p2.shape[0],1]))
	n_label=pd.DataFrame(np.zeros([n2.shape[0],1]))
	
	p_a=pd.concat([p2,p_label],axis=1)
	n_a=pd.concat([n2,n_label],axis=1)

	all=pd.concat([p_a,n_a],axis=0)
	all.columns=['seq','label']
	all1=all.sample(frac=1)
	seq=all1['seq'].tolist()
	label=all1['label']
	if method=='li':
		code=calculate(seq)
		pd.DataFrame(code).to_csv(name,header=None,index=None)
	elif method=='TA':
		code=GetSequenceFromTxt(seq)
		pd.DataFrame(code).to_csv(name,header=None,index=None)
	elif method =='AC':
		code=GetSequenceFromTxt1(seq)
		pd.DataFrame(code).to_csv(name,header=None,index=None)
	elif method=='w2v':
		pre=ol_gram(seq)
		code=embedding(pre,hla_weight)
		pd.DataFrame(code).to_csv(name,header=None,index=None)
	all_norm=pd.DataFrame(norm(code))
	#all_norm1=np.array(pd.concat([pd.DataFrame(all_norm),all1['label']],axis=1))
	
	return all_norm,label,seq

if __name__=='__main__':
	process('AC',1,'AC')
