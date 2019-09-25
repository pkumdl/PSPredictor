import sys


from sklearn.externals import joblib

from sklearn.metrics import classification_report,accuracy_score,f1_score,make_scorer

import sys
from CTD import calculate
import hyper as hp
import numpy as np
from sklearn import metrics,preprocessing
from random import shuffle
import pandas as pd
import re
pat=re.compile(r'No')
from embed1 import ol_gram,embedding,hla
pat1=re.compile(r'(Yes.*)')

hla_weight=hla(hp.weight)


def spe(y_true,y_pred):
	tn=np.sum((y_true==0)&(y_pred==0))
	fp=np.sum((y_true==0)&(y_pred==1))
	s=tn/(tn+fp)
	return s
spe=make_scorer(spe,greater_is_better=True)


def norm(x):
	min_max_scaler = preprocessing.MinMaxScaler()
	x_n= min_max_scaler.fit_transform(x)
	return x_n,min_max_scaler

def process(method,scale):
	data=pd.read_table(hp.test,sep='\t',header=None)
	data.columns=['seq','label']
	seq=data['seq'].tolist()
	label1=data['label'].str.replace(pat,'0')	
	label2=label1.str.replace(pat1,'1')
	label3=label2.astype('int')
	if method=='li':
		code=calculate(seq)
		data=scale.transform(code)
	
	elif method=='w2v':
		pre=ol_gram(seq)
		code=embedding(pre,hla_weight)
		data=scale.transform(code)
	#all_norm1=np.array(pd.concat([pd.DataFrame(all_norm),all1['label']],axis=1))
	
	return data,label3,seq

if __name__=='__main__':
	file1=hp.test
	file2=file1+'predict'
	f=open(file2,'w')
	tn=[]

	data=pd.read_table(hp.code,header=None,sep=',')
	st,scal=norm(data)
	code,label,seq=process(hp.method,scal)
	label_true=np.array(label)
	model=joblib.load(hp.model)

	p1=model.predict_proba(code)
	print(p1)
	pd.DataFrame(p1).to_csv(hp.path+file2,header=None,float_format='%.3f')
