# PSPredict
prediction of phase separation proteins  
written by Tanlin Sun

# how to prepare sequences
Prepare sequences with its label separting by tab("Yes" means a PSP, "No" means not a PSP). If label was not known, add "Yes" or "No" at will.  
For example, prepare file like this:  
MAICTEAISPFALIFLAIFLAFVALFEAFEHKCNCBKHFHALFKHLAK  Yes  
MAICTEAISPFALIFLAIFLAFVALF  Yes

# how to prepare other file
In hyper.py, change "path" to the actual path in your working place, change "test" to where the sequence file is saved.

# predict
python test_model.py

# training samples
pos_seq is the positive training samples for PSPredict training, neg_seq is the negative training samples for PSPredict training
