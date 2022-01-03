#Author EWTDTHK
#The algorithm will use 4 codon "lookbook" before the current codon position to be predicted or generated
#Instead of a algorithm build by hand this alogrithm use decision tree fitting in an attempt to build the logic behind 
#the codonVaccine
#this program takes as input the https://raw.githubusercontent.com/berthubert/bnt162b2/master/side-by-side.csv
#output is CV1.pdf and CV2.pdf, CV3.pdf which correspond to the estimated decision tree build the the codonVaccine  
from sklearn.tree import DecisionTreeClassifier as Classifier
from collections import deque
from sklearn import tree
from sklearn.tree import plot_tree
import graphviz
import time
#PARAMETER for fitting the decision tree
lookback=4
maxtreedepth=16
testpercentage=0.0
#=======================================
xdata=[]
y1data=[]
y2data=[]
y3data=[]
mappingdict={}
mappingdict["A"]=0
mappingdict["C"]=1
mappingdict["G"]=2
mappingdict["T"]=3
#how many colon to lookback
lookbackx=deque(maxlen=lookback)
def readandtestdata(filename):
  for i in range(lookback): 
    #PUTTING [A,T,G] for the first few lookback
    lookbackx.append([2,4,1])
  
  f=open(filename,"r")
  #read and build the data
  cnt=0
  for line in f.readlines():
    if cnt>=1:
      field=line.split(",")
      x=[]
      for i in range(3): x.append(mappingdict[field[1][i]])
      #encode each colon to predict in y1, y2, y3
      y1=mappingdict[field[2][0]]
      y2=mappingdict[field[2][1]]
      y3=mappingdict[field[2][2]]
      y1data.append(y1)
      y2data.append(y2)
      y3data.append(y3)
      xx=[]
      lookbackx.append(x)
      for i in range(lookback):
         xx+=lookbackx[-i]
      xdata.append(xx)
    cnt+=1
  clf1=Classifier(max_depth=maxtreedepth)
  clf2=Classifier(max_depth=maxtreedepth)
  clf3=Classifier(max_depth=maxtreedepth)
  #cal portion of data for training
  trainends=int(testpercentage*len(xdata))
  #training
  if testpercentage!=0:
   clf1.fit(xdata[:-trainends],y1data[:-trainends])
   clf2.fit(xdata[:-trainends],y2data[:-trainends])
   clf3.fit(xdata[:-trainends],y3data[:-trainends])
  else:
   clf1.fit(xdata,y1data)
   clf2.fit(xdata,y2data)
   clf3.fit(xdata,y3data)

  #prediction
  y1hat=clf1.predict(xdata)
  y2hat=clf2.predict(xdata)
  y3hat=clf3.predict(xdata)
  matchc=0
  nucmatch=0
  #validate predict
  for i in range(len(y1data)):
    if (y1hat[i]==y1data[i] and y2hat[i]==y2data[i] and y3hat[i]==y3data[i]): matchc+=1
    if y1hat[i]==y1data[i]: nucmatch+=1
    if y2hat[i]==y2data[i]: nucmatch+=1
    if y3hat[i]==y3data[i]: nucmatch+=1
  #print the codon match percentage
  print("codon matches: ",matchc/len(y1data)*100,"%","nuc match",nucmatch/(len(y1data)*3)*100,"%")
  feature_names=[]
  for i in range(lookback): 
     feature_names.append("P"+str(i)+"_0")
     feature_names.append("P"+str(i)+"_1")
     feature_names.append("P"+str(i)+"_2")
  tree_data = tree.export_graphviz(clf1, out_file=None, 
                      feature_names=feature_names,  
                      class_names=["G","A","C","T"],  
                      filled=True, rounded=True,  
                      special_characters=True)
  graph=graphviz.Source(tree_data)
  graph.render("CV1")
  tree_data = tree.export_graphviz(clf2, out_file=None, 
                      feature_names=feature_names,  
                      class_names=["G","A","C","T"],  
                      filled=True, rounded=True,  
                      special_characters=True)
  graph=graphviz.Source(tree_data)
  graph.render("CV2")
  tree_data = tree.export_graphviz(clf3, out_file=None, 
                      feature_names=feature_names,  
                      class_names=["G","A","C","T"],  
                      filled=True, rounded=True,  
                      special_characters=True)
  graph=graphviz.Source(tree_data)
  graph.render("CV3")
readandtestdata("side-by-side.csv")
