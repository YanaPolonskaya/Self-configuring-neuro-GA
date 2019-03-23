import h5py

from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Dropout
from keras.utils import np_utils
from multiprocessing import Pool 
import math
import numpy 
numpy.random.seed(7)
#dataset = numpy.loadtxt("/home/ipolonskaia/projects/Linux_project/data.txt", delimiter=",")
#X1 = dataset[:,0:8]
#Y1 = dataset[:,8]
Number_of_index=4
What_is_this_function=0;
#import sys
def Get_data(n_val,h5,N,validation_size): #n_val=5  
    i=0
    for id, participant in h5['/dataset/right'].items():  # There is a 'left' and a 'right'
        # Get labels, example: h5['/dataset/left/P20160211_14_m/labels'].value
        labels = participant['labels'].value  # NumPy NdArray
  
  
        # Get features (with a nice access to the elements;) )
        class Features:
            def __init__(self, participant):
                self.respiration = participant['atmung'].value  # Respiration
                self.ecg = participant['ekg'].value  # Heart rate
                self.emg = participant['emg'].value  # Trash
                self.gsr = participant['gsr'].value  # 
                self.gsr_phasic = participant['gsr_phasic'].value  # --"--
                self.gsr_tonic = participant['gsr_tonic'].value  # --"--
                self.geoFeaturesExtended = participant['geoFeaturesExtended'].value  # Video features
                # not added: ekgPQRSTWAVE, headpose, lbptop, phog
  
  
        features = Features(participant)

        # Get some information

        Y=labels[0]
        X = numpy.swapaxes(features.gsr_tonic, 0, 1)
        X = numpy.concatenate((X, numpy.swapaxes(features.gsr_phasic, 0, 1)), axis=1)
        X = numpy.concatenate((X, numpy.delete(numpy.swapaxes(features.gsr, 0, 1), [12,13,52,53,54], 1)), axis=1)

        samples = X.shape[1]  # Count of features
      
      
        indexes=[ind for ind in range(len(Y)) if Y[ind]==0 or Y[ind]==3 ]
        X= numpy.array([sample for j,sample in enumerate(X) if j in indexes]) 
        Y=numpy.array([int(sample/3) for j,sample in enumerate(Y) if j in indexes])

        y = np_utils.to_categorical(Y)

        if n_val==i:   
            test_X=X
            test_Y=y
        elif i==0:
            train_X=X
            train_Y=y
        elif n_val==0 and i==1:
            train_X=X
            train_Y=y
        else:
            train_X=numpy.concatenate((train_X, X), axis=0)
            train_Y=numpy.concatenate((train_Y, y), axis=0)
        i+=1
        if i==validation_size: 
            break
        
        print ("FEARURES" ,train_X.shape[1] )
    return train_X,train_Y,test_X,test_Y

def calcNN (structura): 
    print ("Print structura: ",structura)
    h5 = h5py.File ("/mnt/eshare/DatasetsRW/SENSEmotion/SENSEmotion_Patrick_Features.mat", 'r')
    N = len (h5['/dataset/left'])   # Count of participants
    acc=[]

    validation_size=0
    epohi=0

    if What_is_this_function==0:
        validation_size=5 
        epohi=6
    else:
        validation_size=N
        epohi=25

    length=len(structura)
    for i in range(validation_size):  
        train_X,train_Y,test_X,test_Y=Get_data(i,h5,N, validation_size)
        
        model = Sequential()
        
        for i in range(0,length,Number_of_index):
            if i==0:
                model.add(Dense(structura[i+1], input_dim=train_X.shape[1],activation=structura[i]))
            else:
                model.add(Dense(structura[i+1], activation=structura[i]))
            if structura[i+2]==1:
                model.add(Dropout(structura[i+3]))
        model.add(Dense(train_Y.shape[1], activation='sigmoid'))
        # Compile model
        model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
        # Fit the model
        model.fit(train_X, train_Y, epochs=epohi, batch_size=64)
        scores = model.evaluate(test_X, test_Y)
        acc.append(scores[1])
    del h5
    #print (sum(acc)/5)
    return sum(acc)/validation_size
    #return 1.0/(1.0+sum(acc)/validation_size)
    #print("Accuracy : %.2f%%" % (scores[1]*100))
    #return scores[1] 


def Convertation(strura):


    for i in range(0, len(strura)):
        for j in range(0, len(strura[i]),Number_of_index):
            strura[i][j] = int(strura[i][j])
            strura[i][j+1] = int(strura[i][j+1])
            strura[i][j+2] = int(strura[i][j+2])
            strura[i][j+3]=round(strura[i][j+3],1)
    #for i in range(0, len(strura)):
        #strura[i] = numpy.fromiter(strura[i], dtype = numpy.int)
    print("This returned in Pyth file YANA", strura)
    
    str=[]
    pravila=(0, "softmax", 1, "elu", 2, "selu", 3, "softplus", 4, "softsign", 5, "relu", 6, "tanh", 7, "sigmoid", 8, "hard_sigmoid", 9, "linear" )
    for i in range(0, len(strura)):
        str.append([])
        for j in range(0, len(strura[i]),4):
            if strura[i][j+1]!=0:
                for k in range (0, len(pravila),2):
                    if(strura[i][j]==pravila[k]):
                        str[i].append(pravila[k+1])
                        break
                str[i].append(strura[i][j+1])
                str[i].append(strura[i][j+2])
                str[i].append(strura[i][j+3])
                #str[i].append(strura[i][j])
    return str

def Model (strura):
    
    What_is_this_function=1
    print ("Model strura", strura)

   
    strura1=[]
    for i in range(0, len(strura),4):
        s=strura[i]
        strura1.append(int(s))
        s=strura[i+1]
        strura1.append(int(s))
        s=strura[i+2]
        strura1.append(int(s))   
        s=strura[i+3]
        strura1.append(round(s,1))


    print ("Model strura after ", strura1)


    structura=[]
    
    pravila=(0, "softmax", 1, "elu", 2, "selu", 3, "softplus", 4, "softsign", 5, "relu", 6, "tanh", 7, "sigmoid", 8, "hard_sigmoid", 9, "linear" )
    
    length=len(strura1)
    #for i in range(0, length):
        #if i%2 !=0 and strura[i]!=0:
    for i in range(0, length, Number_of_index):
        if(strura1[i+1]!=0):
            for j in range (0, len(pravila),2):
                    if(strura1[i]==pravila[j]):
                        structura.append(pravila[j+1])
                        for k in range(1, Number_of_index):
                            structura.append(strura1[i+k])
                        break
    print ("Model structura", structura)

    h5 = h5py.File ("/mnt/eshare/DatasetsRW/SENSEmotion/SENSEmotion_Patrick_Features.mat", 'r')
    N = len (h5['/dataset/left'])   # Count of participants
    acc=[]

    validation_size=N
    epohi=20

    length=len(structura)
    for i in range(validation_size):  
        train_X,train_Y,test_X,test_Y=Get_data(i,h5,N, validation_size)
        
        model = Sequential()
        
        for i in range(0,length,Number_of_index):
            if i==0:
                model.add(Dense(structura[i+1], input_dim=train_X.shape[1],activation=structura[i]))
            else:
                model.add(Dense(structura[i+1], activation=structura[i]))
            if structura[i+2]==1:
                model.add(Dropout(structura[i+3]))
        model.add(Dense(train_Y.shape[1], activation='sigmoid'))
        # Compile model
        model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
        # Fit the model
        model.fit(train_X, train_Y, epochs=epohi, batch_size=64)
        scores = model.evaluate(test_X, test_Y)
        acc.append(scores[1])

    del h5
    
    #standard deviation
    my=numpy.sum(acc)/len(acc)
    
    deviation=0
    for i in range(len(acc)):
        deviation+=(acc[i]-my)**2
    deviation=math.sqrt(deviation)/len(acc)

    print("ALL acc:", acc)
    print("DEVIATION:", deviation)
    return sum(acc)/validation_size


    #Fi=calcNN(str)

    #return Fi

def someFunction(strura):
    
    What_is_this_function=0
    print ("Length of list: ", len(strura))
    print ("Strctura in Python", strura)
   

    str=Convertation(strura)

    print("Str convertation: ",str)

    

    pool=Pool(processes=20) 
    Fi = pool.map(calcNN,str) 
    print("Fi ", Fi)
    
    #return scores[1]
    return Fi

if __name__ == '__main__':
    print("It's main function")
    h5 = h5py.File ("/mnt/eshare/DatasetsRW/SENSEmotion/SENSEmotion_Patrick_Features.mat", 'r')

    strura= [1 , 32 , 0 , 0.351562 , 9 , 36 , 0 , 0.34375 , 2 , 13 , 1 , 0.015625 , 3 , 29 , 0 , 0.859375 ]

    strura1=[]
    for i in range(0, len(strura),4):
        s=strura[i]
        strura1.append(int(s))
        s=strura[i+1]
        strura1.append(int(s))
        s=strura[i+2]
        strura1.append(int(s))   
        s=strura[i+3]
        strura1.append(round(s,1))


    print ("Model strura after ", strura1)


    structura=[]
    
    pravila=(0, "softmax", 1, "elu", 2, "selu", 3, "softplus", 4, "softsign", 5, "relu", 6, "tanh", 7, "sigmoid", 8, "hard_sigmoid", 9, "linear" )
    
    length=len(strura1)
    #for i in range(0, length):
        #if i%2 !=0 and strura[i]!=0:
    for i in range(0, length, Number_of_index):
        if(strura1[i+1]!=0):
            for j in range (0, len(pravila),2):
                    if(strura1[i]==pravila[j]):
                        structura.append(pravila[j+1])
                        for k in range(1, Number_of_index):
                            structura.append(strura1[i+k])
                        break
    print ("Model structura", structura)

    h5 = h5py.File ("/mnt/eshare/DatasetsRW/SENSEmotion/SENSEmotion_Patrick_Features.mat", 'r')
    N = len (h5['/dataset/left'])   # Count of participants
    acc=[]

    validation_size=N
    epohi=10

    length=len(structura)
    for i in range(validation_size):  
        train_X,train_Y,test_X,test_Y=Get_data(i,h5,N, validation_size)
        
        model = Sequential()
        
        for i in range(0,length,Number_of_index):
            if i==0:
                model.add(Dense(structura[i+1], input_dim=train_X.shape[1],activation=structura[i]))
            else:
                model.add(Dense(structura[i+1], activation=structura[i]))
            if structura[i+2]==1:
                model.add(Dropout(structura[i+3]))
        model.add(Dense(train_Y.shape[1], activation='sigmoid'))
        # Compile model
        model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
        # Fit the model
        model.fit(train_X, train_Y, epochs=epohi, batch_size=64)
        scores = model.evaluate(test_X, test_Y)
        acc.append(scores[1])

    del h5
    
    #standard deviation
    my=numpy.sum(acc)/len(acc)
    
    deviation=0
    for i in range(len(acc)):
        deviation+=(acc[i]-my)**2
    deviation=math.sqrt(deviation)/len(acc)

    print("ALL acc:", acc)
    print("DEVIATION:", deviation)
    
 




