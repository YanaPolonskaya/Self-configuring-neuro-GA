#!/usr/bin/env python3

#import matplotlib.pyplot as plt
import h5py



# Naive LSTM to learn one-char to one-char mapping

import numpy as np
np.random.seed(7)
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM, Dropout,RepeatVector,SimpleRNN, GRU
from keras.utils import np_utils
import pandas as pn
from numpy import genfromtxt
from sklearn.model_selection import train_test_split

def Get_data(n_val,h5,N):
    
    my_data = genfromtxt('/home/ipolonskaia/Data.csv', delimiter=',')
    my_labels=genfromtxt('/home/ipolonskaia/Labels.csv', delimiter=',')
    my_data = genfromtxt('/home/ipolonskaia/Data.csv', delimiter=',')
    my_labels=genfromtxt('/home/ipolonskaia/Labels.csv', delimiter=',')
    train_X=np.concatenate((my_data, my_labels), axis=1)
    np.random.shuffle(train_X)
    train_1, test_1 = train_X[:198,:], train_X[198:,:]
    test_labels=test_1[:,[77,78]]
    train_labels=train_1[:,[77,78]]
    train=train_1[:,:77]
    test=test_1[:,:77]
    
    #train_X, test_X, train_Y, test_Y = train_test_split(my_data, my_labels, test_size=0.2) 
    #train_y=np.array(train_Y) 
    #test_y=np.array(test_Y)   
    #return train_X,train_Y,test_X,test_Y
    return train,train_labels,test,test_labels

def Model(train_X,train_Y,test_X,test_Y):
    # create and fit the model
    model = Sequential()
    #model.add(RepeatVector(3, input_shape=(train_X.shape[1],)))
    # GP constructed part

    {NN_model_structure}

    #End of GP constructed part

    model.add(Dense(train_Y.shape[1], activation='softmax'))
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    model.fit(train_X, train_Y, epochs=3, batch_size=64, verbose=0)
    #print(model.summary())
    #print("val_acc: %.2f%%" % hist.history['val_acc'][-1])
  
    # summarize performance of the model
    scores = model.evaluate(test_X, test_Y, verbose=0)
    #print("Model Accuracy: %.2f%%" % (scores[1] * 100))
    return scores[1]
  
    # demonstrate some model predictions
    '''print(X[0])
    prediction = model.predict(X, verbose=0)
    print ( [np.argmax(i) for i in prediction] )
    print("and so on ...")'''

def MyRNN(text): 
    h5 = h5py.File ("/mnt/eshare/DatasetsRW/SENSEmotion/SENSEmotion_Patrick_Features.mat", 'r')

    N = len (h5['/dataset/right'])   # Count of participants
    acc=[]

    #X = [0,0,1,1,1,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,1,1,1,1,0,1,1,1,0,1,0,0,0,0,1,1,1,1,1,0,0,1,1,1,1,0,1,1,1,0,0,1,1,0,1,0,1,1,1,1,1,1,0,1,1,1,0,1,1,1,0,1]
    #X = [0,1,0,0,0,1,1,0,0,0,0,1,1,1,1,0,1,0,1,1,0,0,1,1,0,1,0,1,0,1,0,0,0,1,1,1,0,0,0,1,1,0,0,0,1,1,1,1,0,0,0,0,1,0,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,1,1,0]
    for i in range(5):
      train_X,train_Y,test_X,test_Y=Get_data(i,h5,N)

      '''train_X = pn.DataFrame(train_X) 
      test_X = pn.DataFrame(test_X)

      for itr,j in enumerate(X): 
          if (j==0): 
              del train_X[itr] 
              del test_X[itr] '''

      acc.append(Model(train_X,train_Y,test_X,test_Y))
    del h5
    #print (sum(acc)/5)
    return sum(acc)/len(acc) #1.0/(1.0+sum(acc)/5.0)

if __name__ == '__main__':
    print("It's main function")
    
    h5 = h5py.File ("/mnt/eshare/DatasetsRW/SENSEmotion/SENSEmotion_Patrick_Features.mat", 'r')

    N = len (h5['/dataset/right'])   # Count of participants
    acc=[]

    #X = [0,0,1,1,1,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,1,1,1,1,0,1,1,1,0,1,0,0,0,0,1,1,1,1,1,0,0,1,1,1,1,0,1,1,1,0,0,1,1,0,1,0,1,1,1,1,1,1,0,1,1,1,0,1,1,1,0,1]
    #X = [0,1,0,0,0,1,1,0,0,0,0,1,1,1,1,0,1,0,1,1,0,0,1,1,0,1,0,1,0,1,0,0,0,1,1,1,0,0,0,1,1,0,0,0,1,1,1,1,0,0,0,0,1,0,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,1,1,0]
    for i in range(5):
      train_X,train_Y,test_X,test_Y=Get_data(i,h5,N)

      '''train_X = pn.DataFrame(train_X) 
      test_X = pn.DataFrame(test_X)

      for itr,j in enumerate(X): 
          if (j==0): 
              del train_X[itr] 
              del test_X[itr] '''

      acc.append(Model(train_X,train_Y,test_X,test_Y))
    del h5
    #print (sum(acc)/5)
    print("Acc:")
    print( sum(acc)/len(acc)) #1.0/(1.0+sum(acc)/5.0)
  
