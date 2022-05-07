# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 21:25:39 2020

@author: etarmol
following PEP8
"""

from sklearn import svm
from sklearn import feature_extraction
from sklearn.linear_model import SGDClassifier
from sklearn.pipeline import Pipeline
import numpy as np
import json
import pandas as pd
from nltk.corpus import stopwords
from sklearn.datasets import fetch_20newsgroups

class Model:

    #Constructor
    def __init__(self,path_to_traning_file,path_to_testing_file):
        self.traning_file = path_to_traning_file
        self.testing_file = path_to_testing_file
        self.emotions = ['väldigt positivt','positivt','neutral','negativt','väldigt negativt']
        self.training = []
        self.train_target = []
        self.testing = []
        self.test_target = []
    
    
    ##METHODS##
    def read_json(self,path_to_file):
        #INPUT: path to json file
        #OUTPUT: dictionary-like data from json file
        with open(path_to_file) as dataBE:
            data =json.load(dataBE)
        return data
        #print(json.dumps(dataBE, sort_keys = True,indent = 4))

    def building_pipeline(self,classifier):
        #INPUT: specific classifier in sklearn
        #OUTPUT: feature vector to be used by SVM
        
        
        stop_words_swedish = stopwords.words('swedish')

        self.fetch_train_test_data()

        text_clf = Pipeline([('vect' ,  feature_extraction.text.CountVectorizer(stop_words = stop_words_swedish)),
                             ('tfidf',  feature_extraction.text.TfidfTransformer()),
                             ('clf'  ,  classifier(loss='hinge', penalty='l2',
                                                       alpha=1e-3, random_state=42,
                                                       max_iter=5, tol=None))
        ])

        return text_clf
        
    def test_classification(self,predicted):
        print(predicted,len(predicted))

        nr_successful_classification = 0.0
        res_len = len(predicted)
        for i in range(res_len):
            if predicted[i] == self.test_target[i]:
                nr_successful_classification += 1
        

        prediction_success_ratio = nr_successful_classification/res_len
        print(prediction_success_ratio)


    def fetch_train_test_data(self):
        #INPUT: Nothing
        #OUTPUT: Assign training & testing data / target
        #TODO: read input file for larger sets of training/testing
        self.training = ['glada tankar ','jag hatar mitt liv','älskar delfiner','jag vet inte längre','sluta vara jobbig idiot']
        self.train_target = [self.emotions[1],self.emotions[4],self.emotions[0],self.emotions[2],self.emotions[3]]

        self.testing = ['jag är väldigt glad','jag älskar mitt liv','hata allt och alla','längre bort','fan va jobbig han är']
        self.test_target = [self.emotions[1],self.emotions[0],self.emotions[4],self.emotions[2],self.emotions[3]]
        
       
       





    

    
        #categories = ['alt.atheism', 'soc.religion.christian', 'comp.graphics', 'sci.med']
        #twenty_train = fetch_20newsgroups(subset='train',
        #                                  categories=categories, shuffle=True, random_state=42)
       # 
        #twenty_test = fetch_20newsgroups(subset='test',
        #                                  categories=categories, shuffle=True, random_state=42)
        
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
        #print(pd.DataFrame(tf_transform.toarray(), columns=sorted(count_vec.vocabulary_.keys())))
        
        #print(count_vec.vocabulary_.get(u'laers'))
        #print([w for w in sorted(vec.vocabulary_.keys())])
        #print(pd.DataFrame(transform.toarray(), columns=sorted(count_vec.vocabulary_.keys())))


        #counting occurances
        #count_vec = feature_extraction.text.CountVectorizer(stop_words = stop_words_swedish,
        #                                                    analyzer = 'word'
        #                                                    lowercase = False)
        #features = count_vec.fit_transform(texts)
        #features_arr = features.toarray()
        
        #Term Frequency = 'TF'
        #tf_transformer = feature_extraction.text.TfidfTransformer(use_idf=False).fit(transform)
        #use idf
        #tf_transform = tf_transformer.transform(transform)
        #tf_transform = feature_extraction.text.TfidfTransformer().fit_transform(features)