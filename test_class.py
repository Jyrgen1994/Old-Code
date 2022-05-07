# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 22:54:31 2020

@author: etarmol

TESTING CLASS
"""
from model import Model
from sklearn.linear_model import SGDClassifier
def main():
    test = Model(':/',':/')
    # Markus C:/data_MVK/H7091-1.json
    #path_to_file = "/Users/aissata/mySkolfiler2/mvk/python/H7091-1.json"
    #test.read_json(path_to_file)
    
    #Build a pipeline (includes )
    text_clf = test.building_pipeline(SGDClassifier)

    text_clf.fit(test.training,test.train_target)

        
    docs_test = test.testing
    predicted = text_clf.predict(docs_test)
    test.test_classification(predicted)


if __name__ == '__main__':
    main()