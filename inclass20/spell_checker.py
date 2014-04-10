# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 16:55:20 2014

@author: pruvolo
"""

import re
import time
from matplotlib import pyplot
from bisect import bisect_left

def bisect_in(x, t):
    """ Uses bisect to check if x is in t.

        x: element to find
        t: sorted list

        Returns: boolean
    """
    i = bisect_left(t, x)
    return i != len(t) and t[i] == x

def load_words_as_list(skip_factor):
    """ Load a dictionary of words as a list.

        skip_factor: the step size to move through the textfile words.txt when
                     constructing the dictionary
        Returns: a list containing the words in the dictionary
    """
        
    f = open('words.txt','r')
    words = f.readlines()
    selected_words = []
    for i in range(0,len(words),skip_factor):
        selected_words.append(words[i].rstrip())
    f.close()
    return selected_words
    
def load_words_as_set(skip_factor):
    """ Load a dictionary of words as a set.

        skip_factor: the step size to move through the textfile words.txt when
                     constructing the dictionary
        Returns: a set containing the words in the dictionary
    """
    f = open('words.txt','r')
    words = f.readlines()
    selected_words = []
    for i in range(0,len(words),skip_factor):
        selected_words.append(words[i].rstrip())
    f.close()
    return set(selected_words)

def get_words(file_name):
    """ Grabs the words for spellchecking from the specified file.
    
        file_name: the file to spellcheck
        Returns: a list of all words in the input file
    """
    f = open(file_name,'r')
    s = f.read()
    f.close()
    return re.findall('\w+',s)

def spell_check_alg1(skip_factor):
    """ Algorithm 1 for spell checking
    
        skip_factor: the skip factor to use when constructing the dictionary
        Returns: the average time to spell check each word
    """
    all_words = load_words_as_list(skip_factor)
    words_to_spellcheck = get_words('const.txt')
    start_time = time.time()
    for i in range(len(words_to_spellcheck)):
        current_word = words_to_spellcheck[i].lower()
        spelled_correctly = False
        for w in all_words:
            if current_word == w:
                spelled_correctly = True
        if i % 100 == 0 and i > 0:
            current_time = time.time()
            print "time per word: %.10f" % ((current_time - start_time)/(i+1))
            return (current_time - start_time)/(i+1)

def spell_check_alg2(skip_factor):
    """ Algorithm 2 for spell checking
    
        skip_factor: the skip factor to use when constructing the dictionary
        Returns: the average time to spell check each word
    """
    all_words = load_words_as_set(skip_factor)
    words_to_spellcheck = get_words('const.txt')
    start_time = time.time()
    for i in range(len(words_to_spellcheck)):
        current_word = words_to_spellcheck[i].lower()
        spelled_correctly = current_word in all_words
        if i % 1000 == 0 and i > 0:
            current_time = time.time()
            print "time per word: %.10f" % ((current_time - start_time)/(i+1))
            return (current_time - start_time)/(i+1)

def spell_check_alg3(skip_factor):
    """ Algorithm 3 for spell checking
    
        skip_factor: the skip factor to use when constructing the dictionary
        Returns: the average time to spell check each word
    """
    all_words = load_words_as_list(skip_factor)
    words_to_spellcheck = get_words('const.txt')
    start_time = time.time()
    for i in range(len(words_to_spellcheck)):
        current_word = words_to_spellcheck[i].lower()
        spelled_correctly = bisect_in(current_word,all_words)
        if i % 1000 == 0 and i > 0:
            current_time = time.time()
            print "time per word: %.10f" % ((current_time - start_time)/(i+1))
            return (current_time - start_time)/(i+1)

if __name__ == '__main__':
    spell_check = spell_check_alg1
    make_plot = True
    if not(make_plot):
        spell_check(1)
    else:
        average_time = []
        total_words = len(load_words_as_list(1))
        for i in range(1,20):
            average_time.append(spell_check(i))
        pyplot.plot([float(total_words)/i for i in range(1,20)],average_time)
        pyplot.show()