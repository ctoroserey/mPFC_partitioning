#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os
import sys

# get the file name from bash input
if os.environ.get('filename') == None:
    print "No outside input from bash"
    filename = raw_input("Enter file name (without .txt): ")
else:
    print "Input from bash"
    filename = os.environ['filename']
fname = filename + ".txt"

# make target directory
if not os.path.exists(filename):
    os.makedirs(filename)

# open files to read & write
textFile = open(fname,'rU')
textArray = textFile.readlines()
totalRows = len(textArray)

# change dir to store output
os.chdir(filename)

studyNumber = 1

for i in range(totalRows):
    if ('//' in textArray[i]) or (textArray[i] == '\n'):
        # if you want to include all of the coordinates into one file:
        # just use continue and open (coordFile, 'w') outside of the loop
        if ('Subjects=' not in textArray[i]) and ('Reference=' not in textArray[i]) and (textArray[i] != '\n'):
            splitting = textArray[i].split(' ')[1]
            fileName = 'coords_%s_%s.txt' % (studyNumber,splitting[:-1])
            studyNumber += 1
            coordFile = open(fileName,'w')
        else:
            continue
    else:

        print textArray[i]
        num = round(textArray[i])
        print num
        coordFile.write(textArray[i])

coordFile.close
textFile.close
