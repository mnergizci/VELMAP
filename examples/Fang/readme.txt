#!/bin/bash

#To generate the insardir_x lines in the .conf file, run the following bash command in the same directory as where you want to call the velmap.m script.
#This command ensure the folders are ordered according to ascending and descending tracks, numbered as required by readparfile.m and each path is finished with a slash

for i in insar/*A* insar/*D* ; do echo $i >> frames.txt ; done
cat frames.txt | awk '{print "insardir_"NR":   "$0"/"}'
rm frames.txt
