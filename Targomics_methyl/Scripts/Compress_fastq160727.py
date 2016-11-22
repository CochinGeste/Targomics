#!/usr/bin/env python 
import sys
import os
import re

#fichier_fastq=sys.argv[1]

fichier_fastq="/home/hyperion/Documents/Projet/Bisufilte_seq/Projet_batch1_21juillet2016/Input/IonXpress_014_H_R_2015_07_31_15_04_21_user_TOU-12-NGS15_001_GUIA_BISULFITE_Puce_3_.fastq"




name = os.path.splitext(fichier_fastq)[0]
ext = os.path.splitext(fichier_fastq)[1]

name_out=name+"_compress"+ext


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
				
longueur = 	file_len(fichier_fastq)	

fichier_fastq_out = open(name_out,"w")
file1_inp= open(fichier_fastq, "rU")
'''
c=0
for line in file1_inp.readlines() :
	c=c+1
	print c
	print line 
	if c==100:
			break
file1_inp.close()
'''

cmpt=0
tot = 0
for line in file1_inp.readlines() :
	cmpt = cmpt+1
	tot = tot +1
	print "Processing %d / %d"%(tot,longueur)
	# line= file1_inp.readlines(4)
	if (cmpt == 1) :
		fichier_fastq_out.write(line)
	
	if (cmpt == 2) : ## Si il y a une sequence ## 
	#seq = file1_inp.readline()
		seq2= line.strip()
		
		## Phase de protection ##
		seq2 =  re.sub(r'CA',"RR",seq2)
		seq2 =  re.sub(r'CG',"WW",seq2)
		seq2 =  re.sub(r'TG',"XX",seq2)
		
		## Detection des Homopolymere T ##
		p = re.compile("T+")
		pos_T_retires=list()
		for m in p.finditer(seq2):
			pos_T_retires.append( range(    m.start(),  (  m.start()+len(m.group()) )  -1) )
		pos_T_retires = [val for sublist in pos_T_retires for val in sublist]
		
		## Detection des Homopolymere A ##
		p = re.compile("A+")
		pos_A_retires=list()
		for m in p.finditer(seq2):
			pos_A_retires.append( range(    m.start(),  (  m.start()+len(m.group()) )  -1) )
		pos_A_retires = [val for sublist in pos_A_retires for val in sublist]
			
		## Detection des Homopolymere C ##
		p = re.compile("C+")
		pos_C_retires=list()
		for m in p.finditer(seq2):
			pos_C_retires.append( range(    m.start(),  (  m.start()+len(m.group()) )  -1) )
		pos_C_retires = [val for sublist in pos_C_retires for val in sublist]
			
		## Detection des Homopolymere G ##
		p = re.compile("G+")
		pos_G_retires=list()
		for m in p.finditer(seq2):
			pos_G_retires.append( range(    m.start(),  (  m.start()+len(m.group()) )  -1) )
		pos_G_retires = [val for sublist in pos_G_retires for val in sublist]
	
	
		## Concatenation des positions a retirer 
		all_retires =  pos_T_retires + pos_A_retires + pos_C_retires + pos_G_retires
		all_retires = sorted(all_retires, reverse=True)
		
		## On retire les positions ##
		seq2_listed = list(seq2)
		for indice in all_retires :
			seq2_listed.pop(indice)
	
		seq2=''.join(seq2_listed)
	
		fichier_fastq_out.write(seq2+"\n")
	
	
		## Phase de deprotection ##

		
		seq2 =  re.sub(r'RR',"CA",seq2)
		seq2 =  re.sub(r'WW',"CG",seq2)
		seq2 =  re.sub(r'XX',"TG",seq2)
	
	if (cmpt == 3) :
		fichier_fastq_out.write(line)		
		
		
	if (cmpt == 4) : ## Si il y a un score de qualite
		
		seq2_listed = list(line)
		for indice in all_retires :
			seq2_listed.pop(indice)
	
		seq2=''.join(seq2_listed)
	
	
		cmpt= 0
	
		fichier_fastq_out.write(seq2)
		
		



fichier_fastq_out.close()
file1_inp.close()


## Sandbox 


				
		