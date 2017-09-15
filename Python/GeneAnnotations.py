##########################################################################################
# GeneAnnotations.py
# Purpose: script queries OMOP de-id and returns a vector of values for each patient
# version 1.1.0
# date: 08.15.2017
##########################################################################################

# create gene list needed for querying kaBOB
##########################################
######## GENERATING DATA FOR RENO ########
##########################################

# DATE: 03/02/2017

#import needed libraries
import mygene
import csv
import mygene


## STEP 1: PARSE GENE LIST ##

#parse the text file of genes
flist = [x.split('\n')[0].split(' /// ') if ' /// ' in x.split('\n')[0] else x.split('\n')[0] for x in open(
    'R/ignorenet/gene_union.txt').readlines()]
len(flist) #total number of lines in original file (479)

file = open('Python/KaBOB/Input_values.txt', 'w')
for line in set(flist):
    if ' /// ' in line:
        print(line)
    gene_value = 'iaohgnc:HGNC_' + str(line) + '_ICE'
    # print(line, gene_value)
    file.write(str(gene_value) + '\n')

# NOTE: of the 479 symbols we were provided, 49 are not in KaBOB


## STEP 2: CONVERT KABOB RESULTS TO TABLE OF COUNTS - (08/30/17: DID NOT USE THIS)  ##
# Parse results - genes and Reactome Pathways
file = open("results_genes_pathways.csv").read().split('\n')

res_dict = dict()

for line in file[1:]:
    key = str(', '.join(line.split(',')[3:]).replace('"', '').strip())
    key_ICE = str(line.split(',')[2].replace('"', '').split('/')[-1].split('_')[1])
    if (key, key_ICE) in res_dict:
        res_dict[(key, key_ICE)]['ice'].add(str(line.split(',')[0].replace('"', '').split('/')[-1].split('_')[1]))
        res_dict[(key, key_ICE)]['symbol'].add(str(line.split(',')[1].replace('"', '')))
    else:
        res_dict[(key, key_ICE)] = {}
        res_dict[(key, key_ICE)]['ice'] = {}
        res_dict[(key, key_ICE)]['symbol'] = {}
        res_dict[(key, key_ICE)]['ice'] = set([str(line.split(',')[0].replace('"', '').split('/')[-1].split('_')[1])])
        res_dict[(key, key_ICE)]['symbol'] = set([str(line.split(',')[1].replace('"', ''))])


#write out results
for key, value in res_dict.items():
    print key[1], key[0], str('; '.join([x.strip() for x in value['symbol']])), str('; '.join([x.strip() for x in value['ice']]))

f = open('gene_pathway_results.csv', 'w')
try:
    writer = csv.writer(f)
    writer.writerow(('Pathway ID', 'Pathway', 'Gene Symbol', 'Entrez ID'))
    for key, value in res_dict.items():
        writer.writerow((str(key[1]), str(key[0]), str('\015'.join([x.strip() for x in value['symbol']])), str('\015'.join([x.strip() for x in value['ice']]))))
finally:
    f.close()


# Parse results - genes and DrugBank drugs
file = open("results_genes_pathways.csv").read().split('\n')

res_dict = dict()

for line in file[1:]:
    key = str(', '.join(line.split(',')[3:]).replace('"', '').strip())
    key_ICE = str(line.split(',')[2].replace('"', '').split('/')[-1].split('_')[1])
    if (key, key_ICE) in res_dict:
        res_dict[(key, key_ICE)]['ice'].add(str(line.split(',')[0].replace('"', '').split('/')[-1].split('_')[1]))
        res_dict[(key, key_ICE)]['symbol'].add(str(line.split(',')[1].replace('"', '')))
    else:
        res_dict[(key, key_ICE)] = {}
        res_dict[(key, key_ICE)]['ice'] = {}
        res_dict[(key, key_ICE)]['symbol'] = {}
        res_dict[(key, key_ICE)]['ice'] = set([str(line.split(',')[0].replace('"', '').split('/')[-1].split('_')[1])])
        res_dict[(key, key_ICE)]['symbol'] = set([str(line.split(',')[1].replace('"', ''))])


#write out results
for key, value in res_dict.items():
    print key[1], key[0], str('; '.join([x.strip() for x in value['symbol']])), str('; '.join([x.strip() for x in value['ice']]))

f = open('gene_pathway_results.csv', 'w')
try:
    writer = csv.writer(f)
    writer.writerow(('Pathway ID', 'Pathway', 'Gene Symbol', 'Entrez ID'))
    for key, value in res_dict.items():
        writer.writerow((str(key[1]), str(key[0]), str('\015'.join([x.strip() for x in value['symbol']])), str('\015'.join([x.strip() for x in value['ice']]))))
finally:
    f.close()