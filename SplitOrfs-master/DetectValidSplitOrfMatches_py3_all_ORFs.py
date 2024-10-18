#!/usr/bin/python
#

import sys
import pandas as pd

minOrfNum=2  #minimum number ORFs that need to align to a peptide sequence
identityCutoff = 95  #minimum percent identity of the protein alignment to be considered as a valid ORF-peptide match
minLength = 50 #minimum number of amino acids for an ORF to be considered
minAlignmentRate = 0.5 # rate of positions to be covered by the alignment of an ORF (to remove spurious local protein alignments)
#colon=":"

if len(sys.argv) < 3:
    print("usage DetectValidSplitOrfMatches_py3.py BlastpAlign.out Outname.txt")
else :
    #read in blast output
    blast_output = pd.read_csv(sys.argv[1], sep = '\t', header = None, names = 
                               ['source', 'target', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart' ,'send' ,'evalue' ,'bitscore'])

    #only consider alignments with at least 95% identitiy of aligned positions
    blast_output = blast_output[blast_output['pident'] >= identityCutoff]
    
    ###########################################################################################################
    #calculate additional columns
    ###########################################################################################################
    #get additional columns like position alignment in target
    blast_output['protAlignPos'] = blast_output['sstart'].astype(str) + '-' + blast_output['send'].astype(str)
    #ORf positions
    blast_output['OrfPos'] = blast_output['source'].str.split(':').apply(lambda x : x[2]) \
        + '-' + blast_output['source'].str.split(':').apply(lambda x : x[3])

    #calculate alignment length of target positions
    blast_output['alignLength'] = blast_output['send'].astype(int) - blast_output['sstart'].astype(int) + 1
    #calculate the length of the ORF in nucleotides
    blast_output['orfLenNuc'] = blast_output['source'].str.split(':').apply(lambda x : int(x[3])) \
        - blast_output['source'].str.split(':').apply(lambda x : int(x[2])) + 1
    #calculate length of ORF in AA
    blast_output['orfLenAA'] = blast_output['orfLenNuc'] / 3
    

     
    #calculate length of blasted sequence in the source
    blast_output['blastOrfLength'] = blast_output['source'].str.split(':')\
        .apply(lambda x: int(x[4].split('-')[1]) - int(x[4].split('-')[0]))
    #blast_output = blast_output[blast_output['blastOrfLength'] >= 20]
    blast_output['AlignmentRate'] = blast_output['alignLength']/ blast_output['blastOrfLength']
    blast_output['OrfID'] = blast_output['source'].str.split(':').apply(lambda x: x[1])
    
    #get gene ids
    blast_output['gene_id_source'] = blast_output['source'].str.split('|').apply(lambda x: x[0])
    blast_output['gene_id_target'] = blast_output['target'].str.split('|').apply(lambda x: x[0])

    #get transcript ids
    blast_output['OrfTransID'] = blast_output['source'].str.split('|').apply(lambda x: x[1].split(':')[0])
    blast_output['targetTransID'] = blast_output['target'].str.split('|').apply(lambda x: x[1])


    ###########################################################################################################
    #FILTERING
    ###########################################################################################################
    blast_output = blast_output[blast_output['AlignmentRate'] >= minAlignmentRate]
    blast_output = blast_output[blast_output['gene_id_source'] == blast_output['gene_id_target']]


    target_source_grouped = blast_output.groupby(['targetTransID', 'OrfTransID'])['OrfID'].nunique()
    # Find the combinations where n_unique > 1
    at_least_2_ORFs = target_source_grouped >= minOrfNum

    # Combine values of targetTransID and OrfTransID with "|" wherever condition is True
    at_least_2_ORFS_frame = at_least_2_ORFs.reset_index()
    at_least_2_ORFS_frame = at_least_2_ORFS_frame[at_least_2_ORFS_frame['OrfID'] == True]

    print(len(blast_output['source']))
    # Use merge to filter based on common combinations of targetTransID and OrfTransID
    blast_output = blast_output.merge(at_least_2_ORFS_frame[['targetTransID', 'OrfTransID']], 
                                 on=['targetTransID', 'OrfTransID'], 
                                 how='inner')
    print(len(blast_output['source']))
    #only consider the longest match for each ORF
    longest_ORF_alignment = blast_output.groupby(['OrfID', 'targetTransID'])['alignLength'].max().reset_index()
    #print(longest_ORF_alignment)
    blast_output = blast_output.merge(longest_ORF_alignment, 
                                 on=['OrfID', 'targetTransID', 'alignLength'], 
                                 how='inner')
    
    #for the cases where there are several alignments of the same length and sometimes even pid
    sorting = blast_output.sort_values(['OrfTransID', 'OrfID','targetTransID', 'alignLength', 'pident'], 
                               ascending=[True, True, True, False, False])[['OrfTransID', 'OrfID', 'targetTransID', 'alignLength', 'pident']]
    #print(sorting[sorting.duplicated()])
    blast_output = blast_output.sort_values(['OrfTransID', 'OrfID', 'targetTransID', 'alignLength', 'pident'], 
                               ascending=[True, True, True, False, False]).drop_duplicates(['OrfID', 'OrfTransID','targetTransID'])
    #print(blast_output.sort_values(['OrfTransID', 'OrfID', 'targetTransID', 'alignLength', 'pident'], 
     #                          ascending=[True, True, True, False, False])[['OrfTransID','OrfID', 'targetTransID', 'alignLength', 'pident']])
    print(len(blast_output['source']))
    ###########################################################################################################
    #GET OUTPUT FRAME
    ###########################################################################################################
    Valid_SO_frame = blast_output.groupby(['targetTransID', 'OrfTransID'])['OrfID'].nunique().reset_index()
    Valid_SO_frame.rename(columns= {'OrfID': 'NumOrfs'}, inplace = True)
    gene_dict = blast_output.set_index('OrfTransID')['gene_id_source'].to_dict()
    Valid_SO_frame['geneID'] = Valid_SO_frame['OrfTransID'].map(gene_dict)


    blast_output['orfLenNuc'] = blast_output['orfLenNuc'].astype(str)
    blast_output['protAlignPos'] = blast_output['protAlignPos'].astype(str)
    grouped_df = blast_output.groupby(['targetTransID', 'OrfTransID'])
    Valid_SO_frame['OrfIDs'] = grouped_df['OrfID'].agg(lambda x: ','.join(x.dropna())).reset_index()['OrfID']
    Valid_SO_frame['OrfPos'] = grouped_df['OrfPos'].agg(lambda x: ','.join(x.dropna())).reset_index()['OrfPos']
    Valid_SO_frame['OrfLengths'] = grouped_df['orfLenNuc'].agg(lambda x: ','.join(x.dropna())).reset_index()['orfLenNuc']
    Valid_SO_frame['MinSeqIdent'] = grouped_df['pident'].min().reset_index()['pident']
    Valid_SO_frame['MaxSeqIdent'] = grouped_df['pident'].max().reset_index()['pident']
    Valid_SO_frame['protAlignPos'] = grouped_df['protAlignPos'].agg(lambda x: ','.join(x.dropna())).reset_index()['protAlignPos']
    Valid_SO_frame['protCoverage'] = grouped_df['alignLength'].sum().reset_index()['alignLength']

    blast_output['pident'] = blast_output['pident'].astype(str)
    grouped_df = blast_output.groupby(['targetTransID', 'OrfTransID'])
    Valid_SO_frame['OrfSeqIdents'] = grouped_df['pident'].agg(lambda x: ','.join(x.dropna())).reset_index()['pident']
    Valid_SO_frame = Valid_SO_frame[['geneID', 'targetTransID', 'OrfTransID', 'NumOrfs', 'OrfIDs', 'OrfPos',
                                     'OrfLengths', 'OrfSeqIdents', 'MinSeqIdent', 'MaxSeqIdent', 'protAlignPos', 'protCoverage']]
    print(Valid_SO_frame)


    ###########################################################################################################
    #WRITE OUTPUT TO FILE
    ###########################################################################################################
    Valid_SO_frame.sort_values(by = 'OrfTransID').to_csv(sys.argv[2], sep = '\t', index = False)
