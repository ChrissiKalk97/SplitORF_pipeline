import sys
import pandas as pd

transcript_coords_test = pd.read_csv(
    sys.argv[1], sep='\t', header=None, names=['ID', 'gen_start', 'gen_stop', 'transcript_start', 'transcript_stop', 'strand', 'chr'])

transcript_coords_test['transcript_length'] = transcript_coords_test['transcript_stop'] - \
    transcript_coords_test['transcript_start']
transcript_coords_test['genomic_length'] = transcript_coords_test['gen_stop'] - \
    transcript_coords_test['gen_start']

# first check: genomic and transcriptomic regions need to have the same length
assert (transcript_coords_test['transcript_length']
        == transcript_coords_test['genomic_length']).all(), 'transcript and genomic lengths differ'

print(transcript_coords_test.dtypes)


whole_transcript_length = transcript_coords_test.groupby(
    'ID')['transcript_length'].apply(lambda x: x.sum())


# second test: from ENsmebl archives v110, take the transcript lengths and compare
assert whole_transcript_length['ENSG00000163138|ENST00000508952'] == 987, 'transcript length of ENST00000508952 incorrect'
assert whole_transcript_length['ENSG00000126524|ENST00000697861'] == 1339, 'transcript length of ENST00000697861 incorrect'
assert whole_transcript_length['ENSG00000126524|ENST00000697862'] == 1532, 'transcript length of ENST00000697862 incorrect'

print(transcript_coords_test[transcript_coords_test['ID'] ==
                             'ENSG00000126524|ENST00000697862']['gen_start'].min())
# for the minus strand take one transcript and check coordiantes with Ensembl directly
assert transcript_coords_test[transcript_coords_test['ID'] ==
                              'ENSG00000126524|ENST00000697862']['gen_start'].min() == 66987700 - 1
assert transcript_coords_test[transcript_coords_test['ID'] ==
                              'ENSG00000126524|ENST00000697862']['gen_stop'].max() == 66995604
