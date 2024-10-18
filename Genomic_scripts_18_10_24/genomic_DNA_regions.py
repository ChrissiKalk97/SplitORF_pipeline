import sys
import pandas as pd

#helper function definitions
def unique_coordinate_in_exon(coordinate: int, regions: pd.DataFrame):
    for region in regions.itertuples():
        if region.start_trans <= coordinate and region.end_trans >= coordinate:
            return index, regions.loc[region[0], :]
        

def unique_DNA_regions_to_genomic(one_region, exon_trans_coords, f):
    #get corresponding exon coords
    exon_trans_coords_one_region = exon_trans_coords[exon_trans_coords['ID'] == one_region['ID']]
    one_region['ID'] = one_region['ID'] + one_region['ORF']
    #print(exon_trans_coords_one_region)

    start_unique = one_region.iloc[1]
    end_unique = one_region.iloc[2]

    #get start and end exon positions of the transcript corresponding to the unique region        
    start_index, start_exon = unique_coordinate_in_exon(start_unique, exon_trans_coords_one_region)
    end_index, end_exon = unique_coordinate_in_exon(end_unique, exon_trans_coords_one_region)

    strand = int(start_exon.iloc[5])
    chromosome = start_exon.iloc[6]

    start_exon_gen_start = start_exon.iloc[1]
    start_exon_gen_end = start_exon.iloc[2]
    start_exon_trans_start = start_exon.iloc[3]
    end_exon_gen_start = end_exon.iloc[1]

    if start_exon_gen_start != end_exon_gen_start:
        if strand == 1:
            genomic_start = start_exon_gen_start + start_unique - start_exon_trans_start
            genomic_end = start_exon_gen_end
            #note down the first exon as the stop exon comes later, at least this exon is present until the end
            f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + one_region['ID'] + '\n')
            current_exon_index = start_index + 1
            while current_exon_index <= end_index:
                current_exon = exon_trans_coords_one_region.loc[current_exon_index, :]
                genomic_start = current_exon.iloc[1]
                current_exon_trans_start = current_exon.iloc[3]
                current_exon_trans_end = current_exon.iloc[4]
                genomic_end = genomic_start + current_exon_trans_end - current_exon_trans_start
                assert(end_unique >= current_exon_trans_start)
                f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + one_region['ID'] + '\n')
                current_exon_index += 1
        elif strand == -1:
            genomic_start = start_exon_gen_start 
            genomic_end = start_exon_gen_end - (start_unique - start_exon_trans_start)
            #note down the first exon as the stop exon comes later, at least this exon is present until the end
            f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + one_region['ID'] + '\n')
            current_exon_index = start_index + 1
            while current_exon_index <= end_index:
                current_exon = exon_trans_coords_one_region.loc[current_exon_index, :]
                genomic_end = current_exon.iloc[2]
                current_exon_trans_start = current_exon.iloc[3]
                current_exon_trans_end = current_exon.iloc[4]
                genomic_start =  genomic_end - (current_exon_trans_end - current_exon_trans_start)
                assert(end_unique >= current_exon_trans_start)
                f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + one_region['ID'] + '\n')
                current_exon_index += 1

    else:
        if strand == 1:
            genomic_start = start_exon_gen_start + start_unique - start_exon_trans_start
            genomic_end = start_exon_gen_start  + end_unique - start_exon_trans_start
            f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + one_region['ID'] + '\n')
        elif strand == -1:
            genomic_start = start_exon_gen_end - (start_unique - start_exon_trans_start)
            genomic_end = start_exon_gen_end  - (end_unique - start_exon_trans_start)
            f.write(chromosome + '\t' + str(genomic_end) + '\t' + str(genomic_start) + '\t' + one_region['ID'] + '\n')
        

# read in unique regions and exon coordinates
dtypes_unique_regions = {
    'ID': 'string', 
    'start_trans': 'int64', 
    'end_trans': 'int64', 
    'ORF': 'string'        
}
dtypes_exon_coords = {
    'ID': 'string', 
    'start_trans': 'int64', 
    'end_trans': 'int64', 
    'start_gen' : 'int64',
    'end_gen' : 'int64',
    'strand' : 'string', 
    'chr' : 'string'
}


unique_DNA_regions_tr_coords = pd.read_csv(sys.argv[1],
                                            sep = '\t',
                                              header = None,
                                                names = ['ID', 'start_trans', 'end_trans', 'ORF'],
                                                dtype = dtypes_unique_regions)
print(unique_DNA_regions_tr_coords.head())

exon_trans_coords = pd.read_csv(sys.argv[2],
                                 sep = '\t',
                                   header = None,
                                     names = ['ID', 'start_gen', 'end_gen', 'start_trans', 'end_trans', 'strand', 'chr'],
                                     dtype = dtypes_exon_coords)
print(exon_trans_coords.head())

with open(sys.argv[3], 'w') as f:
    for index in unique_DNA_regions_tr_coords.index:
        unique_region = unique_DNA_regions_tr_coords.loc[index,:].copy()
        unique_DNA_regions_to_genomic(unique_region, exon_trans_coords, f)


