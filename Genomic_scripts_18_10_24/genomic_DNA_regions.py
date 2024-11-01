import sys
import pandas as pd

#helper function definitions
def unique_coordinate_in_exon(coordinate: int, regions: pd.DataFrame):
    for region in regions.itertuples():
        if region.start_trans <= coordinate and region.end_trans > coordinate:
            return region[0], regions.loc[region[0], :]
    #print('regions', regions, coordinate)
        

def unique_DNA_regions_to_genomic(one_region, exon_trans_coords, f):
    #get corresponding exon coords
    exon_trans_coords_one_region = exon_trans_coords[exon_trans_coords['ID'] == one_region['ID']]
    #print(exon_trans_coords_one_region)
    one_region['ID'] = one_region['ID'] + ':' + one_region['ORF']

    start_unique = one_region.iloc[1]
    end_unique = one_region.iloc[2]

    #get start and end exon positions of the transcript corresponding to the unique region        
    start_index, start_exon = unique_coordinate_in_exon(start_unique, exon_trans_coords_one_region)
    end_index, end_exon = unique_coordinate_in_exon(end_unique, exon_trans_coords_one_region)
    #print(start_index, start_exon)
    #print(end_index, end_exon)
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
            while current_exon_index < end_index:
                current_exon = exon_trans_coords_one_region.loc[current_exon_index, :]
                genomic_start = current_exon.iloc[1]
                genomic_end = current_exon.iloc[2]
                if not end_unique > current_exon.iloc[4] or not start_unique < current_exon.iloc[3]:
                    print(current_exon)
                    print(start_unique)
                    print(end_unique)
                    print(one_region)
                    print('current', current_exon_index)
                    print('end', end_index)
                f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + one_region['ID'] + '\n')
                current_exon_index += 1
            #if we are in the exon where the unique regions ends, we need to note
            # the region until the point where it ends
            current_exon = exon_trans_coords_one_region.loc[current_exon_index, :]
            current_exon_trans_start = current_exon.iloc[3]
            genomic_start =  current_exon.iloc[1]
            genomic_end = genomic_start + end_unique - current_exon_trans_start
            f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + one_region['ID'] + '\n')
            

        elif strand == -1:
            genomic_start = start_exon_gen_start 
            genomic_end = start_exon_gen_end - (start_unique - start_exon_trans_start)
            #note down the first exon as the stop exon comes later, at least this exon is present until the end
            f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + one_region['ID'] + '\n')
            current_exon_index = start_index + 1
            while current_exon_index < end_index:
                current_exon = exon_trans_coords_one_region.loc[current_exon_index, :]
                genomic_end = current_exon.iloc[2]
                genomic_start =  current_exon.iloc[1]
                assert(end_unique > current_exon.iloc[3])
                f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + one_region['ID'] + '\n')
                current_exon_index += 1
            current_exon = exon_trans_coords_one_region.loc[current_exon_index, :]
            current_exon_trans_start = current_exon.iloc[3]
            genomic_end = current_exon.iloc[2]
            genomic_start =  genomic_end - (end_unique - current_exon_trans_start)
            f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + one_region['ID'] + '\n')


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
        #print(unique_region)
        unique_DNA_regions_to_genomic(unique_region, exon_trans_coords, f)


