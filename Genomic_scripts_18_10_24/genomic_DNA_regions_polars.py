import sys
import polars as pl

def unique_coordinate_in_exon(coordinate: int, regions: pl.DataFrame):
    filtered = regions.filter(
        (regions["start_trans"] <= coordinate) & (regions["end_trans"] >= coordinate)
    )
    if not filtered.is_empty():
        return filtered.head(1) 
    
def unique_trans_coords_to_genomic(exon_trans_coords_one_region, one_region, f):
    # Step 2: Update the 'ID' column in one_region
    ID_with_ORF = one_region[0] + one_region[3]
    start_unique = int(one_region[1])
    end_unique = int(one_region[2])

    #get start and end exon positions of the transcript corresponding to the unique region        
    start_exon = unique_coordinate_in_exon(start_unique, exon_trans_coords_one_region)
    end_exon = unique_coordinate_in_exon(end_unique, exon_trans_coords_one_region)
    start_exon = start_exon.row(0)
    end_exon = end_exon.row(0)
    start_index = int(start_exon[0])
    end_index = int(end_exon[0])

    strand = int(start_exon[6])
    chromosome = start_exon[7]

    start_exon_gen_start = int(start_exon[2])
    start_exon_gen_end = int(start_exon[3])
    start_exon_trans_start = int(start_exon[4])
    end_exon_gen_start = int(end_exon[2])

    if start_exon_gen_start != end_exon_gen_start:
        if strand == 1:
            genomic_start = start_exon_gen_start + start_unique - start_exon_trans_start
            genomic_end = start_exon_gen_end
            #note down the first exon as the stop exon comes later, at least this exon is present until the end
            f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + ID_with_ORF + '\n')
            current_exon_index = start_index + 1
            while current_exon_index <= end_index:
                current_exon = exon_trans_coords_one_region.filter(exon_trans_coords_one_region['index'] == current_exon_index).row(0)
                genomic_start = current_exon[2]
                current_exon_trans_start = current_exon[4]
                current_exon_trans_end = current_exon[5]
                genomic_end = genomic_start + current_exon_trans_end - current_exon_trans_start
                assert(end_unique >= current_exon_trans_start)
                f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + ID_with_ORF + '\n')
                current_exon_index += 1
        elif strand == -1:
            genomic_start = start_exon_gen_start 
            genomic_end = start_exon_gen_end - (start_unique - start_exon_trans_start)
            #note down the first exon as the stop exon comes later, at least this exon is present until the end
            f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + ID_with_ORF + '\n')
            current_exon_index = start_index + 1
            while current_exon_index <= end_index:
                current_exon = exon_trans_coords_one_region.filter(exon_trans_coords_one_region['index'] == current_exon_index).row(0)
                genomic_end = current_exon[3]
                current_exon_trans_start = current_exon[4]
                current_exon_trans_end = current_exon[5]
                genomic_start =  genomic_end - (current_exon_trans_end - current_exon_trans_start)
                assert(end_unique >= current_exon_trans_start)
                f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + ID_with_ORF + '\n')
                current_exon_index += 1

    else:
        if strand == 1:
            genomic_start = start_exon_gen_start + start_unique - start_exon_trans_start
            genomic_end = start_exon_gen_start  + end_unique - start_exon_trans_start
            f.write(chromosome + '\t' + str(genomic_start) + '\t' + str(genomic_end) + '\t' + ID_with_ORF + '\n')
        elif strand == -1:
            genomic_start = start_exon_gen_end - (start_unique - start_exon_trans_start)
            genomic_end = start_exon_gen_end  - (end_unique - start_exon_trans_start)
            f.write(chromosome + '\t' + str(genomic_end) + '\t' + str(genomic_start) + '\t' + ID_with_ORF + '\n')



    

###########################################################################################
# LOAD DATA    
###########################################################################################


# Define the data types for the columns
dtypes_unique_regions = {
    'ID': pl.Utf8,                # Use Utf8 for string columns
    'start_trans': pl.Int32,       # Use Int32 for integer columns
    'end_trans': pl.Int32,         # Use Int32 for integer columns
    'ORF': pl.Utf8              # Use Float32 for floating point numbers
}

dtypes_exon_coords = {
    'ID': pl.Utf8, 
    'start_trans': pl.Int32, 
    'end_trans': pl.Int32, 
    'start_gen' : pl.Int32,
    'end_gen' : pl.Int32,
    'strand' : pl.Utf8, 
    'chr' : pl.Utf8
}

# Read the CSV file
unique_DNA_regions_tr_coords = pl.read_csv(
    sys.argv[1],                   # File path from sys.argv[1]
    separator = '\t',                      # Use tab separator
    has_header=False,              # No header in the file (equivalent to header=None in pandas)
    new_columns=['ID', 'start_trans', 'end_trans', 'ORF'],  # Specify column names
    schema_overrides=dtypes_unique_regions   # Define data types for each column
)

# Display the DataFrame
print(unique_DNA_regions_tr_coords)

exon_trans_coords = pl.read_csv(
    sys.argv[2],
    separator = '\t',
    has_header = False,
    new_columns = ['ID', 'start_gen', 'end_gen', 'start_trans', 'end_trans', 'strand', 'chr'],
    schema_overrides = dtypes_exon_coords
)
exon_trans_coords = exon_trans_coords.with_row_index("index")
print(exon_trans_coords)
print(type(exon_trans_coords))


with open(sys.argv[3], 'w') as f:
    #unique_DNA_regions_tr_coords.with_columns(
    #pl.struct(['ID', 'start_trans', 'end_trans', 'ORF']).apply(lambda row: unique_trans_coords_to_genomic(exon_trans_coords,row, f))
    #)
    for row in unique_DNA_regions_tr_coords.iter_rows():
        # Convert the row (tuple) to a Polars Series
        unique_region = pl.Series(row, strict=False)

        # Filter the exon_trans_coords DataFrame based on the unique_region ID
        exon_trans_coords_one_region = exon_trans_coords.filter(
            pl.col("ID") == unique_region[0]  # Access the ID from the unique_region Series
        ).clone()
        
        # Call your function with the unique_region and other DataFrames
        unique_trans_coords_to_genomic(exon_trans_coords_one_region, unique_region, f)


        


