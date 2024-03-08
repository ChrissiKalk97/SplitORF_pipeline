#filters the UniqueProteinORFPairs for the NMD targets (sources) that have wierd unique regions


import sys
import pandas as pd
from pybedtools import BedTool

def main():
    unique_regions = BedTool(sys.argv[1])

    weird_cases_beginning = pd.read_csv(sys.argv[2], header = 0, engine='python')
    weird_cases_end = pd.read_csv(sys.argv[3], header = 0, engine='python')

    weird_cases_beginning['ORF_filter'] = weird_cases_beginning['ORF_nr'] +':'+\
        weird_cases_beginning['unique_start'].astype(int).astype(str) +':'+weird_cases_beginning['unique_end'].astype(int).astype(str)
    weird_cases_end['ORF_filter'] = weird_cases_end['ORF_nr'] +':'+\
        weird_cases_end['unique_start'].astype(int).astype(str) +':'+weird_cases_end['unique_end'].astype(int).astype(str)
    ORFs_beginning = weird_cases_beginning['ORF_filter'].tolist()
    ORFs_end = weird_cases_end['ORF_filter'].tolist()
    ORFs_to_filter = ORFs_beginning + ORFs_end
    print('ORFs_to_filter', len(ORFs_to_filter))
    print(ORFs_to_filter[0:10])
    file_name = sys.argv[1].split('/')[-1]
    file_name = file_name[:-4]
    unique_bed_filtered = unique_regions.filter(lambda unique_region:\
        unique_region.name.split(':')[0]+':'+str(unique_region.start)+':'+str(unique_region.end) \
            not in set(ORFs_to_filter)).saveas(f'filtered_unique_regions/{file_name}_filtered.bed')
    

    
    
    
    



  



if __name__ == "__main__":
    main()
