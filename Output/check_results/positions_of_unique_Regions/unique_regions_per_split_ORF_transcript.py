#unique_region_locations.py takes as input a unique regions file from a
#run of the SplitORF pipeline and calculates the position of the unique regions 
#with respect to their ORF, the start in stop positions of the unique region
#with respect to the ORF are plotted as density plots
#The type of unique regions needs to be specified as DNA or protein

#usage: python unique_regions_per_split_ORF_transcript.py 
#directory_pipe_output/UniqueProteinMatches.bed
#directory_pipe_output/Unique_DNA_Regions_merged.bed region_type transcript_type

#with region type DNa or protein
#transcript type RNA or NMD


import sys
import pandas as pd
import seaborn as sbn
import matplotlib.pyplot as plt

def main():
    def count_ORFs(unique_regions: pd.DataFrame) -> pd.DataFrame:
        """helper function to count ORFs correctly, accounting for all ORFs found for the source
        and mapped to the respective target and for having several unique regions"""
        unique_regions.groupby(['NMD_transcript'])
        previous_transcript = unique_regions.iloc[0, 3]
        previous_alignment_start = unique_regions.iloc[0, 1]
        for i in range(1, len(unique_regions['NMD_transcript'])):
            current_transcript = unique_regions.iloc[i, 3]
            current_alignment_start = unique_regions.iloc[i, 1]
            if current_transcript == previous_transcript and current_alignment_start  == previous_alignment_start:
                unique_regions.iloc[i, 14] = unique_regions.iloc[i-1, 14]
            elif current_transcript == previous_transcript:
                unique_regions.iloc[i, 14] = unique_regions.iloc[i-1, 14] + 1
            previous_alignment_start = current_alignment_start
            previous_transcript = current_transcript
        return unique_regions


    valid_ORFs = pd.read_csv(sys.argv[1], sep =":|\t", header = None,\
            names = ["target", "protAlignPos_start", "protAlignPos_end", "NMD_transcript", "ORF_nr"], engine='python')
    #UniqueProteinMatches.bed

    unique_regions = pd.read_csv(sys.argv[2], sep =":|\t", \
        names=["gene_transcript", "ORF_nr", "ORF_start", "ORF_end", "unique_start", "unique_end"],\
            engine='python')

    #unique regions were filtered by all valid ORFs not the best protein match filtered ones
    #need to filter for the best matches
    unique_regions = pd.merge(valid_ORFs, unique_regions, on="ORF_nr", how="outer")
    unique_regions = unique_regions[unique_regions['protAlignPos_start'].notnull()]

    unique_regions = unique_regions.sort_values(['NMD_transcript', 'protAlignPos_start'])
    
    transcript_type = sys.argv[4]
    region_type = sys.argv[3]
    if region_type == "protein":
        unique_regions["ORF_length"] = (unique_regions["ORF_end"] - unique_regions["ORF_start"])/3#protein coordinates
    elif region_type == "DNA":
        unique_regions["ORF_length"] = (unique_regions["ORF_end"] - unique_regions["ORF_start"])#DNA coordinates
    else:
        raise ValueError("Please give an appropriate region type as the second argument.\
                         This can be either protein or DNA.")

    unique_regions["start_percent"] = (unique_regions["unique_start"])/unique_regions["ORF_length"]
    unique_regions["stop_percent"] = (unique_regions["unique_end"])/unique_regions["ORF_length"]
    unique_regions["unique_percent"] = unique_regions["stop_percent"] - unique_regions["start_percent"]

    #count number of ORFs per source transcript
    unique_regions['ORFcounter'] = 1
    unique_regions['nr_unique_regon'] = unique_regions.groupby(['NMD_transcript', 'protAlignPos_start']).cumcount()+1

    
    
    unique_regions = count_ORFs(unique_regions)
    print(unique_regions.head(50))
    print(unique_regions.tail(50))
    
    #indicate whether first ORF is True = 1, or if it is last ORF = -1 or intermediate = 0
    unique_regions['firstORF'] = 0
    unique_regions.loc[unique_regions['ORFcounter'] == 1, 'firstORF'] = 1
    indices = unique_regions.groupby(['NMD_transcript'])['ORFcounter'].transform(max) == unique_regions['ORFcounter']
    unique_regions.loc[indices, 'firstORF'] = -1
    
    unique_regions['ORFcounter'] = unique_regions['ORFcounter'].astype('category')
    unique_regions['firstORF'] = unique_regions['firstORF'].astype('category')



    folder = sys.argv[5]
    #get dataframe of cases to inspect more closely
    weird_cases_first_orf = unique_regions[(unique_regions['firstORF'] == 1) & (unique_regions["start_percent"]  < 0.01)]
    weird_cases_last_orf = unique_regions[(unique_regions['firstORF'] == -1) & (unique_regions["stop_percent"]  > 0.99)]
    weird_cases_first_orf.reset_index(inplace =True)
    weird_cases_last_orf.reset_index(inplace =True)
    weird_cases_first_orf.to_csv(f'unique_regions_by_ORF_number/{folder}/first_orf_unique_beginning_{region_type}_{transcript_type}.csv')
    weird_cases_last_orf.to_csv(f'unique_regions_by_ORF_number/{folder}/last_orf_unique_end_{region_type}_{transcript_type}.csv') 



    """ #plotting of all ORFs
    sbn.histplot(data = unique_regions, x = 'stop_percent', hue = 'ORFcounter', multiple="stack",\
                  palette=palette)\
        .set(title=f'Histogram of end positions of unique {transcript_type}_{region_type} regions, by ORF number')
    plt.savefig(f"unique_regions_by_ORF_number/{folder}/unique_regions_stop_positions_hist_by_ORF_number_{transcript_type}_{region_type}_stack.png")
    plt.close()

    sbn.histplot(x = unique_regions["start_percent"], hue = unique_regions['ORFcounter'], multiple="stack", \
                  palette=palette)\
        .set(title=f'Histogram of start positions of unique {transcript_type}_{region_type} regions, by ORF number')
    plt.savefig(f"unique_regions_by_ORF_number/{folder}/unique_regions_start_positions_hist_by_ORF_number_{transcript_type}_{region_type}_stack.png")
    plt.close()"""

   

    palette = sbn.color_palette("tab10") + sbn.color_palette("viridis", n_colors = 15)
    #plotting accroding to first and last ORF
    sbn.histplot(data = unique_regions, x = 'stop_percent', hue = 'firstORF', multiple="stack",\
                  palette=palette)\
        .set(title=f'Histogram of end positions of unique {transcript_type}_{region_type} regions, by first/last ORF')
    plt.savefig(f"unique_regions_by_ORF_number/{folder}/unique_regions_stop_positions_hist__{transcript_type}_{region_type}_first_last_ORF_stack.png")
    plt.close()

    sbn.histplot(x = unique_regions["start_percent"], hue = unique_regions['firstORF'], multiple="stack", \
                  palette=palette)\
        .set(title=f'Histogram of start positions of unique {transcript_type}_{region_type} regions, by first/last ORF')
    plt.savefig(f"unique_regions_by_ORF_number/{folder}/unique_regions_start_positions_hist_{transcript_type}_{region_type}_first_last_ORF_stack.png")
    plt.close()



    #plotting of general info unique regions per ORF
    sbn.histplot(unique_regions.groupby("NMD_transcript").count()["ORF_nr"], binwidth=0.8)\
        .set(title=f'Number of unique regions per transcript',\
              xlabel="Nr of unique regions per transcript")
    plt.savefig(f"unique_regions_by_ORF_number/{folder}/nr_of_unique_{region_type}_regions_per_transcripts_{transcript_type}.png")
    plt.close()

    sbn.histplot(data = unique_regions, x = 'unique_percent', hue = 'firstORF', multiple="stack",\
                  palette=palette)\
        .set(title=f'Histogram of  unique percentage of per ORF {transcript_type}_{region_type} regions, by first/last ORF')
    plt.savefig(f"unique_regions_by_ORF_number/{folder}/unique_region_percentage_{transcript_type}_{region_type}_first_last_ORF_stack.png")
    plt.close()



if __name__ == "__main__":
    main()
