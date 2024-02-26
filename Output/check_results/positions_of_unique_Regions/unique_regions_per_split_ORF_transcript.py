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
    valid_ORFs = pd.read_csv(sys.argv[1], sep =":|\t", header = None,\
            names = ["target", "protAlignPos_start", "protAlignPos_end", "NMD_transcript", "ORF_nr"])
    #UniqueProteinMatches.bed
    print(valid_ORFs.head())
    print(len(valid_ORFs["ORF_nr"]))

    unique_regions = pd.read_csv(sys.argv[2], sep =":|\t", \
        names=["gene_transcript", "ORF_nr", "ORF_start", "ORF_end", "unique_start", "unique_end"],\
            engine='python')
    print(unique_regions.head())
    print(len(unique_regions["ORF_nr"]))

    #print("shit", [r for r in unique_regions["ORF_nr"].unique() if r not in valid_ORFs["ORF_nr"].unique()])
    
    unique_regions = pd.merge(valid_ORFs, unique_regions, on="ORF_nr", how="outer")
    unique_regions = unique_regions[unique_regions['protAlignPos_start'].notnull()]
    #unique_regions["ORF_start"] = unique_regions["ORF_start"].astype(int)
    #unique_regions["protAlignPos_start"] = unique_regions["protAlignPos_start"].astype(int)
    unique_regions = unique_regions.sort_values(['NMD_transcript', 'protAlignPos_start'])
    print(unique_regions.head())
    print(unique_regions.tail())
    
    
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


    unique_regions['ORFcounter'] = unique_regions.groupby(['NMD_transcript']).cumcount()+1
    print(unique_regions.head())
    print(unique_regions.tail())
    unique_regions['ORFcounter'] = unique_regions['ORFcounter'].astype('category')
    

    transcript_type = sys.argv[4]


    palette = sbn.color_palette("tab10") + sbn.color_palette("viridis", n_colors = 15)
    sbn.histplot(data = unique_regions, x = 'stop_percent', hue = 'ORFcounter', multiple="stack",\
                  palette=palette)\
        .set(title=f'Histogram of end positions of unique {transcript_type}_{region_type} regions, by ORF number')
    plt.savefig(f"unique_regions_by_ORF_number/stop_positions_hist_by_ORF_number_{transcript_type}_{region_type}_stack.png")
    plt.close()

    sbn.histplot(x = unique_regions["start_percent"], hue = unique_regions['ORFcounter'], multiple="stack", \
                  palette=palette)\
        .set(title=f'Histogram of start positions of unique {transcript_type}_{region_type} regions, by ORF number')
    plt.savefig(f"unique_regions_by_ORF_number/start_positions_hist_by_ORF_number_{transcript_type}_{region_type}_stack.png")
    plt.close()

    sbn.histplot(data = unique_regions, x = 'stop_percent', hue = 'ORFcounter', multiple="dodge", binwidth=0.8, \
                  palette=palette)\
        .set(title=f'Histogram of end positions of unique {transcript_type}_{region_type} regions, by ORF number')
    plt.savefig(f"unique_regions_by_ORF_number/stop_positions_hist_by_ORF_number_{transcript_type}_{region_type}_dodge.png")
    plt.close()

    sbn.histplot(x = unique_regions["start_percent"], hue = unique_regions['ORFcounter'], multiple="dodge", binwidth=0.8, palette = palette)\
        .set(title=f'Histogram of start positions of unique {transcript_type}_{region_type} regions, by ORF number')
    
    plt.savefig(f"unique_regions_by_ORF_number/start_positions_hist_by_ORF_number_{transcript_type}_{region_type}_dodge.png")
    plt.close()


    sbn.histplot(unique_regions.groupby("NMD_transcript").count()["ORF_nr"], binwidth=0.8)\
        .set(title=f'Number of unique regions per transcript',\
              xlabel="Nr of unique regions per transcript")
    plt.savefig(f"unique_regions_by_ORF_number/nr_of_unique_{region_type}_regions_per_transcripts_{transcript_type}.png")
    plt.close()



if __name__ == "__main__":
    main()
