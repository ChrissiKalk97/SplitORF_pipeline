#unique_region_locations.py takes as input a unique regions file from a
#run of the SplitORF pipeline and calculates the position of the unique regions 
#with respect to their ORF, the start in stop positions of the unique region
#with respect tot he ORf are plotted as density plots
#The type of unique regions needs to be specified as DNA or protein

#usage: python unique_region_locations.py Unique_{type}_Regions_merged_{for type protein: valid}.bed type


import sys
import pandas as pd
import seaborn as sbn
import matplotlib.pyplot as plt

def main():
    valid_ORFs = pd.read_csv(sys.argv[1], sep =":|\t")

    unique_regions = pd.read_csv(sys.argv[2], sep =":|\t", \
        names=["gene_transcript", "ORF_nr", "ORF_start", "ORF_end", "unique_start", "unique_end"],\
            engine='python')
    unique_regions["ORF_start"] = unique_regions["ORF_start"].astype(int)
    unique_regions = unique_regions.sort_values(['gene_transcript', 'ORF_start'])
    
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


    unique_regions['ORFcounter'] = unique_regions.groupby(['gene_transcript','ORF_start']).cumcount()+1
    #unique_regions['ORFcounter'] = unique_regions['ORFcounter'].astype('category')
    print(unique_regions.head())
    print(unique_regions.tail())

    transcript_type = sys.argv[4]

    sbn.histplot(data = unique_regions, x = 'stop_percent', hue = 'ORFcounter', multiple="stack",\
                  palette=["C0", "C1", "C2", "C3", "C4", "C5", "C6"])\
        .set(title=f'Histogram of end positions of unique {transcript_type}_{region_type} regions, by ORF number')
    plt.savefig(f"unique_regions_by_ORF_number/stop_positions_hist_by_ORF_number_{transcript_type}_{region_type}_stack.png")
    plt.close()

    sbn.histplot(x = unique_regions["start_percent"], hue = unique_regions['ORFcounter'], multiple="stack", \
                  palette=["C0", "C1", "C2", "C3", "C4", "C5", "C6"])\
        .set(title=f'Histogram of start positions of unique {transcript_type}_{region_type} regions, by ORF number')
    plt.savefig(f"unique_regions_by_ORF_number/start_positions_hist_by_ORF_number_{transcript_type}_{region_type}_stack.png")
    plt.close()

    sbn.histplot(data = unique_regions, x = 'stop_percent', hue = 'ORFcounter', multiple="dodge",\
                  palette=["C0", "C1", "C2", "C3", "C4", "C5", "C6"])\
        .set(title=f'Histogram of end positions of unique {transcript_type}_{region_type} regions, by ORF number')
    plt.savefig(f"unique_regions_by_ORF_number/stop_positions_hist_by_ORF_number_{transcript_type}_{region_type}_dodge.png")
    plt.close()

    sbn.histplot(x = unique_regions["start_percent"], hue = unique_regions['ORFcounter'], multiple="dodge", \
                  palette=["C0", "C1", "C2", "C3", "C4", "C5", "C6"])\
        .set(title=f'Histogram of start positions of unique {transcript_type}_{region_type} regions, by ORF number')
    plt.savefig(f"unique_regions_by_ORF_number/start_positions_hist_by_ORF_number_{transcript_type}_{region_type}_dodge.png")
    plt.close()


    sbn.histplot(unique_regions.groupby("gene_transcript").count()["ORF_nr"])\
        .set(title=f'Number of unique regions per transcript',\
              xlabel="Nr of unique regions per transcript")
    plt.savefig(f"unique_regions_by_ORF_number/nr_of_unique_{region_type}_regions_per_transcripts_{transcript_type}.png")
    plt.close()

    """fig, ((ax0, ax1)) = plt.subplots(nrows=1, ncols=2)

    colors = ['red', 'tan', 'lime']
    ax0.hist(x = unique_regions["start_percent"], 100, density=True, histtype='bar', color=colors, label=colors)
    ax0.legend(prop={'size': 10})
    ax0.set_title('bars with legend')"""


if __name__ == "__main__":
    main()
