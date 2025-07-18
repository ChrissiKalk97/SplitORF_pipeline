# Make a barplot from a csv file of gProfiler
# with preselected terms


import argparse
import seaborn as sbn
import pandas as pd
import matplotlib.pyplot as plt


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Barplot of selected GO terms")

    parser.add_argument(
        "go_terms",
        type=str,
        help="Path to the selected GO terms gprofiler CSV file"
    )

    parser.add_argument(
        "out_plot",
        type=str,
        help="Path to the output plot file (SVG)"
    )

    return parser.parse_args()


def main(go_terms, out_plot):
    go_terms_df = pd.read_csv(go_terms, sep=',', header=0)
    go_terms_df = go_terms_df.sort_values(
        by='negative_log10_of_adjusted_p_value', ascending=False)
    g = sbn.catplot(data=go_terms_df, x="negative_log10_of_adjusted_p_value",
                    y="term_name", kind="bar", palette="viridis", height=6, aspect=1.5)
    g.set_axis_labels("−log₁₀(adjusted p-value)", "GO Term")

    plt.savefig(f'{out_plot}', format="svg")


if __name__ == "__main__":
    args = parse_arguments()
    main(args.go_terms, args.out_plot)
