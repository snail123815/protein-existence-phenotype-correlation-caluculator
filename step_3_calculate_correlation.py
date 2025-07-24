import gzip
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import pearsonr, pointbiserialr
from tqdm import tqdm

from load_configs import (
    GATHER_MATCH_TSV,
    PRESENCE_TSV,
    STRAINS_PICKLE_FILE,
    TARGET_STRAIN,
)


def load_experimental_data():
    # Load experimental data:
    with open(STRAINS_PICKLE_FILE, "rb") as handle:
        phenotype_strains, all_strains = pickle.load(handle)
    phenotype_strains: dict[str, dict[str, float]]
    all_strains: dict[str, Path]
    with gzip.open(all_strains[TARGET_STRAIN], "rt") as handle:
        all_ref_prots: list[str] = [s.id for s in SeqIO.parse(handle, "fasta")]

    return all_ref_prots, phenotype_strains, all_strains


def gen_presense_absence_table(
    ref_prots: list[str], strains: list[str], match_table_p=GATHER_MATCH_TSV
) -> pd.DataFrame:
    print(f"Reading match data {match_table_p}.")
    match_df = pd.read_csv(match_table_p, sep="\t", index_col=0, header=0)
    df_presence = pd.DataFrame(
        np.zeros((len(ref_prots), len(strains)), dtype=int),
        index=pd.Index(sorted(ref_prots), name="gene"),
        columns=sorted(strains),
    )
    presence_dict: dict[str, set[str]] = {}
    for _, row in tqdm(
        match_df.iterrows(), desc="Reading table", total=match_df.shape[0]
    ):
        gene = row["Query"]
        strain = row["Target strain"]
        if gene not in presence_dict:
            presence_dict[gene] = set()
        presence_dict[gene].add(strain)

    def fill_one_gene(target_row) -> pd.Series:
        gene = target_row.name
        if gene in presence_dict:
            for strain in strains:
                if strain in presence_dict[gene]:
                    target_row[strain] = 1
        return target_row

    print("Making presence/absence table.")
    df_presence = df_presence.apply(fill_one_gene, axis=1)
    return df_presence


def cal_correlation(
    phenotype: str,
    phenotype_strains: dict[str, dict[str, float]],
    all_strains: dict[str, Path],
    presence_df: pd.DataFrame,
) -> pd.DataFrame:
    phenotype_vector = (
        pd.Series(phenotype_strains[phenotype], index=sorted(list(all_strains.keys())))
        .fillna(0.0)
        .to_numpy(dtype=float)
    )

    correlation_dict = {}
    is_binary = False
    print(f"\nCalculating correlation for {phenotype} phenotype.")
    # Columns already sorted but just to be sure
    presence_df = presence_df[sorted(presence_df.columns)]

    for gene, presence_data in tqdm(
        presence_df.iterrows(), desc="Processing genes", total=presence_df.shape[0]
    ):
        unique_presence = set(presence_data)
        if len(unique_presence) == 1:
            # If all strains have the same presence/absence value, skip this gene
            continue
        else:
            if len(unique_presence) == 2 and unique_presence.issubset({0, 1}):
                is_binary = True
        # Calculate point biserial correlation
        if is_binary:
            # If the presence data is binary (0/1), use point biserial correlation
            statistic, pvalue = pointbiserialr(presence_data, phenotype_vector)
        else:
            # If the presence data is continuous, use Pearson correlation
            statistic, pvalue = pearsonr(presence_data, phenotype_vector)
        correlation_dict[gene] = [statistic, pvalue]

    correlation_df = pd.DataFrame(correlation_dict).T
    correlation_df.columns = [
        ("point_biserial_Corr." if is_binary else "pearson_Corr."),
        "p",
    ]
    correlation_df.index.name = "gene"

    return correlation_df


if __name__ == "__main__":
    all_ref_prots, phenotype_strains, all_strains = load_experimental_data()
    if not PRESENCE_TSV.exists():
        presence_df = gen_presense_absence_table(
            all_ref_prots, list(all_strains.keys())
        )
        print(f"Writing presence table {PRESENCE_TSV}.")
        presence_df.to_csv(PRESENCE_TSV, sep="\t")
    else:
        print(f"Found presence table, read from {PRESENCE_TSV}")
        presence_df = pd.read_csv(PRESENCE_TSV, sep="\t", index_col=0)

    correlations = {}
    for phenotype in phenotype_strains:
        correlations[phenotype] = cal_correlation(
            phenotype, phenotype_strains, all_strains, presence_df
        )

    for phenotype, corr_df in correlations.items():
        corr_df.to_csv(GATHER_MATCH_TSV.parent / f"corr_{phenotype}.tsv", sep="\t")
