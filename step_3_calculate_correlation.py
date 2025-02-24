"""
Since we have some strains doing one conjugation, and some are doing two.
It means that we have a "continuous" phenotype. We can calculate Pearson
correlation.

It is also possible that the "continuous" phenotype is just an artifact. Caused
by for example measurement limitation. Then we can treat the phenotype binary,
use point-biserial correlation.
"""

import gzip
import pickle

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import pearsonr, pointbiserialr
from tqdm import tqdm

from project_settings import (GATHER_MATCH_TSV, PRESENCE_TSV, TARGET_STRAIN,
                              STRAINS_PICKLE_FILE)


def load_experimental_data():
    # Load experimental data:
    with open(STRAINS_PICKLE_FILE, "rb") as handle:
        experimental = pickle.load(handle)

    all_strains = {}
    for phenotype in experimental:
        all_strains.update(experimental[phenotype])

    with gzip.open(all_strains[TARGET_STRAIN], "rt") as handle:
        all_ref_prots = [s.id for s in SeqIO.parse(handle, "fasta")]

    return all_ref_prots, experimental, all_strains


def gen_presense_absence_table(
    ref_prots, strains, match_table_p=GATHER_MATCH_TSV
) -> pd.DataFrame:
    print(f"Reading match data {GATHER_MATCH_TSV}")
    match_df = pd.read_csv(GATHER_MATCH_TSV, sep="\t", index_col=0, header=0)
    presence_df = pd.DataFrame(
        np.zeros((len(ref_prots), len(strains)), dtype=int),
        index=pd.Index(sorted(ref_prots), name="gene"),
        columns=sorted(strains),
    )
    presence_dict = {}
    for _, row in tqdm(
        match_df.iterrows(), desc="Reading table", total=match_df.shape[0]
    ):
        gene = row["Query"]
        strain = row["Target strain"]
        if gene not in presence_dict:
            presence_dict[gene] = set()
        presence_dict[gene].add(strain)

    def fill_one_gene(target_row):
        gene = target_row.name
        if gene in presence_dict:
            for strain in strains:
                if strain in presence_dict[gene]:
                    target_row[strain] = 1
        return target_row

    print("Making presence/absence table.")
    presence_df = presence_df.apply(fill_one_gene, axis=1)
    return presence_df


def cal_pointbiserialr(
    phenotypes,  # list of keys in `experimental`
    experimental,  # dict of pheno:{strain:path}
    presence_df,
) -> pd.DataFrame:
    """
    Binary phenotypes
    """
    phenotype_strains = set()
    for phenotype in experimental:
        if phenotype in phenotypes:
            for strain in experimental[phenotype]:
                phenotype_strains.add(strain)
    phenotype_vector = np.zeros((presence_df.shape[1]), dtype=int)
    for i, strain in enumerate(presence_df.columns):
        if strain in phenotype_strains:
            phenotype_vector[i] = 1

    correlation_dict = {}
    for gene, presence_data in presence_df.iterrows():
        unique_presence = set(presence_data)
        if len(unique_presence) == 1:
            continue
        stat = pointbiserialr(presence_data, phenotype_vector)
        correlation_dict[gene] = [stat.correlation, stat.pvalue]

    correlation_df = pd.DataFrame(correlation_dict).T
    correlation_df.columns = ["point_biserial_Corr.", "p"]
    correlation_df.index.name = "gene"

    return correlation_df


def cal_step_pearsonr(
    phenotype_score_setup: dict[str:int],  # {pheno:score} Stepped scoring
    experimental,  # dict of pheno:{strain:path}
    presence_df,
) -> pd.DataFrame:
    """
    Continuous phenotypes
    """
    phenotype_score_dict = {}
    for phenotype in experimental:
        if phenotype in phenotype_score_setup:
            for strain in experimental[phenotype]:
                phenotype_score_dict[strain] = phenotype_score_setup[phenotype]

    phenotype_vector = np.zeros((presence_df.shape[1]), dtype=int)
    for i, strain in enumerate(presence_df.columns):
        if strain in phenotype_score_dict:
            phenotype_vector[i] = phenotype_score_dict[strain]

    correlation_dict = {}
    for gene, presence_data in presence_df.iterrows():
        unique_presence = set(presence_data)
        if len(unique_presence) == 1:
            continue
        stat = pearsonr(presence_data, phenotype_vector)
        correlation_dict[gene] = [stat.correlation, stat.pvalue]

    correlation_df = pd.DataFrame(correlation_dict).T
    correlation_df.columns = ["point_biserial_Corr.", "p"]
    correlation_df.index.name = "gene"

    return correlation_df


if __name__ == "__main__":
    all_ref_prots, experimental, all_strains = load_experimental_data()
    if not PRESENCE_TSV.exists():
        presence_df = gen_presense_absence_table(
            all_ref_prots, list(all_strains.keys())
        )
        print(f"Writing presence table {PRESENCE_TSV}.")
        presence_df.to_csv(PRESENCE_TSV, sep="\t")
    else:
        print(f"Found presence table, read from {PRESENCE_TSV}")
        presence_df = pd.read_csv(PRESENCE_TSV, sep="\t", index_col=0)

    pointbiseral_corr_df_conj = cal_pointbiserialr(
        ["Double conj.", "Single conj."], experimental, presence_df
    )

    pointbiseral_corr_df_double = cal_pointbiserialr(
        ["Double conj."], experimental, presence_df
    )

    pointbiseral_corr_df_single = cal_pointbiserialr(
        ["Single conj."], experimental, presence_df
    )

    pearson_corr_df = cal_step_pearsonr(
        {"Double conj.": 2, "Single conj.": 1}, experimental, presence_df
    )

    for name, df in zip(
        [
            "corr_pointbiseral_conj.tsv",
            "corr_pointbiseral_double.tsv",
            "corr_pointbiseral_single.tsv",
            "corr_pearson.tsv",
        ],
        [
            pointbiseral_corr_df_conj,
            pointbiseral_corr_df_double,
            pointbiseral_corr_df_single,
            pearson_corr_df,
        ],
    ):
        df.to_csv(name, sep="\t")
