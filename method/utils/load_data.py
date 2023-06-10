import os
import pandas as pd

file_location = os.path.dirname(__file__)
DATA_DIR = os.path.join(file_location, "../../data")


def load_METABRIC_survival_data() -> pd.DataFrame:
    assert f"METABRIC-expressions.tsv" in os.listdir(
        DATA_DIR
    ), f"METABRIC data is missing in the data folder."

    expressions = pd.read_csv(
        f"{DATA_DIR}/METABRIC-expressions.tsv",
        sep="\t",
        index_col=[0],
    ).T
    metadata = pd.read_csv(
        f"{DATA_DIR}/METABRIC-metadata.tsv", sep="\t", index_col=[0]
    )

    return pd.concat((expressions, metadata), axis=1)


def load_TCGA_survival_data(tcga_project: str = "CESC") -> pd.DataFrame:
    assert f"TCGA-{tcga_project}-expressions.tsv" in os.listdir(
        DATA_DIR
    ), f"TCGA project '{tcga_project}' data is missing in the data folder."

    expressions = pd.read_csv(
        f"{DATA_DIR}/TCGA-{tcga_project}-expressions.tsv",
        sep="\t",
        index_col=[0],
    ).T
    metadata = pd.read_csv(
        f"{DATA_DIR}/TCGA-{tcga_project}-metadata.tsv", sep="\t", index_col=[0]
    )

    return pd.concat((expressions, metadata), axis=1)