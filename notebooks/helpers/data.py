import os
import pandas as pd

file_location = os.path.dirname(__file__)
DATA_DIR = os.path.join(file_location, "../../data")

def load_TCGA_survival_data(tcga_project: str = "CESC") -> pd.DataFrame:
    """
    Loads TCGA expression and meta data for survival.
    The last two columns of the DataFrame are "time" and "event"

    Arguments
    ---------
    tcga_project: str
        The abbriviation of the TCGA project to load the data.

    Return
    ------
    pd.DataFrame
        A DataFrame with rows as samples and columns as genes. The last two
        columns are survival time ('time') ane event occurance ('event')
    """

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
