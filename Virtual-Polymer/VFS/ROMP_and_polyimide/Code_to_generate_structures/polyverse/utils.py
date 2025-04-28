from rdkit import Chem
import pandas as pd
import multiprocessing as mp
import dask.dataframe as dd


def canonical_mol(molecule):
    """
    Convert a molecule to its canonical form.

    Args:
        molecule (rdkit.Chem.Mol): Input molecule.

    Returns:
        rdkit.Chem.Mol: Canonicalized molecule.
    """
    return Chem.MolFromSmiles(Chem.MolToSmiles(molecule))


def canonical_smiles(smiles):
    """
    Convert a SMILES string to its canonical form.

    Args:
        smiles (str): Input SMILES string.

    Returns:
        str: Canonicalized SMILES string.
    """
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def daskify_df2df(function, input_df, npartitions=-1):
    """
    Parallelize a function's application to a DataFrame using Dask.

    Args:
        function (callable): The function to apply to each partition of the DataFrame.
        input_df (pd.DataFrame): Input DataFrame to apply the function to.
        npartitions (int, optional): Number of partitions for parallel processing.
            Default is -1, which uses 4 times the number of CPU cores.

    Returns:
        pd.DataFrame: Resulting DataFrame after applying the function.
    """
    if npartitions == -1:
        npartitions = 4 * mp.cpu_count()

    if npartitions > 1:
        # Parallelize the function using Dask.
        with ProgressBar(dt=1.0):
            dask_df = dd.from_pandas(input_df, npartitions=npartitions)
            result = (
                dask_df.map_partitions(
                    function,
                    meta=pd.DataFrame,
                )
                .compute(scheduler="processes")
                .values.tolist()
            )
            return pd.DataFrame(result)
    else:
        # No parallelization, apply the function directly.
        result = function(input_df)

        if isinstance(result, pd.DataFrame):
            return result
        elif isinstance(result, pd.Series):
            result = result.tolist()
            return pd.DataFrame(result)
        else:
            raise ValueError(f"Type {type(result)} is not supported.")
