from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
import pandas as pd
from numpy import ceil
import multiprocessing as mp
import itertools as it
from more_itertools import sliced

import polyverse.constants as ks
from polyverse.utils import daskify_df2df, canonical_mol
from polyverse.constraints import (
    DefaultRingSizeConstraint,
    ChargeConstraint,
    BackbiteConstraint,
)

# Disable rdkit warnings.
RDLogger.DisableLog("rdApp.*")

get_monomers_max_chunk_size = 50000  # Hard-coded. Based on experience.
int8_max = 127


class Polymerize:
    def __init__(
        self,
        name: str,
        reaction_class: str,
        rxn_smarts: str,
        polymerization1_smarts: str,
        polymerization2_smarts: str,
        substruct_1_exception: str,
        substruct_2_exception: str,
        npartitions: int = -1,
    ) -> None:
        """
        Initialize the Polymerize object.

        Args:
            name (str): The name of the polymerization process.
            reaction_class (str): The class of the reaction (e.g., "addition" or "condensation").
            rxn_smarts (str): SMARTS pattern for the reaction.
            polymerization1_smarts (str): SMARTS pattern for the first termination step.
            polymerization2_smarts (str): SMARTS pattern for the second termination step.
            substruct_1_exception (str): SMARTS patterns to exclude from molecules of type A.
            substruct_2_exception (str): SMARTS patterns to exclude from molecules of type B.
            npartitions (int, optional): Number of partitions for parallel processing. Default is -1.
        """
        self.name = name
        # Create the RDKit reaction object and initialize it.
        self.rxn_smarts = rxn_smarts
        self.rxn = AllChem.ReactionFromSmarts(self.rxn_smarts)
        self.rxn.Initialize()

        self.reaction_class = reaction_class

        # These are used to terminate (i.e., add star atoms to) the product of the
        # reaction encoded by self.rxn_smarts.
        self.polymerization1_smarts = polymerization1_smarts
        self.polymerization2_smarts = polymerization2_smarts

        # Substructures of the reacting species
        self.template_1 = self.rxn.GetReactantTemplate(0)
        if self.reaction_class == "condensation":
            self.template_2 = self.rxn.GetReactantTemplate(1)
        elif self.reaction_class == "addition":
            self.template_2 = None
        if self.polymerization1_smarts:
            self.P1 = AllChem.ReactionFromSmarts(polymerization1_smarts)
        else:
            self.P1 = None
        if self.polymerization2_smarts:
            self.P2 = AllChem.ReactionFromSmarts(polymerization2_smarts)
        else:
            self.P2 = None

        # Exceptions to exclude from molecules type A
        if substruct_1_exception:
            self.substruct_1_exceptions = [
                Chem.MolFromSmarts(x) for x in substruct_1_exception
            ]
        else:
            self.substruct_1_exceptions = []

        # Exceptions to exclude from molecules type B
        if substruct_2_exception:
            self.substruct_2_exceptions = [
                Chem.MolFromSmarts(x) for x in substruct_2_exception
            ]
        else:
            self.substruct_2_exceptions = []

        # Get the number of partitions to use.
        self.npartitions = npartitions
        if self.npartitions == -1:
            self.npartitions = 4 * mp.cpu_count()

    def discard_molecules(self, molecule_df):
        """
        Discard molecules that do not match the reaction template.

        Args:
            molecule_df (pd.DataFrame): Contains two columns, "id-smiles" and "index".

        Returns:
            pd.DataFrame: Contains three columns, "id-smiles", "template_1matches", "template_2matches".
        """
        print(f"\nStarting with {len(molecule_df)} molecules.", flush=True)

        def get_matches(x):
            """
            Compute the number of times self.template_1 and self.template_2 appear in x["id-smiles"].

            Args:
                x (pd.Series): A row of molecule_df.
            """
            _, smiles = x["id-smiles"]
            x_dict = x.to_dict()
            mol = Chem.MolFromSmiles(smiles)
            x_dict.pop("index", None)  # clear memory
            x_dict["template_1matches"] = min(
                len(mol.GetSubstructMatches(self.template_1)), int8_max
            )
            if self.template_2:
                x_dict["template_2matches"] = min(
                    len(mol.GetSubstructMatches(self.template_2)), int8_max
                )
            else:
                x_dict["template_2matches"] = -1  # replacement for nan

            return x_dict

        apply_fn = lambda df: df.apply(get_matches, axis=1)
        molecule_df = daskify_df2df(apply_fn, molecule_df, self.npartitions)

        # Free up memory by forcing data types.
        molecule_df["template_1matches"] = molecule_df["template_1matches"].astype(
            "int8"
        )
        molecule_df["template_2matches"] = molecule_df["template_2matches"].astype(
            "int8"
        )

        # Discard molecules without matches.
        molecule_df = molecule_df.loc[
            (molecule_df.template_1matches > 0) | (molecule_df.template_2matches > 0)
        ]

        print(
            f"\nDiscard complete. {len(molecule_df)} molecules remaining.", flush=True
        )
        return molecule_df

    def generate_polymers(self, molecule_df):
        """
        Generate polymers by combining the molecules in molecule_df
        according to the rules specified in __init__.
        """
        raise NotImplementedError()


class Polyimide(Polymerize):
    def __init__(
        self,
        npartitions=-1,
        charge_constraint=lambda x: ChargeConstraint(
            "<", 0.0005, False, x[0], x[1]
        ).check(mol=x[2]),
        bb_constraint=BackbiteConstraint(),
    ):
        """
        Initialize a Polyimide object for generating polyimide polymers.

        Args:
            npartitions (int, optional): Number of partitions for parallel processing.
            charge_constraint (function, optional): Charge constraint function.
            bb_constraint (BackbiteConstraint, optional): Backbite constraint object.
        """
        super().__init__(
            name="polyimide",
            reaction_class="condensation",
            rxn_smarts="[C:1](=[O:2])[O:3][C:4](=[O:5]).[NH2:6][#6:7]>>[C:1](=[O:2])[N:6]([#6:7])[C:4](=[O:5])",
            polymerization1_smarts="[C:1](=[O:2])[O:3][C:4](=[O:5])>>[C:1](=[O:2])[N:6]([C:4](=[O:5]))[#0]",
            polymerization2_smarts="[NH2:1][#6:2]>>[#6:2][#0]",
            substruct_1_exception=None,
            substruct_2_exception=None,
            npartitions=npartitions,
        )
        self.charge_constraint = charge_constraint
        self.bb_constraint = bb_constraint

    def check_charge(self, row):
        """
        Check if a molecule meets the charge constraint.

        Args:
            row (pd.Series): A row from the input DataFrame.

        Returns:
            pd.Series: A modified row containing charge constraint results.
        """
        _id, smiles = row["id-smiles"]
        mol = Chem.MolFromSmiles(smiles)
        template_1matches = row["template_1matches"]
        template_2matches = row["template_2matches"]

        # If the molecule does not exist, we cannot check the charge constraint.
        if mol is None:
            return {
                "ID": _id,
                "smiles": smiles,
                "mol": mol,
                "template_1matches": template_1matches,
                "template_2matches": template_2matches,
            }

        # If a molecule has two template1 matches, make sure the matches
        # pass the charge constraint.
        if template_1matches == 2:
            substruct1, substruct2 = mol.GetSubstructMatches(self.template_1)
            pass_constraint = self.charge_constraint((substruct1, substruct2, mol))

        # If a molecule has two template2 matches, make sure the matches
        # pass the charge constraint.
        elif template_2matches == 2:
            substruct1, substruct2 = mol.GetSubstructMatches(self.template_2)
            pass_constraint = self.charge_constraint((substruct1, substruct2, mol))
        else:
            raise ValueError("Unexpected template matches.")

        # If the charge constraint is not passed, set the mol to None.
        if not pass_constraint:
            mol = None
        return {
            "ID": _id,
            "smiles": smiles,
            "mol": mol,
            "template_1matches": template_1matches,
            "template_2matches": template_2matches,
        }

    def link_monomers(self, monomer_data):
        """
        Link monomers to form polymers based on reaction rules.

        This method takes a monomer represented as a pair of a unique identifier (ID)
        and a chemical molecule (Mol) and attempts to create a polymer by applying
        the reaction rules specified during initialization.

        Args:
            monomer_data (tuple): A tuple containing a unique identifier (ID) and a
                chemical molecule (Mol) representing the monomer.

        Returns:
            dict: A dictionary containing information about the generated polymer,
                including its unique identifier (zinc_ids), SMILES representation
                (polymer), and the SMILES representation of the reactants (reactants).

                Example:
                {
                    "zinc_ids": "unique_id",
                    "polymer": "SMILES_representation_of_polymer",
                    "reactants": "SMILES_representation_of_reactants"
                }

                If a valid polymer cannot be formed, the "polymer" value will be None.
        """
        _id1, mol1 = tuple(monomer_data[0][0])
        _id2, mol2 = tuple(monomer_data[0][1])
        reactants = f"{Chem.MolToSmiles(mol1)}.{Chem.MolToSmiles(mol2)}"
        _id = f"{_id1}.{_id2}"

        intermediate = canonical_mol(self.rxn.RunReactants((mol1, mol2))[0][0])
        if self.P1:
            intermediate = canonical_mol(self.P1.RunReactants((intermediate,))[0][0])
        if self.P2:
            product = canonical_mol(self.P2.RunReactants((intermediate,))[0][0])
        else:
            product = intermediate

        # If product is a valid molecule, let's return its SMILES string along
        # with other information
        if product:
            return {
                "zinc_ids": _id,
                "polymer": product,
                "reactants": reactants,
            }

        # If product is not a valid molecule, then make its SMILES string None.
        else:
            return {
                "zinc_ids": _id,
                "polymer": None,
                "reactants": reactants,
            }

    def check_backbiting(self, row):
        """
        Check if a polymer is likely to backbite.

        Args:
            row (pd.Series): A row from the input DataFrame.

        Returns:
            pd.Series: A modified row with information about backbiting.
        """
        polymer = row["polymer"]
        likely_to_backbite = not self.bb_constraint.check(mol=polymer)
        if likely_to_backbite:
            polymer = None
        else:
            polymer = Chem.MolToSmiles(polymer)
        row["polymer"] = polymer
        return row

    def generate_polymers(self, molecule_df):
        """
        Generate polymers by combining monomers according to defined rules.

        This method takes a DataFrame containing molecules, evaluates their
        suitability for polymerization based on predefined criteria, and generates
        polymers by linking monomers together. The process includes discarding
        unsuitable molecules, filtering molecules based on specific criteria, and
        linking monomers to create polymers.

        Args:
            molecule_df (pd.DataFrame): Contains two columns, "id-smiles" and "index".

        Returns:
            pd.DataFrame: A DataFrame containing information about the generated
                polymers, including their unique identifiers ("zinc_ids"),
                SMILES representations of the polymers ("polymer"), and SMILES
                representations of the reactants used to form each polymer
                ("reactants").

        Note:
            The criteria for molecule suitability, reaction rules, and reactants are
            defined during object initialization. The method applies these criteria
            to filter and link monomers to create polymers.
        """
        # Discard molecules without any active sites.
        molecule_df = self.discard_molecules(molecule_df)

        # Keep molecules that only have either 2 amide groups or two 2 anhydride groups
        molecule_df = molecule_df.loc[
            (
                (molecule_df.template_1matches == 2)
                & (molecule_df.template_2matches == 0)
            )
            | (
                (molecule_df.template_1matches == 0)
                & (molecule_df.template_2matches == 2)
            )
        ]

        # Split molecule_df into chunks.
        index_slices = sliced(range(len(molecule_df)), get_monomers_max_chunk_size)
        n_chunks = int(ceil(len(molecule_df) / get_monomers_max_chunk_size))
        chunks = []
        print("\n")
        for ind, index_slice in enumerate(index_slices):
            print(f"Processing chunk {ind+1} of {n_chunks} ...", flush=True)
            chunk = molecule_df.iloc[index_slice]

            # Only keep molecules with certain ring sizes.
            chunk = daskify_df2df(
                lambda df: df.apply(lambda x: self.check_charge(x), axis=1),
                chunk,
                self.npartitions,
            )

            # Drop molecules with mol = None and then append this chunk
            # to a list.
            chunk = chunk.dropna(subset="mol")
            chunks.append(chunk)

        print("\nFinished processing chunks.")

        # Create molecule_df from chunks.
        if not chunks:
            # If there are no chunks, return an empty DataFrame.
            molecule_df = pd.DataFrame(
                {
                    "ID": [],
                    "mol": [],
                    "smiles": [],
                }
            )
        else:
            molecule_df = pd.concat(chunks)
            del chunks  # Free up memory
        print(
            f"\n{len(molecule_df)} {self.name} monomer(s) found.",
            flush=True,
        )

        # Pair monomers and convert them to a pandas Series.
        df1 = molecule_df[molecule_df.template_1matches == 2]
        df2 = molecule_df[molecule_df.template_2matches == 2]
        pair_ls = pd.Series(
            list(
                zip(
                    it.product(
                        df1[["ID", "mol"]].values.tolist(),
                        df2[["ID", "mol"]].values.tolist(),
                    )
                )
            )
        )

        # Link reactants to form polymers
        polymer_df = daskify_df2df(
            lambda df: df.apply((lambda x: self.link_monomers(x))),
            pair_ls,
            self.npartitions,
        )

        # Drop rows that contain NaN values in the "polymer" column.
        if not polymer_df.empty:
            polymer_df = polymer_df.dropna(subset="polymer")

        # Remove polymers with high chance of backbiting
        # TODO: Daskify the function applied below?
        polymer_df = polymer_df.apply(lambda x: self.check_backbiting(x), axis=1)
        polymer_df = polymer_df.dropna(subset="polymer")

        print(f"\n{len(polymer_df)} {self.name} polymer(s) found.", flush=True)
        return polymer_df


class ROMP(Polymerize):
    def __init__(self, npartitions=-1):
        super().__init__(
            name="ROMP",
            reaction_class="addition",
            rxn_smarts=ks.romp_rxn_smarts,
            polymerization1_smarts=None,
            polymerization2_smarts=None,
            substruct_1_exception=None,
            substruct_2_exception=None,
            npartitions=npartitions,
        )
        self.ring_constraint = DefaultRingSizeConstraint(
            debug=False,
            substruct=Chem.MolFromSmarts("[CH1R]=[CH1R]"),
        )

    def check_ring_size(self, row):
        _id, smiles = row["id-smiles"]
        mol = Chem.MolFromSmiles(smiles)
        # If the molecule does not exist, then we cannot check the ring size. Set mol = None.
        # Or, if the molecule fails the ring constraint, set mol = None.
        if mol is None or not self.ring_constraint.check(mol=mol):
            mol = None
        return {
            "ID": _id,
            "smiles": smiles,
            "mol": mol,
        }

    def link_monomers(self, monomer_data):
        """
        Link monomers to form polymers based on reaction rules.

        This method takes a monomer represented as a pair of a unique identifier (ID)
        and a chemical molecule (Mol) and attempts to create a polymer by applying
        the reaction rules specified during initialization.

        Args:
            monomer_data (tuple): A tuple containing a unique identifier (ID) and a
                chemical molecule (Mol) representing the monomer.

        Returns:
            dict: A dictionary containing information about the generated polymer,
                including its unique identifier (zinc_ids), SMILES representation
                (polymer), and the SMILES representation of the reactants (reactants).

                Example:
                {
                    "zinc_ids": "unique_id",
                    "polymer": "SMILES_representation_of_polymer",
                    "reactants": "SMILES_representation_of_reactants"
                }

                If a valid polymer cannot be formed, the "polymer" value will be None.
        """

        _id, mol = tuple(monomer_data)
        reactants = Chem.MolToSmiles(mol)

        intermediate = canonical_mol(self.rxn.RunReactants((mol,))[0][0])
        if self.P1:
            intermediate = canonical_mol(self.P1.RunReactants((intermediate,))[0][0])
        if self.P2:
            product = canonical_mol(self.P2.RunReactants((intermediate,))[0][0])
        else:
            product = intermediate
        # If product is a valid molecule, let's return its SMILES string along
        # with other information
        if product:
            return {
                "zinc_ids": _id,
                "polymer": Chem.MolToSmiles(product),
                "reactants": reactants,
            }
        # If product is not a valid molecule, then make its SMILES string None.
        else:
            return {
                "zinc_ids": _id,
                "polymer": None,
                "reactants": reactants,
            }

    def generate_polymers(self, molecule_df):
        """
        Generate polymers by combining monomers according to defined rules.

        This method takes a DataFrame containing molecules, evaluates their
        suitability for polymerization based on predefined criteria, and generates
        polymers by linking monomers together. The process includes discarding
        unsuitable molecules, filtering molecules based on specific criteria, and
        linking monomers to create polymers.

        Args:
            molecule_df (pd.DataFrame): Contains two columns, "id-smiles" and "index".

        Returns:
            pd.DataFrame: A DataFrame containing information about the generated
                polymers, including their unique identifiers ("zinc_ids"),
                SMILES representations of the polymers ("polymer"), and SMILES
                representations of the reactants used to form each polymer
                ("reactants").

        Note:
            The criteria for molecule suitability, reaction rules, and reactants are
            defined during object initialization. The method applies these criteria
            to filter and link monomers to create polymers.
        """
        # Discard molecules without any active sites.
        molecule_df = self.discard_molecules(molecule_df)

        # Only keep molecules with one active site.
        molecule_df = molecule_df.loc[molecule_df.template_1matches == 1]

        # Split molecule_df into chunks.
        index_slices = sliced(range(len(molecule_df)), get_monomers_max_chunk_size)
        n_chunks = int(ceil(len(molecule_df) / get_monomers_max_chunk_size))
        chunks = []
        print("\n")
        for ind, index_slice in enumerate(index_slices):
            print(f"Processing chunk {ind+1} of {n_chunks} ...", flush=True)
            chunk = molecule_df.iloc[index_slice]

            # Only keep molecules with certain ring sizes.
            chunk = daskify_df2df(
                lambda df: df.apply(lambda x: self.check_ring_size(x), axis=1),
                chunk,
                self.npartitions,
            )

            # Drop molecules with mol = None.
            chunk = chunk.dropna(subset="mol")
            chunks.append(chunk)
        print("\nFinished processing chunks.")

        # Create molecule_df from chunks.
        if not chunks:
            # If there are no chunks, return an empty DataFrame.
            molecule_df = pd.DataFrame(
                {
                    "ID": [],
                    "mol": [],
                    "smiles": [],
                }
            )
        else:
            molecule_df = pd.concat(chunks)
            del chunks  # Free up memory
        print(
            f"\n{len(molecule_df)} {self.name} monomer(s) found.",
            flush=True,
        )

        # Extract monomers as a pandas Series.
        monomer_ls = molecule_df[["ID", "mol"]].values.tolist()
        monomer_ls = pd.Series(monomer_ls)

        # Link reactants to form polymers
        polymer_df = daskify_df2df(
            lambda df: df.apply((lambda x: self.link_monomers(x))),
            monomer_ls,
            self.npartitions,
        )

        # Drop rows that contain NaN values in the "polymer" column.
        if not polymer_df.empty:
            polymer_df = polymer_df.dropna(subset="polymer")

        print(f"\n{len(polymer_df)} {self.name} polymer(s) found.", flush=True)
        return polymer_df
