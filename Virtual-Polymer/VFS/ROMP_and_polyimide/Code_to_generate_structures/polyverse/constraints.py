from rdkit import Chem
import numpy as np
import names


def pretty_print(obj):
    """
    Convert an object to a pretty-printed representation.

    Args:
        obj: The input object.

    Returns:
        str: Pretty-printed representation of the object.
    """
    if isinstance(obj, Chem.rdchem.Mol):
        return Chem.MolToSmarts(obj)
    else:
        return obj


def pretty_print_dict(dictionary):
    """
    Pretty-print a dictionary of values.

    Args:
        dictionary (dict): The input dictionary.

    Returns:
        dict: Dictionary with pretty-printed values.
    """
    return {k: pretty_print(v) for k, v in dictionary.items()}


class Constraint:
    """
    Base class for defining constraints.

    Attributes:
        pass_if (bool): A boolean indicating whether the constraint passes or fails.
        comparator (str): The comparison operator for the constraint.
        threshold: The threshold value for comparison.
        debug (bool): A boolean indicating whether to enable debugging.
        name (str): A random name assigned to the constraint.
    """

    def __init__(self, comparator, threshold, debug=False):
        self.pass_if = True
        self.comparator = comparator
        self.threshold = threshold
        self.debug = debug
        self.name = names.get_full_name()

    def _compute_variable(self, **kwargs):
        """
        Placeholder for computing the constraint variable.

        Args:
            **kwargs: Additional keyword arguments for computation.

        Returns:
            Variable: The computed variable.
        """
        raise NotImplementedError

    def check(self, **kwargs):
        """
        Check if the constraint is satisfied.

        Args:
            **kwargs: Keyword arguments needed for computing the constraint variable.

        Returns:
            bool: True if the constraint is satisfied, otherwise False.
        """
        variable = self._compute_variable(**kwargs)
        if self.comparator == "==":
            return_bool = variable == self.threshold
        elif self.comparator == "<":
            return_bool = variable < self.threshold
        elif self.comparator == ">":
            return_bool = variable > self.threshold
        elif self.comparator == "<=":
            return_bool = variable <= self.threshold
        elif self.comparator == ">=":
            return_bool = variable >= self.threshold
        elif self.comparator == "in":
            return_bool = variable in self.threshold

        if self.debug:
            if not return_bool:
                print(
                    f"{self.__class__.__name__} failed on {pretty_print_dict(kwargs)} with dict {pretty_print_dict(self.__dict__)}"
                )
            return return_bool
        else:
            return return_bool


class ChargeConstraint(Constraint):
    def __init__(self, comparator, threshold, debug, substruct1, substruct2) -> None:
        """
        Charge constraint for comparing partial charges between two substructures.

        Args:
            comparator (str): Comparison operator (e.g., '==', '<', '>').
            threshold (float): Threshold value for the comparison.
            debug (bool): Debug mode flag.
            substruct1 (list): Atom indices of the first substructure to test.
            substruct2 (list): Atom indices of the second substructure to compare with substruct1.
        """
        super().__init__(comparator, threshold, debug)
        self.substruct1 = substruct1
        self.substruct2 = substruct2

    def _compute_variable(self, mol):
        """
        Calculate the maximum absolute difference in partial charges between atoms
        in substruct1 and substruct2.

        Args:
            mol (rdkit.Chem.Mol): RDKit molecule to calculate charges on.

        Returns:
            float: Maximum absolute charge difference.
        """
        mol.ComputeGasteigerCharges()
        return np.max(
            np.abs(
                [
                    mol.GetAtomWithIdx(x).GetDoubleProp("_GasteigerCharge")
                    - mol.GetAtomWithIdx(y).GetDoubleProp("_GasteigerCharge")
                    for (x, y) in zip(self.substruct1, self.substruct2)
                ]
            )
        )


class BackbiteConstraint(Constraint):
    """
    Constraint for backbiting.

    Attributes:
        smartsA (str): SMARTS pattern for atom A.
        smartsB (str): SMARTS pattern for atom B.
        atomA (rdkit.Chem.Mol): RDKit molecule for atom A.
        atomB (rdkit.Chem.Mol): RDKit molecule for atom B.
    """

    def __init__(self, debug=False, min_dist=8):
        """
        Initialize the BackbiteConstraint.

        Args:
            debug (bool): Debug mode flag.
            min_dist (int): Minimum distance between atoms A and B.
        """
        comparator, threshold, smartsA, smartsB = "==", True, "[#0]", "[#0]"
        super().__init__(comparator, threshold, debug)
        self.smartsA = smartsA
        self.smartsB = smartsB
        self.atomA = Chem.MolFromSmarts(smartsA)
        self.atomB = Chem.MolFromSmarts(smartsB)
        self.min_dist = min_dist

    def _compute_variable(self, mol):
        """
        Check if the shortest path distance between atoms A and B is greater than or
        equal to the minimum required distance.

        Args:
            mol (rdkit.Chem.Mol): RDKit molecule to analyze.

        Returns:
            bool: True if the condition is met, False otherwise.
        """
        self.indexA = mol.GetSubstructMatches(self.atomA)
        self.indexB = mol.GetSubstructMatches(self.atomB)

        if self.smartsA == self.smartsB:
            assert len(self.indexA) == 2
            assert len(self.indexB) == 2
            self.indexA = self.indexA[0]
            self.indexB = self.indexB[1]
        else:
            assert len(self.indexA) == 1
            assert len(self.indexB) == 1

        self.indexA, self.indexB = self.indexA[0], self.indexB[0]
        sp = Chem.GetShortestPath(mol, self.indexA, self.indexB)  # shortest path
        dist = len(sp)
        return dist >= self.min_dist


class RingSizeConstraint(Constraint):
    """
    Constraint for evaluating ring size.

    Attributes:
        substruct (rdkit.Chem.Mol): RDKit molecule representing the substructure.
    """

    def __init__(self, threshold, debug, substruct):
        comparator = "in"
        super().__init__(comparator, threshold, debug)
        self.substruct = substruct

    def _compute_variable(self, mol):
        atom_inds = set(mol.GetSubstructMatches(self.substruct)[0])
        ri = mol.GetRingInfo()

        for ring in ri.AtomRings():
            ring = set(ring)
            if not atom_inds.difference(ring):
                return len(ring)


class DefaultRingSizeConstraint(RingSizeConstraint):
    """
    Default constraint for ring size evaluation.

    Attributes:
        substruct (rdkit.Chem.Mol): RDKit molecule representing the substructure.
    """

    def __init__(self, debug, substruct):
        threshold = [3, 4, 5, 7, 8, 9, 10, 11]  # allowable ring sizes
        super().__init__(threshold, debug, substruct)
