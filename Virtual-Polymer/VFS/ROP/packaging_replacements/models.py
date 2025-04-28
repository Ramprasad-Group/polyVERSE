from typing import Optional, Union

from pydantic import BaseModel


class Substructure(BaseModel):
    """
    Attributes:
        smarts: SMARTS substructure for the query
        description: description of the substructure
        reference: Optional reference where the SMARTS was found
        greater_than_or_equal: For necessary substructures, how many of the specific
            substructure must be in the procedure.
        less_than_or_equal: For necessary substructures, how few of the specific
            substructure must be in the procedure. Must be greater than or equal to the
            "greater_than_or_equal" number.
    """

    smarts: str
    description: str
    reference: Optional[str] = None
    greater_than_or_equal: Optional[int] = 1
    less_than_or_equal: Optional[int] = 1


class ReactantsToFind(BaseModel):
    necessary_substructures: list[Substructure]
    unacceptable_substructures: list[Substructure]


class Reaction(BaseModel):
    """
    Attributes:
        smarts: SMARTS substructure for the reaction
        description: Explains the reaction
        reference: Optional reference for where the reaction came from
        reactants: Optional list of additional reactants needed for the reaction
        reactants_to_find: Optional list of reactants to find in a molecule list
            that would be used in this reaction. The reactants are searched to ensure
            they have all necessary substructures and no unacceptable ones.
    """

    smarts: str
    description: str
    reference: Optional[str] = None
    additional_reactants: Optional[list[str]] = []
    reactants_to_find: Optional[list[ReactantsToFind]] = []


class ReactionProcedure(BaseModel):
    """
    Attributes:
        reactions: List of rdkit smarts reactions (e.g.,
            [
                [C:6]1[C:1][O:2][C:3](=[O:4])[C:5]1.[S:7]>>[C:6]1[C:1][O:2][C:3](=[S:7]
                )[C:5]1.[O:4],
                [C:6]1[C:1][O:2][C:3](=[S:7])[C:5]1.[O:4]>>[*][C:3](=[O:2]
                )[C:5][C:6][C:1][S:7][*].[O:4]
            ]
        )
        description: Explains the procedure
        note: Optional notes about the procedure (like if any specific solvents or
            reaction conditions were used)
        reference: Optional doi or description of where the reaction was found
        monomer_generation_step: the step after which the monomer for the polymer is
            generated. If None, then assume the monomer is the starting molecule in the
            reaction.
    """

    reactions: list[Reaction]
    description: str
    note: Optional[str] = None
    reference: Optional[str] = None
    monomer_generation_step: Optional[int] = None


class Homopolymer(BaseModel):
    """
    Attributes:
        smiles: smiles string of the result of the first molecule that
            ran through the reaction procedure and created this polymer
        monomer: smiles of the monomer
    """

    polymer_smiles: str
    monomer_smiles: Optional[str]


class Molecule(BaseModel):
    smiles: str


class ReactionStep(BaseModel):
    additional_reactants: list[str]
    molecules: Optional[list[Molecule]] = []
    smarts: str


class PolymerizationReaction(BaseModel):
    monomer_generation_step: Optional[int] = None
    reactions: list[ReactionStep]
