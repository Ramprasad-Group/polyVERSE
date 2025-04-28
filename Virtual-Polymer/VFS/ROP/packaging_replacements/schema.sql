Table "homopolymers" {
  "hp_id" bigserial [not null, increment]
  "smiles" text [not null]
  "canonical_smiles" text
  "fingerprint" json
  "fingerprint_version" text
  "category" public.categorytype
  "canonicalization_version" text

Indexes {
  hp_id [pk, name: "homopolymers_pkey"]
  canonical_smiles [type: btree, name: "ix_homopolymers_canonical_smiles"]
  category [type: btree, name: "ix_homopolymers_category"]
  smiles [type: btree, name: "ix_homopolymers_smiles"]
}
}

Table "ml_models" {
  "ml_model_id" serial4 [not null]
  "name" text [not null]
  "algorithm" text [not null]

Indexes {
  ml_model_id [pk, name: "ml_models_pkey"]
}
}

Table "molecule_suppliers" {
  "molecule_supplier_id" serial4 [not null]
  "name" text [not null]
  "site" text

Indexes {
  molecule_supplier_id [pk, name: "molecule_suppliers_pkey"]
}
}

Table "molecules" {
  "molecule_id" bigserial [not null, increment]
  "smiles" text [not null]
  "cas" text
  "canonical_smiles" text
  "rdkit_version" text
  "fingerprint" json
  "fingerprint_version" text
  "category" public.categorytype
  "sa_score" float8

Indexes {
  molecule_id [pk, name: "molecules_pkey"]
  canonical_smiles [type: btree, name: "ix_molecules_canonical_smiles"]
  category [type: btree, name: "ix_molecules_category"]
  smiles [type: btree, name: "ix_molecules_smiles"]
}
}

Table "properties" {
  "prop_id" serial4 [not null]
  "name" text [not null]
  "unit" text

Indexes {
  prop_id [pk, name: "properties_pkey"]
}
}

Table "reaction_procedures" {
  "reaction_procedure_id" serial4 [not null]
  "description" text [not null]
  "run_status" public.vfs_run_status
  "monomer_generation_step" int4

Indexes {
  reaction_procedure_id [pk, name: "reaction_procedures_pkey"]
}
}

Table "reaction_reactants" {
  "reactant_id" serial4 [not null]
  "smiles" text [not null]

Indexes {
  reactant_id [pk, name: "reaction_reactants_pkey"]
}
}

Table "reactions" {
  "reaction_id" serial4 [not null]
  "smarts" text [not null]
  "description" text

Indexes {
  reaction_id [pk, name: "reactions_pkey"]
}
}

Table "substructures" {
  "substructure_id" serial4 [not null]
  "smarts" text [not null]
  "description" text

Indexes {
  substructure_id [pk, name: "substructures_pkey"]
  smarts [type: btree, name: "ix_substructures_smarts"]
}
}

Table "ml_model_properties" {
  "ml_model_id" int4 [not null]
  "prop_id" int4 [not null]
  "calculation_method" public.calculationmethod [not null]

Indexes {
  (ml_model_id, prop_id, calculation_method) [pk, name: "ml_model_properties_pkey"]
}
}

Table "ml_model_versions" {
  "ml_model_version_id" serial4 [not null]
  "ml_model_id" int4 [not null]
  "version" text [not null]
  "developer" text [not null]
  "date_uploaded" timestamptz [not null]

Indexes {
  ml_model_version_id [pk, name: "ml_model_versions_pkey"]
}
}

Table "molecule_costs" {
  "molecule_cost_id" serial4 [not null]
  "mol_cost" float8 [not null]
  "cost_unit" text [not null]
  "amount" float8 [not null]
  "amount_unit" text [not null]
  "usd_cost_per_gram" float8 [not null]
  "molecule_id" int8 [not null]
  "molecule_supplier_id" int4 [not null]
  "date_uploaded" timestamptz [not null]

Indexes {
  molecule_cost_id [pk, name: "molecule_costs_pkey"]
  molecule_id [type: btree, name: "ix_molecule_costs_molecule_id"]
}
}

Table "molecule_substructures" {
  "substructure_id" int4 [not null]
  "molecule_id" int4 [not null]
  "count" int4 [not null]

Indexes {
  (substructure_id, molecule_id) [pk, name: "molecule_substructures_pkey"]
  count [type: btree, name: "ix_molecule_substructures_count"]
  molecule_id [type: btree, name: "ix_molecule_substructures_molecule_id"]
}
}

Table "molecule_substructures_searched_to" {
  "substructure_id" int4 [not null]
  "molecule_id" int4 [not null]

Indexes {
  substructure_id [pk, name: "molecule_substructures_searched_to_pkey"]
}
}

Table "molecule_supplier_mol_ids" {
  "molecule_supplier_mol_id" text [not null]
  "molecule_id" int8 [not null]
  "molecule_supplier_id" int4 [not null]

Indexes {
  (molecule_id, molecule_supplier_id) [pk, name: "molecule_supplier_molecule_id_map"]
  molecule_id [type: btree, name: "ix_molecule_supplier_mol_ids_molecule_id"]
}
}

Table "polymers" {
  "pol_id" bigserial [not null, increment]
  "hp_id" int8
  "copol_id" int8

Indexes {
  pol_id [pk, name: "polymers_pkey"]
  copol_id [type: btree, name: "ix_polymers_copol_id"]
  hp_id [type: btree, name: "ix_polymers_hp_id"]
}
}

Table "reaction_polymer_mappings" {
  "reaction_pol_map_id" text [not null]
  "reaction_procedure_id" int4 [not null]
  "pol_id" int8 [not null]
  "monomer_id" int8

Indexes {
  reaction_pol_map_id [pk, name: "reaction_polymer_mappings_pkey"]
  monomer_id [type: btree, name: "ix_reaction_polymer_mappings_monomer_id"]
  pol_id [type: btree, name: "ix_reaction_polymer_mappings_pol_id"]
  reaction_procedure_id [type: btree, name: "ix_reaction_polymer_mappings_reaction_procedure_id"]
}
}

Table "reaction_reactants_to_find" {
  "reactant_to_find_id" serial4 [not null]
  "reaction_procedure_id" int4 [not null]
  "step" int4 [not null]

Indexes {
  reactant_to_find_id [pk, name: "reaction_reactants_to_find_pkey"]
}
}

Table "reaction_step_reactants" {
  "reaction_procedure_id" int4 [not null]
  "reactant_id" int4 [not null]
  "step" int4 [not null]

Indexes {
  (reaction_procedure_id, reactant_id, step) [pk, name: "reaction_step_reactants_pkey"]
}
}

Table "reaction_steps" {
  "reaction_id" int4 [not null]
  "reaction_procedure_id" int4 [not null]
  "step" int4 [not null]

Indexes {
  (reaction_id, reaction_procedure_id, step) [pk, name: "reaction_steps_pkey"]
}
}

Table "ml_model_polymer_molecule_property_predictions" {
  "pol_id" int8 [not null]
  "molecule_id" int8 [not null]
  "prop_id" int4 [not null]
  "ml_model_version_id" int4 [not null]
  "value" float8 [not null]
  "uncertainty_value" float8
  "conditions" json
  "calculation_method" public.calculationmethod [not null, default: `'exp'::calculationmethod`]
  "uncertainty_type" public.propertyuncertaintytype

Indexes {
  (pol_id, molecule_id, prop_id, calculation_method, ml_model_version_id) [pk, name: "ml_model_polymer_molecule_property_predictions_pkey"]
  ml_model_version_id [type: btree, name: "ix_ml_model_polymer_molecule_property_predictions_ml_mo_87a5"]
  molecule_id [type: btree, name: "ix_ml_model_polymer_molecule_property_predictions_molecule_id"]
  pol_id [type: btree, name: "ix_ml_model_polymer_molecule_property_predictions_pol_id"]
  prop_id [type: btree, name: "ix_ml_model_polymer_molecule_property_predictions_prop_id"]
}
}

Table "ml_model_property_predictions" {
  "pol_id" int8 [not null]
  "prop_id" int4 [not null]
  "calculation_method" public.calculationmethod [not null]
  "ml_model_version_id" int4 [not null]
  "value" float8 [not null]
  "uncertainty_value" float8
  "uncertainty_type" public.propertyuncertaintytype
  "conditions" json

Indexes {
  (pol_id, prop_id, calculation_method, ml_model_version_id) [pk, name: "ml_model_property_predictions_pkey"]
  ml_model_version_id [type: btree, name: "ix_ml_model_property_predictions_ml_model_version_id"]
  pol_id [type: btree, name: "ix_ml_model_property_predictions_pol_id"]
  prop_id [type: btree, name: "ix_ml_model_property_predictions_prop_id"]
}
}

Table "polymer_properties" {
  "pol_prop_id" int4 [not null, default: `nextval('homopolymer_properties_hp_prop_id_seq'::regclass)`]
  "pol_id" int8 [not null]
  "prop_id" int4 [not null]
  "value" float8 [not null]
  "uncertainty_value" float8
  "uncertainty_type" public.propertyuncertaintytype
  "calculation_method" public.calculationmethod [not null]
  "reference" text
  "conditions" json
  "note" text
  "date_uploaded" timestamptz [not null]
  "uploader" text

Indexes {
  pol_prop_id [pk, name: "homopolymer_properties_pkey"]
}
}

Table "reaction_found_reactant_map" {
  "molecule_id" int8 [not null]
  "reactant_to_find_id" int4 [not null]
  "reaction_pol_map_id" text [not null]

Indexes {
  (molecule_id, reactant_to_find_id, reaction_pol_map_id) [pk, name: "reaction_found_reactant_map_pkey"]
}
}

Table "reaction_reactant_last_offset" {
  "reactant_to_find_id" int4 [not null]
  "last_offset" int8 [not null]

Indexes {
  reactant_to_find_id [pk, name: "reaction_reactant_last_offset_pkey"]
}
}

Table "reaction_reactant_necessary_substructures" {
  "substructure_id" int4 [not null]
  "reactant_to_find_id" int4 [not null]
  "greater_than_or_equal" int4 [not null]
  "less_than_or_equal" int4

Indexes {
  (substructure_id, reactant_to_find_id) [pk, name: "reaction_reactant_necessary_substructures_pkey"]
}
}

Table "reaction_reactant_unacceptable_substructures" {
  "substructure_id" int4 [not null]
  "reactant_to_find_id" int4 [not null]

Indexes {
  (substructure_id, reactant_to_find_id) [pk, name: "reaction_reactant_unacceptable_substructures_pkey"]
}
}

Ref "ml_model_properties_ml_model_id_fkey":"ml_models"."ml_model_id" < "ml_model_properties"."ml_model_id"

Ref "ml_model_properties_prop_id_fkey":"properties"."prop_id" < "ml_model_properties"."prop_id"

Ref "ml_model_versions_ml_model_id_fkey":"ml_models"."ml_model_id" < "ml_model_versions"."ml_model_id"

Ref "molecule_costs_molecule_id_fkey":"molecules"."molecule_id" < "molecule_costs"."molecule_id"

Ref "molecule_costs_molecule_supplier_id_fkey":"molecule_suppliers"."molecule_supplier_id" < "molecule_costs"."molecule_supplier_id"

Ref "molecule_substructures_molecule_id_fkey":"molecules"."molecule_id" < "molecule_substructures"."molecule_id"

Ref "molecule_substructures_substructure_id_fkey":"substructures"."substructure_id" < "molecule_substructures"."substructure_id"

Ref "molecule_substructures_searched_to_molecule_id_fkey":"molecules"."molecule_id" < "molecule_substructures_searched_to"."molecule_id"

Ref "molecule_substructures_searched_to_substructure_id_fkey":"substructures"."substructure_id" < "molecule_substructures_searched_to"."substructure_id"

Ref "molecule_supplier_mol_ids_molecule_id_fkey":"molecules"."molecule_id" < "molecule_supplier_mol_ids"."molecule_id"

Ref "molecule_supplier_mol_ids_molecule_supplier_id_fkey":"molecule_suppliers"."molecule_supplier_id" < "molecule_supplier_mol_ids"."molecule_supplier_id"

Ref "polymers_hp_id_fkey":"homopolymers"."hp_id" < "polymers"."hp_id"

Ref "reaction_homopolymer_mappings_reaction_procedure_id_fkey":"reaction_procedures"."reaction_procedure_id" < "reaction_polymer_mappings"."reaction_procedure_id"

Ref "reaction_polymer_mappings_monomer_id_fkey":"molecules"."molecule_id" < "reaction_polymer_mappings"."monomer_id"

Ref "reaction_polymer_mappings_pol_id_fkey":"polymers"."pol_id" < "reaction_polymer_mappings"."pol_id"

Ref "reaction_reactants_to_find_reaction_procedure_id_fkey":"reaction_procedures"."reaction_procedure_id" < "reaction_reactants_to_find"."reaction_procedure_id"

Ref "reaction_step_reactants_reaction_id_fkey":"reaction_reactants"."reactant_id" < "reaction_step_reactants"."reactant_id"

Ref "reaction_step_reactants_reaction_procedure_id_fkey":"reaction_procedures"."reaction_procedure_id" < "reaction_step_reactants"."reaction_procedure_id"

Ref "reaction_steps_reaction_id_fkey":"reactions"."reaction_id" < "reaction_steps"."reaction_id"

Ref "reaction_steps_reaction_procedure_id_fkey":"reaction_procedures"."reaction_procedure_id" < "reaction_steps"."reaction_procedure_id"

Ref "ml_model_polymer_molecule_property_pre_ml_model_version_id_fkey":"ml_model_versions"."ml_model_version_id" < "ml_model_polymer_molecule_property_predictions"."ml_model_version_id"

Ref "ml_model_polymer_molecule_property_predictions_molecule_id_fkey":"molecules"."molecule_id" < "ml_model_polymer_molecule_property_predictions"."molecule_id"

Ref "ml_model_polymer_molecule_property_predictions_pol_id_fkey":"polymers"."pol_id" < "ml_model_polymer_molecule_property_predictions"."pol_id"

Ref "ml_model_polymer_molecule_property_predictions_prop_id_fkey":"properties"."prop_id" < "ml_model_polymer_molecule_property_predictions"."prop_id"

Ref "ml_model_property_predictions_ml_model_version_id_fkey":"ml_model_versions"."ml_model_version_id" < "ml_model_property_predictions"."ml_model_version_id"

Ref "ml_model_property_predictions_pol_id_fkey":"polymers"."pol_id" < "ml_model_property_predictions"."pol_id"

Ref "ml_model_property_predictions_prop_id_fkey":"properties"."prop_id" < "ml_model_property_predictions"."prop_id"

Ref "homopolymer_properties_prop_id_fkey":"properties"."prop_id" < "polymer_properties"."prop_id"

Ref "polymer_properties_pol_id_fkey":"polymers"."pol_id" < "polymer_properties"."pol_id"

Ref "reaction_found_reactant_map_molecule_id_fkey":"molecules"."molecule_id" < "reaction_found_reactant_map"."molecule_id"

Ref "reaction_found_reactant_map_reactant_to_find_id_fkey":"reaction_reactants_to_find"."reactant_to_find_id" < "reaction_found_reactant_map"."reactant_to_find_id"

Ref "reaction_found_reactant_map_reaction_pol_map_id_fkey":"reaction_polymer_mappings"."reaction_pol_map_id" < "reaction_found_reactant_map"."reaction_pol_map_id"

Ref "reaction_reactant_last_offset_reactant_to_find_id_fkey":"reaction_reactants_to_find"."reactant_to_find_id" < "reaction_reactant_last_offset"."reactant_to_find_id"

Ref "reaction_reactant_necessary_substructu_reactant_to_find_id_fkey":"reaction_reactants_to_find"."reactant_to_find_id" < "reaction_reactant_necessary_substructures"."reactant_to_find_id"

Ref "reaction_reactant_necessary_substructures_substructure_id_fkey":"substructures"."substructure_id" < "reaction_reactant_necessary_substructures"."substructure_id"

Ref "reaction_reactant_unacceptable_substru_reactant_to_find_id_fkey":"reaction_reactants_to_find"."reactant_to_find_id" < "reaction_reactant_unacceptable_substructures"."reactant_to_find_id"

Ref "reaction_reactant_unacceptable_substructur_substructure_id_fkey":"substructures"."substructure_id" < "reaction_reactant_unacceptable_substructures"."substructure_id"
