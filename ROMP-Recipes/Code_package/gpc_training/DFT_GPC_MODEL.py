#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import os
from pathlib import Path
from typing import Dict, List, Tuple

import joblib
import numpy as np
import pandas as pd
from sklearn.feature_selection import VarianceThreshold
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import ConstantKernel, RBF, WhiteKernel
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


DEFAULT_EXCLUDE_COLUMNS = {
    "ID", "id", "name", "Name", "SampleName",
    "target_property", "reactivity", "reactivity.1",
    "monomer_smiles", "Monomer_smiles", "smiles", "SMILES", "smiles1", "PSMILES",
    "Polymer_smiles", "polymer_smiles", "structure", "Structure",
}

LABEL_CANDIDATES = ["target_property", "reactivity", "reactivity.1"]


def to_jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    return obj


def read_table(path: str, sheet_name=0) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(path, sheet_name=sheet_name)
    if ext == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"Unsupported input extension: {ext}. Use .xlsx, .xls, or .csv")


def load_feature_columns(path: str) -> List[str]:
    feature_path = Path(path)
    if not feature_path.exists():
        raise FileNotFoundError(f"Feature-definition file not found: {feature_path}")
    with feature_path.open("r", encoding="utf-8") as handle:
        values = json.load(handle)
    if not isinstance(
        values,
        list
    ) or not values or not all(isinstance(x, str) and x.strip() for x in values):
        raise ValueError("The feature-definition file must contain a non-empty JSON list of feature names.")
    features = [x.strip() for x in values]
    duplicates = sorted({x for x in features if features.count(x) > 1})
    if duplicates:
        raise ValueError(f"Duplicate feature names were found: {duplicates}")
    return features


def infer_label_column(df: pd.DataFrame, user_label_col: str | None = None) -> str:
    if user_label_col:
        if user_label_col not in df.columns:
            raise ValueError(f"User-specified label column not found: {user_label_col}")
        return user_label_col
    for candidate in LABEL_CANDIDATES:
        if candidate in df.columns:
            return candidate
    raise ValueError(f"Cannot find label column. Tried: {LABEL_CANDIDATES}")


def validate_and_clean_labels(df: pd.DataFrame, label_col: str) -> pd.DataFrame:
    out = df.copy()
    labels = pd.to_numeric(out[label_col], errors="coerce")
    n_missing = int(labels.isna().sum())
    if n_missing:
        print(f"[WARN] Dropping {n_missing} rows with missing or non-numeric labels in {label_col}")
    out = out.loc[labels.notna()].copy()
    out[label_col] = pd.to_numeric(out[label_col], errors="coerce").astype(int)
    values = sorted(out[label_col].unique().tolist())
    if not set(values).issubset({0, 1}):
        raise ValueError(f"{label_col} must contain binary labels 0/1. Found: {values}")
    if len(values) < 2:
        raise ValueError(f"{label_col} contains only one class: {values}")
    return out


def select_feature_columns(
    df: pd.DataFrame,
    label_col: str,
    feature_json: str,
    feature_cols: str | None = None,
    use_all_numeric_features: bool = False,


) -> Tuple[List[str], str]:
    if feature_cols:
        features = [x.strip() for x in feature_cols.split(",") if x.strip()]
        if not features:
            raise ValueError("--feature_cols did not contain any feature names.")
        missing = [x for x in features if x not in df.columns]
        if missing:
            raise ValueError(f"User-specified feature columns not found: {missing}")
        return features, "command_line_feature_cols"

    if use_all_numeric_features:
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        excluded = set(DEFAULT_EXCLUDE_COLUMNS)
        excluded.add(label_col)
        features = [x for x in numeric_cols if x not in excluded]
        if not features:
            raise RuntimeError("No numeric descriptor features remained after exclusions.")
        return features, "all_numeric_non_id_non_label_columns"

    features = load_feature_columns(feature_json)
    missing = [x for x in features if x not in df.columns]
    if missing:
        raise ValueError(
            "Input table is missing model-required DFT features:\n"
            + "\n".join(missing)
        )
    return features, str(Path(feature_json))


def build_dft_matrix(
    df: pd.DataFrame,
    feature_columns: List[str],


) -> Tuple[np.ndarray, Dict, Dict[str, float]]:
    missing = [x for x in feature_columns if x not in df.columns]
    if missing:
        raise ValueError(f"Input is missing feature columns: {missing}")

    matrix = df[feature_columns].apply(pd.to_numeric, errors="coerce")
    nan_counts = {x: int(matrix[x].isna().sum()) for x in feature_columns}
    n_nan_rows = int(matrix.isna().any(axis=1).sum())

    medians = matrix.median(numeric_only=True)
    all_nan_features = [x for x in feature_columns if pd.isna(medians.get(x))]
    if all_nan_features:
        raise ValueError(
            "The following required features contain no numeric values:\n"
            + "\n".join(all_nan_features)
        )

    training_medians = {x: float(medians[x]) for x in feature_columns}
    matrix = matrix.fillna(pd.Series(training_medians))

    if matrix.isna().any().any():
        unresolved = [x for x in feature_columns if matrix[x].isna().any()]
        raise ValueError(f"Unresolved missing values remain in: {unresolved}")

    diagnostics = {
        "n_nan_descriptor_rows_before_imputation": n_nan_rows,
        "nan_counts_by_feature": nan_counts,
        "imputation_source": "training_set_medians",
    }
    return matrix.to_numpy(dtype=np.float64), diagnostics, training_medians


def make_gpc(seed: int, args) -> Pipeline:
    kernel = (
        ConstantKernel(
            args.kernel_constant,
            constant_value_bounds=(args.constant_lower_bound, args.constant_upper_bound),
        )
        * RBF(
            length_scale=args.rbf_length_scale,
            length_scale_bounds=(args.length_scale_lower_bound, args.length_scale_upper_bound),
        )
        + WhiteKernel(
            noise_level=args.noise_level,
            noise_level_bounds=(args.noise_lower_bound, args.noise_upper_bound),
        )
    )

    classifier = GaussianProcessClassifier(
        kernel=kernel,
        optimizer=None if args.disable_kernel_optimizer else "fmin_l_bfgs_b",
        n_restarts_optimizer=args.n_restarts_optimizer,
        max_iter_predict=args.max_iter_predict,
        warm_start=False,
        copy_X_train=True,
        random_state=seed,
        multi_class="one_vs_rest",
        n_jobs=args.n_jobs,
    )

    return Pipeline([
        ("variance_filter", VarianceThreshold(threshold=args.variance_threshold)),
        ("scaler", StandardScaler()),
        ("gpc", classifier),
    ])


def evaluate_training_fit(model: Pipeline, X: np.ndarray, y: np.ndarray) -> Dict:
    probability = model.predict_proba(X)[:, 1]
    prediction = model.predict(X).astype(int)

    metrics = {
        "training_accuracy_apparent": float(accuracy_score(y, prediction)),
        "training_balanced_accuracy_apparent": float(balanced_accuracy_score(y, prediction)),
        "training_precision_pos_apparent": float(precision_score(y, prediction, zero_division=0)),
        "training_recall_pos_apparent": float(recall_score(y, prediction, zero_division=0)),
        "training_f1_pos_apparent": float(f1_score(y, prediction, zero_division=0)),
        "training_macro_f1_apparent": float(f1_score(
            y,
            prediction,
            average="macro",
            zero_division=0
        )),
        "training_weighted_f1_apparent": float(
            f1_score(y, prediction, average="weighted", zero_division=0)
        ),
        "training_confusion_matrix_labels_0_1": confusion_matrix(
            y, prediction, labels=[0, 1]
        ).tolist(),
        "training_classification_report_apparent": classification_report(
            y,
            prediction,
            labels=[0, 1],
            output_dict=True,
            zero_division=0,
        ),
        "training_roc_auc_apparent": float(roc_auc_score(y, probability)),
        "training_pr_auc_apparent": float(average_precision_score(y, probability)),
    }
    return metrics


def main():
    default_feature_json = Path(__file__).with_name("dft_feature_columns.json")

    parser = argparse.ArgumentParser(
        description="Train a binary ROMP DFT-GPC model on a labeled descriptor table."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Labeled DFT descriptor table in CSV or Excel format."
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory for model and training outputs."
    )
    parser.add_argument(
        "--sheet_name",
        default=0,
        help="Excel sheet name or index; default is the first sheet."
    )
    parser.add_argument("--label_col", default=None, help="Label column; otherwise auto-detected.")
    parser.add_argument(
        "--feature_json",
        default=str(default_feature_json),
        help="JSON list defining the default model feature columns.",
    )
    parser.add_argument(
        "--feature_cols",
        default=None,
        help="Optional comma-separated feature names overriding --feature_json.",
    )
    parser.add_argument(
        "--use_all_numeric_features",
        action="store_true",
        help="Use all numeric non-ID/non-label columns instead of the JSON feature definition.",
    )
    parser.add_argument(
        "--drop_duplicate_ids",
        action="store_true",
        help="Drop duplicate IDs and keep the last occurrence when the ID column exists.",
    )
    parser.add_argument(
        "--id_col",
        default="ID",
        help="ID column used for optional duplicate removal."
    )

    parser.add_argument("--random_state", type=int, default=40)
    parser.add_argument("--variance_threshold", type=float, default=0.0)
    parser.add_argument("--kernel_constant", type=float, default=1.0)
    parser.add_argument("--rbf_length_scale", type=float, default=1.0)
    parser.add_argument("--noise_level", type=float, default=1e-5)
    parser.add_argument("--constant_lower_bound", type=float, default=1e-3)
    parser.add_argument("--constant_upper_bound", type=float, default=1e3)
    parser.add_argument("--length_scale_lower_bound", type=float, default=1e-2)
    parser.add_argument("--length_scale_upper_bound", type=float, default=1e2)
    parser.add_argument("--noise_lower_bound", type=float, default=1e-8)
    parser.add_argument("--noise_upper_bound", type=float, default=1e-1)
    parser.add_argument("--disable_kernel_optimizer", action="store_true")
    parser.add_argument("--n_restarts_optimizer", type=int, default=5)
    parser.add_argument("--max_iter_predict", type=int, default=100)
    parser.add_argument("--n_jobs", type=int, default=-1)
    args = parser.parse_args()

    if args.constant_lower_bound <= 0 or args.constant_upper_bound <= args.constant_lower_bound:
        raise ValueError("ConstantKernel bounds must satisfy 0 < lower < upper.")
    if args.length_scale_lower_bound <= 0 or args.length_scale_upper_bound <= args.length_scale_lower_bound:
        raise ValueError("RBF length-scale bounds must satisfy 0 < lower < upper.")
    if args.noise_lower_bound <= 0 or args.noise_upper_bound <= args.noise_lower_bound:
        raise ValueError("WhiteKernel noise bounds must satisfy 0 < lower < upper.")
    if args.variance_threshold < 0:
        raise ValueError("variance_threshold must be non-negative.")
    if args.n_restarts_optimizer < 0:
        raise ValueError("n_restarts_optimizer must be non-negative.")
    if args.max_iter_predict < 1:
        raise ValueError("max_iter_predict must be at least 1.")

    sheet_name = int(args.sheet_name) if isinstance(
        args.sheet_name,
        str
    ) and args.sheet_name.isdigit() else args.sheet_name
    os.makedirs(args.output_dir, exist_ok=True)

    df = read_table(args.input, sheet_name=sheet_name)
    raw_shape = df.shape

    duplicate_count = None
    if args.drop_duplicate_ids and args.id_col in df.columns:
        duplicate_count = int(df[args.id_col].duplicated().sum())
        if duplicate_count:
            print(f"[WARN] Dropping {duplicate_count} duplicate IDs and keeping the last occurrence.")
            df = df.drop_duplicates(subset=[args.id_col], keep="last").copy()

    label_col = infer_label_column(df, args.label_col)
    df = validate_and_clean_labels(df, label_col)

    feature_columns, feature_source = select_feature_columns(
        df=df,
        label_col=label_col,
        feature_json=args.feature_json,
        feature_cols=args.feature_cols,
        use_all_numeric_features=args.use_all_numeric_features,
    )

    X, matrix_diagnostics, training_medians = build_dft_matrix(df, feature_columns)
    y = df[label_col].astype(int).to_numpy()

    model = make_gpc(args.random_state, args)
    model.fit(X, y)
    apparent_metrics = evaluate_training_fit(model, X, y)

    package = {
        "model": model,
        "feature_columns": feature_columns,
        "training_medians": training_medians,
        "args": vars(args),
        "label_col": label_col,
        "model_type": "dft_gpc_binary_full_labeled_dataset",
    }

    model_path = os.path.join(args.output_dir, "model_dft_gpc.joblib")
    joblib.dump(package, model_path)

    feature_json_path = os.path.join(args.output_dir, "dft_feature_columns.json")
    with open(feature_json_path, "w", encoding="utf-8") as handle:
        json.dump(feature_columns, handle, indent=2, ensure_ascii=False)

    training_medians_path = os.path.join(args.output_dir, "training_medians.json")
    with open(training_medians_path, "w", encoding="utf-8") as handle:
        json.dump(training_medians, handle, indent=2, ensure_ascii=False)

    summary = {
        "input": args.input,
        "raw_shape": [int(raw_shape[0]), int(raw_shape[1])],
        "training_shape_after_label_cleaning": [int(df.shape[0]), int(df.shape[1])],
        "label_col": label_col,
        "label_counts": {
            str(k): int(v)
            for k, v in df[label_col].value_counts().sort_index().items()
        },
        "duplicate_ids_removed_keep_last": duplicate_count,
        "feature_source": feature_source,
        "feature_count": int(len(feature_columns)),
        "feature_columns": feature_columns,
        "matrix_diagnostics": matrix_diagnostics,
        "kernel_after_fit": str(model.named_steps["gpc"].kernel_),
        "note": "The reported training metrics are apparent in-sample fit metrics, not cross-validation or held-out test performance.",
        **apparent_metrics,
        "model_path": model_path,
        "feature_json_path": feature_json_path,
        "training_medians_path": training_medians_path,
    }

    summary_path = os.path.join(args.output_dir, "production_training_summary.json")
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(to_jsonable(summary), handle, indent=2, ensure_ascii=False)

    config_path = os.path.join(args.output_dir, "run_config.json")
    with open(config_path, "w", encoding="utf-8") as handle:
        json.dump(to_jsonable(vars(args)), handle, indent=2, ensure_ascii=False)

    print("=== ROMP DFT-GPC training complete ===")
    print(f"Rows used: {len(df)}")
    print(f"Label counts: {df[label_col].value_counts().sort_index().to_dict()}")
    print(f"Feature count: {len(feature_columns)}")
    print(f"Rows with missing descriptors before imputation: {matrix_diagnostics['n_nan_descriptor_rows_before_imputation']}")
    print(f"Saved model: {model_path}")
    print(f"Saved feature columns: {feature_json_path}")
    print(f"Saved training medians: {training_medians_path}")
    print(f"Saved summary: {summary_path}")


if __name__ == "__main__":
    main()
