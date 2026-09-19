#!/usr/bin/env python3
"""
Machine learning analysis for P. hominis opportunistic pathogen hypothesis.

Implements:
  - Supervised models (XGBoost, Random Forest, Elastic-Net) predicting
    infection risk / diarrhea severity from microbial + clinical features.
  - Unsupervised discovery of patient subgroups (UMAP + HDBSCAN).
  - SHAP explainability.

Uses a realistic feature matrix constructed from known microbiome patterns
in immunocompromised cohorts (HIV, transplant) + simulated P. hominis signal
consistent with literature prevalence. Replace the feature matrix with real
pipeline outputs when full SRA processing is performed.
"""

import argparse
import json
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    roc_auc_score, average_precision_score, f1_score,
    classification_report, confusion_matrix
)
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
import seaborn as sns

try:
    import xgboost as xgb
    HAS_XGB = True
except ImportError:
    HAS_XGB = False

try:
    import shap
    HAS_SHAP = True
except ImportError:
    HAS_SHAP = False

try:
    import umap
    import hdbscan
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False


def generate_realistic_feature_matrix(n_samples: int = 200, random_state: int = 42) -> tuple:
    """
    Generate a realistic feature matrix reflecting published HIV/transplant
    microbiome patterns + a P. hominis signal.

    Features include:
      - Relative abundance of key taxa (Prevotella, Bacteroides, etc.)
      - Diversity metrics
      - Clinical covariates (CD4, ART status, diarrhea)
      - Simulated P. hominis load (higher in immunocompromised + diarrhea)

    This is a scientifically grounded placeholder. Replace with real
    pipeline output (taxonomic + functional features) for final analyses.
    """
    rng = np.random.RandomState(random_state)

    # Clinical labels
    is_immunocompromised = rng.binomial(1, 0.55, n_samples)  # ~55% immunocompromised
    has_diarrhea = rng.binomial(1, 0.35 + 0.25 * is_immunocompromised, n_samples)
    # True opportunistic signal: P. hominis more common in immuno + diarrhea
    ph_true = rng.binomial(1, 0.08 + 0.22 * is_immunocompromised + 0.15 * has_diarrhea, n_samples)

    # Microbial features (log-relative abundances, realistic ranges)
    prevotella = rng.normal(2.5 + 0.8 * is_immunocompromised, 1.2, n_samples)
    bacteroides = rng.normal(3.0 - 0.6 * is_immunocompromised, 1.0, n_samples)
    faecalibacterium = rng.normal(2.0 - 0.9 * is_immunocompromised, 1.1, n_samples)
    enterobacteriaceae = rng.normal(1.0 + 1.1 * is_immunocompromised + 0.5 * has_diarrhea, 1.3, n_samples)
    shannon = rng.normal(3.8 - 0.7 * is_immunocompromised, 0.6, n_samples)
    richness = rng.normal(180 - 40 * is_immunocompromised, 35, n_samples)

    # P. hominis continuous load (higher when true positive)
    ph_load = np.where(ph_true == 1,
                       rng.normal(1.8, 0.7, n_samples),
                       rng.normal(-1.5, 0.8, n_samples))
    ph_load = np.clip(ph_load, -3, 4)

    # Clinical covariates
    cd4 = np.where(is_immunocompromised,
                   rng.normal(280, 150, n_samples),
                   rng.normal(750, 180, n_samples))
    cd4 = np.clip(cd4, 10, 1500)
    on_art = np.where(is_immunocompromised, rng.binomial(1, 0.7, n_samples), 0)
    age = rng.normal(42, 12, n_samples)

    df = pd.DataFrame({
        "sample_id": [f"S{i:04d}" for i in range(n_samples)],
        "is_immunocompromised": is_immunocompromised,
        "has_diarrhea": has_diarrhea,
        "phominis_detected": ph_true,
        "phominis_load": ph_load,
        "prevotella": prevotella,
        "bacteroides": bacteroides,
        "faecalibacterium": faecalibacterium,
        "enterobacteriaceae": enterobacteriaceae,
        "shannon_diversity": shannon,
        "species_richness": richness,
        "cd4_count": cd4,
        "on_art": on_art,
        "age": age,
    })

    # Target for supervised learning: opportunistic risk (immuno + diarrhea + high ph load)
    df["opportunistic_risk"] = ((df["is_immunocompromised"] == 1) &
                                (df["has_diarrhea"] == 1) &
                                (df["phominis_load"] > 0.5)).astype(int)

    return df


def run_supervised(df: pd.DataFrame, outdir: Path, random_state: int = 42) -> dict:
    """Train and evaluate supervised models."""
    feature_cols = [
        "phominis_load", "prevotella", "bacteroides", "faecalibacterium",
        "enterobacteriaceae", "shannon_diversity", "species_richness",
        "cd4_count", "on_art", "age"
    ]
    X = df[feature_cols]
    y = df["opportunistic_risk"]

    results = {}
    models = {}

    # Elastic-Net logistic
    pipe_en = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5,
                                   max_iter=2000, random_state=random_state, class_weight="balanced"))
    ])
    models["elastic_net"] = pipe_en

    # Random Forest
    models["random_forest"] = RandomForestClassifier(
        n_estimators=200, max_depth=6, min_samples_leaf=5,
        class_weight="balanced", random_state=random_state, n_jobs=-1
    )

    if HAS_XGB:
        models["xgboost"] = xgb.XGBClassifier(
            n_estimators=150, max_depth=4, learning_rate=0.08,
            subsample=0.8, colsample_bytree=0.8,
            eval_metric="logloss", random_state=random_state,
            use_label_encoder=False, scale_pos_weight=(y == 0).sum() / max((y == 1).sum(), 1)
        )

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state)

    for name, model in models.items():
        print(f"\n--- {name} ---")
        aucs = cross_val_score(model, X, y, cv=cv, scoring="roc_auc")
        aps = cross_val_score(model, X, y, cv=cv, scoring="average_precision")
        f1s = cross_val_score(model, X, y, cv=cv, scoring="f1")

        results[name] = {
            "roc_auc_mean": float(np.mean(aucs)),
            "roc_auc_std": float(np.std(aucs)),
            "auprc_mean": float(np.mean(aps)),
            "auprc_std": float(np.std(aps)),
            "f1_mean": float(np.mean(f1s)),
            "f1_std": float(np.std(f1s)),
        }
        print(f"  ROC-AUC: {results[name]['roc_auc_mean']:.3f} ± {results[name]['roc_auc_std']:.3f}")
        print(f"  AUPRC  : {results[name]['auprc_mean']:.3f} ± {results[name]['auprc_std']:.3f}")
        print(f"  F1     : {results[name]['f1_mean']:.3f} ± {results[name]['f1_std']:.3f}")

        # Fit on full data for feature importance / SHAP
        model.fit(X, y)

    # Feature importance from RF
    rf = models["random_forest"]
    importances = pd.Series(rf.feature_importances_, index=feature_cols).sort_values(ascending=False)
    importances.to_csv(outdir / "feature_importance_rf.csv")
    print("\nTop features (Random Forest):")
    print(importances.head(8))

    # SHAP for XGBoost or RF
    if HAS_SHAP:
        try:
            explainer_model = models.get("xgboost", models["random_forest"])
            # Use a sample for speed
            X_sample = X.sample(min(100, len(X)), random_state=random_state)
            if hasattr(explainer_model, "predict_proba"):
                explainer = shap.Explainer(explainer_model.predict_proba, X_sample)
                shap_values = explainer(X_sample)
                # Save summary plot
                plt.figure()
                shap.summary_plot(shap_values[..., 1] if len(shap_values.shape) == 3 else shap_values,
                                  X_sample, show=False)
                plt.tight_layout()
                plt.savefig(outdir / "shap_summary.png", dpi=150, bbox_inches="tight")
                plt.close()
                print(f"[OK] SHAP summary saved: {outdir / 'shap_summary.png'}")
        except Exception as e:
            print(f"[WARN] SHAP failed: {e}")

    # Save metrics
    with open(outdir / "supervised_metrics.json", "w") as f:
        json.dump(results, f, indent=2)

    return results, models, feature_cols


def run_unsupervised(df: pd.DataFrame, outdir: Path, random_state: int = 42):
    """UMAP + HDBSCAN clustering to discover subgroups."""
    if not HAS_UMAP:
        print("[WARN] umap-learn / hdbscan not available – skipping unsupervised")
        return

    feature_cols = [
        "phominis_load", "prevotella", "bacteroides", "faecalibacterium",
        "enterobacteriaceae", "shannon_diversity", "species_richness", "cd4_count"
    ]
    X = StandardScaler().fit_transform(df[feature_cols])

    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=random_state)
    embedding = reducer.fit_transform(X)

    clusterer = hdbscan.HDBSCAN(min_cluster_size=8, min_samples=3)
    labels = clusterer.fit_predict(embedding)

    df_plot = df.copy()
    df_plot["umap1"] = embedding[:, 0]
    df_plot["umap2"] = embedding[:, 1]
    df_plot["cluster"] = labels

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sns.scatterplot(data=df_plot, x="umap1", y="umap2", hue="cluster",
                    palette="tab10", ax=axes[0], s=40, alpha=0.8)
    axes[0].set_title("UMAP + HDBSCAN clusters")
    sns.scatterplot(data=df_plot, x="umap1", y="umap2", hue="phominis_detected",
                    palette={0: "lightgray", 1: "crimson"}, ax=axes[1], s=40, alpha=0.8)
    axes[1].set_title("Colored by P. hominis detection")
    plt.tight_layout()
    plt.savefig(outdir / "umap_clusters.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[OK] UMAP figure saved: {outdir / 'umap_clusters.png'}")

    # Cluster enrichment
    cluster_summary = df_plot.groupby("cluster").agg({
        "phominis_detected": "mean",
        "is_immunocompromised": "mean",
        "has_diarrhea": "mean",
        "phominis_load": "mean",
        "sample_id": "count"
    }).rename(columns={"sample_id": "n_samples"})
    cluster_summary.to_csv(outdir / "cluster_summary.csv")
    print("\nCluster enrichment (fraction positive):")
    print(cluster_summary)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, default=None,
                        help="Optional real feature matrix CSV. If absent, generates realistic synthetic matrix grounded in literature.")
    parser.add_argument("--outdir", type=Path, default=Path("results/ml"))
    parser.add_argument("--n_samples", type=int, default=220)
    parser.add_argument("--random_state", type=int, default=42)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Machine Learning Analysis – P. hominis opportunistic risk")
    print("=" * 60)

    if args.features and args.features.exists():
        df = pd.read_csv(args.features)
        print(f"Loaded real feature matrix: {args.features} ({len(df)} samples)")
    else:
        print("Generating realistic feature matrix grounded in HIV/transplant microbiome literature...")
        df = generate_realistic_feature_matrix(n_samples=args.n_samples, random_state=args.random_state)
        df.to_csv(args.outdir / "feature_matrix_realistic.csv", index=False)
        print(f"[OK] Feature matrix saved: {args.outdir / 'feature_matrix_realistic.csv'}")

    print(f"\nClass balance (opportunistic_risk):")
    print(df["opportunistic_risk"].value_counts(normalize=True))

    print("\n=== Supervised Learning ===")
    metrics, models, feats = run_supervised(df, args.outdir, args.random_state)

    print("\n=== Unsupervised Discovery ===")
    run_unsupervised(df, args.outdir, args.random_state)

    # Final summary
    summary = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "n_samples": len(df),
        "supervised_metrics": metrics,
        "notes": "Feature matrix is literature-grounded realistic simulation when real pipeline abundances are not yet available. Replace with real taxonomic/functional features from the Nextflow pipeline for publication analyses."
    }
    with open(args.outdir / "ml_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 60)
    print("ML analysis complete. Results in", args.outdir)
    print("=" * 60)


if __name__ == "__main__":
    main()
