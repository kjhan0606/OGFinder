"""Train SVM classifier for LSBG candidate filtering.

Usage:
    python3 -m lsbg.training.train_svm \
        --data lsbg/data/training/labeled.tsv \
        --output lsbg/data/checkpoints/svm_lsbg_best.joblib

Pipeline: SimpleImputer(median) -> StandardScaler -> SVC(rbf, balanced)
Tuning: GridSearchCV with StratifiedKFold(5), scoring='f1'
"""

import argparse
import os
import sys
import warnings

import numpy as np

# Suppress convergence warnings during grid search
warnings.filterwarnings('ignore', category=UserWarning)


def load_labeled_data(path):
    """Load labeled TSV into list of dicts + labels array."""
    with open(path, 'r') as f:
        lines = f.read().strip().split('\n')
    if len(lines) < 2:
        return [], np.array([])

    headers = [h.strip() for h in lines[0].split('\t')]
    rows = []
    labels = []

    label_idx = headers.index('LABEL') if 'LABEL' in headers else -1
    if label_idx < 0:
        print("ERROR: No LABEL column found", file=sys.stderr)
        sys.exit(1)

    for line in lines[1:]:
        fields = line.split('\t')
        if len(fields) < len(headers):
            continue
        row = {}
        for h, v in zip(headers, fields):
            row[h] = v.strip()
        try:
            label = int(row['LABEL'])
        except (ValueError, KeyError):
            continue
        rows.append(row)
        labels.append(label)

    return rows, np.array(labels, dtype=int)


def main():
    parser = argparse.ArgumentParser(
        description='Train LSBG SVM classifier')
    parser.add_argument('--data', required=True,
                        help='Labeled training data TSV (from make_labels)')
    parser.add_argument('--output', required=True,
                        help='Output model path (.joblib)')
    parser.add_argument('--n-splits', type=int, default=5,
                        help='Number of CV folds')
    parser.add_argument('--scoring', default='f1',
                        help='Scoring metric for GridSearchCV')
    args = parser.parse_args()

    from sklearn.impute import SimpleImputer
    from sklearn.model_selection import GridSearchCV, StratifiedKFold
    from sklearn.metrics import (classification_report, confusion_matrix,
                                 roc_auc_score)
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC
    import joblib

    # Add project root for lsbg import
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
    if project_root not in sys.path:
        sys.path.insert(0, project_root)

    from lsbg.classify import extract_features, FEATURE_NAMES

    # Load data
    print(f"Loading data: {args.data}", file=sys.stderr)
    rows, labels = load_labeled_data(args.data)
    if len(rows) == 0:
        print("ERROR: No valid training data", file=sys.stderr)
        sys.exit(1)

    n_pos = int(labels.sum())
    n_neg = len(labels) - n_pos
    print(f"  {len(rows)} samples: {n_pos} positive, {n_neg} negative",
          file=sys.stderr)

    if n_pos < 2 or n_neg < 2:
        print("ERROR: Need at least 2 samples per class", file=sys.stderr)
        sys.exit(1)

    # Extract features
    print(f"Extracting {len(FEATURE_NAMES)} features...", file=sys.stderr)
    X = extract_features(rows)
    y = labels

    # Report NaN statistics
    nan_counts = np.isnan(X).sum(axis=0)
    for fname, nc in zip(FEATURE_NAMES, nan_counts):
        if nc > 0:
            print(f"  {fname}: {nc}/{len(X)} NaN ({100*nc/len(X):.1f}%)",
                  file=sys.stderr)

    # Build pipeline
    pipeline = Pipeline([
        ('imputer', SimpleImputer(strategy='median')),
        ('scaler', StandardScaler()),
        ('svc', SVC(kernel='rbf', probability=True, class_weight='balanced')),
    ])

    # Hyperparameter grid
    param_grid = {
        'svc__C': [0.1, 1, 10, 100],
        'svc__gamma': ['scale', 'auto', 0.01, 0.1],
    }

    # Adjust n_splits if too few positive samples
    n_splits = min(args.n_splits, n_pos)
    if n_splits < 2:
        n_splits = 2

    print(f"GridSearchCV: {n_splits}-fold, scoring={args.scoring}",
          file=sys.stderr)

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    grid = GridSearchCV(
        pipeline, param_grid, cv=cv, scoring=args.scoring,
        refit=True, n_jobs=-1, verbose=0,
    )
    grid.fit(X, y)

    # Best parameters
    print(f"\nBest parameters:", file=sys.stderr)
    for k, v in grid.best_params_.items():
        print(f"  {k}: {v}", file=sys.stderr)
    print(f"Best CV {args.scoring}: {grid.best_score_:.4f}", file=sys.stderr)

    # Evaluate on full training set (for reporting, not validation)
    best_model = grid.best_estimator_
    y_pred = best_model.predict(X)
    y_proba = best_model.predict_proba(X)[:, 1]

    print(f"\n=== Training Set Evaluation ===", file=sys.stderr)
    print(classification_report(y, y_pred, target_names=['contaminant', 'LSBG'],
                                zero_division=0), file=sys.stderr)

    cm = confusion_matrix(y, y_pred)
    print(f"Confusion matrix:", file=sys.stderr)
    print(f"  TN={cm[0,0]}  FP={cm[0,1]}", file=sys.stderr)
    print(f"  FN={cm[1,0]}  TP={cm[1,1]}", file=sys.stderr)

    try:
        auc = roc_auc_score(y, y_proba)
        print(f"ROC AUC: {auc:.4f}", file=sys.stderr)
    except ValueError:
        pass

    # Cross-validation detailed results
    print(f"\n=== Cross-Validation Results ===", file=sys.stderr)
    cv_results = grid.cv_results_
    best_idx = grid.best_index_
    mean_score = cv_results['mean_test_score'][best_idx]
    std_score = cv_results['std_test_score'][best_idx]
    print(f"  CV {args.scoring}: {mean_score:.4f} +/- {std_score:.4f}",
          file=sys.stderr)

    # Save model
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    joblib.dump(best_model, args.output)
    print(f"\nModel saved: {args.output}", file=sys.stderr)

    # Print summary for downstream use
    print(f"\n=== Summary ===", file=sys.stderr)
    print(f"  Features: {len(FEATURE_NAMES)}", file=sys.stderr)
    print(f"  Samples: {len(rows)} ({n_pos}+/{n_neg}-)", file=sys.stderr)
    print(f"  Best C={grid.best_params_['svc__C']}, "
          f"gamma={grid.best_params_['svc__gamma']}", file=sys.stderr)
    print(f"  CV F1={mean_score:.4f}", file=sys.stderr)


if __name__ == '__main__':
    main()
