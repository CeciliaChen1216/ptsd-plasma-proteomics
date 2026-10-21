"""Read-only coverage and factor/MDD checks against supplied result snapshots.

Run with Python, pandas and scipy installed. No file or directory is written.
"""
from pathlib import Path
import math
import pandas as pd


def main():
    package = Path(__file__).resolve().parents[1]
    snapshots = package / 'results_snapshot'
    sem = snapshots / 'sem_expansion_20261019'
    components = snapshots / 'current_manuscript_components'
    pqtl = pd.read_csv(sem / 'pQTL_harmonization_audit.csv')
    raw_audit = pd.read_csv(sem / 'PTSD_filter_audit.csv')
    primary = pd.read_csv(components / 'TableS3a.csv')
    assert len(pqtl) == 1860 and len(primary) == 1258
    coverage = pqtl.merge(
        primary[['protein', 'SNP']].assign(in_primary_screen=True),
        on=['protein', 'SNP'], how='left', validate='one_to_one',
    ).merge(
        raw_audit[['SNP', 'raw_present', 'reason']],
        on='SNP', how='left', validate='many_to_one',
    )
    coverage['in_primary_screen'] = coverage['in_primary_screen'].eq(True)
    excluded = coverage[~coverage['in_primary_screen']]
    assert len(excluded) == 602
    assert excluded['raw_present'].eq(False).sum() == 571
    assert excluded['raw_present'].eq(True).sum() == 31

    factor = pd.read_csv(sem / 'primary' / 'factor_MR_expanded.csv')
    mdd = pd.read_csv(components / 'TableS3b.csv')
    comparison = factor.merge(
        mdd, on=['protein', 'SNP'], suffixes=('_factor', '_MDD'),
        validate='one_to_one',
    )
    assert len(comparison) == 1162
    factor_sig = comparison['mr_fdr'] < 0.05
    mdd_sig = comparison['fdr'] < 0.05
    summary = {
        'matched_proteins': len(comparison),
        'Pearson_reported_Wald_effects': comparison['mr_beta_factor'].corr(comparison['mr_beta_MDD']),
        'Spearman_reported_Wald_effects': comparison['mr_beta_factor'].corr(comparison['mr_beta_MDD'], method='spearman'),
        'same_effect_sign': int((comparison['mr_beta_factor'] * comparison['mr_beta_MDD'] > 0).sum()),
        'both_full_screen_FDR_significant': int((factor_sig & mdd_sig).sum()),
        'factor_only_full_screen_FDR': int((factor_sig & ~mdd_sig).sum()),
        'MDD_only_full_screen_FDR': int((~factor_sig & mdd_sig).sum()),
    }
    for metric, value in {
        'same_effect_sign': 1010,
        'both_full_screen_FDR_significant': 29,
        'factor_only_full_screen_FDR': 22,
        'MDD_only_full_screen_FDR': 6,
    }.items():
        assert summary[metric] == value, (metric, summary[metric], value)
    recorded = pd.read_csv(
        snapshots / 'reviewer_audit_20261020' / 'factor_MDD_summary.csv'
    ).set_index('metric')['value']
    for metric, value in summary.items():
        assert math.isclose(float(recorded[metric]), float(value), rel_tol=1e-12, abs_tol=1e-12), (metric, recorded[metric], value)
    print('PASS: 602 excluded proteins; 571 absent and 31 present by raw rsID lookup.')
    print('PASS: 1,162 factor/MDD comparisons; sign/FDR counts 1010/29/22/6.')
    print('PASS: supplied summary values agree; no snapshot files were written.')


if __name__ == '__main__':
    main()
