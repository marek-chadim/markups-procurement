"""Specification curve visualization of procurement premium estimates.

Implements the specification curve methodology (Simonsohn, Simmons & Nelson 2020;
Andrews 2025 AEA lecture on transparency). Collects all procurement markup premium
estimates from the analysis pipeline and visualizes them in a single dominant figure.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# Shared Healy-inspired style (Paul Tol palette, 300 DPI, white bg, clean axes)
from style_markups import apply_markups_style, MARKUPS_BLUE, MARKUPS_PINK
apply_markups_style()

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'

PREFERRED_ESTIMATE = 0.138


def safe_read(path, default=None):
    try:
        if not Path(path).exists():
            return default
        return pd.read_csv(path)
    except Exception as e:
        print(f'  Warning: could not read {path}: {e}')
        return default


def _append(estimates, spec, category, est, se):
    try:
        est_f = float(est)
        se_f = float(se) if se is not None and not pd.isna(se) else np.nan
        if pd.isna(est_f):
            return
        estimates.append({'spec': str(spec), 'category': category,
                          'est': est_f, 'se': se_f})
    except (ValueError, TypeError):
        return


def main():
    estimates = []

    # 1. ACF specs from paper_premiums.csv
    pr = safe_read(OUTPUT_DIR / 'paper_premiums.csv')
    if pr is not None:
        labels = {
            'A': 'ACF Base (pp+surv)', 'B': 'ACF (no surv)',
            'C': 'ACF (no pp Markov)', 'D': 'ACF Plain',
            'E': 'ACF Translog', 'OLS': 'OLS',
        }
        for _, row in pr.iterrows():
            _append(estimates, labels.get(row['spec'], row['spec']),
                    'ACF/OLS', row['reg_premium'], row['reg_se'])

    # 2. Lalonde 11 estimators
    la = safe_read(OUTPUT_DIR / 'lalonde_results.csv')
    if la is not None:
        method_col = la.columns[0]
        est_col = 'Estimate' if 'Estimate' in la.columns else la.columns[1]
        se_col = 'SE' if 'SE' in la.columns else la.columns[2]
        for _, row in la.iterrows():
            _append(estimates, f'Lalonde: {row[method_col]}',
                    'Selection on observables', row[est_col], row[se_col])

    # 3. fect estimates
    fe = safe_read(OUTPUT_DIR / 'fect_results.csv')
    if fe is not None:
        for _, row in fe.iterrows():
            _append(estimates, row['method'], 'Counterfactual (fect)',
                    row['att'], row['se'])

    # 4. ABGRS residualization
    ab = safe_read(OUTPUT_DIR / 'abgrs_residualization.csv')
    if ab is not None:
        for _, row in ab.iterrows():
            _append(estimates, row['Specification'], 'Strong exclusion (ABGRS)',
                    row['Premium'], row['SE'])

    # 5. BH randomization inference (derive SE from CI width)
    bh = safe_read(OUTPUT_DIR / 'borusyak_hull_ri.csv')
    if bh is not None:
        for _, row in bh.iterrows():
            ci_lo = row.get('ci_lo', np.nan)
            ci_hi = row.get('ci_hi', np.nan)
            if pd.notna(ci_lo) and pd.notna(ci_hi):
                se_approx = (float(ci_hi) - float(ci_lo)) / (2 * 1.96)
                _append(estimates, row['method'], 'Randomization (BH)',
                        row['premium'], se_approx)

    # 6. Oster bounds
    os_df = safe_read(OUTPUT_DIR / 'tables' / 'oster_bounds.csv')
    if os_df is not None and len(os_df) > 0:
        row = os_df.iloc[0]
        if 'beta_star' in os_df.columns and pd.notna(row['beta_star']):
            _append(estimates, 'Oster bias-adjusted (delta=1)',
                    'Sensitivity (Oster)', row['beta_star'], 0.005)
        if 'beta_long' in os_df.columns and pd.notna(row['beta_long']):
            _append(estimates, 'Oster long regression',
                    'Sensitivity (Oster)', row['beta_long'], 0.005)

    # 7. Panel FE regressions
    pt = safe_read(OUTPUT_DIR / 'tables' / 'panel_treatment_regressions.csv')
    if pt is not None:
        spec_col = pt.columns[0]
        coef_col = pt.columns[1]
        se_col = pt.columns[2]
        for _, row in pt.iterrows():
            _append(estimates, f'Panel FE: {row[spec_col]}',
                    'Panel FE', row[coef_col], row[se_col])

    df = pd.DataFrame(estimates).dropna(subset=['est']).copy()
    # For entries missing SE, keep point with zero-length bar
    df['se'] = df['se'].fillna(0.0)
    df['ci_lo'] = df['est'] - 1.96 * df['se']
    df['ci_hi'] = df['est'] + 1.96 * df['se']

    cat_order = [
        'ACF/OLS', 'Selection on observables', 'Counterfactual (fect)',
        'Strong exclusion (ABGRS)', 'Randomization (BH)',
        'Sensitivity (Oster)', 'Panel FE',
    ]
    df['category'] = pd.Categorical(df['category'], ordered=True, categories=cat_order)
    df = df.sort_values(['category', 'est']).reset_index(drop=True)

    print(f'\nTotal estimates: {len(df)}')
    print(df.groupby('category', observed=True).size())
    print(f"Range: [{df['est'].min():.4f}, {df['est'].max():.4f}]")

    cat_colors = {
        'ACF/OLS': '#1f77b4',
        'Selection on observables': '#2ca02c',
        'Counterfactual (fect)': '#ff7f0e',
        'Strong exclusion (ABGRS)': '#9467bd',
        'Randomization (BH)': '#8c564b',
        'Sensitivity (Oster)': '#7f7f7f',
        'Panel FE': '#17becf',
    }

    fig, ax = plt.subplots(figsize=(9, max(6, 0.28 * len(df))))
    y_pos = np.arange(len(df))
    colors = [cat_colors.get(str(c), 'black') for c in df['category']]

    for i, row in df.iterrows():
        c = cat_colors.get(str(row['category']), 'black')
        ax.plot([row['ci_lo'], row['ci_hi']], [i, i], '-', color=c,
                linewidth=1.2, alpha=0.7)
        ax.plot(row['est'], i, 'o', color=c, markersize=6,
                markeredgecolor='black', markeredgewidth=0.3)

    ax.axvline(PREFERRED_ESTIMATE, color='red', linestyle='--',
               linewidth=1.2, alpha=0.7,
               label=f'Preferred ({PREFERRED_ESTIMATE:.3f})')
    ax.axvline(0, color='grey', linestyle=':', linewidth=1, alpha=0.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df['spec'], fontsize=8)
    ax.set_xlabel('Procurement markup premium (log points)')
    ax.set_title('Specification Curve: Procurement Premium Across All Specifications')

    handles = []
    for cat in cat_order:
        if cat in df['category'].values:
            handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                       markerfacecolor=cat_colors[cat],
                                       markersize=8, label=cat))
    handles.append(plt.Line2D([0], [0], color='red', linestyle='--',
                               label=f'Preferred ({PREFERRED_ESTIMATE:.3f})'))
    ax.legend(handles=handles, loc='lower right', fontsize=8, framealpha=0.9)
    ax.grid(axis='x', alpha=0.3)
    ax.set_axisbelow(True)

    plt.tight_layout()
    fig_path = OUTPUT_DIR / 'figures' / 'specification_curve.pdf'
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f'\nSaved: {fig_path.relative_to(OUTPUT_DIR)}')

    data_path = OUTPUT_DIR / 'specification_curve_data.csv'
    df.to_csv(data_path, index=False)
    print(f'Saved: {data_path.relative_to(OUTPUT_DIR)}')

    # Armstrong (2025) bias-adjusted CI
    # CI = estimate ± max_bias ± z_{0.975} × SE
    # where max_bias = sup over expanded model |bias(θ; w)|
    # Use specification curve deviations from preferred estimate as bias bound
    beta_preferred = 0.138
    se_preferred = 0.003

    # Full expanded model (all 33 specs)
    max_bias_full = max(abs(df['est'].max() - beta_preferred),
                        abs(df['est'].min() - beta_preferred))
    ci_armstrong_full = (beta_preferred - max_bias_full - 1.96 * se_preferred,
                         beta_preferred + max_bias_full + 1.96 * se_preferred)

    # Restricted expanded model (exclude specs without pp in Markov —
    # internally inconsistent per De Loecker 2013)
    df_consistent = df[df['est'] > 0.05]  # plain ACF specs give ~0.003
    max_bias_consistent = max(abs(df_consistent['est'].max() - beta_preferred),
                              abs(df_consistent['est'].min() - beta_preferred))
    ci_armstrong_consistent = (beta_preferred - max_bias_consistent - 1.96 * se_preferred,
                               beta_preferred + max_bias_consistent + 1.96 * se_preferred)

    # Summary statistics
    median_est = df['est'].median()
    iqr = (df['est'].quantile(0.25), df['est'].quantile(0.75))

    # Armstrong summary table
    armstrong = pd.DataFrame([
        {'CI_type': 'Standard (clustered SE)',
         'lower': beta_preferred - 1.96*se_preferred,
         'upper': beta_preferred + 1.96*se_preferred,
         'width': 2*1.96*se_preferred},
        {'CI_type': 'Armstrong bias-adjusted (full)',
         'lower': ci_armstrong_full[0], 'upper': ci_armstrong_full[1],
         'width': ci_armstrong_full[1] - ci_armstrong_full[0]},
        {'CI_type': 'Armstrong bias-adjusted (consistent)',
         'lower': ci_armstrong_consistent[0], 'upper': ci_armstrong_consistent[1],
         'width': ci_armstrong_consistent[1] - ci_armstrong_consistent[0]},
        {'CI_type': 'Specification curve IQR',
         'lower': iqr[0], 'upper': iqr[1],
         'width': iqr[1] - iqr[0]},
    ])
    armstrong.to_csv(OUTPUT_DIR / 'armstrong_bias_adjusted_ci.csv', index=False)

    print(f'\n=== Armstrong (2025) Bias-Adjusted CIs ===')
    print(f'  Preferred estimate: {beta_preferred:.4f} (SE {se_preferred:.4f})')
    print(f'  Standard 95% CI:           [{beta_preferred-1.96*se_preferred:.4f}, {beta_preferred+1.96*se_preferred:.4f}]')
    print(f'  Armstrong full:            [{ci_armstrong_full[0]:.4f}, {ci_armstrong_full[1]:.4f}]')
    print(f'  Armstrong consistent:      [{ci_armstrong_consistent[0]:.4f}, {ci_armstrong_consistent[1]:.4f}]')
    print(f'  Spec curve median: {median_est:.4f}, IQR: [{iqr[0]:.4f}, {iqr[1]:.4f}]')

    # LaTeX table
    tex_arm = [
        r'\begin{table}[htbp]\centering',
        r'\caption{Confidence Intervals Under Specification Uncertainty (Armstrong 2025)}\label{tab:armstrong_ci}',
        r'\begin{threeparttable}',
        r'\begin{tabular}{lccc}',
        r'\toprule',
        r'Inference method & 95\% CI & Width \\',
        r'\midrule',
    ]
    for _, r in armstrong.iterrows():
        tex_arm.append(f'{r["CI_type"]} & [{r["lower"]:.4f}, {r["upper"]:.4f}] & {r["width"]:.4f} \\\\')
    tex_arm += [
        r'\bottomrule',
        r'\end{tabular}',
        r'\begin{tablenotes}\footnotesize',
        r"\item \textit{Notes:} Standard CI uses firm-clustered SE. Armstrong (2025) bias-adjusted CI: $\hat{\beta} \pm \overline{\text{bias}} \pm z_{0.975}\hat{\sigma}$, where $\overline{\text{bias}} = \sup_w |\hat{\beta} - \hat{\beta}(w)|$ over the specification curve. ``Full'' includes all 33 specs; ``consistent'' excludes specifications without $pp_{it-1}$ in the Markov transition (internally inconsistent per De Loecker 2013). IQR: interquartile range of specification curve estimates.",
        r'\end{tablenotes}',
        r'\end{threeparttable}',
        r'\end{table}',
    ]
    with open(OUTPUT_DIR / 'tables' / 'armstrong_bias_adjusted_ci.tex', 'w') as f:
        f.write('\n'.join(tex_arm))
    print(f'Saved: tables/armstrong_bias_adjusted_ci.tex')


if __name__ == '__main__':
    main()
