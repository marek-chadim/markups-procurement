"""
DLEU Replication: The Rise of Market Power — Czech Construction.

Full implementation of De Loecker, Eeckhout, and Unger (2020, QJE)
"The Rise of Market Power and the Macroeconomic Implications" methodology,
adapted to Czech construction firms (CZ-NACE F, 2006-2021).

Replicates (adapted):
  - Figure I:   Aggregate markup time series
  - Figure II:  Sensitivity to θ and weighting scheme
  - Figure III: Distribution (kernel density + percentiles)
  - Figure IV:  Cumulative decomposition with counterfactuals (eq. 9)
  - Figure V:   Micro aggregation vs industry averages
  - Figure VII: Cost share evolution
  - Figure VIII: Profit rate and distribution
  - Figure XII: Cost-share-based markup and output elasticity
  - Table I:    Sectoral decomposition (eq. 10) by NACE 41/42/43
  - Appendix 12: Industry-specific trends

Extensions unique to this paper:
  - Procurement vs non-procurement decomposition
  - Czech transparency reform markers (2012, 2016)
  - Markup-procurement premium by percentile

References:
  De Loecker, Eeckhout & Unger (2020). QJE 135(2), 561-644.
  Haltiwanger (1997). Measuring and interpreting productivity growth.
  De Loecker, Eeckhout & Unger (2025). Reply to BMY. QJE.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Shared Healy-inspired style (Paul Tol palette, 300 DPI, white bg, clean axes)
from style_markups import apply_markups_style, MARKUPS_BLUE, MARKUPS_PINK
apply_markups_style()
from matplotlib.ticker import MaxNLocator
from scipy.stats import gaussian_kde
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR.parent / 'input'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'
FIG_DIR = OUTPUT_DIR / 'figures'
TAB_DIR = OUTPUT_DIR / 'tables'
DAT_DIR = OUTPUT_DIR / 'data'

USER_COST_RATE = 0.12  # DLEU: r = (I - Π) + Δ ≈ 12%
THETA_CALIBRATED = 0.85  # DLEU benchmark calibration

NACE_LABELS = {41: 'Buildings', 42: 'Civil Eng.', 43: 'Specialized'}

# Consistent style
plt.rcParams.update({
    'figure.figsize': (8, 5),
    'axes.grid': True,
    'grid.alpha': 0.3,
    'font.size': 11,
    'axes.labelsize': 12,
    'legend.fontsize': 9,
})

REFORM_LINES = [
    (2012, 'Act 55/2012\n(single-bid ban)', 'red'),
    (2016, 'Act 134/2016\n(MEAT criteria)', 'green'),
]


# ---------------------------------------------------------------------------
# Data loading and markup construction
# ---------------------------------------------------------------------------

def load_data():
    """Load panel data and merge markups."""
    df = pd.read_stata(str(INPUT_DIR / 'data.dta'))
    mk = pd.read_stata(str(DAT_DIR / 'paper_markups.dta'))

    for d in [df, mk]:
        if hasattr(d['year'].iloc[0], 'year'):
            d['year'] = d['year'].dt.year
        d['year'] = d['year'].astype(int)

    panel = df.merge(mk[['id', 'year', 'markup_A', 'markup_B', 'alphahat']],
                     on=['id', 'year'], how='inner')
    panel = panel[panel['markup_A'].notna() & (panel['markup_A'] > 0)].copy()

    # Recover levels from logs
    panel['Sales'] = np.exp(panel['go'])
    panel['COGS'] = np.exp(panel['cogs'])
    panel['K'] = np.exp(panel['k'])

    # Capital expenditure (user cost × capital stock) — DLEU convention
    panel['kexp'] = USER_COST_RATE * panel['K']

    # Employment midpoint (if available)
    if 'empl_mid' in panel.columns:
        panel['Emp'] = panel['empl_mid']
    elif 'le' in panel.columns:
        panel['Emp'] = np.exp(panel['le'])
    else:
        panel['Emp'] = np.nan

    # Wages
    if 'rW' in panel.columns:
        panel['Wages'] = panel['rW']
    elif 'w' in panel.columns:
        panel['Wages'] = np.exp(panel['w'])
    else:
        panel['Wages'] = np.nan

    # Procurement dummy
    if 'pp_dummy' not in panel.columns:
        panel['pp_dummy'] = 0

    return panel


def construct_markups(panel):
    """
    Construct multiple markup measures following DLEU.

    mu_cs085:  calibrated θ=0.85 (DLEU mu_0)
    mu_cs1:    firm-level cost share COGS/(COGS+rK) (DLEU mu_1)
    mu_cs_ind: industry-median cost share by nace2×year (DLEU mu_3)
    mu_acf:    ACF-estimated (markup_A from estimator) (DLEU mu_10)
    """
    # Cost shares
    panel['cs_firm'] = panel['COGS'] / (panel['COGS'] + panel['kexp'])
    panel['cs_ind'] = panel.groupby(['nace2', 'year'])['cs_firm'].transform('median')

    # Markup measures: θ × (Sales / COGS)
    ratio = panel['Sales'] / panel['COGS']
    panel['mu_cs085'] = THETA_CALIBRATED * ratio
    panel['mu_cs1'] = panel['cs_firm'] * ratio
    panel['mu_cs_ind'] = panel['cs_ind'] * ratio
    panel['mu_acf'] = panel['markup_A']

    # Total costs
    panel['totcost'] = panel['COGS'] + panel['kexp']

    # Profit rate: π = (Sales - COGS - rK) / Sales
    panel['profit_rate'] = (panel['Sales'] - panel['COGS'] - panel['kexp']) / panel['Sales']

    # Model-based profit rate: π = 1 - θ/μ - rK/S (DLEU eq. 13, simplified)
    theta_acf = panel['mu_acf'] * panel['alphahat']  # recover θ = μ × α
    panel['profit_model'] = 1 - theta_acf / panel['mu_acf'] - panel['kexp'] / panel['Sales']
    # Simplifies to: 1 - alphahat - kexp/Sales

    # Labor share (if wages available)
    if panel['Wages'].notna().any():
        panel['labor_share'] = panel['Wages'] / panel['Sales']
    panel['capital_share'] = panel['kexp'] / panel['Sales']
    panel['cogs_share'] = panel['COGS'] / panel['Sales']

    return panel


# ---------------------------------------------------------------------------
# Aggregation functions
# ---------------------------------------------------------------------------

def aggregate_markups(panel, mu_col='mu_acf'):
    """Compute aggregate markups with multiple weighting schemes."""
    results = []
    for year, grp in panel.groupby('year'):
        n = len(grp)
        total_sales = grp['Sales'].sum()
        total_cogs = grp['COGS'].sum()
        total_emp = grp['Emp'].sum() if grp['Emp'].notna().all() else np.nan
        total_cost = grp['totcost'].sum()

        sw = grp['Sales'] / total_sales  # sales weights
        cw = grp['COGS'] / total_cogs    # COGS weights
        tcw = grp['totcost'] / total_cost  # total cost weights

        mu = grp[mu_col]

        m_sw = (sw * mu).sum()
        m_cw = (cw * mu).sum()
        m_tcw = (tcw * mu).sum()
        m_mean = mu.mean()
        m_harmonic = 1.0 / (sw / mu).sum()

        # Employment-weighted
        if pd.notna(total_emp) and total_emp > 0:
            ew = grp['Emp'] / total_emp
            m_ew = (ew * mu).sum()
        else:
            m_ew = np.nan

        # Sales-weighted percentiles
        sorted_idx = mu.sort_values().index
        sw_sorted = sw.loc[sorted_idx]
        mu_sorted = mu.loc[sorted_idx]
        cum_w = sw_sorted.cumsum()

        pcts = {}
        for p in [10, 25, 50, 75, 90]:
            idx = (cum_w >= p / 100).idxmax()
            pcts[f'P{p}'] = mu_sorted.loc[idx]

        results.append({
            'year': year, 'N': n,
            'M_sw': m_sw, 'M_cogs': m_cw, 'M_totcost': m_tcw,
            'M_mean': m_mean, 'M_harmonic': m_harmonic, 'M_emp': m_ew,
            'total_sales': total_sales,
            **pcts,
        })

    return pd.DataFrame(results).sort_values('year').reset_index(drop=True)


# ---------------------------------------------------------------------------
# Decomposition: DLEU equation (9) — Haltiwanger (1997)
# ---------------------------------------------------------------------------

def decompose_dleu(panel, mu_col='mu_acf'):
    """
    Firm-level decomposition following DLEU eq. (9) with Haltiwanger demeaning.

    ΔM_t = Σ m_{i,t-1} Δμ_it                          [within]
          + Σ (μ_{i,t-1} - M_{t-1}) Δm_it              [market share]
          + Σ Δμ_it Δm_it                               [cross]
          + Σ_{entry} (μ_it - M_{t-1}) m_it
            - Σ_{exit} (μ_{i,t-1} - M_{t-1}) m_{i,t-1}  [net entry]
    """
    # Pre-compute aggregate and shares by year
    yearly = {}
    for year, grp in panel.groupby('year'):
        grp = grp.set_index('id')
        total_s = grp['Sales'].sum()
        grp = grp.copy()
        grp['share'] = grp['Sales'] / total_s
        grp['mu'] = grp[mu_col]
        M = (grp['share'] * grp['mu']).sum()
        yearly[year] = (grp[['share', 'mu']], M)

    years = sorted(yearly.keys())
    decomp = []

    for i in range(1, len(years)):
        t = years[i]
        t1 = years[i - 1]
        curr_df, M_t = yearly[t]
        prev_df, M_t1 = yearly[t1]

        incumbents = curr_df.index.intersection(prev_df.index)
        entrants = curr_df.index.difference(prev_df.index)
        exiters = prev_df.index.difference(curr_df.index)

        if len(incumbents) == 0:
            continue

        c = curr_df.loc[incumbents]
        p = prev_df.loc[incumbents]

        d_mu = c['mu'] - p['mu']
        d_share = c['share'] - p['share']
        mu_tilde_prev = p['mu'] - M_t1  # Haltiwanger demeaning

        within = (p['share'] * d_mu).sum()
        mkt_share = (mu_tilde_prev * d_share).sum()
        cross = (d_mu * d_share).sum()

        # Net entry (demeaned by M_{t-1})
        if len(entrants) > 0:
            e = curr_df.loc[entrants]
            entry = ((e['mu'] - M_t1) * e['share']).sum()
        else:
            entry = 0.0

        if len(exiters) > 0:
            x = prev_df.loc[exiters]
            exit_ = ((x['mu'] - M_t1) * x['share']).sum()
        else:
            exit_ = 0.0

        net_entry = entry - exit_
        reallocation = mkt_share + cross
        total = M_t - M_t1

        decomp.append({
            'year': t, 'dM': total,
            'within': within, 'mkt_share': mkt_share,
            'cross': cross, 'reallocation': reallocation,
            'entry': entry, 'exit': exit_, 'net_entry': net_entry,
            'N_inc': len(incumbents), 'N_ent': len(entrants), 'N_exit': len(exiters),
        })

    return pd.DataFrame(decomp)


def decompose_sectoral(panel, mu_col='mu_acf', period_length=5):
    """
    Sectoral decomposition following DLEU eq. (10).

    ΔM_t = Σ_s m_{s,t-1} Δμ_st  [within]
          + Σ_s μ_{s,t-1} Δm_st  [between]
          + Σ_s Δμ_st Δm_st      [cross]

    where s indexes NACE 2-digit sub-industries.
    """
    # Compute sector-level aggregates
    sector_yr = []
    for (nace, year), grp in panel.groupby(['nace2', 'year']):
        total_sales = grp['Sales'].sum()
        sw = grp['Sales'] / total_sales
        mu_s = (sw * grp[mu_col]).sum()
        sector_yr.append({
            'nace2': nace, 'year': year,
            'mu_s': mu_s, 'sales_s': total_sales,
        })
    sdf = pd.DataFrame(sector_yr)

    # Add sector shares in total economy
    total_by_year = sdf.groupby('year')['sales_s'].sum().rename('total_sales')
    sdf = sdf.merge(total_by_year, on='year')
    sdf['share_s'] = sdf['sales_s'] / sdf['total_sales']

    # Economy-wide aggregate
    agg_by_year = sdf.groupby('year').apply(
        lambda g: (g['share_s'] * g['mu_s']).sum()).rename('M_agg')
    sdf = sdf.merge(agg_by_year, on='year')

    years = sorted(sdf['year'].unique())

    # Use period_length-year changes
    results = []
    for t in years:
        t_prev = t - period_length
        if t_prev not in years:
            continue

        curr = sdf[sdf['year'] == t].set_index('nace2')
        prev = sdf[sdf['year'] == t_prev].set_index('nace2')
        common = curr.index.intersection(prev.index)
        if len(common) == 0:
            continue

        c = curr.loc[common]
        p = prev.loc[common]

        d_mu = c['mu_s'] - p['mu_s']
        d_share = c['share_s'] - p['share_s']

        within = (p['share_s'] * d_mu).sum()
        between = (p['mu_s'] * d_share).sum()
        cross_term = (d_mu * d_share).sum()
        d_markup = c['M_agg'].iloc[0] - p['M_agg'].iloc[0]

        results.append({
            'year': t, 'period': f'{t_prev}-{t}',
            'M_agg': c['M_agg'].iloc[0], 'dM': d_markup,
            'within': within, 'between': between, 'cross': cross_term,
        })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def add_reform_lines(ax, ymin=None, ymax=None):
    """Add Czech reform vertical lines."""
    for yr, label, color in REFORM_LINES:
        ax.axvline(yr, color=color, linestyle='--', alpha=0.4, linewidth=0.8)


def fig1_aggregate_markup(agg, agg_specs):
    """Figure I: Aggregate markup time series (benchmark = ACF)."""
    fig, ax = plt.subplots()
    ax.plot(agg['year'], agg['M_sw'], 'o-', color='C0', linewidth=2,
            markersize=4, label='ACF (benchmark)')
    ax.plot(agg_specs['cs085']['year'], agg_specs['cs085']['M_sw'],
            's--', color='C1', linewidth=1.5, markersize=3,
            alpha=0.7, label=r'Calibrated $\theta=0.85$')
    ax.plot(agg_specs['cs1']['year'], agg_specs['cs1']['M_sw'],
            '^:', color='C2', linewidth=1.5, markersize=3,
            alpha=0.7, label='Firm cost share')
    ax.plot(agg_specs['cs_ind']['year'], agg_specs['cs_ind']['M_sw'],
            'v-.', color='C3', linewidth=1.5, markersize=3,
            alpha=0.7, label='Industry cost share')
    add_reform_lines(ax)
    ax.set_xlabel('Year')
    ax.set_ylabel('Sales-weighted aggregate markup')
    ax.set_title('Aggregate Markups — Czech Construction (CZ-NACE F)')
    ax.legend(loc='best')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'dleu_fig1_aggregate_markup.pdf', dpi=300)
    plt.close()


def fig2_sensitivity(agg, panel):
    """Figure II: Sensitivity to θ source (Panel A) and weighting (Panel B)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: benchmark vs constant θ
    agg_cal = aggregate_markups(panel, 'mu_cs085')
    ax1.plot(agg['year'], agg['M_sw'], '-', color='red', linewidth=2,
             label='Benchmark (ACF)')
    ax1.plot(agg_cal['year'], agg_cal['M_sw'], '--', color='green',
             linewidth=2, label=r'Constant $\theta = 0.85$')
    add_reform_lines(ax1)
    ax1.set_xlabel('Year')
    ax1.set_ylabel('Aggregate markup')
    ax1.set_title('(A) Constant elasticity')
    ax1.legend(loc='best')
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

    # Panel B: sales-weighted vs input-weighted
    ax2.plot(agg['year'], agg['M_sw'], '-', color='red', linewidth=2,
             label='Sales-weighted (benchmark)')
    ax2.plot(agg['year'], agg['M_cogs'], '--', color='green', linewidth=2,
             label='Input-weighted (COGS)')
    ax2.plot(agg['year'], agg['M_totcost'], ':', color='blue', linewidth=2,
             label='Input-weighted (total cost)')
    add_reform_lines(ax2)
    ax2.set_xlabel('Year')
    ax2.set_title('(B) Input-weighted')
    ax2.legend(loc='best')
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.tight_layout()
    plt.savefig(FIG_DIR / 'dleu_fig2_sensitivity.pdf', dpi=300)
    plt.close()


def fig3_distribution(panel, agg):
    """Figure III: Markup distribution — density + sales-weighted percentiles."""
    years = sorted(panel['year'].unique())
    y_first, y_last = years[0], years[-1]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: kernel density (unweighted)
    for yr, ls, lbl in [(y_first, '--', str(y_first)), (y_last, '-', str(y_last))]:
        mu = panel.loc[panel['year'] == yr, 'mu_acf'].dropna()
        mu = mu[(mu > 0.3) & (mu < 5)]  # trim for plotting
        kde = gaussian_kde(mu, bw_method=0.3)
        x = np.linspace(0.3, 5, 300)
        ax1.plot(x, kde(x), ls, linewidth=2, label=lbl)

    ax1.set_xlabel(r'Markup $\mu_{it}$')
    ax1.set_ylabel('Density')
    ax1.set_title('(A) Kernel density (unweighted)')
    ax1.legend()

    # Panel B: sales-weighted percentiles
    for p, ls, lbl in [('P90', '--', 'P90'), ('P75', ':', 'P75'),
                        ('P50', '-.', 'P50')]:
        ax2.plot(agg['year'], agg[p], ls, linewidth=2, label=lbl)
    ax2.plot(agg['year'], agg['M_sw'], '-', linewidth=2,
             color='black', label='Average')
    add_reform_lines(ax2)
    ax2.set_xlabel('Year')
    ax2.set_ylabel(r'Markup $\mu_{it}$')
    ax2.set_title('(B) Percentiles (sales-weighted)')
    ax2.legend(loc='best')
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.tight_layout()
    plt.savefig(FIG_DIR / 'dleu_fig3_distribution.pdf', dpi=300)
    plt.close()


def fig4_decomposition(decomp, agg):
    """
    Figure IV: Cumulative decomposition from base year with counterfactuals.

    Three counterfactual paths starting from base-year markup:
      - Within only:       base + cumΣ(within)
      - Reallocation only: base + cumΣ(mkt_share + cross)
      - Net entry only:    base + cumΣ(net_entry)
    """
    if decomp.empty:
        return

    base_year = agg['year'].iloc[0]
    base_mu = agg.loc[agg['year'] == base_year, 'M_sw'].iloc[0]

    cum_within = decomp['within'].cumsum()
    cum_realloc = decomp['reallocation'].cumsum()
    cum_net_entry = decomp['net_entry'].cumsum()

    fig, ax = plt.subplots()
    ax.plot(agg['year'], agg['M_sw'], '-', color='red', linewidth=2.5,
            label='Markup (benchmark)')
    ax.plot(decomp['year'], base_mu + cum_within, '--', color='blue',
            linewidth=2, label='Within')
    ax.plot(decomp['year'], base_mu + cum_realloc, ':', color='black',
            linewidth=2, label='Reallocation')
    ax.plot(decomp['year'], base_mu + cum_net_entry, '-.', color='green',
            linewidth=2, label='Net entry')
    add_reform_lines(ax)
    ax.set_xlabel('Year')
    ax.set_ylabel('Aggregate markup')
    ax.set_title('Decomposition of Markup Growth at the Firm Level')
    ax.legend(loc='best')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'dleu_fig4_decomposition.pdf', dpi=300)
    plt.close()


def fig5_micro_vs_agg(panel, agg):
    """
    Figure V: Micro aggregation vs industry averages.

    Demonstrates Jensen's inequality: firm-level sales-weighted aggregate >
    industry-average aggregate when markup distribution is skewed.
    """
    # Industry-level: average markup by nace2, then weight by sector sales share
    ind_agg = []
    for year, grp in panel.groupby('year'):
        total_sales = grp['Sales'].sum()
        ind_mu = []
        for nace, sg in grp.groupby('nace2'):
            ind_sales = sg['Sales'].sum()
            ind_share = ind_sales / total_sales
            # Simple within-industry average (unweighted)
            mu_ind_avg = sg['mu_acf'].mean()
            ind_mu.append({'nace2': nace, 'share': ind_share, 'mu': mu_ind_avg})
        idf = pd.DataFrame(ind_mu)
        m_ind = (idf['share'] * idf['mu']).sum()
        # Economy-wide simple average (one industry)
        m_econ = grp['mu_acf'].mean()
        ind_agg.append({
            'year': year, 'M_ind_avg': m_ind, 'M_econ_avg': m_econ,
        })
    ind_df = pd.DataFrame(ind_agg)

    fig, ax = plt.subplots()
    ax.plot(agg['year'], agg['M_sw'], '-', color='red', linewidth=2.5,
            label='Firm-level (sales-weighted)')
    ax.plot(ind_df['year'], ind_df['M_ind_avg'], '--', color='blue',
            linewidth=2, label='Industry averages')
    ax.plot(ind_df['year'], ind_df['M_econ_avg'], ':', color='black',
            linewidth=2, label='Economy-wide average')
    add_reform_lines(ax)
    ax.set_xlabel('Year')
    ax.set_ylabel('Aggregate markup')
    ax.set_title('Micro Aggregation vs Industry Averages')
    ax.legend(loc='best')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'dleu_fig5_micro_vs_agg.pdf', dpi=300)
    plt.close()


def fig7_cost_shares(panel):
    """Figure VII: Cost share evolution over time."""
    shares = panel.groupby('year').agg(
        cs_cogs=('cogs_share', 'mean'),
        cs_capital=('capital_share', 'mean'),
    ).reset_index()

    fig, ax = plt.subplots()
    ax.plot(shares['year'], shares['cs_cogs'], '-', color='blue',
            linewidth=2, label='COGS share (COGS/Sales)')
    ax.plot(shares['year'], shares['cs_capital'], '--', color='orange',
            linewidth=2, label='Capital share (rK/Sales)')
    add_reform_lines(ax)
    ax.set_xlabel('Year')
    ax.set_ylabel('Cost share of sales')
    ax.set_title('Cost Shares of Total Sales')
    ax.legend(loc='best')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'dleu_fig7_cost_shares.pdf', dpi=300)
    plt.close()


def fig8_profit_rate(panel):
    """Figure VIII: Profit rate — time series and distribution."""
    years = sorted(panel['year'].unique())
    y_first, y_last = years[0], years[-1]

    # Sales-weighted profit rate
    pr_agg = []
    for year, grp in panel.groupby('year'):
        sw = grp['Sales'] / grp['Sales'].sum()
        pr = (sw * grp['profit_rate']).sum()
        pr_agg.append({'year': year, 'profit_rate': pr})
    pr_df = pd.DataFrame(pr_agg)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: time series
    ax1.plot(pr_df['year'], pr_df['profit_rate'], 'o-', color='red',
             linewidth=2, markersize=4)
    add_reform_lines(ax1)
    ax1.set_xlabel('Year')
    ax1.set_ylabel('Profit rate')
    ax1.set_title('(A) Average profit rate (sales-weighted)')
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

    # Panel B: kernel density
    for yr, ls, lbl in [(y_first, '--', str(y_first)), (y_last, '-', str(y_last))]:
        pr = panel.loc[panel['year'] == yr, 'profit_rate'].dropna()
        pr = pr[(pr > -0.5) & (pr < 0.8)]
        if len(pr) > 10:
            kde = gaussian_kde(pr, bw_method=0.3)
            x = np.linspace(-0.5, 0.8, 300)
            ax2.plot(x, kde(x), ls, linewidth=2, label=lbl)
    ax2.set_xlabel('Profit rate')
    ax2.set_ylabel('Density')
    ax2.set_title('(B) Profit rate distribution')
    ax2.legend()

    plt.tight_layout()
    plt.savefig(FIG_DIR / 'dleu_fig8_profit_rate.pdf', dpi=300)
    plt.close()

    return pr_df


def fig12_cost_share_markup(panel, agg):
    """Figure XII: Cost-share-based markup and output elasticity comparison."""
    # Compute cost-share-based aggregate markup (CS approach)
    cs_agg = aggregate_markups(panel, 'mu_cs1')

    # Sector-weighted average cost share vs ACF θ
    theta_cs = panel.groupby('year')['cs_firm'].median().reset_index()
    theta_cs.columns = ['year', 'cs_median']

    # Recover ACF θ from markup_A * alphahat
    theta_acf = panel.copy()
    theta_acf['theta'] = theta_acf['mu_acf'] * theta_acf['alphahat']
    theta_yr = []
    for year, grp in theta_acf.groupby('year'):
        sw = grp['Sales'] / grp['Sales'].sum()
        theta_yr.append({'year': year, 'theta_acf': (sw * grp['theta']).sum()})
    theta_df = pd.DataFrame(theta_yr)
    theta_df = theta_df.merge(theta_cs, on='year')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: CS-based aggregate markup
    ax1.plot(cs_agg['year'], cs_agg['M_sw'], '-', color='red', linewidth=2,
             label='Cost-share markup')
    ax1.plot(agg['year'], agg['M_sw'], '--', color='blue', linewidth=2,
             label='ACF markup')
    add_reform_lines(ax1)
    ax1.set_xlabel('Year')
    ax1.set_ylabel('Aggregate markup')
    ax1.set_title('(A) Aggregate markup: CS vs ACF')
    ax1.legend()
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

    # Panel B: θ comparison
    ax2.plot(theta_df['year'], theta_df['theta_acf'], '-', color='red',
             linewidth=2, label=r'ACF $\hat{\theta}$')
    ax2.plot(theta_df['year'], theta_df['cs_median'], '--', color='green',
             linewidth=2, label='Cost share (median)')
    add_reform_lines(ax2)
    ax2.set_xlabel('Year')
    ax2.set_ylabel(r'Output elasticity $\theta^V$')
    ax2.set_title(r'(B) Output elasticity $\theta^V$: ACF vs cost share')
    ax2.legend()
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.tight_layout()
    plt.savefig(FIG_DIR / 'dleu_fig12_cs_vs_pf.pdf', dpi=300)
    plt.close()


def fig_industry_trends(panel):
    """Appendix 12: Industry-specific markup trends by NACE 41/42/43."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)

    for idx, (nace, label) in enumerate(NACE_LABELS.items()):
        ax = axes[idx]
        sub = panel[panel['nace2'] == nace]
        if sub.empty:
            continue

        agg_n = []
        for year, grp in sub.groupby('year'):
            sw = grp['Sales'] / grp['Sales'].sum()
            agg_n.append({
                'year': year,
                'M_sw': (sw * grp['mu_acf']).sum(),
                'M_cs': (sw * grp['mu_cs1']).sum(),
                'N': len(grp),
            })
        ndf = pd.DataFrame(agg_n)

        ax.plot(ndf['year'], ndf['M_sw'], 'o-', color='red', linewidth=2,
                markersize=3, label='ACF')
        ax.plot(ndf['year'], ndf['M_cs'], 's--', color='blue', linewidth=1.5,
                markersize=3, label='Cost share')
        add_reform_lines(ax)
        ax.set_title(f'NACE {nace}: {label}')
        ax.set_xlabel('Year')
        if idx == 0:
            ax.set_ylabel('Sales-weighted markup')
        ax.legend(fontsize=8)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.suptitle('Industry-Specific Markups (CZ-NACE F)', fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'dleu_fig_industry_trends.pdf', dpi=300,
                bbox_inches='tight')
    plt.close()


def fig_procurement_decomposition(panel):
    """
    Extension: Decomposition by procurement status.

    Compare aggregate markup trajectories for procurement vs non-procurement
    firms, and show within vs reallocation for each group.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    for pp_val, label, color in [(1, 'Procurement', 'red'), (0, 'Non-procurement', 'blue')]:
        sub = panel[panel['pp_dummy'] == pp_val]
        if sub.empty or len(sub) < 50:
            continue

        agg = aggregate_markups(sub, 'mu_acf')
        decomp = decompose_dleu(sub, 'mu_acf')

        ax1.plot(agg['year'], agg['M_sw'], '-' if pp_val else '--',
                 color=color, linewidth=2, label=label)

        if not decomp.empty:
            base = agg['M_sw'].iloc[0]
            ax2.plot(decomp['year'], base + decomp['within'].cumsum(),
                     '-' if pp_val else '--', color=color, linewidth=2,
                     label=f'{label}: within')
            ax2.plot(decomp['year'], base + decomp['reallocation'].cumsum(),
                     ':' if pp_val else '-.', color=color, linewidth=1.5,
                     alpha=0.7, label=f'{label}: realloc.')

    add_reform_lines(ax1)
    add_reform_lines(ax2)
    ax1.set_xlabel('Year')
    ax1.set_ylabel('Sales-weighted markup')
    ax1.set_title('(A) Aggregate markup by procurement status')
    ax1.legend(loc='best')
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

    ax2.set_xlabel('Year')
    ax2.set_ylabel('Cumulative change from base')
    ax2.set_title('(B) Within vs reallocation by status')
    ax2.legend(fontsize=8, loc='best')
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.tight_layout()
    plt.savefig(FIG_DIR / 'dleu_fig_procurement_split.pdf', dpi=300)
    plt.close()


def fig_labor_capital_shares(panel):
    """Factor share evolution — DLEU Section V analog."""
    shares = []
    for year, grp in panel.groupby('year'):
        sw = grp['Sales'] / grp['Sales'].sum()
        row = {'year': year, 'capital_share': (sw * grp['capital_share']).sum()}
        if grp['Wages'].notna().any():
            wage_sub = grp[grp['Wages'].notna()]
            if len(wage_sub) > 10:
                sw_w = wage_sub['Sales'] / wage_sub['Sales'].sum()
                row['labor_share'] = (sw_w * wage_sub['labor_share']).sum()
        shares.append(row)
    sdf = pd.DataFrame(shares)

    fig, ax = plt.subplots()
    ax.plot(sdf['year'], sdf['capital_share'], '--', color='orange',
            linewidth=2, label='Capital share (rK/Sales)')
    if 'labor_share' in sdf.columns and sdf['labor_share'].notna().any():
        ax.plot(sdf['year'], sdf['labor_share'], '-', color='green',
                linewidth=2, label='Labor share (W/Sales)')
    add_reform_lines(ax)
    ax.set_xlabel('Year')
    ax.set_ylabel('Factor share of sales')
    ax.set_title('Factor Shares of Sales — Czech Construction')
    ax.legend(loc='best')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'dleu_fig_factor_shares.pdf', dpi=300)
    plt.close()


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

def table_sectoral_decomposition(sec_decomp):
    """Table I equivalent: Sectoral decomposition by NACE."""
    if sec_decomp.empty:
        return

    tex = [
        r'\begin{table}[htbp]\centering',
        r'\caption{Sectoral Decomposition of Markup Change (DLEU Eq.~10)}\label{tab:dleu_sectoral}',
        r'\begin{threeparttable}',
        r'\begin{tabular}{lccccc}',
        r'\toprule',
        r'Period & Markup & $\Delta$Markup & $\Delta$Within & $\Delta$Between & $\Delta$Cross \\',
        r'\midrule',
    ]
    for _, r in sec_decomp.iterrows():
        tex.append(
            f'{r["period"]} & {r["M_agg"]:.3f} & {r["dM"]:+.3f} & '
            f'{r["within"]:+.3f} & {r["between"]:+.3f} & {r["cross"]:+.3f} \\\\'
        )
    tex += [
        r'\bottomrule',
        r'\end{tabular}',
        r'\begin{tablenotes}\footnotesize',
        r'\item \textit{Notes:} Sectoral decomposition of changes in aggregate markup following De Loecker, Eeckhout, and Unger (2020, eq.~10). Sectors are CZ-NACE 41 (buildings), 42 (civil engineering), 43 (specialized construction). $\Delta$Within: change in markup at the industry level weighted by lagged share. $\Delta$Between: change in industry share weighted by lagged markup. $\Delta$Cross: joint change.',
        r'\end{tablenotes}',
        r'\end{threeparttable}',
        r'\end{table}',
    ]
    with open(TAB_DIR / 'dleu_table1_sectoral_decomp.tex', 'w') as f:
        f.write('\n'.join(tex))


def table_firm_decomposition(decomp):
    """Appendix 4 equivalent: Year-by-year firm-level decomposition (eq. 9)."""
    if decomp.empty:
        return

    tex = [
        r'\begin{table}[htbp]\centering',
        r'\caption{Decomposition of Aggregate Markups at the Firm Level (DLEU Eq.~9)}\label{tab:dleu_firm_decomp}',
        r'\begin{threeparttable}',
        r'\begin{tabular}{lcccccc}',
        r'\toprule',
        r'Year & $\Delta M$ & $\Delta$Within & $\Delta$Mkt.\ Share & $\Delta$Cross & Net Entry & $N$ \\',
        r'\midrule',
    ]
    for _, r in decomp.iterrows():
        tex.append(
            f'{int(r["year"])} & {r["dM"]:+.3f} & {r["within"]:+.3f} & '
            f'{r["mkt_share"]:+.3f} & {r["cross"]:+.3f} & '
            f'{r["net_entry"]:+.3f} & {int(r["N_inc"]):,} \\\\'
        )
    tex += [
        r'\bottomrule',
        r'\end{tabular}',
        r'\begin{tablenotes}\footnotesize',
        r'\item \textit{Notes:} Firm-level decomposition following De Loecker, Eeckhout, and Unger (2020, eq.~9) with Haltiwanger (1997) demeaning. $\Delta$Within: $\sum_i m_{i,t-1}\Delta\mu_{it}$. $\Delta$Mkt.\ Share: $\sum_i \tilde{\mu}_{i,t-1}\Delta m_{it}$ where $\tilde{\mu}_{i,t-1} = \mu_{i,t-1} - M_{t-1}$. $\Delta$Cross: $\sum_i \Delta\mu_{it}\Delta m_{it}$. Net Entry: entrant contribution minus exiter contribution. $N$: number of incumbent firms.',
        r'\end{tablenotes}',
        r'\end{threeparttable}',
        r'\end{table}',
    ]
    with open(TAB_DIR / 'dleu_table_firm_decomp.tex', 'w') as f:
        f.write('\n'.join(tex))


def table_summary_by_spec(agg_specs, agg_acf):
    """Summary table: aggregate markup under different specifications."""
    years_show = [2007, 2010, 2014, 2018, 2021]

    tex = [
        r'\begin{table}[htbp]\centering',
        r'\caption{Aggregate Markup by Specification and Year}\label{tab:dleu_markup_specs}',
        r'\begin{threeparttable}',
        r'\begin{tabular}{lccccc}',
        r'\toprule',
        r'Specification & ' + ' & '.join(str(y) for y in years_show) + r' \\',
        r'\midrule',
    ]

    specs = [
        ('ACF (benchmark)', agg_acf),
        (r'Calibrated $\theta=0.85$', agg_specs['cs085']),
        ('Firm cost share', agg_specs['cs1']),
        ('Industry cost share', agg_specs['cs_ind']),
    ]
    for label, df in specs:
        vals = []
        for y in years_show:
            row = df[df['year'] == y]
            vals.append(f'{row["M_sw"].iloc[0]:.3f}' if len(row) > 0 else '--')
        tex.append(f'{label} & ' + ' & '.join(vals) + r' \\')

    tex += [
        r'\bottomrule',
        r'\end{tabular}',
        r'\begin{tablenotes}\footnotesize',
        r'\item \textit{Notes:} Sales-weighted aggregate markup under alternative specifications. ACF uses production-function-estimated output elasticities $\hat{\theta}_{st}$. Calibrated sets $\theta = 0.85$ (constant). Firm cost share uses $\alpha^V_{it} = \text{COGS}_{it}/(\text{COGS}_{it} + r_t K_{it})$. Industry cost share uses $\text{median}_{i \in s}(\alpha^V_{it})$ within NACE 2-digit $\times$ year cells.',
        r'\end{tablenotes}',
        r'\end{threeparttable}',
        r'\end{table}',
    ]
    with open(TAB_DIR / 'dleu_table_markup_specs.tex', 'w') as f:
        f.write('\n'.join(tex))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print('=' * 60)
    print('DLEU Replication — Czech Construction')
    print('=' * 60)

    # Load and construct
    panel = load_data()
    panel = construct_markups(panel)
    print(f'Panel: {len(panel):,} obs, {panel["id"].nunique():,} firms, '
          f'{panel["year"].min()}-{panel["year"].max()}')

    # 1. Aggregate markups — multiple specifications
    agg_acf = aggregate_markups(panel, 'mu_acf')
    agg_specs = {
        'cs085': aggregate_markups(panel, 'mu_cs085'),
        'cs1': aggregate_markups(panel, 'mu_cs1'),
        'cs_ind': aggregate_markups(panel, 'mu_cs_ind'),
    }

    print('\nAggregate Markups (sales-weighted):')
    summary = agg_acf[['year', 'M_sw', 'M_mean', 'M_cogs', 'P50', 'P75', 'P90', 'N']]
    print(summary.to_string(index=False, float_format='%.3f'))

    # 2. Firm-level decomposition (DLEU eq. 9)
    decomp = decompose_dleu(panel, 'mu_acf')
    print('\nFirm-level Decomposition (DLEU eq. 9, Haltiwanger demeaning):')
    print(decomp[['year', 'dM', 'within', 'mkt_share', 'cross', 'net_entry']
                 ].to_string(index=False, float_format='%+.4f'))

    # Cumulative totals
    if not decomp.empty:
        cum = decomp[['within', 'reallocation', 'net_entry']].sum()
        total_change = agg_acf['M_sw'].iloc[-1] - agg_acf['M_sw'].iloc[0]
        print(f'\nCumulative ({panel["year"].min()}-{panel["year"].max()}):')
        print(f'  Total ΔM:     {total_change:+.4f}')
        print(f'  Within:       {cum["within"]:+.4f} ({cum["within"]/total_change*100:+.1f}%)')
        print(f'  Reallocation: {cum["reallocation"]:+.4f} ({cum["reallocation"]/total_change*100:+.1f}%)')
        print(f'  Net entry:    {cum["net_entry"]:+.4f} ({cum["net_entry"]/total_change*100:+.1f}%)')

    # 3. Sectoral decomposition (DLEU eq. 10)
    sec_decomp = decompose_sectoral(panel, 'mu_acf', period_length=5)
    if not sec_decomp.empty:
        print('\nSectoral Decomposition (5-year changes, DLEU eq. 10):')
        print(sec_decomp[['period', 'M_agg', 'dM', 'within', 'between', 'cross']
                         ].to_string(index=False, float_format='%+.4f'))

    # 4. Generate all figures
    print('\nGenerating figures...')
    fig1_aggregate_markup(agg_acf, agg_specs)
    fig2_sensitivity(agg_acf, panel)
    fig3_distribution(panel, agg_acf)
    fig4_decomposition(decomp, agg_acf)
    fig5_micro_vs_agg(panel, agg_acf)
    fig7_cost_shares(panel)
    pr_df = fig8_profit_rate(panel)
    fig12_cost_share_markup(panel, agg_acf)
    fig_industry_trends(panel)
    fig_procurement_decomposition(panel)
    fig_labor_capital_shares(panel)
    print(f'  Saved 11 figures to {FIG_DIR}/')

    # 5. Generate tables
    print('Generating tables...')
    table_sectoral_decomposition(sec_decomp)
    table_firm_decomposition(decomp)
    table_summary_by_spec(agg_specs, agg_acf)
    print(f'  Saved 3 tables to {TAB_DIR}/')

    # 6. Save data outputs
    agg_acf.to_csv(DAT_DIR / 'dleu_aggregate_markups.csv', index=False)
    decomp.to_csv(DAT_DIR / 'dleu_firm_decomposition.csv', index=False)
    if not sec_decomp.empty:
        sec_decomp.to_csv(DAT_DIR / 'dleu_sectoral_decomposition.csv', index=False)

    # Save Stata-compatible
    agg_acf.to_stata(str(DAT_DIR / 'dleu_aggregate_markups.dta'),
                     write_index=False, version=118)
    print(f'  Saved data to {DAT_DIR}/')

    # Key comparison with DLEU
    print('\n' + '=' * 60)
    print('KEY FINDINGS — Czech Construction vs DLEU (US)')
    print('=' * 60)
    first_yr = agg_acf['year'].iloc[0]
    last_yr = agg_acf['year'].iloc[-1]
    first_mu = agg_acf['M_sw'].iloc[0]
    last_mu = agg_acf['M_sw'].iloc[-1]
    print(f'  Czech aggregate markup: {first_mu:.3f} ({first_yr}) → '
          f'{last_mu:.3f} ({last_yr})  [Δ = {last_mu - first_mu:+.3f}]')
    print(f'  DLEU US benchmark:      1.21 (1980) → 1.61 (2016)  [Δ = +0.40]')
    gap = agg_acf['M_sw'].max() - agg_acf['M_mean'].max()
    print(f'  SW-mean gap (max):      {gap:.3f}  [DLEU: ~0.20-0.40]')
    if not decomp.empty:
        cum = decomp[['within', 'reallocation', 'net_entry']].sum()
        print(f'  Within share:           {cum["within"]/total_change*100:.1f}%  '
              f'[DLEU: ~33%]')
        print(f'  Reallocation share:     {cum["reallocation"]/total_change*100:.1f}%  '
              f'[DLEU: ~67%]')

    print('\nDone.')


if __name__ == '__main__':
    main()
