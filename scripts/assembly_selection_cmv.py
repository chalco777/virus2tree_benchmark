# Select 34 assemblies stratified by CONTINENT × time window, boosting 2020–2025,
# and preferring known COLLECTION DATE. No per-country limits.
# Output files + a displayed table. Ends with criteria summary to provide in prose.

# This script was created with Chatgpt

import pandas as pd, numpy as np, re
from pathlib import Path
from ace_tools import display_dataframe_to_user

runs_path = Path("/mnt/data/cmv_runs.xlsx")
assem_path = Path("/mnt/data/cmv_assem.xlsx")

# ---------- Helpers ----------
def clean_cols(df):
    df = df.copy()
    df.columns = (
        df.columns.str.strip()
        .str.lower()
        .str.replace(r"\s+", "_", regex=True)
        .str.replace(r"[^0-9a-zA-Z_]", "", regex=True)
    )
    return df

def read_best_sheet(p):
    xl = pd.ExcelFile(p)
    cand = [s for s in xl.sheet_names if re.search(r"(main|runs|assembl|sheet1|data)", s, re.I)]
    sheet = cand[0] if cand else xl.sheet_names[0]
    return clean_cols(xl.parse(sheet))

def clean_country(x):
    if pd.isna(x): return np.nan
    s = str(x)
    s = re.split(r"[:;/|,]", s)[0].strip()
    repl = {"United States":"USA","United States of America":"USA","U.S.A.":"USA","U.S.A":"USA",
            "UK":"United Kingdom","England":"United Kingdom","Scotland":"United Kingdom",
            "Wales":"United Kingdom","Northern Ireland":"United Kingdom",
            "Korea, Republic of":"South Korea","Korea, South":"South Korea",
            "Viet Nam":"Vietnam","Czech Republic":"Czechia","Russian Federation":"Russia",
            "PR China":"China","People's Republic of China":"China"}
    return repl.get(s, s)

# Minimal country->continent map for observed values (extend as needed)
CONTINENT = {
    # Europe
    "United Kingdom":"Europe","Germany":"Europe","Italy":"Europe","Austria":"Europe","Czechia":"Europe",
    "France":"Europe","Belgium":"Europe","Netherlands":"Europe","Greece":"Europe","Spain":"Europe",
    "Portugal":"Europe","Poland":"Europe","Switzerland":"Europe","Denmark":"Europe","Sweden":"Europe",
    "Norway":"Europe","Finland":"Europe","Ireland":"Europe","Hungary":"Europe","Romania":"Europe",
    "Bulgaria":"Europe","Serbia":"Europe","Slovakia":"Europe","Slovenia":"Europe","Croatia":"Europe",
    # North America
    "USA":"North America","Canada":"North America","Mexico":"North America",
    # South America
    "Brazil":"South America","Argentina":"South America","Peru":"South America","Chile":"South America",
    "Colombia":"South America","Uruguay":"South America","Paraguay":"South America","Ecuador":"South America",
    # Asia
    "Israel":"Asia","China":"Asia","South Korea":"Asia","Japan":"Asia","India":"Asia","Singapore":"Asia",
    "Taiwan":"Asia","Thailand":"Asia","Vietnam":"Asia","Bangladesh":"Asia","Saudi Arabia":"Asia","Qatar":"Asia",
    "United Arab Emirates":"Asia","Turkey":"Asia",
    # Africa
    "South Africa":"Africa","Zambia":"Africa","Nigeria":"Africa","Uganda":"Africa","Kenya":"Africa","Ghana":"Africa",
    "Ethiopia":"Africa","Egypt":"Africa","Morocco":"Africa","Tunisia":"Africa",
    # Oceania
    "Australia":"Oceania","New Zealand":"Oceania",
}

def to_year(x):
    if pd.isna(x): return np.nan
    m = re.search(r"(19|20)\d{2}", str(x))
    return int(m.group(0)) if m else np.nan

def year_bin(y):
    if pd.isna(y): return "unknown"
    y = int(y)
    if y < 2010: return "<2010"
    if y < 2015: return "2010-2014"
    if y < 2020: return "2015-2019"
    if y <= 2025: return "2020-2025"
    return ">2025"

def continent_of(country):
    if pd.isna(country): return "Unknown"
    return CONTINENT.get(str(country), "Unknown")

def completeness_score(row):
    # Strongly reward collection date known
    score = 0
    if pd.notna(row.get('coll_year')): score += 3  # collection known
    if pd.notna(row.get('country')): score += 1
    if pd.notna(row.get('bioproject')): score += 1
    if pd.notna(row.get('biosample')): score += 1
    if pd.notna(row.get('isolate')): score += 1
    if pd.notna(row.get('tissue_specimen_source')): score += 1
    return score

def pick_diverse_scored(df, k):
    if df.empty or k <= 0: return df.iloc[0:0]
    df = df.copy()
    # Score
    df['_meta'] = df.apply(completeness_score, axis=1)
    # Within-year preference; boost newer bins
    y_order = {"2020-2025":4,"2015-2019":3,"2010-2014":2,"<2010":1,"unknown":0}
    df['_ybin'] = df['year_bin'].map(y_order).fillna(0)
    comp = pd.to_numeric(df.get('nuc_completeness'), errors='coerce')
    df['_comp'] = comp.fillna(comp.median(skipna=True) if not np.isnan(comp).all() else 0)
    df['_len'] = pd.to_numeric(df.get('length'), errors='coerce').fillna(0)
    df['tissue_f'] = df.get('tissue_specimen_source').fillna('unknown')
    # Sort by (metadata score desc, newer bin desc, completeness desc, length desc)
    df = df.sort_values(by=['_meta','_ybin','_comp','_len'], ascending=[False,False,False,False])
    # Interleave by tissue to avoid same-source monotony
    picks = []
    for t, g in df.groupby('tissue_f', dropna=False):
        for _, row in g.iterrows():
            picks.append(row)
    seen, out = set(), []
    for row in picks:
        acc = row.get('accession')
        if acc not in seen:
            out.append(row)
            seen.add(acc)
        if len(out) >= k: break
    return pd.DataFrame(out)

# ---------- Load & derive ----------
runs = read_best_sheet(runs_path)
assem = read_best_sheet(assem_path)

# RUNS: country -> continent, year bins
runs['country'] = pd.Series(np.nan, index=runs.index)
for col in ['geo_loc_name_country','geographic_location_country_andor_sea','geo_loc_name','geographic_location']:
    if col in runs.columns:
        runs['country'] = runs['country'].fillna(runs[col].map(clean_country))
runs['continent'] = runs['country'].map(continent_of)
runs['coll_year'] = runs.get('collection_date', pd.Series(np.nan, index=runs.index)).map(to_year)
runs['rel_year'] = runs.get('releasedate', pd.Series(np.nan, index=runs.index)).map(to_year)
runs['year'] = runs['coll_year'].where(~runs['coll_year'].isna(), runs['rel_year'])
runs['year_bin'] = runs['year'].map(year_bin)

# Assemblies: country -> continent, year from collection (preferred) else release
assem['country'] = assem.get('geo_location').map(clean_country)
assem['continent'] = assem['country'].map(continent_of)
assem['coll_year'] = assem.get('collection_date', pd.Series(np.nan, index=assem.index)).map(to_year)
assem['rel_year'] = assem.get('release_date', pd.Series(np.nan, index=assem.index)).map(to_year)
assem['year'] = assem['coll_year'].where(~assem['coll_year'].isna(), assem['rel_year'])
assem['year_bin'] = assem['year'].map(year_bin)

# ---------- Build continent×time weights from RUNS ----------
runs_ct = runs[runs['continent']!="Unknown"].pivot_table(index='continent', columns='year_bin', values='run', aggfunc='count', fill_value=0)
avail_ct = assem[assem['continent']!="Unknown"].groupby(['continent','year_bin']).size().unstack(fill_value=0)

# Align columns
all_cols = sorted(set(runs_ct.columns) | set(avail_ct.columns))
runs_ct = runs_ct.reindex(columns=all_cols, fill_value=0)
avail_ct = avail_ct.reindex(columns=all_cols, fill_value=0)

# Boost 2020-2025
weights = runs_ct.copy().astype(float)
if '2020-2025' in weights.columns:
    weights['2020-2025'] = weights['2020-2025'] * 1.8  # stronger bias to recent
# Normalize
total_w = weights.values.sum()
if total_w == 0:
    weights[:] = 1.0
    total_w = weights.values.sum()
weights = weights / total_w

# ---------- Quotas ----------
TOTAL = 34
frac_q = weights * TOTAL
quota = frac_q.round().astype(int)

# Cap by availability and adjust to sum to TOTAL
import numpy as np
def adjust(q, target, availability, pref):
    q = q.copy()
    availability = availability.reindex_like(q).fillna(0).astype(int)
    # cap
    q = pd.DataFrame(np.minimum(q.values, availability.values), index=q.index, columns=q.columns).astype(int)
    # fill up
    while q.values.sum() < target:
        diff = (pref - q).stack()
        resid = (availability - q).stack()
        resid = resid[resid > 0]
        if resid.empty: break
        cand = diff[resid.index].idxmax()
        q.loc[cand] += 1
    # trim if needed
    while q.values.sum() > target:
        nonz = q.stack()
        nonz = nonz[nonz > 0]
        if nonz.empty: break
        w = pref.stack()[nonz.index]
        cand = w.idxmin()
        q.loc[cand] -= 1
    return q

quota = adjust(quota, TOTAL, avail_ct, frac_q)

# ---------- Selection ----------
selected_parts = []
for cont in quota.index:
    for yb in quota.columns:
        k = int(quota.loc[cont, yb])
        if k <= 0: continue
        cand = assem[(assem['continent']==cont) & (assem['year_bin']==yb)].copy()
        # Prefer collection-known: split candidates
        cand_coll = cand[cand['coll_year'].notna()]
        cand_rel  = cand[cand['coll_year'].isna()]
        take1 = pick_diverse_scored(cand_coll, k)
        remaining = k - len(take1)
        if remaining > 0 and not cand_rel.empty:
            take2 = pick_diverse_scored(cand_rel[~cand_rel['accession'].isin(take1['accession'])], remaining)
            take = pd.concat([take1, take2], ignore_index=True)
        else:
            take = take1
        selected_parts.append(take)

selected = pd.concat(selected_parts, ignore_index=True) if selected_parts else assem.iloc[0:0]

# If underfilled (scarcity in some cells), top-up within same continent first, preferring 2020–2025 & coll known
def top_up(selection, target, pool, q, availability, pref):
    if len(selection) >= target: return selection
    chosen = set(selection['accession'])
    # (continent, year_bin) order by pref desc
    order = pref.stack().sort_values(ascending=False).index
    picks = []
    for cont, yb in order:
        if len(selection) + len(picks) >= target: break
        # favor collection-known within continent, any year_bin
        cont_pool = pool[(pool['continent']==cont) & (~pool['accession'].isin(chosen))].copy()
        # sort by collection-known first and 2020–2025 next
        cont_pool['coll_known'] = cont_pool['coll_year'].notna()
        y_order = {"2020-2025":4,"2015-2019":3,"2010-2014":2,"<2010":1,"unknown":0}
        cont_pool['_ybin'] = cont_pool['year_bin'].map(y_order).fillna(0)
        comp = pd.to_numeric(cont_pool.get('nuc_completeness'), errors='coerce')
        cont_pool['_comp'] = comp.fillna(comp.median(skipna=True) if not np.isnan(comp).all() else 0)
        cont_pool = cont_pool.sort_values(by=['coll_known','_ybin','_comp','length'], ascending=[False,False,False,False])
        for _, row in cont_pool.iterrows():
            if len(selection) + len(picks) >= target: break
            acc = row['accession']
            if acc in chosen: continue
            picks.append(row)
            chosen.add(acc)
    if picks:
        selection = pd.concat([selection, pd.DataFrame(picks)], ignore_index=True).head(target)
    return selection

selected = top_up(selected, TOTAL, assem, quota, avail_ct.reindex_like(quota).fillna(0).astype(int), frac_q)

# ---------- Output ----------
cols_out = [
    'accession','organism_name','species','continent','country',
    'tissue_specimen_source','coll_year','rel_year','year','year_bin',
    'length','nuc_completeness','bioproject','biosample','isolate','assembly','release_date'
]
cols_out = [c for c in cols_out if c in selected.columns]
selected = selected[cols_out].sort_values(by=['continent','year','country','accession'], ascending=[True, True, True, True])

xlsx = "/mnt/data/cmv_selected_34_continent_time_collpref.xlsx"
csv  = "/mnt/data/cmv_selected_34_continent_time_collpref.csv"
txt  = "/mnt/data/cmv_selected_34_accessions_continent_time.txt"
selected.to_excel(xlsx, index=False)
selected.to_csv(csv, index=False)
with open(txt, "w") as f:
    f.write("\n".join(selected['accession'].astype(str).tolist()))

display_dataframe_to_user("CMV seleccionados (34) por continente×ventana temporal, priorizando fecha de colección", selected)

print(f"Selected: {len(selected)}")
print(xlsx)
print(csv)
print(txt)

# Quick check: distribution by continent/year_bin and collection-known fraction
print("\nBy continent × year_bin:")
print(selected.groupby(['continent','year_bin']).size().unstack(fill_value=0))

print("\nCollection known counts:")
print(selected['coll_year'].notna().value_counts())
