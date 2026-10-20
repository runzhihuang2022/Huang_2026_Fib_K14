from pathlib import Path
import csv
import json
import math
import random
import time
import urllib.parse
import urllib.request
from collections import defaultdict

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


INPUT_XLSX = Path(r"data/external/local_47a98e1a880c")
OUT_DIR = Path(__file__).resolve().parent
CACHE_DIR = OUT_DIR / "cache"
CACHE_DIR.mkdir(exist_ok=True)

CONTROL = "Control"
CASE = "IGFL1"
ORGANISM = "human"
GO_LIBRARY = "GO_Biological_Process_2023"
KEGG_LIBRARY = "KEGG_2021_Human"
RANDOM_SEED = 20260624


def font(size=18, bold=False):
    candidates = [
        r"data/external/local_28ac8473f586" if bold else r"data/external/local_16a706108b0a",
        r"data/external/local_bb3719575541" if bold else r"data/external/local_474970ae2dcd",
    ]
    for p in candidates:
        if Path(p).exists():
            return ImageFont.truetype(p, size)
    return ImageFont.load_default()


F12, F14, F16, F18, F20, F24, F30 = font(12), font(14), font(16), font(18), font(20), font(24), font(30, True)
FB16, FB18, FB22 = font(16, True), font(18, True), font(22, True)


COLORS = {
    "ink": (40, 44, 52),
    "muted": (105, 113, 125),
    "grid": (222, 226, 232),
    "panel": (250, 250, 249),
    "up": (190, 45, 55),
    "down": (45, 105, 175),
    "neutral": (178, 183, 190),
    "accent": (48, 132, 112),
    "orange": (242, 156, 18),
    "green": (33, 140, 90),
}


def draw_text(draw, xy, text, fill=None, fnt=None, anchor=None):
    draw.text(xy, str(text), fill=fill or COLORS["ink"], font=fnt or F16, anchor=anchor)


def save_figure(img, out):
    out = Path(out)
    if out.suffix.lower() != ".png":
        out = out.with_suffix(".png")
    img.save(out)
    img.save(out.with_suffix(".pdf"), "PDF", resolution=300.0)


def ellipsize(text, max_chars):
    text = str(text)
    return text if len(text) <= max_chars else text[:max_chars - 3] + "..."


def wrap_full_text(text, max_chars):
    text = str(text)
    words = text.replace("_", " ").split()
    lines, current = [], ""
    for word in words:
        candidate = word if not current else current + " " + word
        if len(candidate) <= max_chars:
            current = candidate
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    return lines or [text]


def safe_name(name):
    keep = []
    for ch in name:
        keep.append(ch if ch.isalnum() or ch in "-_." else "_")
    return "".join(keep)[:120]


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        val = ranked[i] * n / (i + 1)
        prev = min(prev, val)
        adj[order[i]] = min(prev, 1.0)
    return adj


def normal_p_from_z(z):
    return float(math.erfc(abs(float(z)) / math.sqrt(2.0)))


def log_comb(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k_minus_1, M, n, N):
    # P[X >= k] for population M, successes n, draws N.
    k = k_minus_1 + 1
    hi = min(n, N)
    denom = log_comb(M, N)
    vals = [math.exp(log_comb(n, x) + log_comb(M - n, N - x) - denom) for x in range(k, hi + 1)]
    return min(1.0, max(0.0, sum(vals)))


def urlopen_retry(req, timeout=30, tries=3):
    last = None
    for i in range(tries):
        try:
            return urllib.request.urlopen(req, timeout=timeout)
        except Exception as exc:
            last = exc
            time.sleep(1.5 * (i + 1))
    raise last


def fetch_mygene_mapping(ensembl_ids):
    cache = CACHE_DIR / "mygene_mapping.json"
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    ids = [g.split(".")[0] for g in sorted(set(ensembl_ids))]
    mapping = {}
    url = "https://mygene.info/v3/query"
    for start in range(0, len(ids), 1000):
        chunk = ids[start:start + 1000]
        data = urllib.parse.urlencode({
            "q": ",".join(chunk),
            "scopes": "ensembl.gene",
            "fields": "symbol,entrezgene,name",
            "species": ORGANISM,
            "size": 1000,
        }).encode("utf-8")
        req = urllib.request.Request(url, data=data, headers={"User-Agent": "Codex bulk RNA-seq analysis"})
        with urlopen_retry(req, timeout=45) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
        hits = payload if isinstance(payload, list) else payload.get("hits", [])
        for hit in hits:
            ens = hit.get("query")
            if ens:
                mapping[ens] = {
                    "symbol": hit.get("symbol", ens),
                    "entrezgene": hit.get("entrezgene", ""),
                    "name": hit.get("name", ""),
                }
    cache.write_text(json.dumps(mapping, ensure_ascii=False, indent=2), encoding="utf-8")
    return mapping


def fetch_enrichr_library(library_name):
    cache = CACHE_DIR / f"{safe_name(library_name)}.gmt.json"
    if cache.exists():
        data = json.loads(cache.read_text(encoding="utf-8"))
        return {k: set(v) for k, v in data.items()}
    url = "https://maayanlab.cloud/Enrichr/geneSetLibrary?" + urllib.parse.urlencode({
        "mode": "text",
        "libraryName": library_name,
    })
    req = urllib.request.Request(url, headers={"User-Agent": "Codex bulk RNA-seq analysis"})
    with urlopen_retry(req, timeout=90) as resp:
        text = resp.read().decode("utf-8", errors="replace")
    gene_sets = {}
    for line in text.splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            term = parts[0]
            genes = {g.upper() for g in parts[2:] if g}
            if genes:
                gene_sets[term] = genes
    cache.write_text(json.dumps({k: sorted(v) for k, v in gene_sets.items()}, ensure_ascii=False), encoding="utf-8")
    return gene_sets


def read_inputs():
    counts = pd.read_excel(INPUT_XLSX, sheet_name="counts")
    if "gene_id" not in counts.columns:
        counts = counts.rename(columns={counts.columns[0]: "gene_id"})
    sample_cols = [c for c in counts.columns if c not in {"gene_id", "gene_name"}]
    meta = pd.DataFrame({
        "specimen": sample_cols,
        "condition": [CONTROL if str(c).upper().startswith("PBS") else CASE for c in sample_cols],
    })
    counts[sample_cols] = counts[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0).astype(float)
    return counts, meta, sample_cols


def size_factors(count_mat):
    gm = np.exp(np.mean(np.log(np.where(count_mat > 0, count_mat, np.nan)), axis=1))
    valid = np.isfinite(gm) & (gm > 0)
    ratios = count_mat[valid, :] / gm[valid, None]
    sf = np.nanmedian(ratios, axis=0)
    sf = sf / np.exp(np.mean(np.log(sf)))
    return sf


def differential_expression(counts, meta, sample_cols):
    mat = counts[sample_cols].to_numpy(dtype=float)
    keep = mat.sum(axis=1) >= 10
    counts_f = counts.loc[keep].copy().reset_index(drop=True)
    mat = counts_f[sample_cols].to_numpy(dtype=float)
    sf = size_factors(mat)
    norm = mat / sf[None, :]
    lib = norm.sum(axis=0)
    logcpm = np.log2((norm / lib[None, :] * 1e6) + 1.0)

    cond = meta.set_index("specimen").loc[sample_cols, "condition"].to_numpy()
    c_idx = np.where(cond == CONTROL)[0]
    k_idx = np.where(cond == CASE)[0]
    c = logcpm[:, c_idx]
    k = logcpm[:, k_idx]
    log2fc = k.mean(axis=1) - c.mean(axis=1)
    base_mean = norm.mean(axis=1)
    var_c = c.var(axis=1, ddof=1)
    var_k = k.var(axis=1, ddof=1)
    se = np.sqrt(var_c / len(c_idx) + var_k / len(k_idx) + 1e-6)
    z = log2fc / se
    pvalue = np.array([normal_p_from_z(v) for v in z])
    padj = bh_adjust(pvalue)
    res = pd.DataFrame({
        "gene_id": counts_f["gene_id"],
        "input_gene_name": counts_f["gene_name"] if "gene_name" in counts_f.columns else "",
        "baseMean": base_mean,
        "log2FoldChange": log2fc,
        "stat": z,
        "pvalue": pvalue,
        "padj": padj,
    })
    return res, pd.DataFrame(logcpm, columns=sample_cols), sf, keep


def annotate(res):
    mapping = fetch_mygene_mapping(res["gene_id"].tolist())
    symbols, entrez, names = [], [], []
    input_names = res["input_gene_name"].astype(str).tolist() if "input_gene_name" in res.columns else [""] * len(res)
    for gid, original in zip(res["gene_id"], input_names):
        key = gid.split(".")[0]
        m = mapping.get(key, {})
        symbols.append(m.get("symbol", original if original and original.lower() != "nan" else key))
        entrez.append(m.get("entrezgene", ""))
        names.append(m.get("name", ""))
    out = res.copy()
    out.insert(1, "ensembl_base", out["gene_id"].str.split(".").str[0])
    out.insert(2, "symbol", symbols)
    out.insert(3, "entrezgene", entrez)
    out.insert(4, "description", names)
    return out


def ora(de_genes, universe_genes, gene_sets, label):
    de = {g.upper() for g in de_genes if isinstance(g, str) and g}
    universe = {g.upper() for g in universe_genes if isinstance(g, str) and g}
    M = len(universe)
    N = len(de)
    rows = []
    for term, genes in gene_sets.items():
        gs = set(genes) & universe
        if len(gs) < 5:
            continue
        overlap = sorted(de & gs)
        if not overlap:
            continue
        p = hypergeom_sf(len(overlap) - 1, M, len(gs), N)
        rows.append({
            "category": label,
            "term": term,
            "overlap_n": len(overlap),
            "term_size": len(gs),
            "query_size": N,
            "pvalue": p,
            "genes": ";".join(overlap[:50]),
            "top8_genes": ", ".join(overlap[:8]),
        })
    df = pd.DataFrame(rows)
    if len(df):
        df["padj"] = bh_adjust(df["pvalue"].values)
        df = df.sort_values(["padj", "pvalue", "overlap_n"], ascending=[True, True, False])
    return df


def enrichment_score(ranked_genes, scores, gene_set):
    genes = np.array([g.upper() for g in ranked_genes])
    hits = np.isin(genes, list(gene_set))
    Nh = hits.sum()
    N = len(genes)
    if Nh == 0 or Nh == N:
        return 0.0
    weights = np.abs(scores)
    hit_weights = weights[hits]
    hit_norm = hit_weights.sum()
    miss_norm = N - Nh
    running = np.cumsum(np.where(hits, weights / hit_norm, -1.0 / miss_norm))
    max_es = running.max()
    min_es = running.min()
    return float(max_es if abs(max_es) >= abs(min_es) else min_es)


def running_enrichment_profile(ranked_genes, scores, gene_set):
    genes = np.array([g.upper() for g in ranked_genes])
    hits = np.isin(genes, list(gene_set))
    Nh = hits.sum()
    N = len(genes)
    if Nh == 0 or Nh == N:
        return np.zeros(N), np.where(hits)[0]
    weights = np.abs(scores)
    hit_norm = weights[hits].sum()
    miss_norm = N - Nh
    running = np.cumsum(np.where(hits, weights / hit_norm, -1.0 / miss_norm))
    return running, np.where(hits)[0]


def ranked_gene_table(de):
    ranked = de.dropna(subset=["symbol", "stat", "log2FoldChange"]).copy()
    ranked["rank_score"] = (
        np.sign(ranked["log2FoldChange"])
        * (-np.log10(ranked["pvalue"].clip(lower=1e-300)))
        * np.abs(ranked["log2FoldChange"])
    )
    return ranked.sort_values("rank_score", ascending=False)


def run_gsea(de, gene_sets, nperm=200):
    rng = random.Random(RANDOM_SEED)
    ranked = ranked_gene_table(de)
    genes = ranked["symbol"].str.upper().to_numpy()
    scores = ranked["rank_score"].to_numpy(dtype=float)
    universe = set(genes)
    top_rank = set(genes[:1000])
    candidates = []
    for term, gs in gene_sets.items():
        gs2 = set(gs) & universe
        if 15 <= len(gs2) <= 500:
            candidates.append((len(gs2 & top_rank), term, gs2))
    candidates = sorted(candidates, reverse=True)[:350]
    rows = []
    universe_list = sorted(universe)
    for _, term, gs2 in candidates:
        es = enrichment_score(genes, scores, gs2)
        null = []
        size = len(gs2)
        for _ in range(nperm):
            sampled = set(rng.sample(universe_list, size))
            null.append(enrichment_score(genes, scores, sampled))
        null = np.array(null)
        if es >= 0:
            denom = np.mean(null[null >= 0]) if np.any(null >= 0) else 1.0
            p = (np.sum(null >= es) + 1) / (np.sum(null >= 0) + 1)
        else:
            denom = abs(np.mean(null[null < 0])) if np.any(null < 0) else 1.0
            p = (np.sum(null <= es) + 1) / (np.sum(null < 0) + 1)
        nes = es / denom if denom else es
        leading = [g for g in genes[:max(1, int(len(genes) * 0.25))] if g in gs2][:20]
        rows.append({"term": term, "size": len(gs2), "ES": es, "NES": nes, "pvalue": p, "leading_edge": ";".join(leading)})
    df = pd.DataFrame(rows)
    if len(df):
        df["padj"] = bh_adjust(df["pvalue"].values)
        df = df.sort_values(["padj", "NES"], ascending=[True, False])
    return df


def ssGSEA_like(logcpm, de, gene_sets, terms):
    expr = logcpm.copy()
    expr.index = de["symbol"].str.upper().values
    expr = expr.groupby(expr.index).mean()
    z = (expr.sub(expr.mean(axis=1), axis=0)).div(expr.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    rows = []
    for term in terms:
        genes = [g for g in gene_sets.get(term, set()) if g in z.index]
        if len(genes) < 5:
            continue
        vals = z.loc[genes].mean(axis=0)
        rows.append(pd.Series(vals, name=term))
    return pd.DataFrame(rows)


def axes(draw, x0, y0, x1, y1, xlabel="", ylabel=""):
    draw.line((x0, y1, x1, y1), fill=COLORS["ink"], width=3)
    draw.line((x0, y0, x0, y1), fill=COLORS["ink"], width=3)
    draw_text(draw, ((x0 + x1) / 2, y1 + 38), xlabel, COLORS["ink"], F18, anchor="mm")
    draw_text(draw, (x0, y0 - 25), ylabel, COLORS["ink"], F18)


def draw_pca(logcpm, meta, out):
    X = logcpm.T.to_numpy()
    X = X - X.mean(axis=0)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    pcs = U[:, :2] * S[:2]
    var = (S ** 2) / np.sum(S ** 2)
    w, h = 720, 520
    img = Image.new("RGB", (w, h), "white")
    d = ImageDraw.Draw(img)
    draw_text(d, (42, 26), "PCA", COLORS["ink"], F30)
    x0, y0, x1, y1 = 105, 82, 640, 410
    axes(d, x0, y0, x1, y1, f"PC1 ({var[0] * 100:.1f}%)", f"PC2 ({var[1] * 100:.1f}%)")
    xs, ys = pcs[:, 0], pcs[:, 1]
    padx = max(1e-6, (xs.max() - xs.min()) * 0.15)
    pady = max(1e-6, (ys.max() - ys.min()) * 0.15)
    def sx(v): return x0 + (v - xs.min() + padx) / (xs.max() - xs.min() + 2 * padx) * (x1 - x0)
    def sy(v): return y1 - (v - ys.min() + pady) / (ys.max() - ys.min() + 2 * pady) * (y1 - y0)
    meta2 = meta.set_index("specimen")
    palette = {CONTROL: COLORS["down"], CASE: COLORS["up"]}
    for i, sample in enumerate(logcpm.columns):
        cond = meta2.loc[sample, "condition"]
        x, y = sx(xs[i]), sy(ys[i])
        d.ellipse((x - 12, y - 12, x + 12, y + 12), fill=palette.get(cond, COLORS["neutral"]), outline=COLORS["ink"], width=2)
        draw_text(d, (x + 16, y - 10), sample, COLORS["ink"], FB16)
    draw_text(d, (500, 460), f"{CONTROL}", COLORS["down"], FB18)
    draw_text(d, (610, 460), f"{CASE}", COLORS["up"], FB18)
    save_figure(img, out)


def draw_heatmap(matrix, row_labels, col_labels, out, title):
    rows, cols = matrix.shape
    long_rows = max([len(str(r)) for r in row_labels] or [0]) > 24
    cell_w, cell_h = 72, (54 if long_rows else 34)
    left, top = (560 if long_rows else 210), 82
    w = left + cols * cell_w + 70
    h = top + rows * cell_h + 80
    img = Image.new("RGB", (w, h), "white")
    d = ImageDraw.Draw(img)
    draw_text(d, (28, 24), title, COLORS["ink"], FB22)
    v = np.nanmax(np.abs(matrix)) or 1
    def col(x):
        x = max(-v, min(v, x)) / v
        if x >= 0:
            return (245, int(245 - 130 * x), int(245 - 125 * x))
        return (int(245 + 150 * x), int(245 + 100 * x), 245)
    for j, c in enumerate(col_labels):
        draw_text(d, (left + j * cell_w + cell_w / 2, top - 24), c, COLORS["ink"], FB16, anchor="mm")
    for i, r in enumerate(row_labels):
        y = top + i * cell_h
        if long_rows:
            lines = wrap_full_text(r, 46)
            for k, line in enumerate(lines[:2]):
                draw_text(d, (left - 8, y + 8 + k * 18), line, COLORS["ink"], F14, anchor="rm")
        else:
            draw_text(d, (left - 8, y + cell_h / 2), str(r), COLORS["ink"], FB16, anchor="rm")
        for j in range(cols):
            x = left + j * cell_w
            d.rectangle((x, y, x + cell_w, y + cell_h), fill=col(matrix[i, j]), outline="white")
    draw_text(d, (left, h - 34), "row z-score", COLORS["muted"], F16)
    save_figure(img, out)


def draw_volcano(de, out):
    df = de.copy()
    df["raw_neglog10padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    nonzero = df.loc[df["padj"] > 0, "raw_neglog10padj"].replace([np.inf, -np.inf], np.nan).dropna()
    high_ref = float(nonzero.quantile(0.995)) if len(nonzero) else 80.0
    y_cap = min(320.0, max(60.0, high_ref * 1.35))
    y_threshold = y_cap * 0.72
    max_raw = max(float(df["raw_neglog10padj"].replace([np.inf, -np.inf], np.nan).max()), y_threshold + 1)
    compressed = df["raw_neglog10padj"].copy()
    high = compressed > y_threshold
    compressed.loc[high] = y_threshold + (
        np.log1p(compressed.loc[high] - y_threshold)
        / max(1e-9, np.log1p(max_raw - y_threshold))
    ) * (y_cap - y_threshold) * 0.86
    zero_mask = df["padj"] <= 0
    if zero_mask.any():
        order = np.linspace(0, 1, int(zero_mask.sum()), endpoint=True)
        df.loc[zero_mask, "plot_neglog10padj"] = y_cap * (0.78 + 0.17 * order)
        df.loc[~zero_mask, "plot_neglog10padj"] = compressed.loc[~zero_mask]
    else:
        df["plot_neglog10padj"] = compressed
    top = df.reindex(df["log2FoldChange"].abs().sort_values(ascending=False).index).head(10)
    w, h = 900, 650
    img = Image.new("RGB", (w, h), "white")
    d = ImageDraw.Draw(img)
    draw_text(d, (38, 26), "Volcano", COLORS["ink"], F30)
    x0, y0, x1, y1 = 105, 92, 780, 510
    axes(d, x0, y0, x1, y1, f"log2 fold change ({CASE} / {CONTROL})", "-log10 adjusted P")
    xs = df["log2FoldChange"].clip(-8, 8).to_numpy()
    xmin, xmax = -max(abs(xs.min()), abs(xs.max())), max(abs(xs.min()), abs(xs.max()))
    ymax = max(2, y_cap)
    def sx(v): return x0 + (v - xmin) / (xmax - xmin + 1e-9) * (x1 - x0)
    def sy(v): return y1 - v / ymax * (y1 - y0)
    d.line((sx(-1), y0, sx(-1), y1), fill=COLORS["grid"], width=1)
    d.line((sx(1), y0, sx(1), y1), fill=COLORS["grid"], width=1)
    d.line((x0, sy(-math.log10(0.05)), x1, sy(-math.log10(0.05))), fill=COLORS["grid"], width=1)
    for _, r in df.sample(min(8000, len(df)), random_state=RANDOM_SEED).iterrows():
        sig_up = r["padj"] < 0.05 and r["log2FoldChange"] > 1
        sig_dn = r["padj"] < 0.05 and r["log2FoldChange"] < -1
        color = COLORS["up"] if sig_up else COLORS["down"] if sig_dn else COLORS["neutral"]
        x, y = sx(max(xmin, min(xmax, r["log2FoldChange"]))), sy(min(ymax, r["plot_neglog10padj"]))
        d.ellipse((x - 2, y - 2, x + 2, y + 2), fill=color)
    left_i, right_i = 0, 0
    for _, r in top.iterrows():
        x, y = sx(max(xmin, min(xmax, r["log2FoldChange"]))), sy(min(ymax, r["plot_neglog10padj"]))
        d.ellipse((x - 6, y - 6, x + 6, y + 6), fill=(20, 20, 20))
        if r["log2FoldChange"] >= 0:
            label_x = x1 - 104
            label_y = y0 + 30 + right_i * 24
            right_i += 1
        else:
            label_x = 18
            label_y = y0 + 30 + left_i * 24
            left_i += 1
        d.line((x, y, label_x, label_y + 7), fill=COLORS["ink"], width=1)
        draw_text(d, (label_x, label_y), r["symbol"], COLORS["ink"], F12)
    save_figure(img, out)


def draw_barplot(enr, out, title):
    df = enr.head(10).copy()
    if df.empty:
        return
    df["score"] = -np.log10(df["padj"].clip(lower=1e-300))
    w, h = 1400, 980
    img = Image.new("RGB", (w, h), "white")
    d = ImageDraw.Draw(img)
    draw_text(d, (34, 24), title, COLORS["ink"], F30)
    x0, y0, x1 = 560, 88, 900
    row_h = 84
    maxv = max(df["score"].max(), 1)
    for i, (_, r) in enumerate(df.iloc[::-1].iterrows()):
        y = y0 + i * row_h
        term_lines = wrap_full_text(r["term"], 42)
        for j, line in enumerate(term_lines[:3]):
            draw_text(d, (x0 - 14, y + 3 + j * 22), line, COLORS["ink"], F16, anchor="rm")
        bw = (r["score"] / maxv) * (x1 - x0)
        d.rectangle((x0, y + 8, x0 + bw, y + 34), fill=COLORS["green"])
        draw_text(d, (x0 + bw + 8, y + 11), f"{r['score']:.1f}", COLORS["ink"], FB16)
        for j, line in enumerate(wrap_full_text(str(r["top8_genes"]), 44)[:2]):
            draw_text(d, (950, y + 4 + j * 20), line, COLORS["ink"], F14)
    draw_text(d, (x0, h - 48), "-log10(FDR)", COLORS["muted"], F16)
    draw_text(d, (950, 58), "Top genes", COLORS["muted"], FB16)
    save_figure(img, out)


def draw_gsea(gsea, out):
    df = gsea.reindex(gsea["NES"].abs().sort_values(ascending=False).index).head(10).copy()
    if df.empty:
        return
    df["score"] = df["NES"]
    w, h = 1300, 940
    img = Image.new("RGB", (w, h), "white")
    d = ImageDraw.Draw(img)
    draw_text(d, (36, 26), "GSEA NES", COLORS["ink"], F30)
    x0, y0, xmid = 540, 80, 650
    row_h = 82
    maxv = max(abs(df["score"]).max(), 1)
    d.line((xmid, y0 - 8, xmid, y0 + row_h * len(df) + 5), fill=COLORS["ink"], width=2)
    for i, (_, r) in enumerate(df.iloc[::-1].iterrows()):
        y = y0 + i * row_h
        for j, line in enumerate(wrap_full_text(r["term"], 44)[:3]):
            draw_text(d, (x0 - 12, y + 2 + j * 22), line, COLORS["ink"], F16, anchor="rm")
        bw = abs(r["score"]) / maxv * 260
        if r["score"] >= 0:
            d.rectangle((xmid, y + 8, xmid + bw, y + 28), fill=COLORS["up"])
        else:
            d.rectangle((xmid - bw, y + 8, xmid, y + 28), fill=COLORS["down"])
        draw_text(d, (1030, y + 9), f"NES={r['NES']:.2f}", COLORS["ink"], FB16)
    save_figure(img, out)


def draw_classic_gsea(term, row, de, gene_sets, out):
    ranked = ranked_gene_table(de)
    genes = ranked["symbol"].str.upper().to_numpy()
    scores = ranked["rank_score"].to_numpy(dtype=float)
    gene_set = set(gene_sets.get(term, set())) & set(genes)
    running, hit_pos = running_enrichment_profile(genes, scores, gene_set)
    n = len(genes)
    if n == 0 or len(hit_pos) == 0:
        return

    w, h = 760, 560
    img = Image.new("RGB", (w, h), "white")
    d = ImageDraw.Draw(img)
    for j, line in enumerate(wrap_full_text(term, 64)[:2]):
        draw_text(d, (34, 18 + j * 25), line, COLORS["ink"], FB18)
    draw_text(d, (34, 70), f"NES={row['NES']:.2f}  FDR={row['padj']:.3g}", COLORS["muted"], FB16)

    x0, x1 = 78, 720
    y0, y1 = 105, 290
    ymin, ymax = float(running.min()), float(running.max())
    pad = max(0.05, (ymax - ymin) * 0.12)
    ymin -= pad
    ymax += pad
    def sx(i): return x0 + i / max(1, n - 1) * (x1 - x0)
    def sy(v): return y1 - (v - ymin) / (ymax - ymin + 1e-9) * (y1 - y0)

    zero_y = sy(0)
    d.line((x0, zero_y, x1, zero_y), fill=(170, 170, 170), width=2)
    d.line((x0, y1, x1, y1), fill=COLORS["ink"], width=3)
    d.line((x0, y0, x0, y1), fill=COLORS["ink"], width=3)
    pts = [(sx(i), sy(v)) for i, v in enumerate(running)]
    if len(pts) > 1:
        d.line(pts, fill=COLORS["green"], width=4)
    draw_text(d, (16, (y0 + y1) / 2), "ES", COLORS["ink"], FB18)

    tick_y0, tick_y1 = 310, 360
    d.rectangle((x0, tick_y0, x1, tick_y1), outline=COLORS["ink"], width=2)
    for pos in hit_pos:
        x = sx(int(pos))
        d.line((x, tick_y0, x, tick_y1), fill=(20, 20, 20), width=1)

    metric_y0, metric_y1 = 385, 505
    d.line((x0, metric_y1, x1, metric_y1), fill=COLORS["ink"], width=3)
    d.line((x0, metric_y0, x0, metric_y1), fill=COLORS["ink"], width=3)
    max_abs = max(1e-9, float(np.max(np.abs(scores))))
    step = max(1, n // 420)
    for i in range(0, n, step):
        x = sx(i)
        v = max(-1, min(1, scores[i] / max_abs))
        y = metric_y1 - abs(v) * (metric_y1 - metric_y0)
        color = COLORS["up"] if scores[i] >= 0 else COLORS["down"]
        d.line((x, metric_y1, x, y), fill=color, width=2)
    draw_text(d, ((x0 + x1) / 2, 535), "Rank in ordered dataset", COLORS["ink"], FB16, anchor="mm")
    draw_text(d, (18, 444), "Rank", COLORS["ink"], FB16)
    save_figure(img, out)


def draw_top10_classic_gsea(gsea, de, gene_sets, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for old in list(out_dir.glob("GSEA_*.png")) + list(out_dir.glob("GSEA_*.pdf")):
        old.unlink()
    top = gsea.reindex(gsea["NES"].abs().sort_values(ascending=False).index).head(10).copy()
    rows = []
    for i, (_, row) in enumerate(top.iterrows(), start=1):
        term = row["term"]
        fname = f"GSEA_{i:02d}_{safe_name(term)}.png"
        draw_classic_gsea(term, row, de, gene_sets, out_dir / fname)
        rows.append({
            "rank": i,
            "category": row.get("category", ""),
            "term": term,
            "NES": row["NES"],
            "padj": row["padj"],
            "png": fname,
            "pdf": fname.replace(".png", ".pdf"),
        })
    pd.DataFrame(rows).to_csv(out_dir / "GSEA_top10_classic_index.csv", index=False, encoding="utf-8-sig")


def main():
    counts, meta, sample_cols = read_inputs()
    de_raw, logcpm, size_factor, keep = differential_expression(counts, meta, sample_cols)
    de = annotate(de_raw)
    de["direction"] = np.where((de["padj"] < 0.05) & (de["log2FoldChange"] > 1), "Up",
                        np.where((de["padj"] < 0.05) & (de["log2FoldChange"] < -1), "Down", "NS"))
    de = de.sort_values(["padj", "pvalue", "log2FoldChange"], ascending=[True, True, False])
    de.to_csv(OUT_DIR / "DEG_results_all.csv", index=False, encoding="utf-8-sig")

    universe = de["symbol"].dropna().astype(str).tolist()
    sig = de[de["direction"].isin(["Up", "Down"])]
    if len(sig) < 20:
        sig = pd.concat([
            de.sort_values("log2FoldChange", ascending=False).head(150),
            de.sort_values("log2FoldChange", ascending=True).head(150),
        ]).drop_duplicates("gene_id")
    sig_genes = sig["symbol"].dropna().astype(str).tolist()

    go_sets = fetch_enrichr_library(GO_LIBRARY)
    kegg_sets = fetch_enrichr_library(KEGG_LIBRARY)
    go_enr = ora(sig_genes, universe, go_sets, "GO_BP")
    kegg_enr = ora(sig_genes, universe, kegg_sets, "KEGG")
    go_enr.to_csv(OUT_DIR / "GO_BP_enrichment.csv", index=False, encoding="utf-8-sig")
    kegg_enr.to_csv(OUT_DIR / "KEGG_enrichment.csv", index=False, encoding="utf-8-sig")

    gsea_go = run_gsea(de, go_sets, nperm=50)
    gsea_kegg = run_gsea(de, kegg_sets, nperm=50)
    gsea = pd.concat([gsea_go.assign(category="GO_BP"), gsea_kegg.assign(category="KEGG")], ignore_index=True)
    gsea = gsea.sort_values(["padj", "NES"], ascending=[True, False])
    gsea.to_csv(OUT_DIR / "GSEA_preranked_results.csv", index=False, encoding="utf-8-sig")

    gsva_terms = list(pd.concat([go_enr.head(6), kegg_enr.head(6)])["term"].drop_duplicates())
    all_sets = dict(go_sets)
    all_sets.update(kegg_sets)
    gsva = ssGSEA_like(logcpm, de_raw.assign(symbol=de["symbol"].values), all_sets, gsva_terms)
    gsva.to_csv(OUT_DIR / "GSVA_ssGSEA_like_scores.csv", encoding="utf-8-sig")

    top10 = de.head(10)
    lc = logcpm.copy()
    lc.index = de_raw["gene_id"].values
    heat = lc.loc[top10["gene_id"]].to_numpy()
    heat_z = (heat - heat.mean(axis=1, keepdims=True)) / (heat.std(axis=1, keepdims=True) + 1e-9)

    draw_pca(logcpm, meta, OUT_DIR / "PCA_normalized_expression.png")
    draw_heatmap(heat_z, top10["symbol"].tolist(), sample_cols, OUT_DIR / "heatmap_top10_DEG.png", "Top 10 DEGs heatmap")
    draw_volcano(de, OUT_DIR / "volcano_top10_labeled.png")
    draw_barplot(go_enr, OUT_DIR / "GO_BP_barplot_top_terms_top8_genes.png", "GO BP enrichment")
    draw_barplot(kegg_enr, OUT_DIR / "KEGG_barplot_top_terms_top8_genes.png", "KEGG enrichment")
    draw_gsea(gsea, OUT_DIR / "GSEA_top_pathways.png")
    draw_top10_classic_gsea(gsea, de, all_sets, OUT_DIR / "GSEA_classic_top10")
    if not gsva.empty:
        vals = gsva.to_numpy()
        vals_z = (vals - vals.mean(axis=1, keepdims=True)) / (vals.std(axis=1, keepdims=True) + 1e-9)
        draw_heatmap(vals_z, list(gsva.index), list(gsva.columns), OUT_DIR / "GSVA_ssGSEA_like_heatmap.png", "GSVA / ssGSEA-like pathway activity")

    summary = {
        "input": str(INPUT_XLSX),
        "organism_assumption": "Input gene IDs are ENSG human Ensembl IDs; analysis used human annotation and human GO/KEGG gene sets.",
        "samples": sample_cols,
        "size_factors": dict(zip(sample_cols, [float(x) for x in size_factor])),
        "tested_genes": int(len(de)),
        "significant_up": int((de["direction"] == "Up").sum()),
        "significant_down": int((de["direction"] == "Down").sum()),
        "go_terms": int(len(go_enr)),
        "kegg_terms": int(len(kegg_enr)),
        "gsea_terms": int(len(gsea)),
        "gsva_terms": int(len(gsva)),
        "de_note": "DE statistics are logCPM Welch/Wald-style z tests with BH FDR because DESeq2/edgeR were not available in this runtime; sample sizes are Control n=3 and IGFL1 n=2, so dispersion/FDR estimates are exploratory.",
    }
    (OUT_DIR / "analysis_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
