# routes/textmining.py
import io
import os
import time
import base64
from typing import List, Dict
from flask import Blueprint, request, jsonify
import pandas as pd

# Matplotlib headless backend
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from wordcloud import WordCloud
from lxml import etree
from Bio import Entrez

textmining_blueprint = Blueprint("textmining", __name__)

# -------------------- NCBI / PubMed config --------------------
ENTREZ_EMAIL = os.getenv("ENTREZ_EMAIL", "")
ENTREZ_KEY = os.getenv("ENTREZ_KEY")
if ENTREZ_EMAIL:
    Entrez.email = ENTREZ_EMAIL
else:
    print("[textmining] WARNING: ENTREZ_EMAIL not set; set it for NCBI compliance.")
if ENTREZ_KEY:
    Entrez.api_key = ENTREZ_KEY

RATE_SLEEP = float(os.getenv("PUBMED_SLEEP_SEC", "0.34"))
TOP_N_PER_GENE = int(os.getenv("PUBMED_TOPN_PER_GENE", "3"))
TOP_N_GENES_DEFAULT = int(os.getenv("PUBMED_TOP_N_GENES", "20"))  # barplot only

# -------------------- Helpers --------------------
def _clean_genes(text: str) -> List[str]:
    """Split by newlines/commas/tabs/semicolons, trim, uppercase, unique (preserve order)."""
    text = text.replace("\r", "\n")
    for ch in [",", ";", "\t"]:
        text = text.replace(ch, "\n")
    raw = [s.strip() for s in text.split("\n") if s.strip()]
    seen = set()
    out = []
    for g in raw:
        u = g.upper()
        if u not in seen:
            seen.add(u)
            out.append(u)
    return out

def _read_genes_file(fs) -> List[str]:
    blob = fs.read()
    try:
        txt = blob.decode("utf-8", errors="ignore")
    except Exception:
        txt = blob.decode("latin1", errors="ignore")
    genes = _clean_genes(txt)
    if not genes:
        raise ValueError("Gene list is empty after cleaning.")
    return genes

def _build_query(gene: str, mode: str, term: str) -> str:
    if mode == "gene_term" and term.strip():
        return f"{gene} AND ({term.strip()})"
    return gene

def _esearch_count(term: str) -> int:
    h = Entrez.esearch(db="pubmed", term=term, retmode="xml", retmax=0)
    rec = Entrez.read(h)
    h.close()
    return int(rec["Count"])

def _esearch_pmids(term: str, retmax: int) -> List[str]:
    h = Entrez.esearch(db="pubmed", term=term, retmode="xml", retmax=retmax)
    rec = Entrez.read(h)
    h.close()
    return list(rec.get("IdList", []))

def _esummary(pmids: List[str]) -> List[Dict]:
    if not pmids:
        return []
    h = Entrez.esummary(db="pubmed", id=",".join(pmids), retmode="xml")
    rec = Entrez.read(h)
    h.close()
    out = []
    for item in rec:
        title = item.get("Title")
        journal = item.get("FullJournalName") or item.get("Source")
        pubdate = str(item.get("PubDate") or "")
        year = pubdate  # keep raw as requested
        out.append({"Title": title, "Journal": journal, "Year": year})
    return out

def _efetch_snippets(pmids: List[str]) -> List[str]:
    """(Not used; kept for future abstract snippets if desired.)"""
    if not pmids:
        return []
    h = Entrez.efetch(db="pubmed", id=",".join(pmids), rettype="abstract", retmode="xml")
    xml = h.read()
    h.close()
    root = etree.fromstring(xml.encode() if isinstance(xml, str) else xml)
    arts = root.findall(".//PubmedArticle")
    out = []
    for art in arts[:len(pmids)]:
        texts = ["".join(t.itertext()).strip() for t in art.findall(".//Abstract/AbstractText")]
        txt = " ".join([t for t in texts if t]).strip()
        txt = " ".join(txt.split())
        out.append(txt[:300] if txt else None)
    if len(out) < len(pmids):
        out += [None] * (len(pmids) - len(out))
    return out[:len(pmids)]

def _pmid_url(pmid: str) -> str:
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

def _df_to_csv_b64(df: pd.DataFrame) -> str:
    buf = io.StringIO()
    df.to_csv(buf, index=False)
    return base64.b64encode(buf.getvalue().encode("utf-8")).decode("utf-8")

def _fig_to_b64(fig, fmt="png") -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, bbox_inches="tight", dpi=300)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")

def _wordcloud_from_freq(freq_map: Dict[str, int]):
    wc = WordCloud(
        width=1200, height=800, background_color="white",
        prefer_horizontal=0.85, random_state=42, max_words=200,
        relative_scaling=0.5, colormap="tab10"
    )
    wc.generate_from_frequencies(freq_map)
    return wc.to_image()

def _pil_png_b64(img) -> str:
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")

def _image_pdf_b64(img) -> str:
    fig = plt.figure(figsize=(11, 8.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img)
    ax.axis("off")
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")

# -------------------- Route --------------------
@textmining_blueprint.route("", methods=["POST"])
def textmining():
    try:
        # --- Inputs ---
        if "genes" not in request.files:
            return jsonify({"success": False, "error": "Missing 'genes' file."}), 400
        f = request.files["genes"]

        mode = (request.form.get("mode") or "gene_only").strip()
        if mode not in ("gene_only", "gene_term"):
            return jsonify({"success": False, "error": "Invalid mode. Use 'gene_only' or 'gene_term'."}), 400
        term = (request.form.get("term") or "").strip()

        # Optional: allow frontend to pass top_n_genes for BARPLOT; default to env 20
        try:
            top_n_genes = int(request.form.get("top_n_genes", TOP_N_GENES_DEFAULT))
        except ValueError:
            top_n_genes = TOP_N_GENES_DEFAULT

        genes = _read_genes_file(f)

        # --- 1) Hit counts for ALL genes ---
        rows = []
        for g in genes:
            q = _build_query(g, mode, term)
            try:
                n = _esearch_count(q)
            except Exception:
                n = 0
            rows.append({"Gene": g, "PubMed_Hits": int(n)})
            time.sleep(RATE_SLEEP)

        hit_df = pd.DataFrame(rows)  # ALL genes

        # --- 2) BARPLOT uses TOP-N only ---
        hit_df_bar = hit_df.sort_values("PubMed_Hits", ascending=False).head(top_n_genes)

        # Titles
        bar_title = (
            f"Top {top_n_genes} PubMed Hits: Gene AND {term}"
            if (mode == "gene_term" and term) else
            f"Top {top_n_genes} Genes by PubMed Hits"
        )

        # Dynamic height for readable labels
        plot_height = max(6, min(40, len(hit_df_bar) * 0.33))

        def _barplot_base64(fmt: str):
            fig, ax = plt.subplots(figsize=(9, plot_height))
            # reverse to put largest at top
            ax.barh(list(hit_df_bar["Gene"][::-1]), list(hit_df_bar["PubMed_Hits"][::-1]))
            ax.set_xlabel("Number of Articles")
            ax.set_title(bar_title)
            ax.margins(y=0.02)
            plt.tight_layout()
            return _fig_to_b64(fig, fmt=fmt)

        barplot_png_base64 = _barplot_base64("png")
        barplot_pdf_base64 = _barplot_base64("pdf")

        # --- 3) WORD CLOUD uses ALL genes ---
        freq_map_wc = {r.Gene: int(r.PubMed_Hits) for r in hit_df.itertuples(index=False) if r.PubMed_Hits > 0}
        if not freq_map_wc:
            img_wc = WordCloud(width=1200, height=800, background_color="white").generate("NO_DATA").to_image()
        else:
            img_wc = _wordcloud_from_freq(freq_map_wc)
        wordcloud_png_base64 = _pil_png_b64(img_wc)
        wordcloud_pdf_base64 = _image_pdf_b64(img_wc)

        # --- 4) Per-gene summaries (ALL genes; your requested fields) ---
        summaries_rows = []
        for g in genes:
            q = _build_query(g, mode, term)
            try:
                pmids = _esearch_pmids(q, retmax=min(TOP_N_PER_GENE, 50))
            except Exception:
                pmids = []

            meta = []
            if pmids:
                try:
                    meta = _esummary(pmids)  # Title, Journal, Year
                except Exception:
                    meta = []

            for i, pmid in enumerate(pmids):
                m = meta[i] if i < len(meta) else {}
                summaries_rows.append({
                    "Gene": g,
                    "PMID": pmid,
                    "Title": m.get("Title"),
                    "Journal": m.get("Journal"),
                    "Year": m.get("Year"),
                    "Score": (TOP_N_PER_GENE - i),   # rank-based score
                    "URL": _pmid_url(pmid)
                })
            time.sleep(RATE_SLEEP)

        summaries_df = pd.DataFrame(
            summaries_rows,
            columns=["Gene", "PMID", "Title", "Journal", "Year", "Score", "URL"]
        )

        # --- 5) CSVs (ALL genes) ---
        hitcount_csv_base64 = _df_to_csv_b64(hit_df)
        summary_csv_base64 = _df_to_csv_b64(
            summaries_df if not summaries_df.empty else pd.DataFrame(
                columns=["Gene", "PMID", "Title", "Journal", "Year", "Score", "URL"]
            )
        )

        # --- 6) Payload ---
        payload = {
            "success": True,
            "params": {
                "mode": mode,
                "term": term,
                "top_n_per_gene": TOP_N_PER_GENE,
                "top_n_genes": top_n_genes
            },

            # plots
            "barplot_png_base64": barplot_png_base64,   # Top-N only
            "barplot_pdf_base64": barplot_pdf_base64,   # Top-N only
            "wordcloud_png_base64": wordcloud_png_base64,  # ALL genes
            "wordcloud_pdf_base64": wordcloud_pdf_base64,  # ALL genes

            # CSVs (ALL genes)
            "hitcount_csv_base64": hitcount_csv_base64,
            "summary_csv_base64": summary_csv_base64,

            # Inline tables (ALL genes)
            "table": hit_df.to_dict(orient="records"),
            "summaries": summaries_df.to_dict(orient="records"),

            # Downloads object kept for React fallback logic
            "downloads": {
                "barplot_png": None,
                "barplot_pdf": None,
                "wordcloud_png": None,
                "wordcloud_pdf": None,
                "hitcount_csv": None,
                "summary_csv": None
            }
        }

        return jsonify(payload), 200

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500
