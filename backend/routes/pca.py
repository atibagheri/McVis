import io
import base64
import numpy as np
import pandas as pd
import plotly.express as px
from flask import Blueprint, request, jsonify
from sklearn.decomposition import PCA

pca_blueprint = Blueprint('pca', __name__)

@pca_blueprint.route('/', methods=['POST'])  # ✅ correct
def run_pca():
    expr_file = request.files['expression_file']
    sample_file = request.files['sample_file']

    expr_df = pd.read_csv(expr_file, sep="\t", index_col=0)
    sample_df = pd.read_csv(sample_file, sep="\t", index_col=0)

    # Align columns with sample names
    expr_df = expr_df.loc[:, sample_df.index]

    # Convert to CPM
    cpm = expr_df.div(expr_df.sum(axis=0), axis=1) * 1e6
    log_cpm = np.log2(cpm + 0.5).T

    # PCA
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(log_cpm)
    pca_df = pd.DataFrame(pcs, columns=["PC1", "PC2"], index=log_cpm.index)
    final_df = pca_df.join(sample_df)

    var_exp = pca.explained_variance_ratio_ * 100
    pc1_var, pc2_var = round(var_exp[0], 1), round(var_exp[1], 1)

    fig = px.scatter(
        final_df,
        x="PC1", y="PC2",
        color="Condition", symbol="Condition", size="Day",
        labels={"PC1": f"PC1 ({pc1_var}%)", "PC2": f"PC2 ({pc2_var}%)"},
        title="PCA Plot"
    )

    # Encode PNG
    buf_png = io.BytesIO()
    fig.write_image(buf_png, format="png")
    buf_png.seek(0)
    img_base64 = base64.b64encode(buf_png.read()).decode('utf-8')

    # Encode PDF
    buf_pdf = io.BytesIO()
    fig.write_image(buf_pdf, format="pdf")
    buf_pdf.seek(0)
    pdf_base64 = base64.b64encode(buf_pdf.read()).decode('utf-8')

    return jsonify({"png": img_base64, "pdf": pdf_base64})
