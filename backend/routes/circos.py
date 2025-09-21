# backend/routes/circos.py
from flask import Blueprint, request, jsonify
import tempfile, os, json, requests, shutil
from werkzeug.utils import secure_filename

circos_blueprint = Blueprint("circos", __name__)

# Talk to the R container by service name
R_BASE = os.environ.get("R_SERVICE_URL") or os.environ.get("R_API_BASE") or "http://rapi:8000"
PLUMBER_URL = f"{R_BASE.rstrip('/')}/circos"

# Save uploads where BOTH containers can see them
SHARED_DIR = os.environ.get("SHARED_DIR", "/shared")
os.makedirs(SHARED_DIR, exist_ok=True)

REQUEST_TIMEOUT = (10, 180)

def _extract_labels_from_form(form, n_files: int):
    labels = form.getlist("labels")
    if not labels and form.get("labels_json"):
        try:
            labels = json.loads(form.get("labels_json"))
        except Exception:
            labels = []
    if labels and len(labels) == n_files:
        return labels
    return None  # let R derive

@circos_blueprint.route("", methods=["POST"])
def circos():
    try:
        files = request.files.getlist("files")
        if len(files) < 2:
            return jsonify({"error": "Please upload at least 2 files"}), 400

        gene_col_user = request.form.get("gene_col_user") or None
        try:
            min_shared = int(request.form.get("min_shared", "1"))
        except ValueError:
            return jsonify({"error": "min_shared must be an integer"}), 400
        labels = _extract_labels_from_form(request.form, len(files))

        # Save under shared mount so R can read the same absolute paths
        tmpdir = tempfile.mkdtemp(prefix="circos_", dir=SHARED_DIR)
        saved_paths = []
        try:
            for f in files:
                fname = secure_filename(f.filename) or "upload.txt"
                path = os.path.join(tmpdir, fname)
                f.save(path)
                saved_paths.append(os.path.abspath(path))

            payload = {"paths": saved_paths, "min_shared": min_shared}
            if gene_col_user:
                payload["gene_col_user"] = gene_col_user
            if labels:
                payload["labels"] = labels

            r = requests.post(PLUMBER_URL, json=payload, timeout=REQUEST_TIMEOUT)
        finally:
            # remove temp dir after R finishes (comment out while debugging)
            shutil.rmtree(tmpdir, ignore_errors=True)

        if r.status_code != 200:
            # pass through plumber error details if JSON, else raw text
            try:
                return jsonify(r.json()), r.status_code
            except Exception:
                return jsonify({"error": "Plumber error", "details": r.text}), r.status_code

        data = r.json()
        if not data.get("success", False):
            # R uses 422 for “no links to plot”; propagate with diagnostics
            return jsonify(data), 422

        return jsonify({
            "png_base64": data.get("png_base64"),
            "pdf_base64": data.get("pdf_base64")
        }), 200

    except Exception as e:
        return jsonify({"error": str(e)}), 500
