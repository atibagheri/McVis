from flask import Blueprint, request, jsonify
import os, json, tempfile, shutil, requests
from werkzeug.utils import secure_filename

venn_upset_blueprint = Blueprint("venn_upset", __name__)

R_BASE = os.getenv("R_SERVICE_URL") or os.getenv("R_API_BASE") or "http://rapi:8000"
PLUMBER_URL = f"{R_BASE.rstrip('/')}/venn-upset"
REQUEST_TIMEOUT = (10, 180)

SHARED_DIR = os.getenv("SHARED_DIR", "/shared")
os.makedirs(SHARED_DIR, exist_ok=True)

@venn_upset_blueprint.route("", methods=["POST"])
def venn_upset():
    files = request.files.getlist("files")
    if len(files) < 2:
        return jsonify({"error": "Please upload at least 2 files"}), 400

    tmpdir = tempfile.mkdtemp(prefix="venn_", dir=SHARED_DIR)
    paths = []
    try:
        for fs in files:
            fname = secure_filename(fs.filename) or "upload.txt"
            p = os.path.join(tmpdir, fname)
            fs.save(p)
            paths.append(os.path.abspath(p))

        payload = {"paths": paths}
        r = requests.post(PLUMBER_URL, json=payload, timeout=REQUEST_TIMEOUT)
        try:
            data = r.json()
        except ValueError:
            return jsonify({"error": "Non-JSON from R", "detail": r.text}), 502

        if r.status_code != 200 or not data.get("success", False):
            return jsonify({"error": data.get("error") or "R service error", "detail": data}), (
                500 if r.status_code >= 500 else 400
            )

        return jsonify({
            "png": data.get("png_base64"),
            "pdf": data.get("pdf_base64"),
        }), 200

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
