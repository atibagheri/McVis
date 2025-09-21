# backend/routes/keggmap.py
from flask import Blueprint, request, jsonify
import os, tempfile, requests, shutil
from werkzeug.utils import secure_filename

keggmap_blueprint = Blueprint("keggmap", __name__)

# Talk to the R container by service name
R_BASE = os.getenv("R_SERVICE_URL") or os.getenv("R_API_BASE") or "http://rapi:8000"
PLUMBER_URL = f"{R_BASE.rstrip('/')}/keggmap"

REQUEST_TIMEOUT = (10, 300)  # (connect, read) seconds
ALLOWED_MODES = {"genelist", "genelist_fc"}
ALLOWED_SPECIES = {"mouse", "human"}

# Save uploads where BOTH containers can see them
SHARED_DIR = os.getenv("SHARED_DIR", "/shared")
os.makedirs(SHARED_DIR, exist_ok=True)

@keggmap_blueprint.route("", methods=["POST"])
def run_keggmap():
    print("✅ Incoming /api/keggmap")
    mode = (request.form.get("mode") or "").strip()
    species = (request.form.get("species") or "").strip()
    file = request.files.get("file")

    # Basic validation
    if not file:
        return jsonify(success=False, error="Missing file"), 400
    if mode not in ALLOWED_MODES:
        return jsonify(success=False, error=f"Invalid mode '{mode}'. Use one of {sorted(ALLOWED_MODES)}"), 400
    if species not in ALLOWED_SPECIES:
        return jsonify(success=False, error=f"Invalid species '{species}'. Use one of {sorted(ALLOWED_SPECIES)}"), 400

    # Optional: tiny guard on empty uploads
    file.stream.seek(0, os.SEEK_END)
    size = file.stream.tell()
    file.stream.seek(0)
    if size == 0:
        return jsonify(success=False, error="Uploaded file is empty"), 400

    # Save to shared volume (NOT /tmp)
    tmpdir = tempfile.mkdtemp(prefix="kegg_", dir=SHARED_DIR)
    try:
        fname = secure_filename(file.filename) or "upload.txt"
        input_path = os.path.abspath(os.path.join(tmpdir, fname))
        file.save(input_path)
        print(f"📂 Saved to shared: {input_path}")

        payload = {"mode": mode, "species": species, "file_path": input_path}
        print(f"📤 POST → {PLUMBER_URL} payload={payload}")

        try:
            r = requests.post(PLUMBER_URL, json=payload, timeout=REQUEST_TIMEOUT)
        except requests.exceptions.Timeout:
            return jsonify(success=False, error="KEGG service timed out"), 504
        except requests.exceptions.ConnectionError as ce:
            return jsonify(success=False, error=f"Cannot reach KEGG service: {ce}"), 502
        except Exception as e:
            return jsonify(success=False, error=f"Unexpected error calling KEGG: {e}"), 500

        print(f"📥 R status: {r.status_code}")

        # Try to parse JSON either way to expose plumber errors
        try:
            data = r.json()
        except ValueError:
            # Non-JSON from R
            return jsonify(success=False, error="Non-JSON from KEGG service", detail=r.text), 502

        # Propagate error details clearly
        if r.status_code != 200 or not data.get("success", False):
            status = 500 if r.status_code >= 500 else 400
            return jsonify(success=False, error=data.get("error", "KEGG error"), details=data), status

        # Pass through all base64 outputs, filenames, etc.
        return jsonify(data), 200

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
        print(f"🧹 Cleaned: {tmpdir}")
