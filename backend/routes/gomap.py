# backend/routes/gomap.py
from flask import Blueprint, request, jsonify
import os, tempfile, requests, shutil
from werkzeug.utils import secure_filename

gomap_blueprint = Blueprint("gomap", __name__)

# ✅ Build plumber URL from R base (container name), not localhost
R_BASE = os.getenv("R_SERVICE_URL") or os.getenv("R_API_BASE") or "http://rapi:8000"
PLUMBER_URL = f"{R_BASE.rstrip('/')}/gomap"

REQUEST_TIMEOUT = (10, 180)  # (connect, read)
ALLOWED_MODES = {"genelist", "genelist_fc"}
ALLOWED_SPECIES = {"mouse", "human"}

# ✅ Shared dir path (mounted in compose in both backend and rapi)
SHARED_DIR = os.getenv("SHARED_DIR", "/shared")
os.makedirs(SHARED_DIR, exist_ok=True)

@gomap_blueprint.route("", methods=["POST"])
def run_gomap():
    print("✅ Incoming request to /api/gomap")
    mode = (request.form.get("mode") or "").strip()
    species = (request.form.get("species") or "").strip()
    file = request.files.get("file")

    if not file or not mode or not species:
        return jsonify({"success": False, "error": "Missing parameters (file, mode, species)"}), 400
    if mode not in ALLOWED_MODES:
        return jsonify({"success": False, "error": f"Invalid mode '{mode}'. Use one of {sorted(ALLOWED_MODES)}"}), 400
    if species not in ALLOWED_SPECIES:
        return jsonify({"success": False, "error": f"Invalid species '{species}'. Use one of {sorted(ALLOWED_SPECIES)}"}), 400

    # ✅ Save into the shared mount, not /tmp
    tmpdir = tempfile.mkdtemp(prefix="gomap_", dir=SHARED_DIR)
    try:
        filename = secure_filename(file.filename) or "upload.txt"
        input_path = os.path.abspath(os.path.join(tmpdir, filename))
        file.save(input_path)
        print(f"📂 Saved uploaded file to (shared): {input_path}")

        payload = {"mode": mode, "species": species, "file_path": input_path}
        print("📦 Payload to R Plumber:", payload)

        try:
            r = requests.post(PLUMBER_URL, json=payload, timeout=REQUEST_TIMEOUT)
        except requests.exceptions.Timeout:
            return jsonify({"success": False, "error": "GOMap service timed out"}), 504
        except requests.exceptions.ConnectionError as ce:
            return jsonify({"success": False, "error": f"Cannot reach GOMap service: {ce}"}), 502
        except Exception as e:
            return jsonify({"success": False, "error": f"Unexpected error calling GOMap: {e}"}), 500

        print("🌐 Status from R Plumber:", r.status_code)

        # Try to parse JSON (surface plumber errors nicely)
        try:
            result_json = r.json()
            print("Keys from plumber:", list(result_json.keys()))
            print("zip_base64 length:", len(result_json.get("zip_base64","")))
        except ValueError:
            if r.status_code != 200:
                return jsonify({"success": False, "error": r.text or "Non-JSON error from GOMap"}), 502
            return jsonify({"success": False, "error": "GOMap returned non-JSON response"}), 502

        if r.status_code != 200 or not result_json.get("success", False):
            status = 500 if r.status_code >= 500 else 400
            return jsonify({
                "success": False,
                "error": result_json.get("error", "GOMap reported an error"),
                "details": result_json
            }), status

        return jsonify(result_json), 200

    finally:
        # You may comment this during debugging to inspect files
        shutil.rmtree(tmpdir, ignore_errors=True)
        print(f"🧹 Cleaned up temp dir: {tmpdir}")
