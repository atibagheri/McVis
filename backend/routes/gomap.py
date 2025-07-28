# backend/routes/gomap.py
from flask import Blueprint, request, jsonify
import os, tempfile, requests
from werkzeug.utils import secure_filename
from pprint import pprint

# Create blueprint
gomap_blueprint = Blueprint("gomap", __name__)

# ✅ Change this to your R plumber GOMap endpoint
PLUMBER_URL = "http://localhost:8000/gomap"

@gomap_blueprint.route("", methods=["POST"])
def run_gomap():
    print("✅ Incoming request to /api/gomap")
    print("🔎 DEBUG: request.form:", request.form)
    print("🔎 DEBUG: request.files:", request.files)

    mode = request.form.get("mode")
    species = request.form.get("species")
    file = request.files.get("file")

    # --- Validate ---
    if not file or not mode or not species:
        return jsonify({"success": False, "error": "Missing parameters"}), 400

    # --- Save file to a temporary directory ---
    temp_dir = tempfile.mkdtemp()
    filename = secure_filename(file.filename)
    input_path = os.path.join(temp_dir, filename)
    file.save(input_path)

    print(f"📂 Saved uploaded file to: {input_path}")

    # --- Build payload for R Plumber ---
    payload = {
        "mode": mode,
        "species": species,
        "file_path": input_path
    }
    print("📦 Payload to R plumber:", payload)

    # --- Call the R plumber endpoint ---
    try:
        r = requests.post(PLUMBER_URL, json=payload)
        print("🌐 Status from R plumber:", r.status_code)

        if r.status_code != 200:
            print("⚠️ Non-200 response from R plumber:", r.text)
            return jsonify({"success": False, "error": r.text}), 500

        result_json = r.json()
        print("✅ Received JSON from plumber:")
        pprint(result_json)

        # ✅ Return whatever plumber returned (contains barplot_base64, zip_file, etc.)
        return jsonify(result_json)

    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({"success": False, "error": str(e)}), 500
