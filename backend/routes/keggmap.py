from flask import Blueprint, request, jsonify, send_from_directory
import os, tempfile, requests, subprocess
from werkzeug.utils import secure_filename
from pprint import pprint

keggmap_blueprint = Blueprint("keggmap", __name__)
PLUMBER_URL = "http://localhost:8000/keggmap"  # your R plumber endpoint

@keggmap_blueprint.route("", methods=["POST"])
def run_keggmap():
    print("✅ Incoming request to /api/keggmap")
    print("🔎 DEBUG: request.form:", request.form)
    print("🔎 DEBUG: request.files:", request.files)

    mode = request.form.get("mode")
    species = request.form.get("species")
    file = request.files.get("file")

    if not file or not mode or not species:
        return jsonify({"success": False, "error": "Missing parameters"}), 400

    # Save the uploaded file    
    
    # Save file to a temporary directory
    temp_dir = tempfile.mkdtemp()
    filename = secure_filename(file.filename)
    input_path = os.path.join(temp_dir, filename)
    file.save(input_path)


    payload = {
        "mode": mode,
        "species": species,
        "file_path": input_path,
    }
    print("📦 Payload to R:", payload)



    print("DEBUG mode:", mode)
    print("DEBUG species:", species)
    print("DEBUG input_path:", input_path)
    try:
        r = requests.post(PLUMBER_URL, json=payload)
        print("DEBUG plumber status:", r.status_code)
        print("DEBUG plumber text:")
        pprint(r.json())
        if r.status_code != 200:
            return jsonify({"success": False, "error": r.text}), 500
        return jsonify(r.json())
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({"success": False, "error": str(e)}), 500

    

    