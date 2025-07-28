
from flask import Blueprint, request, jsonify, send_from_directory
import os
import tempfile
import base64
import json
import requests

enrichment_blueprint = Blueprint("enrichment", __name__)

@enrichment_blueprint.route("", methods=["POST"])
def enrichment():
    file_storages = request.files.getlist("file_path")
    saved_paths = []

    # Save each uploaded file to a temp directory
    for file in file_storages:
        temp_dir = tempfile.gettempdir()
        save_path = os.path.join(temp_dir, file.filename)
        file.save(save_path)
        saved_paths.append(save_path)

    print("Enrichment saved_paths:", saved_paths)

    species = request.form.get("species")
    print("Enrichment species:", species)
    
    # URL of your R Plumber enrichment endpoint
    url = "http://localhost:8000/enrichment"
    headers = {"Content-Type": "application/json"}
    data = {"file_path": saved_paths[0], "species": species}

    # Call the R API
    try:
        # response = requests.post(url, headers=headers, data=json.dumps(data))
        # response_json = response.json()
        response_json = {'success': [True], 
            'go_pdf': ['/Users/student/Desktop/AppDEG/backend/degviz_api/Enrichment/go_enrichment_plot.pdf'], 
            'kegg_pdf': ['/Users/student/Desktop/AppDEG/backend/degviz_api/Enrichment/kegg_enrichment_plot.pdf'],
            'go_txt': ['/Users/student/Desktop/AppDEG/backend/degviz_api/Enrichment/go_enrichment_results.txt'],
            'kegg_txt': ['/Users/student/Desktop/AppDEG/backend/degviz_api/Enrichment/kegg_enrichment_results.txt']}
        print("Enrichment response_json:", response_json)
        return jsonify({"go_pdf": response_json.get("go_pdf"), 
                                "kegg_pdf": response_json.get("kegg_pdf"),
                                "go_txt": response_json.get("go_txt"),
                                "kegg_txt": response_json.get("kegg_txt")})


        # if response.status_code == 200:
        #     return jsonify({"go_pdf": response_json.get("go_pdf"), 
        #                         "kegg_pdf": response_json.get("kegg_pdf"),
        #                         "go_txt": response_json.get("go_txt"),
        #                         "kegg_txt": response_json.get("kegg_txt")})
        # else:
        #     print(f"Error: {response.status_code}")
        #     print(response.text)
        #     return jsonify({"error": response.text}), response.status_code

    except requests.exceptions.RequestException as e:
        print("Error communicating with R API:", str(e))
        return jsonify({"error": str(e)}), 500
