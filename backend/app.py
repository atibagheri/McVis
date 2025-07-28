from flask import Flask, request, jsonify
from flask_cors import CORS
from routes.pca import pca_blueprint
from routes.venn_upset import venn_upset_blueprint
from routes.enrichment import enrichment_blueprint
from routes.gomap import gomap_blueprint
from routes.keggmap import keggmap_blueprint
from flask import send_from_directory
import os
import requests
app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})

# Register blueprints
app.register_blueprint(pca_blueprint, url_prefix="/api/pca")
app.register_blueprint(venn_upset_blueprint, url_prefix="/api/venn-upset")
app.register_blueprint(enrichment_blueprint, url_prefix="/api/enrichment")
app.register_blueprint(gomap_blueprint, url_prefix="/api/gomap")
app.register_blueprint(keggmap_blueprint, url_prefix="/api/keggmap")


# Serve static files (barplot images, zip files, etc.)
@app.route("/downloads/gomap/<path:filename>")
def download_gomap_file(filename):
    output_dir = os.path.join(os.getcwd(), "./degviz_api/gomap_output")
    return send_from_directory(output_dir, filename, as_attachment=True)

@app.route("/downloads/kegg/<path:filename>")
def download_kegg_file(filename):
    output_dir = os.path.join(os.getcwd(), "./degviz_api/Keggmap_output")
    return send_from_directory(output_dir, filename, as_attachment=True)

@app.after_request
def apply_cors_headers(response):
    response.headers["Access-Control-Allow-Origin"] = "*"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type,Authorization"
    response.headers["Access-Control-Allow-Methods"] = "GET,POST,PUT,DELETE,OPTIONS"
    return response

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5050, debug=True)
