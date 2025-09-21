from flask import Flask, jsonify
from flask_cors import CORS
import os

# blueprints
from routes.transpca import transpca_blueprint
from routes.pca import pca_blueprint
from routes.venn_upset import venn_upset_blueprint
from routes.gomap import gomap_blueprint
from routes.keggmap import keggmap_blueprint
from routes.textmining import textmining_blueprint
from routes.circos import circos_blueprint

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})

app.config.update(
    R_SERVICE_URL=os.getenv("R_SERVICE_URL", "http://127.0.0.1:8000"),
    MOUSE_TPM_PATH=os.getenv(
        "MOUSE_TPM_PATH",
        "/Users/student/Desktop/AppDEG/backend/degviz_api/data/MedianNorm_TPMs_order_Mouse.txt"
    ),
)

# --- register blueprints (NO extra indentation) ---
app.register_blueprint(transpca_blueprint, url_prefix="/api/transpca")
app.register_blueprint(pca_blueprint, url_prefix="/api/pca")
app.register_blueprint(venn_upset_blueprint, url_prefix="/api/venn-upset")
app.register_blueprint(gomap_blueprint, url_prefix="/api/gomap")
app.register_blueprint(keggmap_blueprint, url_prefix="/api/keggmap")
app.register_blueprint(textmining_blueprint, url_prefix="/api/textmining")
app.register_blueprint(circos_blueprint, url_prefix="/api/circos")

@app.route("/healthz")
def healthz():
    return jsonify({"status": "ok"}), 200

# optional: extra CORS headers (CORS() already handles most cases)
@app.after_request
def apply_cors_headers(response):
    response.headers["Access-Control-Allow-Origin"] = "*"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type,Authorization"
    response.headers["Access-Control-Allow-Methods"] = "GET,POST,PUT,DELETE,OPTIONS"
    return response

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=int(os.getenv("FLASK_PORT", 5050)), debug=True)
