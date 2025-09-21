# routes/transpca.py
import base64, requests
from flask import Blueprint, current_app, request, jsonify, Response

transpca_blueprint = Blueprint("transpca_api", __name__)

def _b64(fs):
    if fs is None: return None
    fs.stream.seek(0)
    return base64.b64encode(fs.read()).decode("ascii")

@transpca_blueprint.post("/")
def transpca_run():
    f = request.files.get("file")
    if not f:
        return jsonify(status="error", message="Missing 'file'"), 400

    payload = {"file": {"name": f.filename or "user.txt", "b64": _b64(f)}}

    meta = request.files.get("meta")
    if meta:
        payload["meta"] = {"name": meta.filename or "meta.txt", "b64": _b64(meta)}

    # pass-thru optional aesthetics
    for k in ("label_col", "color_col", "size_col", "shape_col"):
        v = request.form.get(k)
        if v: payload[k] = v

    r = requests.post(
        f"{current_app.config['R_SERVICE_URL']}/transpca_json",
        json=payload, timeout=300
    )
    return (r.text, r.status_code,
            {"Content-Type": r.headers.get("Content-Type", "application/json")})

@transpca_blueprint.get("/mouse-tpm-preview")
def mouse_tpm_preview():
    r = requests.get(f"{current_app.config['R_SERVICE_URL']}/mouse-tpm-preview",
                     params=request.args, timeout=60)
    return (r.text, r.status_code, {"Content-Type": "application/json"})

@transpca_blueprint.get("/mouse-tpm-download")
def mouse_tpm_download():
    r = requests.get(f"{current_app.config['R_SERVICE_URL']}/mouse-tpm-download",
                     stream=True, timeout=300)
    resp = Response(r.iter_content(8192), status=r.status_code)
    ct = r.headers.get("Content-Type")
    cd = r.headers.get("Content-Disposition")
    if ct: resp.headers["Content-Type"] = ct
    if cd: resp.headers["Content-Disposition"] = cd
    return resp
