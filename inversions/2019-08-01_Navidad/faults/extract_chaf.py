"""Extrae de CHAF (Melnick, Maldonado & Contreras 2020, PANGAEA doi:10.1594/PANGAEA.922241,
CC-BY 4.0) las fallas de la zona de Navidad a chaf_navidad.geojson.

Uso: python extract_chaf.py <doc.kml del CHAF_Pangaea_v1.kmz>
"""
import json
import re
import sys
from pathlib import Path

BOX = (-73.3, -70.9, -35.3, -33.4)  # lon min, lon max, lat min, lat max
FIELDS = ("F_id", "F_system", "F_name", "FT_name", "type", "strike", "dip", "dipdir", "sense",
          "length_km", "max_z_km", "width_km", "activity", "recent_act", "ass_seism", "refs")


def main(kml):
    s = Path(kml).read_text(encoding="utf-8", errors="ignore")
    feats = []
    for p in re.findall(r"<Placemark>(.*?)</Placemark>", s, re.S):
        c = re.search(r"<coordinates>(.*?)</coordinates>", p, re.S)
        if not c:
            continue
        pts = [[float(v) for v in t.split(",")[:2]] for t in c.group(1).split()]
        if not any(BOX[0] < lo < BOX[1] and BOX[2] < la < BOX[3] for lo, la in pts):
            continue
        props = {k: v.strip() for k, v in re.findall(r'<SimpleData name="(\w+)">(.*?)</SimpleData>', p, re.S) if k in FIELDS}
        feats.append({"type": "Feature", "properties": props, "geometry": {"type": "LineString", "coordinates": pts}})
    out = {"type": "FeatureCollection", "name": "CHAF_Navidad",
           "attribution": "CHAF v1: Melnick, Maldonado & Contreras (2020), PANGAEA doi:10.1594/PANGAEA.922241, CC-BY 4.0",
           "features": feats}
    Path(__file__).with_name("chaf_navidad.geojson").write_text(json.dumps(out, ensure_ascii=False, indent=1))
    for f in feats:
        p = f["properties"]
        print(f"{p.get('F_name') or '-':<22} {p.get('FT_name',''):<12} {p.get('activity',''):<9} {p.get('sense',''):<8} "
              f"strike {p.get('strike','')}/{p.get('dip','')} {p.get('dipdir','')}  n={len(f['geometry']['coordinates'])}")


if __name__ == "__main__":
    main(sys.argv[1])
