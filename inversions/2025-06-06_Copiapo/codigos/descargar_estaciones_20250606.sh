#!/bin/bash

# ============================================================
# Descarga de registros sísmicos - Evento 2025-06-06 17:15:05
# Magnitud 6.4 | Lat: -26.87 | Lon: -70.14 | Prof: 65 km
# ============================================================

EVENT_ID="5628aae4a0c998ed562f667cc202569b"
BASE_URL="https://evtdb.csn.uchile.cl/write/${EVENT_ID}"
OUTPUT_DIR="./registros_20250606"
ZIP_DIR="${OUTPUT_DIR}/zips"

ESTACIONES=(
    GO03
    A02C
    AC01
    AC04
    AC05
    C16O
    C28O
    PB19
    A28C
    A29C
    A30C
    A32C
    A24F
    A07C
    A08C
    A09C
    A10C
    A18C
    A19C
    A06C
    C15O
    C14O
    C03O
    A28F
    A26F
    A21F
    C18O
    C01O
    A15C
    A16C
    A17C
    A22C
    A25C
    A23F
    A05F
)

mkdir -p "$OUTPUT_DIR"
mkdir -p "$ZIP_DIR"

echo "============================================"
echo " CSN - Descarga de estaciones sísmicas"
echo " Evento: 2025-06-06 17:15:05 | Mag 6.4"
echo " Total estaciones: ${#ESTACIONES[@]}"
echo "============================================"
echo ""

TOTAL=${#ESTACIONES[@]}
COUNT=0
ERRORES=()

for ESTACION in "${ESTACIONES[@]}"; do
    COUNT=$((COUNT + 1))
    URL="${BASE_URL}/${ESTACION}"
    ZIP_FILE="${ZIP_DIR}/${ESTACION}.zip"

    echo "[${COUNT}/${TOTAL}] Descargando ${ESTACION}..."

    wget -q --show-progress -O "$ZIP_FILE" "$URL"

    if [ $? -eq 0 ] && [ -s "$ZIP_FILE" ]; then
        ESTACION_DIR="${OUTPUT_DIR}/${ESTACION}"
        mkdir -p "$ESTACION_DIR"

        unzip -q -o "$ZIP_FILE" -d "$ESTACION_DIR"

        if [ $? -eq 0 ]; then
            echo "    ✓ ${ESTACION} OK"
        else
            echo "    ✗ Error al descomprimir ${ESTACION}"
            ERRORES+=("$ESTACION (unzip falló)")
        fi
    else
        echo "    ✗ Error al descargar ${ESTACION}"
        ERRORES+=("$ESTACION (wget falló)")
        rm -f "$ZIP_FILE"
    fi
done

echo ""
echo "============================================"
echo " Resumen"
echo "============================================"
echo " Total: ${TOTAL} | OK: $((TOTAL - ${#ERRORES[@]})) | Errores: ${#ERRORES[@]}"

if [ ${#ERRORES[@]} -gt 0 ]; then
    echo " Estaciones con error:"
    for E in "${ERRORES[@]}"; do
        echo "   - $E"
    done
fi

echo " Archivos en: ${OUTPUT_DIR}/"
echo "============================================"
