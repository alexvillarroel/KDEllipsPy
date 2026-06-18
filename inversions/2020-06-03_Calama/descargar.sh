#!/bin/bash

# ============================================================
# Descarga de registros sísmicos - Evento 2020-06-03 07:35:34
# Magnitud 6.9 | Lat: -23.25 | Lon: -68.53 | Prof: 123 km
# ============================================================

EVENT_ID="887487c5d255a4433574d2476200088b"
BASE_URL="https://evtdb.csn.uchile.cl/write/${EVENT_ID}"
OUTPUT_DIR="./registros_20200603"
ZIP_DIR="${OUTPUT_DIR}/zips"

ESTACIONES=(
    AF01
    PB01
    PB02
    PB03
    PB05
    PB06
    PB09
    PB19
    A14F
    A07F
    A11F
    T15A
    A05C
    A07C
    A14C
    A27C
    A17C
    A28C
    A19C
    T24A
    T10A
    T07A
    A16F
    A29F
    A06F
    A21F
    A01F
    A04F
    A13F
    A19F
    T13A
    T11A
    T20A
    T19A
    A28F
    A12F
    A03F
    A15F
    A22F
    A08F
    A23F
    A05F
    A09F
    A17F
    A18F
    A02F
    A24F
    A10F
)

mkdir -p "$OUTPUT_DIR"
mkdir -p "$ZIP_DIR"

echo "============================================"
echo " CSN - Descarga de estaciones sísmicas"
echo " Evento: 2020-06-03 07:35:34 | Mag 6.9"
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
