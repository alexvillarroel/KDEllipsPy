#!/bin/bash

# ============================================================
# Descarga de registros sísmicos - Evento 2026-05-25 21:52:18
# Magnitud 6.9 | Lat: -22.38 | Lon: -68.76 | Prof: 114 km
# ============================================================

EVENT_ID="b13d6fbb0c255df0ab8335948603c06e"
BASE_URL="https://evtdb.csn.uchile.cl/write/${EVENT_ID}"
OUTPUT_DIR="./registros_sismicos"
ZIP_DIR="${OUTPUT_DIR}/zips"

# Lista completa de estaciones (47 estaciones)
ESTACIONES=(
    A16P
    AF01
    T07A
    T09A
    T16A
    T24A
    TA01
    HMBCX
    PATCX
    PB01
    PB02
    PB03
    PB06
    PB08
    PB09
    PB15
    PB19
    PB23
    A14P
    A16P
    T09A
    T24A
    T16A
    T13A
    T11A
    T10A
    T07A
    T14A
    A09P
    A13F
    A14F
    A06F
    A08P
    A02F
    A03F
    A19F
    A20F
    A12P
    A22F
    A28F
    A18P
    A04F
)

# Eliminar duplicados manteniendo orden
ESTACIONES_UNICAS=($(echo "${ESTACIONES[@]}" | tr ' ' '\n' | awk '!seen[$0]++'))

# Crear directorios
mkdir -p "$OUTPUT_DIR"
mkdir -p "$ZIP_DIR"

echo "============================================"
echo " CSN - Descarga de estaciones sísmicas"
echo " Evento: ${EVENT_ID}"
echo " Total estaciones únicas: ${#ESTACIONES_UNICAS[@]}"
echo "============================================"
echo ""

TOTAL=${#ESTACIONES_UNICAS[@]}
COUNT=0
ERRORES=()

for ESTACION in "${ESTACIONES_UNICAS[@]}"; do
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
            echo "    ✓ ${ESTACION} descomprimido en ${ESTACION_DIR}/"
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
echo " Total estaciones: ${TOTAL}"
echo " Exitosas: $((TOTAL - ${#ERRORES[@]}))"
echo " Errores: ${#ERRORES[@]}"

if [ ${#ERRORES[@]} -gt 0 ]; then
    echo " Estaciones con error:"
    for E in "${ERRORES[@]}"; do
        echo "   - $E"
    done
fi

echo ""
echo " Archivos guardados en: ${OUTPUT_DIR}/"
echo "============================================"
