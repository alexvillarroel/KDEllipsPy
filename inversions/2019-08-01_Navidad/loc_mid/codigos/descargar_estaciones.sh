#!/bin/bash

# ============================================================
# Descarga de registros sismicos - Evento 2019-08-01 18:28:03 UTC
# Magnitud 6.8 (Mww, USGS) | Lat: -34.28 | Lon: -72.51 | Prof: 13-25 km
# 96 km SW de San Antonio, Chile
# ============================================================

EVENT_ID="168a11de801928c2da4f589b6d00e63a"
BASE_URL="https://evtdb.csn.uchile.cl/write/${EVENT_ID}"
OUTPUT_DIR="./registros_sismicos"
ZIP_DIR="${OUTPUT_DIR}/zips"

# Lista completa de estaciones (76 estaciones, desde evtdb)
ESTACIONES=(
    B01I B02I B02O B03O B04O B05I B05O B06I B06O B07O B08O B09O B10O B13I B20I B21I
    BO01 BO02 BO03
    LMEL
    M01L M03L M05L M09L M10L M12L M13L M14L M17L M19L
    ML02
    MT01 MT05 MT09 MT10 MT14 MT15 MT18
    R01M R02M R03M R04M R05M R06M R07M R08M R09M R10M R11M R12M R13M R14M R16M R17M R18M R19M R20M R23M R26M
    V02A V03A V04A V05A V06A V07A V08A V12A V14A V15A V16A V18A V19A V21A V23A
    VA03 VA05
)

# Eliminar duplicados manteniendo orden
ESTACIONES_UNICAS=($(echo "${ESTACIONES[@]}" | tr ' ' '\n' | awk '!seen[$0]++'))

# Crear directorios
mkdir -p "$OUTPUT_DIR"
mkdir -p "$ZIP_DIR"

echo "============================================"
echo " CSN - Descarga de estaciones sismicas"
echo " Evento: ${EVENT_ID}"
echo " Total estaciones unicas: ${#ESTACIONES_UNICAS[@]}"
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

    wget -q -O "$ZIP_FILE" "$URL"

    if [ $? -eq 0 ] && [ -s "$ZIP_FILE" ]; then
        ESTACION_DIR="${OUTPUT_DIR}/${ESTACION}"
        mkdir -p "$ESTACION_DIR"

        unzip -q -o "$ZIP_FILE" -d "$ESTACION_DIR"

        if [ $? -eq 0 ]; then
            echo "    OK ${ESTACION} descomprimido en ${ESTACION_DIR}/"
        else
            echo "    ERROR al descomprimir ${ESTACION}"
            ERRORES+=("$ESTACION (unzip fallo)")
        fi
    else
        echo "    ERROR al descargar ${ESTACION}"
        ERRORES+=("$ESTACION (wget fallo)")
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
