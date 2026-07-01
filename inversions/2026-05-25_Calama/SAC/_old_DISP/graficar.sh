#!/bin/bash

# 1. Definimos las estaciones 
STATIONS=("A20F" "PB09" "A06F" "PB03" "PB15" "PB01" "PB19")

# 2. Coordenadas Y para las 7 filas
YMIN=(0.82 0.70 0.58 0.46 0.34 0.22 0.10)
YMAX=(0.92 0.80 0.68 0.56 0.44 0.32 0.20)

# 3. Coordenadas X para las 3 columnas
XMIN=(0.10 0.40 0.70)
XMAX=(0.35 0.65 0.95)

# Archivo temporal
MACRO="macro_temp.m"

# Iniciamos el archivo macro
echo "bg x" > $MACRO
echo "beginframe" >> $MACRO
echo "qdp off" >> $MACRO
echo "title off" >> $MACRO

# 4. Bucle para recorrer cada estación
for i in "${!STATIONS[@]}"; do
    sta="${STATIONS[$i]}"
    ymin="${YMIN[$i]}"
    ymax="${YMAX[$i]}"

    # --- COMPONENTE ESTE (E) ---
    echo "read $sta.H*E.disp.sac" >> $MACRO
    echo "rmean" >> $MACRO
    echo "div &1,depmax" >> $MACRO
    echo "xvport ${XMIN[0]} ${XMAX[0]}" >> $MACRO
    echo "yvport $ymin $ymax" >> $MACRO
    echo "p" >> $MACRO

    # --- COMPONENTE NORTE (N) ---
    echo "read $sta.H*N.disp.sac" >> $MACRO
    echo "rmean" >> $MACRO
    echo "div &1,depmax" >> $MACRO   
    echo "xvport ${XMIN[1]} ${XMAX[1]}" >> $MACRO
    echo "yvport $ymin $ymax" >> $MACRO
    echo "p" >> $MACRO

    # --- COMPONENTE VERTICAL (Z) ---
    echo "read $sta.H*Z.disp.sac" >> $MACRO
    echo "rmean" >> $MACRO
    echo "div &1,depmax" >> $MACRO   
    echo "xvport ${XMIN[2]} ${XMAX[2]}" >> $MACRO
    echo "yvport $ymin $ymax" >> $MACRO
    echo "p" >> $MACRO
done

# Terminamos la composición y SALIMOS de SAC
#echo "endframe" >> $MACRO
echo "saveimg grilla_estaciones.pdf" >> $MACRO   # Guarda la imagen en alta calidad
echo "quit" >> $MACRO

sac < $MACRO 2> /dev/null
echo "endframe" >> $MACRO
cat $MACRO - | sac 2> /dev/null
# 5. Ejecutamos SAC silenciando las advertencias de libxml2 (2> /dev/null)
sac < $MACRO 2> /dev/null

# 6. Borramos el archivo temporal
rm $MACRO
