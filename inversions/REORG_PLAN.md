# Plan de reorganización de `inversions/` (NO ejecutado)

Decisiones de Alex:
- Alcance: **solo plan, revisar antes de ejecutar**.
- Carpetas legacy duplicadas (`calama2020`, `copiapo-06-06-2025`) → **a `archive/`** (siguen
  en disco, fuera de git).

Principios:
- `git mv` para lo que se **conserva en git** (preserva historial).
- `mv` plano + `git rm -r --cached` para lo que pasa a **`archive/`** (queda en disco, sale de git).
- No se toca ningún `output/` canónico, `input.ctl`, `event.ctl`, ni código en `codigos/`.
- **Sin reescritura de historial** (el `.git` no adelgaza, pero deja de crecer con scratch).

## Estructura objetivo

```
inversions/
  2020-06-03_Calama/      DATA SAC codigos output egf robustez legacy legacy_equiv
                          archive/   <- output_*_killed/_prev/_sel/_t34, latex scratch
  2025-06-06_Copiapo/     DATA SAC codigos output
                          archive/   <- legacy_copiapo (viejo), output_dc_z80, *.log
  2026-05-25_Calama/      DATA SAC codigos output Dynamic_inversion
                          archive/   <- output_*_backup, output_prev_42k
  _codigo_legacy/         <- calama_legacy/ (Fortran Kinematic_inversion; NO es un evento)
  comparacion_2020_vs_2026/
```

## Paso 0 — gitignore (primero, para que archive/ no se re-trackee)

Agregar a `/home/alex/KDEllipsPy/.gitignore`:

```gitignore
# Scratch de inversiones (resultados regenerables, backups, logs)
inversions/**/archive/
inversions/**/*.log
inversions/**/*.bak
# Latex scratch
*.aux
*.fdb_latexmk
*.fls
*.synctex.gz
inversions/**/build/
```

## Paso 1 — renombrar a `YYYY-MM-DD_Lugar` (git mv, conserva historial)

```bash
cd /home/alex/KDEllipsPy
git mv inversions/2020-06-03 inversions/2020-06-03_Calama
git mv inversions/calama_legacy inversions/_codigo_legacy
# 2025-06-06_Copiapo y 2026-05-25_Calama ya cumplen el formato.
```

## Paso 2 — mover legacy duplicadas a archive/ (a disco, fuera de git)

```bash
cd /home/alex/KDEllipsPy/inversions
# Calama 2020 legacy AXITRA
mkdir -p 2020-06-03_Calama/archive
mv calama2020 2020-06-03_Calama/archive/legacy_calama2020
git rm -r --cached calama2020            # quita del índice; archivos quedan en disco

# Copiapo viejo (scripts test_*.py, ipynb pesados)
mkdir -p 2025-06-06_Copiapo/archive
mv copiapo-06-06-2025 2025-06-06_Copiapo/archive/legacy_copiapo
git rm -r --cached copiapo-06-06-2025
```

## Paso 3 — mover scratch output_* a archive/ por evento

```bash
cd /home/alex/KDEllipsPy/inversions

# --- Calama 2020 ---
cd 2020-06-03_Calama && mkdir -p archive
mv output_band002-01_killed output_band006-015_sel output_disp_t34_full \
   output_plano64_sel output_prev_002-01 output_prev_9sta_plano32 \
   output_prueba_t34_300 output_vel_t34_006-015 archive/
# latex scratch (deja .tex y .pdf en su sitio)
mv resumen_trabajo.aux resumen_trabajo.log resumen_trabajo.out \
   resumen_trabajo.fdb_latexmk resumen_trabajo.fls resumen_trabajo.synctex.gz archive/ 2>/dev/null
cd ..
git -C /home/alex/KDEllipsPy rm -r --cached \
   inversions/2020-06-03_Calama/output_band002-01_killed \
   inversions/2020-06-03_Calama/output_band006-015_sel \
   inversions/2020-06-03_Calama/output_disp_t34_full \
   inversions/2020-06-03_Calama/output_plano64_sel \
   inversions/2020-06-03_Calama/output_prev_002-01 \
   inversions/2020-06-03_Calama/output_prev_9sta_plano32 \
   inversions/2020-06-03_Calama/output_prueba_t34_300 \
   inversions/2020-06-03_Calama/output_vel_t34_006-015 2>/dev/null

# --- Calama 2026 ---
cd 2026-05-25_Calama && mkdir -p archive
mv output_A22F_backup output_AF01_backup output_noAF01_42k_backup \
   output_noAF01_7k_backup output_prev_42k output_test_64km_backup archive/
cd ..
git -C /home/alex/KDEllipsPy rm -r --cached \
   inversions/2026-05-25_Calama/output_A22F_backup \
   inversions/2026-05-25_Calama/output_AF01_backup \
   inversions/2026-05-25_Calama/output_noAF01_42k_backup \
   inversions/2026-05-25_Calama/output_noAF01_7k_backup \
   inversions/2026-05-25_Calama/output_prev_42k \
   inversions/2026-05-25_Calama/output_test_64km_backup 2>/dev/null

# --- Copiapo 2025 ---
cd 2025-06-06_Copiapo && mkdir -p archive
mv output_dc_z80 output_dc_z80_launch.log archive/ 2>/dev/null
cd ..
```

## Paso 4 — verificación (antes de commitear)

```bash
cd /home/alex/KDEllipsPy
git status                     # revisar renames (R) y deletions cacheadas (D)
ls inversions/                 # debe verse: 3 eventos ISO_Lugar + _codigo_legacy + comparacion
du -sh inversions/*/archive    # el scratch quedó en disco
git check-ignore inversions/2020-06-03_Calama/archive   # debe imprimir la ruta (ignorada)
```

## Notas / riesgos
- **Única rotura real del renombrado** `2020-06-03` → `2020-06-03_Calama`: 3 referencias en
  `comparacion_2020_vs_2026/codigos/`. Hay que editarlas (cambiar `"2020-06-03"` → `"2020-06-03_Calama"`):
  - `comparar_modelos.py:39`            `"Calama 2020": INVERSIONS / "2020-06-03"`
  - `comparar_modelo_velocidad.py:29`   `INPUT_CTL = INVERSIONS / "2020-06-03" / "input.ctl"`
  - `mapa_mecanismos.py:28`             `dict(dir=INVERSIONS / "2020-06-03", ...)`
  El resto de coincidencias (`calama2020`, `2020-06-03T07:...`, títulos) son comentarios o
  strings de fecha → **no** se tocan.
- Los scripts dentro de cada evento (`real_disp.py`, `espectros_snr.py`, etc.) usan rutas
  relativas a su propia carpeta → el renombrado del contenedor no los afecta.
- Si algún script abre un `output_*` que movimos a `archive/`, ajustar esa ruta (no se detectó
  ninguno en el grep, pero verificar tras mover).
- El `.conda`/notebooks que apunten a rutas viejas también.
- Nada de esto recupera espacio del `.git` (eso requeriría `git filter-repo`, decisión aparte).
