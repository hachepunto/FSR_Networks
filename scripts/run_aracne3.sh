#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# ARACNe3
#
# Construye y consolida una red regulatoria a partir de una
# matriz de expresión y una lista de reguladores (factores de
# transcripción), usando el fork corregido de ARACNe3:
#   https://github.com/hachepunto/ARACNe3
#
# A diferencia de ARACNe-AP, ARACNe3 no tiene un paso separado
# de cálculo de umbral de información mutua: una sola invocación
# genera las subredes y las consolida en consolidated-net_<runid>.tsv.
# ============================================================

# ============================================================
# ARGUMENTOS
#
# 1. Directorio raíz de ARACNe3 (donde está build/src/app/ARACNe3_app_release)
# 2. Matriz de expresión (tumor, ya normalizada — ver preprocess_expression.R)
# 3. Lista de reguladores (factores de transcripción)
# 4. Directorio de salida
# 5. Identificador de corrida (runid) — usar el analysis_id (p.ej. BRCA_basal)
# 6. Núcleos
# 7. Subredes -x (opcional, default: 100)
# 8. Semilla (opcional, default: 42 — coincide con network_inference.aracne3.seed en config.yaml)
# 9. Alpha / FDR (opcional, default: 0.05, el propio default de ARACNe3)
# ============================================================

if [[ $# -lt 6 || $# -gt 9 ]]; then
    echo "Uso: bash scripts/run_aracne3.sh <aracne3_dir> <expression_file> <regulators_file> <output_dir> <runid> <threads> [subnets] [seed] [alpha]"
    exit 1
fi

# ============================================================
# SALIDA
#
#   <output_dir>/consolidated-net_<runid>.tsv
#   <output_dir>/subnets/subnet<N>_<runid>.tsv
#   <output_dir>/subnets_log/log_subnet<N>_<runid>.txt
#   <output_dir>/log_<runid>.txt
# ============================================================

ARACNE3_DIR="$1"
EXPFILE="$2"
REGULATORS="$3"
OUTDIR="$4"
RUNID="$5"
THREADS="$6"
SUBNETS="${7:-100}"
SEED="${8:-42}"
ALPHA="${9:-0.05}"

[[ -d "$ARACNE3_DIR" ]] || { echo "ERROR: No se encontró el directorio de ARACNe3: $ARACNE3_DIR"; exit 1; }
[[ -f "$EXPFILE" ]] || { echo "ERROR: No se encontró la matriz de expresión: $EXPFILE"; exit 1; }
[[ -f "$REGULATORS" ]] || { echo "ERROR: No se encontró la lista de reguladores: $REGULATORS"; exit 1; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: El número de núcleos debe ser un entero positivo"; exit 1; }
[[ "$SUBNETS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: El número de subredes (-x) debe ser un entero positivo"; exit 1; }
[[ "$SEED" =~ ^[0-9]+$ ]] || { echo "ERROR: La semilla debe ser un entero no negativo"; exit 1; }

command -v realpath >/dev/null || { echo "ERROR: No se encontró realpath"; exit 1; }

ARACNE3_BIN="$(realpath "$ARACNE3_DIR")/build/src/app/ARACNe3_app_release"

[[ -x "$ARACNE3_BIN" ]] || {
    echo "ERROR: No se encontró el ejecutable compilado de ARACNe3: $ARACNE3_BIN"
    echo "Compílalo primero con CMake en modo Release (ver README del fork)."
    exit 1
}

EXPFILE="$(realpath "$EXPFILE")"
REGULATORS="$(realpath "$REGULATORS")"
OUTDIR="$(realpath -m "$OUTDIR")"

mkdir -p "$OUTDIR"

CORES_AVAILABLE="$(nproc)"
CORES="$THREADS"

if (( CORES > CORES_AVAILABLE )); then
    echo ">> El servidor tiene $CORES_AVAILABLE núcleos disponibles. Se utilizarán esos."
    CORES="$CORES_AVAILABLE"
fi

echo ">> ARACNe3: runid=$RUNID | núcleos=$CORES | subredes=$SUBNETS | semilla=$SEED | alpha=$ALPHA"
echo ">> Matriz: $EXPFILE"
echo ">> Reguladores: $REGULATORS"
echo ">> Salida: $OUTDIR"

# Se corre desde OUTDIR para que el directorio de caché interno de
# ARACNe3 (./ARACNe3_cached/, usado para los modelos nulos de MI)
# quede confinado a esta corrida y no se mezcle con el de otro
# analysis_id ni con el de otro usuario del servidor.
cd "$OUTDIR"

"$ARACNE3_BIN" \
    -e "$EXPFILE" \
    -r "$REGULATORS" \
    -o "$OUTDIR" \
    -x "$SUBNETS" \
    --threads "$CORES" \
    --seed "$SEED" \
    --alpha "$ALPHA" \
    --runid "$RUNID"

NETWORK_FILE="$OUTDIR/consolidated-net_${RUNID}.tsv"

if [[ ! -s "$NETWORK_FILE" ]]; then
    echo "ERROR: ARACNe3 no generó una red consolidada válida: $NETWORK_FILE"
    exit 1
fi

echo ">> Red consolidada generada correctamente: $NETWORK_FILE"
