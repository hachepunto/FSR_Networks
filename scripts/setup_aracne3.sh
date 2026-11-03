#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# SETUP DE ARACNe3
#
# Clona (o actualiza) ARACNe3 en network_inference.aracne3.dir,
# fija el commit network_inference.aracne3.commit de config.yaml
# y compila en modo Release. Al terminar escribe
# build/ARACNe3_BUILT_COMMIT, que run_aracne3.sh verifica antes
# de correr — así una actualización del commit en config.yaml
# obliga a recompilar en vez de usar un binario viejo en silencio.
#
# Requiere en PATH: git, cmake, un compilador C++17 con OpenMP y
# Rscript con el paquete yaml (p.ej. el entorno conda tcga_env).
# ============================================================

# ============================================================
# ARGUMENTOS
#
# 1. Hilos para compilar (opcional, default: nproc/2 - 1)
# ============================================================

if [[ $# -gt 1 ]]; then
    echo "Uso: bash scripts/setup_aracne3.sh [build_threads]"
    exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="$REPO_ROOT/config/config.yaml"

for cmd in git cmake Rscript; do
    command -v "$cmd" >/dev/null || { echo "ERROR: No se encontró $cmd en PATH"; exit 1; }
done

# readLines + yaml.load en vez de read_yaml: read_yaml trunca el archivo
# en el primer carácter no ASCII si el locale no es UTF-8.
aracne3_cfg() {
    Rscript -e 'a <- commandArgs(TRUE); cat(yaml::yaml.load(readLines(a[1], warn = FALSE, encoding = "UTF-8"))$network_inference$aracne3[[a[2]]])' \
        "$CONFIG" "$1"
}

ARACNE3_DIR="$(realpath -m "$(aracne3_cfg dir)")"   # absoluta: más abajo se hace cd a este directorio
ARACNE3_REPO="$(aracne3_cfg repo)"
ARACNE3_COMMIT="$(aracne3_cfg commit)"

[[ -n "$ARACNE3_DIR" && -n "$ARACNE3_REPO" && -n "$ARACNE3_COMMIT" ]] || {
    echo "ERROR: Faltan dir/repo/commit en network_inference.aracne3 de $CONFIG"
    exit 1
}

DEFAULT_THREADS=$(( $(nproc) / 2 - 1 ))
(( DEFAULT_THREADS >= 1 )) || DEFAULT_THREADS=1
BUILD_THREADS="${1:-$DEFAULT_THREADS}"
[[ "$BUILD_THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: build_threads debe ser un entero positivo"; exit 1; }

echo ">> ARACNe3: $ARACNE3_REPO @ $ARACNE3_COMMIT"
echo ">> Directorio: $ARACNE3_DIR"

# ── 1. Clonar o validar el clon existente ──

if [[ ! -e "$ARACNE3_DIR" ]]; then
    git clone "$ARACNE3_REPO" "$ARACNE3_DIR"
else
    [[ -d "$ARACNE3_DIR/.git" ]] || { echo "ERROR: $ARACNE3_DIR existe pero no es un repositorio git"; exit 1; }

    ORIGIN="$(git -C "$ARACNE3_DIR" remote get-url origin)"
    if [[ "${ORIGIN%.git}" != "${ARACNE3_REPO%.git}" ]]; then
        echo "ERROR: el remoto origin de $ARACNE3_DIR ($ORIGIN) no coincide con config.yaml ($ARACNE3_REPO)"
        exit 1
    fi
fi

cd "$ARACNE3_DIR"

if ! git cat-file -e "${ARACNE3_COMMIT}^{commit}" 2>/dev/null; then
    git fetch origin
fi
git cat-file -e "${ARACNE3_COMMIT}^{commit}" 2>/dev/null || {
    echo "ERROR: el commit $ARACNE3_COMMIT no existe en $ARACNE3_REPO"
    exit 1
}

# ── 2. Clon limpio ──
# Se ignora el contenido sucio dentro de submódulos (no afecta el
# checkout); un submódulo en otro commit sí se reporta.

if [[ -n "$(git status --porcelain --untracked-files=no --ignore-submodules=dirty)" ]]; then
    echo "ERROR: $ARACNE3_DIR tiene cambios locales; respáldalos o descártalos antes de volver a correr este script."
    git status --short --untracked-files=no --ignore-submodules=dirty
    exit 1
fi

# ── 3. Checkout y submódulos ──

git checkout --quiet --detach "$ARACNE3_COMMIT"
git submodule update --init --recursive

# ── 4. Compilación ──

rm -f build/ARACNe3_BUILT_COMMIT

echo ">> Compilando con $BUILD_THREADS hilos..."
# POLICY_VERSION_MINIMUM: CMake >= 4 rechaza subproyectos de deps/boost
# que declaran compatibilidad < 3.5; en CMake < 4 no tiene efecto.
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build -j "$BUILD_THREADS"

ARACNE3_BIN="$ARACNE3_DIR/build/src/app/ARACNe3_app_release"
[[ -x "$ARACNE3_BIN" ]] || { echo "ERROR: la compilación no generó $ARACNE3_BIN"; exit 1; }

echo "$ARACNE3_COMMIT" > build/ARACNe3_BUILT_COMMIT

echo ">> ARACNe3 listo: $ARACNE3_BIN ($ARACNE3_COMMIT)"
