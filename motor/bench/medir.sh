#!/usr/bin/env bash
# Mide los filtros nuevos y la cadena completa con 1, 2, 4 y 8 hilos sobre una
# imagen de 4096x4096 y deja la tabla en mediciones/motor.csv.
#
#   make arnes motor generar && bash bench/medir.sh
set -euo pipefail

cd "$(dirname "$0")/.."

# Un hilo por núcleo físico mientras haya núcleos libres. Sin esto el
# sistema puede poner dos hilos en hermanos SMT del mismo núcleo y la
# medición con 2 hilos sale pesimista.
export OMP_PLACES="${OMP_PLACES:-cores}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"

repeticiones=5
salida="../mediciones/motor.csv"
trabajo="${TMPDIR:-/tmp}/medicion-motor"

for binario in motor arnes generar; do
    if [ ! -x "./$binario" ]; then
        echo "falta ./$binario: corre 'make arnes motor generar' primero" >&2
        exit 2
    fi
done

mkdir -p "$trabajo" "$(dirname "$salida")"
color="$trabajo/medicion.ppm"
grises="$trabajo/medicion.pgm"

# La imagen no se versiona: se genera cada vez, siempre igual.
./generar "$color" 4096 4096 > /dev/null
./motor gris "$color" "$grises"

echo "filtro,hilos,repeticion,milisegundos" > "$salida"

medir() {
    local imagen="$1"
    shift
    local filtro hilos
    for filtro in "$@"; do
        for hilos in 1 2 4 8; do
            OMP_NUM_THREADS="$hilos" ./arnes "$filtro" "$imagen" "$repeticiones" >> "$salida"
        done
    done
}

# Reciben el PPM de color.
medir "$color" gris cadena

# Reciben el PGM de un canal.
medir "$grises" gaussiano gaussiano-directo sobel otsu ecualizar \
    histograma-compartido histograma-privado histograma-relleno

# Línea base secuencial del histograma: una sola configuración, un hilo.
OMP_NUM_THREADS=1 ./arnes histograma-secuencial "$grises" "$repeticiones" --secuencial >> "$salida"

echo "escrito $salida con $(( $(wc -l < "$salida") - 1 )) filas"