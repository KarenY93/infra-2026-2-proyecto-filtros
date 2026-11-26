#!/usr/bin/env python3
"""Calcula T(p), S(p) y E(p) a partir de las sesiones de medición.

Metodología: en cada sesión, de las cinco repeticiones de cada configuración
se descarta la primera y se toma la mediana de las cuatro restantes. Con
varias sesiones, el valor reportado es la mediana de esas medianas.

Uso:  python3 motor/bench/analizar.py [sesion1.csv sesion2.csv ...]
Sin argumentos usa mediciones/motor.csv.
Deja mediciones/tablas.md y, si hay matplotlib, mediciones/speedup.png.
"""
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[2]
ENTRADA = RAIZ / "mediciones" / "motor.csv"
TABLAS = RAIZ / "mediciones" / "tablas.md"
GRAFICA = RAIZ / "mediciones" / "speedup.png"

HILOS = [1, 2, 4, 8]
FILTROS = ["gris", "gaussiano", "gaussiano-directo", "sobel", "ecualizar", "otsu", "cadena"]
HISTOGRAMAS = ["histograma-compartido", "histograma-privado", "histograma-relleno"]
CURVAS = ["gaussiano", "sobel", "otsu", "cadena"]
DISPERSION_MAXIMA = 0.15   # (máx - mín) / mediana entre las medianas de las sesiones


def medianas_de_sesion(ruta):
    corridas = defaultdict(list)
    with open(ruta, newline="") as f:
        for fila in csv.DictReader(f):
            if int(fila["repeticion"]) > 1:   # se descarta la primera
                corridas[(fila["filtro"], int(fila["hilos"]))].append(float(fila["milisegundos"]))
    return {clave: statistics.median(v) for clave, v in corridas.items() if v}


def leer(rutas):
    sesiones = [medianas_de_sesion(r) for r in rutas]
    medianas, dispersion = {}, {}
    for clave in set().union(*sesiones):
        valores = [s[clave] for s in sesiones if clave in s]
        mediana = statistics.median(valores)
        medianas[clave] = mediana
        dispersion[clave] = (max(valores) - min(valores)) / mediana if len(valores) > 1 else 0.0
    return medianas, dispersion


def tabla_filtro(filtro, T):
    if (filtro, 1) not in T:
        return None
    base = T[(filtro, 1)]
    lineas = [f"### {filtro}", "", "| p | T(p) ms | S(p) | E(p) |", "|--:|--:|--:|--:|"]
    for p in HILOS:
        t = T.get((filtro, p))
        if t is None:
            continue
        s = base / t
        lineas.append(f"| {p} | {t:.2f} | {s:.2f} | {s / p:.2f} |")
    return "\n".join(lineas)


def tabla_histogramas(T):
    secuencial = T.get(("histograma-secuencial", 1))
    lineas = ["### Histograma: variantes", "", "T(p) en ms:", "",
              "| variante | p=1 | p=2 | p=4 | p=8 |", "|---|--:|--:|--:|--:|"]
    if secuencial is not None:
        lineas.append(f"| secuencial | {secuencial:.2f} | - | - | - |")
    for h in HISTOGRAMAS:
        celdas = [f"{T[(h, p)]:.2f}" if (h, p) in T else "-" for p in HILOS]
        lineas.append(f"| {h.replace('histograma-', '')} | " + " | ".join(celdas) + " |")
    if secuencial is not None:
        lineas += ["", "Aceleración respecto al histograma secuencial, T_sec / T(p):", "",
                   "| variante | p=1 | p=2 | p=4 | p=8 |", "|---|--:|--:|--:|--:|"]
        for h in HISTOGRAMAS:
            celdas = [f"{secuencial / T[(h, p)]:.2f}" if (h, p) in T else "-" for p in HILOS]
            lineas.append(f"| {h.replace('histograma-', '')} | " + " | ".join(celdas) + " |")
    return "\n".join(lineas)


def grafica(T):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib no está instalado: no se dibujó la gráfica")
        return False
    plt.figure(figsize=(6, 4.5))
    for filtro in CURVAS:
        if (filtro, 1) not in T:
            continue
        ps = [p for p in HILOS if (filtro, p) in T]
        plt.plot(ps, [T[(filtro, 1)] / T[(filtro, p)] for p in ps], marker="o", label=filtro)
    plt.plot(HILOS, HILOS, "k--", label="ideal")
    plt.xticks(HILOS)
    plt.xlabel("hilos")
    plt.ylabel("speedup S(p)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(GRAFICA, dpi=150)
    return True


def main():
    rutas = [Path(a) for a in sys.argv[1:]] or [ENTRADA]
    T, dispersion = leer(rutas)
    encabezado = f"Mediana de las medianas de {len(rutas)} sesión(es).\n\n"
    bloques = [b for b in (tabla_filtro(f, T) for f in FILTROS) if b]
    bloques.append(tabla_histogramas(T))
    texto = encabezado + "\n\n".join(bloques) + "\n"
    TABLAS.write_text(texto, encoding="utf-8")
    print(texto)

    if len(rutas) > 1:
        ruidosas = sorted((d, k) for k, d in dispersion.items() if d > DISPERSION_MAXIMA)
        if ruidosas:
            print(f"Configuraciones con más de {DISPERSION_MAXIMA * 100:.0f} % de diferencia entre sesiones:")
            for d, (filtro, hilos) in reversed(ruidosas):
                print(f"  {filtro} con {hilos} hilo(s): {d * 100:.0f} %")
    if grafica(T):
        print(f"gráfica en {GRAFICA}")
    print(f"tablas en {TABLAS}")


main()