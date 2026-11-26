#!/usr/bin/env python3
"""Referencia independiente: implementa las definiciones del enunciado en
Python puro, sin mirar el C++.
Uso: referencia.py gaussiano|sobel|ecualizar|otsu entrada salida"""
import sys, math
from fractions import Fraction

def leer(ruta):
    d = open(ruta, 'rb').read(); pos = 0; campos = []
    while len(campos) < 4:
        while d[pos:pos+1].isspace(): pos += 1
        if d[pos:pos+1] == b'#':
            while d[pos:pos+1] != b'\n': pos += 1
            continue
        ini = pos
        while not d[pos:pos+1].isspace(): pos += 1
        campos.append(d[ini:pos].decode())
    pos += 1
    mag, w, h, _ = campos; w, h = int(w), int(h); c = 3 if mag == 'P6' else 1
    return w, h, c, list(d[pos:pos + w*h*c])

def escribir(ruta, w, h, c, px):
    open(ruta, 'wb').write(('P6' if c == 3 else 'P5').encode() + b'\n%d %d\n255\n' % (w, h) + bytes(px))

def redondeo(v):  # al más cercano, empates hacia arriba, saturado
    return max(0, min(255, math.floor(v + Fraction(1, 2))))

def gaussiano(w, h, c, px):
    k = [Fraction(n, 16) for n in (1, 4, 6, 4, 1)]
    cl = lambda v, n: max(0, min(n - 1, v))
    T = {(x, y, ch): sum(k[i+2] * px[(y*w + cl(x+i, w))*c + ch] for i in range(-2, 3))
         for y in range(h) for x in range(w) for ch in range(c)}
    return [redondeo(sum(k[j+2] * T[(x, cl(y+j, h), ch)] for j in range(-2, 3)))
            for y in range(h) for x in range(w) for ch in range(c)]

def sobel(w, h, c, px):
    cl = lambda v, n: max(0, min(n - 1, v))
    I = lambda x, y: px[cl(y, h)*w + cl(x, w)]
    GX = [[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]]; GY = [[-1, -2, -1], [0, 0, 0], [1, 2, 1]]
    out = []
    for y in range(h):
        for x in range(w):
            gx = sum(GX[j+1][i+1] * I(x+i, y+j) for j in (-1, 0, 1) for i in (-1, 0, 1))
            gy = sum(GY[j+1][i+1] * I(x+i, y+j) for j in (-1, 0, 1) for i in (-1, 0, 1))
            out.append(min(255, redondeo(Fraction(math.sqrt(gx*gx + gy*gy)))))
    return out

def histograma(px):
    h = [0] * 256
    for v in px:
        h[v] += 1
    return h

def ecualizar(w, h, c, px):
    hist = histograma(px)
    N = w * h
    acum, cmin, visto, total = [], 0, False, 0
    for v in range(256):
        total += hist[v]
        acum.append(total)
        if not visto and total > 0:
            cmin, visto = total, True
    if N == cmin:                       # imagen constante: se deja igual
        return list(px)
    return [redondeo(Fraction((acum[v] - cmin) * 255, N - cmin)) for v in px]

def otsu(w, h, c, px):
    hist = histograma(px)
    N = w * h
    suma_total = sum(v * hist[v] for v in range(256))
    n0, suma0 = 0, 0
    mejor, umbral = None, 0
    for t in range(255):                # 0 <= t < 255
        n0 += hist[t]
        suma0 += t * hist[t]
        n1 = N - n0
        if n0 == 0 or n1 == 0:
            continue
        mu0 = Fraction(suma0, n0)
        mu1 = Fraction(suma_total - suma0, n1)
        varianza = Fraction(n0, N) * Fraction(n1, N) * (mu0 - mu1) ** 2
        if mejor is None or varianza > mejor:   # con empate gana el más bajo
            mejor, umbral = varianza, t
    return [0 if v <= umbral else 255 for v in px]

f, ent, sal = sys.argv[1:4]
w, h, c, px = leer(ent)
filtros = {'gaussiano': gaussiano, 'sobel': sobel, 'ecualizar': ecualizar, 'otsu': otsu}
if f not in filtros:
    sys.exit('filtro desconocido: ' + f)
escribir(sal, w, h, c, filtros[f](w, h, c, px))