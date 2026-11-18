#include "filtros.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

Imagen gris_secuencial(const Imagen& e) {
    if (e.canales != 3) throw std::runtime_error("gris espera una imagen de tres canales");
    Imagen s;
    s.ancho = e.ancho;
    s.alto = e.alto;
    s.canales = 1;
    s.datos.resize(static_cast<std::size_t>(e.ancho) * e.alto);
    for (int y = 0; y < e.alto; ++y) {
        for (int x = 0; x < e.ancho; ++x) {
            const double v = 0.299 * e.en(x, y, 0) + 0.587 * e.en(x, y, 1) + 0.114 * e.en(x, y, 2);
            s.poner(x, y, 0, redondear(v));
        }
    }
    return s;
}

Imagen gris_paralelo(const Imagen& e) {
    if (e.canales != 3) throw std::runtime_error("gris espera una imagen de tres canales");
    Imagen s;
    s.ancho = e.ancho;
    s.alto = e.alto;
    s.canales = 1;
    s.datos.resize(static_cast<std::size_t>(e.ancho) * e.alto);
    // Cada fila es independiente de las demás: no hay nada que compartir
    // entre hilos ni dependencias entre iteraciones.
    #pragma omp parallel for schedule(static)
    for (int y = 0; y < e.alto; ++y) {
        for (int x = 0; x < e.ancho; ++x) {
            const double v = 0.299 * e.en(x, y, 0) + 0.587 * e.en(x, y, 1) + 0.114 * e.en(x, y, 2);
            s.poner(x, y, 0, redondear(v));
        }
    }
    return s;
}

namespace {

constexpr int NUCLEO_GAUSSIANO[5] = {1, 4, 6, 4, 1};

void validar_gaussiano(const Imagen& e) {
    if (e.canales != 1 && e.canales != 3) {
        throw std::runtime_error("gaussiano espera una imagen de uno o tres canales");
    }
}

Imagen gaussiano_separable(const Imagen& e, bool paralelo) {
    validar_gaussiano(e);
    Imagen s = e;
    std::vector<double> intermedio(e.datos.size());

    #pragma omp parallel for schedule(static) if(paralelo)
    for (int y = 0; y < e.alto; ++y) {
        for (int x = 0; x < e.ancho; ++x) {
            for (int c = 0; c < e.canales; ++c) {
                double suma = 0.0;
                for (int k = -2; k <= 2; ++k) {
                    suma += NUCLEO_GAUSSIANO[k + 2] *
                            e.en(recortar(x + k, e.ancho), y, c);
                }
                const std::size_t indice =
                    (static_cast<std::size_t>(y) * e.ancho + x) * e.canales + c;
                intermedio[indice] = suma;
            }
        }
    }

    #pragma omp parallel for schedule(static) if(paralelo)
    for (int y = 0; y < e.alto; ++y) {
        for (int x = 0; x < e.ancho; ++x) {
            for (int c = 0; c < e.canales; ++c) {
                double suma = 0.0;
                for (int k = -2; k <= 2; ++k) {
                    const int fila = recortar(y + k, e.alto);
                    const std::size_t indice =
                        (static_cast<std::size_t>(fila) * e.ancho + x) * e.canales + c;
                    suma += NUCLEO_GAUSSIANO[k + 2] * intermedio[indice];
                }
                s.poner(x, y, c, redondear(suma / 256.0));
            }
        }
    }
    return s;
}

Imagen gaussiano_directo(const Imagen& e, bool paralelo) {
    validar_gaussiano(e);
    Imagen s = e;

    #pragma omp parallel for schedule(static) if(paralelo)
    for (int y = 0; y < e.alto; ++y) {
        for (int x = 0; x < e.ancho; ++x) {
            for (int c = 0; c < e.canales; ++c) {
                double suma = 0.0;
                for (int ky = -2; ky <= 2; ++ky) {
                    const int fila = recortar(y + ky, e.alto);
                    for (int kx = -2; kx <= 2; ++kx) {
                        suma += NUCLEO_GAUSSIANO[ky + 2] *
                                NUCLEO_GAUSSIANO[kx + 2] *
                                e.en(recortar(x + kx, e.ancho), fila, c);
                    }
                }
                s.poner(x, y, c, redondear(suma / 256.0));
            }
        }
    }
    return s;
}

void validar_sobel(const Imagen& e) {
    if (e.canales != 1) {
        throw std::runtime_error("sobel espera una imagen de un canal");
    }
}

Imagen sobel(const Imagen& e, bool paralelo) {
    validar_sobel(e);
    Imagen s;
    s.ancho = e.ancho;
    s.alto = e.alto;
    s.canales = 1;
    s.datos.resize(static_cast<std::size_t>(e.ancho) * e.alto);
    constexpr int KERNEL_X[3][3] = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1},
    };
    constexpr int KERNEL_Y[3][3] = {
        {-1, -2, -1},
        {0, 0, 0},
        {1, 2, 1},
    };

    #pragma omp parallel for schedule(static) if(paralelo)
    for (int y = 0; y < e.alto; ++y) {
        for (int x = 0; x < e.ancho; ++x) {
            int gx = 0;
            int gy = 0;
            for (int ky = -1; ky <= 1; ++ky) {
                const int fila = recortar(y + ky, e.alto);
                for (int kx = -1; kx <= 1; ++kx) {
                    const int valor = e.en(recortar(x + kx, e.ancho), fila, 0);
                    gx += KERNEL_X[ky + 1][kx + 1] * valor;
                    gy += KERNEL_Y[ky + 1][kx + 1] * valor;
                }
            }
            const double magnitud = std::sqrt(static_cast<double>(gx) * gx +
                                              static_cast<double>(gy) * gy);
            s.poner(x, y, 0, redondear(magnitud));
        }
    }
    return s;
}

// El promedio de la vecindad de 3x3 de un canal, con las coordenadas
// recortadas al rectángulo de la imagen.
double promedio_vecindad(const Imagen& e, int x, int y, int c) {
    double suma = 0.0;
    for (int j = -1; j <= 1; ++j) {
        const int fila = recortar(y + j, e.alto);
        for (int i = -1; i <= 1; ++i) {
            suma += e.en(recortar(x + i, e.ancho), fila, c);
        }
    }
    return suma / 9.0;
}

}  // namespace

Imagen caja_secuencial(const Imagen& e) {
    Imagen s = e;
    for (int y = 0; y < e.alto; ++y) {
        for (int x = 0; x < e.ancho; ++x) {
            for (int c = 0; c < e.canales; ++c) {
                s.poner(x, y, c, redondear(promedio_vecindad(e, x, y, c)));
            }
        }
    }
    return s;
}

Imagen caja_paralelo(const Imagen& e) {
    Imagen s = e;
    // El ciclo externo recorre filas y el interno columnas, que es el orden
    // en que la imagen está en memoria. Repartir por filas deja a cada hilo
    // con una banda contigua.
    #pragma omp parallel for schedule(static)
    for (int y = 0; y < e.alto; ++y) {
        for (int x = 0; x < e.ancho; ++x) {
            for (int c = 0; c < e.canales; ++c) {
                s.poner(x, y, c, redondear(promedio_vecindad(e, x, y, c)));
            }
        }
    }
    return s;
}

Imagen gaussiano_secuencial(const Imagen& e) {
    return gaussiano_separable(e, false);
}

Imagen gaussiano_paralelo(const Imagen& e) {
    return gaussiano_separable(e, true);
}

Imagen gaussiano_directo_secuencial(const Imagen& e) {
    return gaussiano_directo(e, false);
}

Imagen gaussiano_directo_paralelo(const Imagen& e) {
    return gaussiano_directo(e, true);
}

Imagen sobel_secuencial(const Imagen& e) {
    return sobel(e, false);
}

Imagen sobel_paralelo(const Imagen& e) {
    return sobel(e, true);
}

Imagen otsu_secuencial(const Imagen&) {
    throw std::runtime_error("otsu: sin implementar");
}

Imagen otsu_paralelo(const Imagen&) {
    throw std::runtime_error("otsu: sin implementar");
}
