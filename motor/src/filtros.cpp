#include "filtros.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>
#include <array>
#include <omp.h>
#include <string>



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

namespace {

void exigir_un_canal(const Imagen& e, const char* filtro) {
    if (e.canales != 1) {
        throw std::runtime_error(std::string(filtro) + " espera una imagen de un canal");
    }
}

// Un histograma por hilo, alineado a línea de caché para que dos hilos
// nunca escriban sobre la misma línea.
struct alignas(64) HistogramaAlineado {
    Histograma h{};
};

using Tabla = std::array<unsigned char, 256>;

// Ecualización: la tabla sale del histograma acumulado. Si todos los píxeles
// tienen el mismo valor (N == cmin) la fórmula divide por cero, y la imagen
// se deja igual.
Tabla tabla_ecualizacion(const Histograma& h, long long total) {
    std::array<long long, 256> c{};
    long long acumulado = 0;
    long long cmin = 0;
    bool visto = false;
    for (int v = 0; v < 256; ++v) {
        acumulado += h[v];
        c[v] = acumulado;
        if (!visto && acumulado > 0) {
            cmin = acumulado;
            visto = true;
        }
    }
    Tabla tabla{};
    for (int v = 0; v < 256; ++v) {
        if (total == cmin) {
            tabla[v] = static_cast<unsigned char>(v);
        } else {
            // Una sola división, con el numerador entero: así el cociente
            // exacto no pierde precisión entre dos operaciones en double.
            tabla[v] = redondear(static_cast<double>((c[v] - cmin) * 255) /
                                 static_cast<double>(total - cmin));
        }
    }
    return tabla;
}

// Barrido de Otsu sobre el histograma: es siempre secuencial, son 255
// iteraciones. Con empate en la varianza gana el umbral más bajo.
int umbral_otsu(const Histograma& h, long long total) {
    long long suma_total = 0;
    for (int v = 0; v < 256; ++v) suma_total += v * h[v];

    long long n0 = 0;
    long long suma0 = 0;
    double mejor = -1.0;
    int umbral = 0;
    for (int t = 0; t < 255; ++t) {
        n0 += h[t];
        suma0 += t * h[t];
        const long long n1 = total - n0;
        if (n0 == 0 || n1 == 0) continue;
        const double mu0 = static_cast<double>(suma0) / static_cast<double>(n0);
        const double mu1 = static_cast<double>(suma_total - suma0) / static_cast<double>(n1);
        const double w0 = static_cast<double>(n0) / static_cast<double>(total);
        const double w1 = static_cast<double>(n1) / static_cast<double>(total);
        const double varianza = w0 * w1 * (mu0 - mu1) * (mu0 - mu1);
        if (varianza > mejor) {
            mejor = varianza;
            umbral = t;
        }
    }
    return umbral;
}

// Salida = tabla[entrada]. Es un map: el reparto va por filas.
Imagen aplicar_tabla(const Imagen& e, const Tabla& tabla, bool paralelo) {
    Imagen s;
    s.ancho = e.ancho;
    s.alto = e.alto;
    s.canales = 1;
    s.datos.resize(static_cast<std::size_t>(e.ancho) * e.alto);
    #pragma omp parallel for schedule(static) if(paralelo)
    for (int y = 0; y < e.alto; ++y) {
        for (int x = 0; x < e.ancho; ++x) {
            s.poner(x, y, 0, tabla[e.en(x, y, 0)]);
        }
    }
    return s;
}

long long total_pixeles(const Imagen& e) {
    return static_cast<long long>(e.ancho) * e.alto;
}

}  // namespace

Histograma histograma_secuencial(const Imagen& e) {
    exigir_un_canal(e, "histograma");
    Histograma h{};
    for (int y = 0; y < e.alto; ++y) {
        for (int x = 0; x < e.ancho; ++x) ++h[e.en(x, y, 0)];
    }
    return h;
}

// Un solo arreglo para todos los hilos: cada incremento es atómico, así que
// el resultado es correcto, pero los hilos compiten por los mismos 256
// contadores.
Histograma histograma_compartido(const Imagen& e) {
    exigir_un_canal(e, "histograma");
    Histograma h{};
    long long* cuenta = h.data();
    #pragma omp parallel for schedule(static)
    for (int y = 0; y < e.alto; ++y) {
        for (int x = 0; x < e.ancho; ++x) {
            const unsigned char v = e.en(x, y, 0);
            #pragma omp atomic
            ++cuenta[v];
        }
    }
    return h;
}

// Un histograma por hilo en un vector contiguo, sin alineación: los
// extremos de histogramas vecinos pueden caer en la misma línea de caché.
Histograma histograma_privado(const Imagen& e) {
    exigir_un_canal(e, "histograma");
    std::vector<Histograma> locales(static_cast<std::size_t>(omp_get_max_threads()));
    #pragma omp parallel
    {
        Histograma& local = locales[static_cast<std::size_t>(omp_get_thread_num())];
        #pragma omp for schedule(static)
        for (int y = 0; y < e.alto; ++y) {
            for (int x = 0; x < e.ancho; ++x) ++local[e.en(x, y, 0)];
        }
    }
    Histograma h{};
    for (const Histograma& local : locales) {
        for (int v = 0; v < 256; ++v) h[v] += local[v];
    }
    return h;
}

// Igual que el anterior, pero cada histograma empieza en un múltiplo de 64
// bytes: ninguna línea de caché es compartida entre hilos.
Histograma histograma_privado_relleno(const Imagen& e) {
    exigir_un_canal(e, "histograma");
    std::vector<HistogramaAlineado> locales(static_cast<std::size_t>(omp_get_max_threads()));
    #pragma omp parallel
    {
        Histograma& local = locales[static_cast<std::size_t>(omp_get_thread_num())].h;
        #pragma omp for schedule(static)
        for (int y = 0; y < e.alto; ++y) {
            for (int x = 0; x < e.ancho; ++x) ++local[e.en(x, y, 0)];
        }
    }
    Histograma h{};
    for (const HistogramaAlineado& local : locales) {
        for (int v = 0; v < 256; ++v) h[v] += local.h[v];
    }
    return h;
}

Imagen ecualizar_secuencial(const Imagen& e) {
    const Histograma h = histograma_secuencial(e);
    return aplicar_tabla(e, tabla_ecualizacion(h, total_pixeles(e)), false);
}

Imagen ecualizar_paralelo(const Imagen& e) {
    const Histograma h = histograma_privado_relleno(e);
    return aplicar_tabla(e, tabla_ecualizacion(h, total_pixeles(e)), true);
}

Imagen otsu_secuencial(const Imagen& e) {
    const Histograma h = histograma_secuencial(e);
    const int t = umbral_otsu(h, total_pixeles(e));
    Tabla tabla{};
    for (int v = 0; v < 256; ++v) tabla[v] = v <= t ? 0 : 255;
    return aplicar_tabla(e, tabla, false);
}

Imagen otsu_paralelo(const Imagen& e) {
    const Histograma h = histograma_privado_relleno(e);
    const int t = umbral_otsu(h, total_pixeles(e));   // barrido secuencial
    Tabla tabla{};
    for (int v = 0; v < 256; ++v) tabla[v] = v <= t ? 0 : 255;
    return aplicar_tabla(e, tabla, true);
}
