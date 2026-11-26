#include "filtros.hpp"

#include <omp.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

Imagen generar(int ancho, int alto, int canales) {
    Imagen e;
    e.ancho = ancho;
    e.alto = alto;
    e.canales = canales;
    e.datos.resize(static_cast<std::size_t>(ancho) * alto * canales);
    for (int y = 0; y < alto; ++y) {
        for (int x = 0; x < ancho; ++x) {
            for (int c = 0; c < canales; ++c) {
                e.poner(x, y, c, static_cast<unsigned char>(
                                      (x * 37 + y * 61 + c * 89 + x * y * 17) % 256));
            }
        }
    }
    return e;
}

Imagen constante(int ancho, int alto, unsigned char valor) {
    Imagen e = generar(ancho, alto, 1);
    for (unsigned char& p : e.datos) p = valor;
    return e;
}

Imagen desde(int ancho, int alto, const std::vector<unsigned char>& datos) {
    Imagen e;
    e.ancho = ancho;
    e.alto = alto;
    e.canales = 1;
    e.datos = datos;
    return e;
}

void exigir_iguales(const Imagen& a, const Imagen& b, const std::string& caso) {
    if (a.ancho != b.ancho || a.alto != b.alto || a.canales != b.canales ||
        a.datos != b.datos) {
        throw std::runtime_error("salidas distintas: " + caso);
    }
}

template <typename F>
bool lanza(F f) {
    try {
        f();
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

const int DIMENSIONES[][2] = {{1, 1}, {7, 4}, {64, 37}, {200, 150}};
const int HILOS[] = {1, 2, 4, 8};

void probar_histograma() {
    const Histograma h = histograma_secuencial(desde(4, 1, {0, 0, 5, 255}));
    if (h[0] != 2 || h[5] != 1 || h[255] != 1) {
        throw std::runtime_error("histograma: conteo a mano incorrecto");
    }
    long long total = 0;
    for (long long cuenta : h) total += cuenta;
    if (total != 4) throw std::runtime_error("histograma: la suma debe ser N");

    const int original = omp_get_max_threads();
    for (int hilos : HILOS) {
        omp_set_num_threads(hilos);
        for (const auto& d : DIMENSIONES) {
            const Imagen e = generar(d[0], d[1], 1);
            const Histograma ref = histograma_secuencial(e);
            const std::string caso = " con " + std::to_string(hilos) + " hilos";
            if (histograma_compartido(e) != ref) {
                throw std::runtime_error("histograma compartido" + caso);
            }
            if (histograma_privado(e) != ref) {
                throw std::runtime_error("histograma privado" + caso);
            }
            if (histograma_privado_relleno(e) != ref) {
                throw std::runtime_error("histograma con relleno" + caso);
            }
        }
    }
    omp_set_num_threads(original);
}

void probar_ecualizar() {
    // N=4, cmin=1: los acumulados 1,2,3,4 dan 0, 85, 170 y 255.
    const Imagen salida = ecualizar_secuencial(desde(4, 1, {10, 20, 30, 40}));
    const std::vector<unsigned char> esperado = {0, 85, 170, 255};
    if (salida.datos != esperado) {
        throw std::runtime_error("ecualizar: rampa a mano incorrecta");
    }

    // Imagen constante: N == cmin, la fórmula divide por cero y la imagen
    // se deja igual.
    const Imagen plana = constante(5, 3, 100);
    exigir_iguales(ecualizar_secuencial(plana), plana, "ecualizar imagen constante");
    exigir_iguales(ecualizar_paralelo(plana), plana, "ecualizar paralelo imagen constante");
}

void probar_otsu() {
    // Dos clases bien separadas: lo claro a 255 y lo oscuro a 0.
    const Imagen salida = otsu_secuencial(desde(8, 1, {10, 10, 10, 10, 200, 200, 200, 200}));
    const std::vector<unsigned char> esperado = {0, 0, 0, 0, 255, 255, 255, 255};
    if (salida.datos != esperado) {
        throw std::runtime_error("otsu: bimodal a mano incorrecto");
    }

    // Imagen constante: el umbral queda en 0, así que 0 sale 0 y cualquier
    // otro valor sale 255.
    for (unsigned char p : otsu_secuencial(constante(5, 3, 100)).datos) {
        if (p != 255) throw std::runtime_error("otsu: constante 100 debe dar 255");
    }
    for (unsigned char p : otsu_secuencial(constante(5, 3, 0)).datos) {
        if (p != 0) throw std::runtime_error("otsu: constante 0 debe dar 0");
    }
}

void probar_paralelo_igual_secuencial() {
    const int original = omp_get_max_threads();
    for (int hilos : HILOS) {
        omp_set_num_threads(hilos);
        for (const auto& d : DIMENSIONES) {
            const Imagen e = generar(d[0], d[1], 1);
            const std::string caso = " con " + std::to_string(hilos) + " hilos";
            exigir_iguales(ecualizar_paralelo(e), ecualizar_secuencial(e), "ecualizar" + caso);
            exigir_iguales(otsu_paralelo(e), otsu_secuencial(e), "otsu" + caso);
        }
    }
    omp_set_num_threads(original);
}

void probar_validacion() {
    const Imagen color = generar(2, 2, 3);
    if (!lanza([&] { histograma_secuencial(color); })) {
        throw std::runtime_error("histograma debe rechazar color");
    }
    if (!lanza([&] { histograma_privado_relleno(color); })) {
        throw std::runtime_error("histograma con relleno debe rechazar color");
    }
    if (!lanza([&] { ecualizar_secuencial(color); })) {
        throw std::runtime_error("ecualizar debe rechazar color");
    }
    if (!lanza([&] { otsu_paralelo(color); })) {
        throw std::runtime_error("otsu debe rechazar color");
    }
}

}  // namespace

int main() {
    try {
        probar_histograma();
        std::cout << "  ok    histograma: las cuatro variantes cuentan igual con 1, 2, 4 y 8 hilos\n";
        probar_ecualizar();
        std::cout << "  ok    ecualización: rampa a mano e imagen constante\n";
        probar_otsu();
        std::cout << "  ok    Otsu: imagen bimodal e imagen constante\n";
        probar_paralelo_igual_secuencial();
        std::cout << "  ok    ecualización y Otsu paralelos idénticos a los secuenciales\n";
        probar_validacion();
        std::cout << "  ok    rechazo de imágenes de color\n";
    } catch (const std::exception& error) {
        std::cerr << "prueba_reduccion: " << error.what() << "\n";
        return 1;
    }
    return 0;
}