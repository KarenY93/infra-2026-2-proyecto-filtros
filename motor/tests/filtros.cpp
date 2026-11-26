#include "filtros.hpp"

#include <iostream>
#include <stdexcept>

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

void exigir_iguales(const Imagen& a, const Imagen& b, const char* caso) {
    if (a.ancho != b.ancho || a.alto != b.alto || a.canales != b.canales ||
        a.datos != b.datos) {
        throw std::runtime_error(std::string("salidas distintas: ") + caso);
    }
}

void probar_gaussiano() {
    const int dimensiones[][2] = {{1, 1}, {2, 1}, {7, 4}};
    for (int canales : {1, 3}) {
        for (const auto& dimension : dimensiones) {
            const Imagen entrada = generar(dimension[0], dimension[1], canales);
            const Imagen directa = gaussiano_directo_secuencial(entrada);
            exigir_iguales(gaussiano_secuencial(entrada), directa, "gaussiano separable/directo");
            exigir_iguales(gaussiano_paralelo(entrada), directa, "gaussiano paralelo/directo");
            exigir_iguales(gaussiano_directo_paralelo(entrada), directa,
                           "gaussiano directo paralelo/secuencial");
        }
    }
}

void probar_sobel() {
    const Imagen entrada = generar(9, 7, 1);
    exigir_iguales(sobel_paralelo(entrada), sobel_secuencial(entrada),
                   "sobel paralelo/secuencial");

    Imagen rampa_x;
    rampa_x.ancho = 3;
    rampa_x.alto = 3;
    rampa_x.canales = 1;
    rampa_x.datos = {10, 20, 30, 10, 20, 30, 10, 20, 30};
    if (sobel_secuencial(rampa_x).en(1, 1, 0) != 80) {
        throw std::runtime_error("Sobel debe calcular gx en una rampa horizontal");
    }

    Imagen rampa_diagonal;
    rampa_diagonal.ancho = 3;
    rampa_diagonal.alto = 3;
    rampa_diagonal.canales = 1;
    rampa_diagonal.datos = {0, 10, 20, 10, 20, 30, 20, 30, 40};
    if (sobel_secuencial(rampa_diagonal).en(1, 1, 0) != 113) {
        throw std::runtime_error("Sobel debe combinar gx y gy para la magnitud");
    }

    Imagen plana = generar(3, 3, 1);
    for (unsigned char& pixel : plana.datos) pixel = 120;
    const Imagen salida = sobel_secuencial(plana);
    for (unsigned char pixel : salida.datos) {
        if (pixel != 0) throw std::runtime_error("Sobel de imagen plana debe ser cero");
    }

    bool rechazo = false;
    try {
        sobel_secuencial(generar(2, 2, 3));
    } catch (const std::runtime_error&) {
        rechazo = true;
    }
    if (!rechazo) throw std::runtime_error("Sobel debería rechazar imágenes de color");
}

}  // namespace

int main() {
    try {
        probar_gaussiano();
        std::cout << "  ok    gaussiano separable igual a la convolución directa 5x5\n";
        probar_sobel();
        std::cout << "  ok    Sobel paralelo y secuencial, imagen plana y validación\n";
    } catch (const std::exception& error) {
        std::cerr << "prueba_filtros: " << error.what() << "\n";
        return 1;
    }
    return 0;
}
