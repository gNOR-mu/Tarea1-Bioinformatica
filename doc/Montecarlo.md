# Probability Density Analysis (Monte Carlo)

# Descripción general

El algoritmo **Probability Density Analysis** (Análisis de Densidad de Probabilidad, comúnmente llamado estrategia de **Monte Carlo**) es una de las estrategias de resolución de Battleship más eficientes en términos de cantidad de disparos. A diferencia de los métodos heurísticos rígidos o puramente aleatorios, esta estrategia utiliza un enfoque probabilístico continuo para determinar el siguiente disparo óptimo.

En cada turno, el resolvedor evalúa todas las formas posibles en que los barcos restantes podrían estar posicionados en el tablero, respetando la información acumulada (disparos fallidos e impactos). Con base en esas evaluaciones genera un **Mapa de Calor** (*Heatmap*) o matriz de densidad de probabilidad, y dispara a la coordenada libre con mayor probabilidad de contener un barco.

---

# Funcionamiento detallado

### 1. Mapa de Calor Inicial Precalculado

Al cargar la clase, se calcula **una única vez** un `INITIAL_HEAT_MAP` estático que representa las posiciones válidas de todos los barcos sobre un tablero limpio. Para cada celda se cuentan las ventanas horizontales y verticales de cada barco que pasan por ella. Esta plantilla se copia con `System.arraycopy` al inicio de cada partida, eliminando el coste del cálculo inicial en las 500 000 simulaciones.

### 2. Actualización Incremental del Mapa de Calor

En lugar de recalcular todo el mapa en cada turno, la implementación aplica **actualizaciones mínimas** según el resultado del último disparo:

**Disparo fallido (`onMiss`)**
* Solo se restan las contribuciones de las ventanas que pasaban por la celda disparada y que eran válidas antes del disparo.
* En modo caza (sin barco objetivo activo) se restan las ventanas de todos los barcos no hundidos; en modo objetivo solo las del barco activo.
* Complejidad: O(barcos × longitud²) por miss.

**Impacto (`onHit`)**
* Se reinicia el heatmap (`Arrays.fill`) y se reconstruye exclusivamente para el barco relevante:
  * Si el barco fue **hundido**: se reconstruye para todos los barcos restantes (modo caza). Ocurre como máximo 5 veces por partida.
  * Si el barco **no fue hundido**: se reconstruye solo para ese barco, restringiendo el barrido a las ventanas que pasan por el primer hit conocido (región acotada a ≤ 2 × longitud celdas).

### 3. Modos de Operación

El resolvedor alterna dinámicamente entre dos modos según el estado del tablero:

**Modo Caza (Hunt Mode)**
* No hay impactos activos pendientes.
* El heatmap contiene contribuciones de todos los barcos no hundidos evaluadas sobre el tablero completo.
* Las celdas centrales acumulan naturalmente mayor densidad porque los barcos largos disponen de más espacio para acomodarse sin salir de los límites.

**Modo Objetivo (Target Mode)**
* Hay al menos un impacto activo (barco tocado pero no hundido).
* El heatmap contiene solo las contribuciones del barco objetivo, exigiendo que las ventanas contengan exactamente `N` impactos activos (donde `N` = hits acumulados en el barco).
* Esto genera un pico de probabilidad en las celdas adyacentes a los impactos, dirigiendo los disparos a lo largo de la orientación del barco.

### 4. Lógica de Hunt-Target Auxiliar

Tras actualizar el heatmap, se ejecuta `huntTarget` para reforzar las celdas vecinas del barco objetivo:
* **1 hit activo:** se elevan al valor máximo (`Byte.MAX_VALUE`) las celdas adyacentes (izquierda, derecha, arriba, abajo) que sean transitables y donde el barco pueda caber.
* **2+ hits activos:** se determina la orientación (horizontal si `minHit % DIM ≠ maxHit % DIM`, vertical en caso contrario) y se eleva al máximo solo el extremo libre en esa dirección.

### 5. Selección del Disparo

* Se selecciona la coordenada con el valor más alto en el heatmap (`getBestCoordFrom`).
* La celda seleccionada se marca con `Byte.MIN_VALUE` para no ser elegida de nuevo.

---

# Optimizaciones de Rendimiento

| Optimización | Detalle |
|---|---|
| `INITIAL_HEAT_MAP` estático | Calculado una vez, copiado con `System.arraycopy` en cada `reset()` |
| Actualización incremental en misses | Solo se restan las ventanas afectadas por la celda disparada, no todo el tablero |
| Reconstrucción acotada en hits | En modo objetivo se barren solo las ventanas que cruzan el primer hit conocido |
| Arreglos planos 1D | `int[]` y `byte[]` sin objetos, sin presión sobre el GC |
| Poda temprana | Los bucles de validación abortan en cuanto encuentran una celda bloqueada |
| `hMin`/`hMax` en un solo escaneo | `huntTarget` calcula el primer y último hit del barco en un único recorrido del tablero |

---

# Resultados (500 000 partidas, tablero 10×10)

```
| Estrategia   | Turnos Prom. | Mejor Juego | Peor Juego |
|--------------|:------------:|:-----------:|:----------:|
| Montecarlo   |    44,74     |     18      |     72     |
```
