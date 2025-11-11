# Proyecto Individual – Comparación de Pérdidas Térmicas (Cilindro vs. Esfera)
 
Pontificia Universidad Católica de Chile  
IIQ3843 – Procesamiento de Hidrógeno para Energías Sostenibles

Emilia Violic Montalba

## Introducción

Este proyecto compara las pérdidas térmicas en estado transitorio y estado estacionario, entre dos cuerpos de igual volumen:  
un cilindro circular recto y una esfera para distintas condiciones de intercambio térmico con el ambiente. 

El trabajo se basa en la ecuación de conducción de calor transitoria en coordenadas radiales, utilizando el Método de Líneas (MOL)
para resolver la ecuación de conducción unidimensional (1D) sujeta a convección en la superficie externa.

El estudio busca aportar criterios de diseño para estanques térmicos más eficientes y sostenibles en aplicaciones energéticas 
en Chile.

## 1. Ecuación de Conducción de Calor Transitoria

Modelo 1D radial para estudiar la evolución temporal de la temperatura en un estanque.

Para el estanque cilindro, se resuelve la ecuación:

$$
\rho c_p \frac{\partial T}{\partial t} = k \left[ \frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial T}{\partial r} \right) \right]
$$

Para el estanque esférico, se resuelve la ecuación:
$$
\rho c_p \frac{\partial T}{\partial t} = k \left[ \frac{1}{r^2} \frac{\partial}{\partial r} \left( r^2 \frac{\partial T}{\partial r} \right) \right]
$$


y los parámetros son:
- $k$ : conductividad térmica [W/m·K]  
- $\rho$ : densidad [kg/m³]  
- $c_p$ : calor específico [J/kg·K]  
- $T(r,t)$ : temperatura dependiente de radio y tiempo [K]

Estas ecuaciones están sujeta a condiciones de borde y una condición inicial.

A continuación, se detallan las condiciones aplicadas en el modelo:

**CB1:** Eje de Simetría (r = 0)

$$
\frac{\partial T}{\partial r}(r=0, t) = 0
$$

Esta condición establece que el gradiente de temperatura es nulo en el centro, lo que corresponde a una condición de simetría. 

**CB2:** Superficie Externa (r = rₒ)

$$
-k \frac{\partial T}{\partial r}(r=r_o, t) = h \, [T_s - T_\infty]
$$

Esta condición de borde (tipo Neumann) representa el intercambio de calor por convección entre la superficie del estanque y el aire del ambiente. 
El término del lado derecho corresponde a la Ley de Enfriamiento de Newton, donde el flujo de calor hacia el ambiente es proporcional a la diferencia de temperaturas entre la pared del estanque $T_s$ y la temperatura del aire $T_\infty$.  

**Condición Inicial (t = 0)**

$$
T(r, 0) = T_i
$$

Esta condición define que en el instante inicial todo el cuerpo se encuentra a una temperatura uniforme $T_i$.  
Esta condición es fundamental para establecer el punto de partida de la simulación y observar la evolución temporal.

## 2. Discretización de la Ecuación de Conducción de Calor  
### Método de Diferencias Finitas

La ecuación de conducción de calor transitoria en coordenadas cilíndricas se expresa como:

$$
\rho c_p \frac{\partial T}{\partial t} = k \left[ \frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial T}{\partial r} \right) \right]
$$

Al desarrollar la derivada espacial y definir el parámetro de **difusividad térmica** $\alpha = \frac{k}{\rho c_p}$, se obtiene:

$$
\frac{\partial T}{\partial t} = \alpha \left( \frac{\partial^2 T}{\partial r^2} + \frac{1}{r}\frac{\partial T}{\partial r} \right)
$$

Debido a que la ecuación es una EDP parabólica (hay una derivada temporal y una segunda derivada espacial de la variable de interés), se utiliza para resolver el **método de líneas**, que consiste en discretizar espacialmente las ecuaciones diferenciales parciales (EDPs) para convertirlas en un sistema de ecuaciones diferenciales ordinarias (EDOs), y resolverlo con un integrador numérico. Para discretizar se usa el método de diferencias finitas. 


#### a) Discretización espacial

El dominio radial se divide en N nodos uniformemente espaciados:

$$
r_i = i \, \Delta r, \quad \text{con } i = 0, 1, 2, \dots, N-1
$$

donde el tamaño de paso espacial es:

$$
\Delta r = \frac{R}{N-1}
$$

Esto permite reemplazar las derivadas espaciales por aproximaciones por diferencias finitas, lo que transforma la ecuación diferencial parcial en un conjunto de ecuaciones algebraicas.

#### b) Aproximaciones por Diferencias Finitas

Se utilizan aproximaciones de segundo orden en el espacio.

Segunda derivada central:
$$
\left( \frac{\partial^2 T}{\partial r^2} \right)_i \approx \frac{T_{i+1} - 2T_i + T_{i-1}}{(\Delta r)^2}
$$

Primera derivada central:
$$
\left( \frac{\partial T}{\partial r} \right)_i \approx \frac{T_{i+1} - T_{i-1}}{2\Delta r}
$$


#### c) Sustitución en la ecuación de conducción

Reemplazando las expresiones discretizadas en la ecuación diferencial original **cilíndrica** se obtiene:

$$
\frac{\partial T_i}{\partial t} = \alpha \left[ 
\frac{T_{i+1} - 2T_i + T_{i-1}}{(\Delta r)^2}
+ \frac{1}{r_i} \frac{T_{i+1} - T_{i-1}}{2\Delta r}
\right]
$$

Esta forma discretizada representa la evolución temporal de la temperatura en el nodo i, considerando la conducción radial entre sus nodos vecinos i-1 y i+1.

De manera análoga para la geometría **esférica** se obtiene la ecuación discretizada, 
$$
\frac{\partial T_i}{\partial t} = \alpha \left[ 
\frac{T_{i+1} - 2T_i + T_{i-1}}{(\Delta r)^2}
+ \frac{1}{r_i} \frac{T_{i+1} - T_{i-1}}{\Delta r}
\right]
$$



Este enfoque permite resolver el problema en forma **numérica** mediante un integrador temporal, como el `solve_ivp` de SciPy, dentro del esquema del **Método de Líneas**, donde el espacio se discretiza y el tiempo se integra de forma continua.

## 3. Método Numérico – Método de Líneas
### 3.1 Ecuación de Conducción de Calor
Una vez discretizada la ecuación de conducción de calor en el espacio, se obtiene una ecuación diferencial ordinaria (EDO) para cada nodo i:

$$
\frac{dT_i}{dt} = \alpha \left[ 
\frac{T_{i+1} - 2T_i + T_{i-1}}{(\Delta r)^2} 
+ \frac{1}{r_i} \frac{T_{i+1} - T_{i-1}}{2\Delta r}
\right]
$$

El sistema completo puede expresarse en forma matricial como:

$$
\frac{d\mathbf{T}}{dt} = A \mathbf{T} + \mathbf{b}
$$

donde:

- $ \mathbf{T} = [T_1, T_2, \ldots, T_N]^T $ es el vector de temperaturas en los nodos.  
- $ A $ es la matriz  que contiene los coeficientes de conducción en dirección radial.  
- $ \mathbf{b} $ es el vector que incluye los efectos de las condiciones de borde.

Agrupando los términos para cada nodo interior se obtiene:

$$
A_{i,i-1} = \alpha \left[\frac{1}{(\Delta r)^2} - \frac{1}{2r_i\Delta r}\right]
$$

$$
A_{i,i} = -2\alpha / (\Delta r)^2
$$

$$
A_{i,i+1} = \alpha \left[\frac{1}{(\Delta r)^2} + \frac{1}{2r_i\Delta r}\right]
$$

Estos coeficientes conforman la matriz A, que es tridiagonal porque cada nodo depende solo de sus vecinos inmediatos. 

Se puede expresar de manera análoga para la geometría esférica. 


### 3.2 Incorporación de las Condiciones de Borde e Iniciales


#### CB1: Eje de Simetría
En el eje del dominio (\(r = 0\)), se impone que el flujo de calor sea nulo, es decir, que el gradiente térmico desaparezca:

$$
\frac{\partial T}{\partial r}\bigg|_{r=0} = 0
$$

Partiendo de la ecuación de conducción en coordenadas cilíndricas:

$$
\frac{\partial T}{\partial t} = \alpha \left[ \frac{\partial^2 T}{\partial r^2} + \frac{1}{r}\frac{\partial T}{\partial r} \right]
$$

el término $ \frac{1}{r}\frac{\partial T}{\partial r} $ no puede evaluarse directamente en $r=0$.  
Aplicando el límite cuando $r \to 0$ y usando la condición de simetría, se obtiene:

$$
\frac{\partial T}{\partial t}\bigg|_{r=0} = 2\alpha \frac{\partial^2 T}{\partial r^2}\bigg|_{r=0}
$$

Utilizando diferencias finitas de segundo orden centradas:

$$
\frac{\partial^2 T}{\partial r^2}\bigg|_{r=0} \approx \frac{T_{-1} - 2T_0 + T_1}{(\Delta r)^2}
$$

y considerando la simetría $T_{-1} = T_1$, se obtiene:

$$
\frac{\partial T_0}{\partial t} = \alpha \left( \frac{-4T_0 + 4T_1}{(\Delta r)^2} \right)
$$

Por lo tanto, los coeficientes en la primera fila de la matriz $$A$ son:

$$
A_{0,0} = -\frac{4\alpha}{(\Delta r)^2}, \quad
A_{0,1} = \frac{4\alpha}{(\Delta r)^2}, \quad
b_0 = 0
$$

#### CB2: Convección en Borde Externo (Neumann)

Para el último nodo N, se aplica un balance de energía entre la conducción que llega desde el nodo interior y la convección que sale hacia el aire:

$$
\underbrace{\rho c_p \frac{dT_N}{dt} \Delta r}_{\text{Acumulación}} 
= \underbrace{k \frac{T_{N-1} - T_N}{\Delta r}}_{\text{Conducción entrante}} 
- \underbrace{h (T_N - T_\infty)}_{\text{Convección saliente}}
$$

Dividiendo por $\rho c_p \, \Delta r$ y reemplazando $\alpha = \frac{k}{\rho c_p}$, se obtiene la ecuación diferencial para el nodo externo:

$$
\frac{dT_N}{dt} = 
\frac{\alpha}{(\Delta r)^2} (T_{N-1} - T_N)
- \frac{h}{\rho c_p \Delta r} (T_N - T_\infty)
$$

Desarrollando los términos, puede escribirse como:

$$
\frac{dT_N}{dt} =
\left(\frac{\alpha}{(\Delta r)^2}\right) T_{N-1} 
+ \left(-\frac{\alpha}{(\Delta r)^2} - \frac{h}{\rho c_p \Delta r}\right) T_N
+ \frac{h T_\infty}{\rho c_p \Delta r}
$$

Por lo tanto, los coeficientes asociados en la última fila de la matriz $A$ y el vector fuente $\mathbf{b}$ son:

$$
A_{N,N-1} = \frac{\alpha}{(\Delta r)^2}, \quad
A_{N,N} = -\frac{\alpha}{(\Delta r)^2} - \frac{h}{\rho c_p \Delta r}, \quad
b_N = \frac{h T_\infty}{\rho c_p \Delta r}
$$

#### Condición Inicial
Se considera una temperatura uniforme en todo el dominio al inicio del proceso:

$$
T(r, 0) = T_i
$$

## 4. Implementación Numérica y Resolución del Sistema

### 4.1 Parámetros Geométricos y Físicos Utilizados


Los parámetros seleccionados para la simulación se basan en valores reportados en la literatura reciente sobre almacenamiento térmico a gran escala.  
Particularmente, Keçebaş et al. (*Journal of Energy Storage*, 2023) evaluaron un tanque esférico subterráneo de acero inoxidable destinado a almacenamiento estacional de energía térmica, considerando rangos de radios entre 0.25 y 1.5 m, conductividad térmica del acero de 16.2 W/m·K y densidad de 7.99 g/cm³. 

Los valores utilizados en el presente modelo $ k = 16.2 [W/m·K] $, $ \rho = 7990  [kg/m³] $, $ c_p = 500  [J/kg·K] $, y $ h = 10  [W/m²·K] $ son coherentes con estos rangos y corresponden a condiciones representativas de estanques metálicos industriales en contacto con aire o suelo.  

El coeficiente de convección ($ h = 10 [W/m²·K]$) se encuentra dentro del rango reportado por Dahash et al. (2021) para sistemas TES de baja temperatura (5–15 W/m²·K). 

Las temperaturas de operación se definieron considerando un escenario de almacenamiento térmico de media temperatura, típico de aplicaciones industriales donde el fluido caloportador (aceite térmico o sales fundidas) se encuentra a altas temperaturas.  Se estableció una temperatura inicial de 150 °C (423 K) y una temperatura ambiente de 20 °C (293 K) como entorno de enfriamiento natural. Estos valores son coherentes con los rangos utilizados en estudios industriales de almacenamiento de energía térmica a media temperatura, como los reportados por Dahash et al. (2021) y Keçebaş et al. (2023).

Los parámetros geométricos utilizados en la simulación se seleccionaron para modelar un estanque a una escala industrialmente relevante ($\mathbf{V = 20.0 \text{ m}^3}$) y con proporciones de diseño optimizadas para el almacenamiento de energía térmica (TES). 

La literatura de ingeniería térmica para estanques cilindricos de TES enfatiza la necesidad de una alta estratificación térmica para maximizar la eficiencia. Esto se logra mediante una relación de Largo a Diámetro ($L/D$) superior a 1. Se fijó la relación de $L/D = 1.5$ (lo que implica $L=3R$), resultando en un radio de $\mathbf{R_{\text{cilindro}} \approx 1.284 \text{ m}}$ y un largo de $\mathbf{L_{\text{cilindro}} \approx 3.852 \text{ m}}$. Este diseño representa el esfuerzo por optimizar la geometría frente a las pérdidas. El diseño de la esfera ($\mathbf{R_{\text{esfera}} \approx 1.684 \text{ m}}$) se calculó manteniendo el mismo volumen ($V=20.0 \text{ m}^3$) para la comparación de eficiencia. 


### 4.2 Método de Integración Temporal

Para resolver el sistema  se utilizó la función solve_ivp del paquete SciPy, la cual permite integrar sistemas de EDOs con distintos métodos numéricos. 

En este caso, se seleccionó el integrador BDF (Backward Differentiation Formula), apropiado para sistemas rígidos como el de la conducción de calor.

El integrador recibe la función derivada $d\mathbf{T}/dt$, los parámetros físicos y las condiciones iniciales. 
Se obtienen las temperaturas en todos los nodos del dominio para distintos instantes de tiempo.


### 4.3 Cálculo de la Pérdida de Calor Instantánea

A partir de la temperatura superficial $T_s(t)$, se calcula la tasa de pérdida de calor instantánea:

$$
\dot{Q}(t) = A_{\text{superficie}} \, h \, [T_s(t) - T_\infty]
$$

donde el área superficial depende de la geometría:

- Cilindro: $ A_c = 2 \pi R_c L_c $
- Esfera: $ A_e = 4 \pi R_e^2 $


### 4.4 Visualización y Análisis

Se generaron tres tipos de gráficos para analizar el comportamiento térmico:

1. **Perfil radial de temperatura** en un instante intermedio, mostrando cómo varía $T(r)$ de cada geometría.  
2. **Curva de enfriamiento superficial** $T_s(t)$, que compara como disminuye la temperatura superficial. 
3. **Pérdida de calor instantánea** $\dot{Q}(t)$, que permite evaluar la tasa de transferencia de energía hacia el ambiente.

Estos resultados permiten cuantificar cómo la geometría influye en la rapidez del enfriamiento y en la cantidad de calor disipada.

## 5. Resultados y Análisis Comparativo

El **Gráfico 1** muestra la distribución radial de temperatura en el tiempo (t ≈ 28 h) para ambas geometrías con el mismo volumen.  
Se observa que la **esfera mantiene mayores temperaturas internas** que el cilindro, evidenciando una **menor pérdida térmica global** debido a su menor relación área/volumen.  
El gradiente de temperatura radial es más pronunciado en el cilindro, indicando una conducción más significativa hacia el ambiente.


El **Gráfico 2** presenta la evolución temporal de la temperatura en la superficie exterior.  
Ambas geometrías tienden asintóticamente a la temperatura ambiente (T∞ = 293 K), pero la **esfera presenta un enfriamiento ligeramente más rápido al inicio** debido a su mayor área de exposición instantánea, alcanzando luego un equilibrio más estable que el cilindro.


El **Gráfico 3** compara la **tasa de pérdida de calor instantánea (Q)** para ambas geometrías.  
Inicialmente, la pérdida de calor es más alta en la esfera por su área superficial mayor, pero **disminuye más rápidamente**, reflejando una **mejor eficiencia térmica acumulada**.  El cilindro, en cambio, conserva pérdidas más sostenidas en el tiempo, lo que implica una disipación energética menos eficiente para un mismo volumen almacenado.

## 6. Conclusiones

- El método de líneas permite resolver de manera estable y precisa la ecuación de conducción transitoria en geometrías radiales.  
- La comparación entre cilindro y esfera evidencia que **la forma geométrica tiene un efecto directo en la eficiencia térmica**.  
- La **esfera** resulta más eficiente para conservar energía térmica, mientras que el **cilindro** facilita un enfriamiento más rápido.  
- Los resultados son coherentes con los principios físicos de transferencia de calor y proporcionan una base para el diseño de estanques térmicos más sostenibles.

