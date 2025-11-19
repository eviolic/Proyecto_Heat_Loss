# Proyecto Individual – Comparación de Pérdidas Térmicas (Cilindro vs. Esfera)
 
Pontificia Universidad Católica de Chile  
IIQ3843 – Procesamiento de Hidrógeno para Energías Sostenibles

Emilia Violic Montalba

## Introducción

El almacenamiento de energía corresponde al proceso mediante el cual se guarda energía, en forma térmica, eléctrica, química o mecánica, para utilizarla posteriormente según la demanda. Este principio es similar al funcionamiento de una batería recargable o la pila de un automóvil: la energía se almacena en un momento y se libera cuando es necesaria. En los últimos años, este tipo de sistemas ha adquirido un rol estratégico en el sector energético, ya que permite complementar fuentes renovables variables como la solar y la eólica.

Chile posee uno de los mayores potenciales solares del mundo, con más de 1.100 GW estimados en energía solar y cerca de 190 GW en energía eólica. Sin embargo, ambas dependen de condiciones climáticas que no siempre son constantes. Esto hace que el almacenamiento de energía sea fundamental para asegurar un suministro eléctrico continuo, estable y compatible con los compromisos de descarbonización del país.

Entre las distintas tecnologías de almacenamiento, el almacenamiento térmico destaca por su simplicidad, bajo costo y versatilidad. Consiste en acumular energía en forma de calor mediante la elevación de temperatura de un medio de almacenamiento (como agua, sales fundidas o materiales con cambio de fase) y liberarla posteriormente para generar calor o electricidad. Este es el principio utilizado por las plantas solares de concentración (CSP), donde tanques de sales fundidas permiten producir energía incluso durante la noche.

No obstante, uno de los desafíos más relevantes en la operación de estos sistemas es minimizar las pérdidas térmicas hacia el ambiente. En varias plantas termosolares, más del 30% de la energía capturada se pierde como calor no aprovechado, principalmente debido a estanques con diseños térmicamente subóptimos. Optimizar la geometría y las condiciones de intercambio térmico de estos estanques permite aumentar la eficiencia global, reducir costos operacionales y mejorar la competitividad de las tecnologías solares en Chile.

En este contexto, el presente proyecto analiza cómo la **geometría de los estanques** influye en las pérdidas de calor. En particular, se comparan dos cuerpos de igual volumen:  
1. un **cilindro circular recto**, y  
2. una **esfera**,  
evaluando su comportamiento tanto en estado estacionario como en régimen transitorio, bajo distintas condiciones de transferencia de calor con el ambiente.

El análisis se basa en las ecuaciones de conducción de calor en coordenadas radiales, derivadas a partir del balance de energía y los principios de transferencia de calor en fenómenos de transporte. Estas ecuaciones se implementan utilizando el **Método de Líneas (MOL)**, discretizando el dominio espacial y resolviendo la evolución temporal mediante un sistema de ecuaciones diferenciales ordinarias.

El objetivo final es identificar qué geometría presenta un desempeño térmico más eficiente y bajo qué condiciones ambientales, aportando criterios de diseño aplicables al desarrollo de tecnologías de almacenamiento térmico más sostenibles y adecuadas al contexto energético chileno.


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

  
El **método de diferencias finitas** permite aproximar las derivadas espaciales de la ecuación de conducción mediante relaciones algebraicas discretas entre nodos, transformando la EDP original en un conjunto de EDOs en función del tiempo.  

A este procedimiento se le conoce como el **método de líneas**, ya que discretiza solo la variable espacial mientras mantiene continua la variable temporal.  De esta manera, el problema físico se reduce a integrar numéricamente la evolución temporal de las temperaturas en cada nodo. En este caso, como la ecuación es una EDP parabólica (hay una derivada temporal y una segunda derivada espacial de la variable de interés), se puede utilizar el método de líneas. 

Para resolver este sistema de EDOs se utiliza el integrador `solve_ivp` de SciPy, el cual aplica métodos adaptativos (como BDF o Runge–Kutta) para obtener una solución estable y precisa.  El esquema **BDF** fue elegido por su robustez frente a ecuaciones **rígidas**, típicas en problemas de transferencia de calor transitorio.

Además, el tamaño del paso espacial y del paso temporal  determinan la **resolución y estabilidad numérica** del método: pasos más pequeños mejoran la precisión de la solución, aunque aumentan el costo computacional.


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

La literatura de ingeniería térmica para estanques cilindricos de TES enfatiza la necesidad de una alta estratificación térmica para maximizar la eficiencia. Esto se logra mediante una relación de Largo a Diámetro ($L/D$) superior a 1. Se fijó la relación de $L/D = 1.5$ (lo que implica $L=3R$), resultando en un radio de $\mathbf{R_{\text{cilindro}} \approx 1.284 \text{ m}}$ y un largo de $\mathbf{L_{\text{cilindro}} \approx 3.855 \text{ m}}$. Este diseño representa el esfuerzo por optimizar la geometría frente a las pérdidas. El diseño de la esfera ($\mathbf{R_{\text{esfera}} \approx 1.684 \text{ m}}$) se calculó manteniendo el mismo volumen ($V=20.0 \text{ m}^3$) para la comparación de eficiencia. 


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

Se generaron distintos tipos de gráficos para analizar el comportamiento térmico:


**Gráfico 1 – Perfiles radiales de temperatura a t ≈ 28,1 h**

El Gráfico 1 muestra la distribución radial de temperatura para el cilindro y la esfera, ambos con igual volumen, en un instante intermedio del enfriamiento (t ≈ 28,1 h). En ambos casos la temperatura disminuye desde el centro hacia la pared, pero se observa que **la esfera mantiene temperaturas internas más altas que el cilindro en todo el dominio radial**, lo que indica que ha perdido menos calor hacia el ambiente.

El **gradiente de temperatura** es más pronunciado en el cilindro (curva más inclinada), evidenciando una conducción de calor más intensa hacia la superficie expuesta. En cambio, la esfera presenta un perfil más suave, consistente con su **menor relación área/volumen**, que reduce las pérdidas de calor por convección. Esto sugiere que, para un mismo volumen de almacenamiento, la esfera es térmicamente más eficiente que el cilindro en retener energía.

**Gráfico 2 – Curva de enfriamiento superficial en el tiempo**

El Gráfico 2 muestra la evolución temporal de la temperatura en la superficie exterior para ambas geometrías, comparadas con la temperatura ambiente (T∞ = 293 K). Ambas curvas presentan un comportamiento típico de enfriamiento transitorio: una caída rápida al inicio seguida de una aproximación asintótica hacia el ambiente.

Se observa que **la esfera se enfría ligeramente más rápido en las primeras horas**, lo cual se explica por su mayor área superficial inicial expuesta a convección. Sin embargo, conforme avanza el tiempo, la tasa de enfriamiento se reduce y la esfera **alcanza temperaturas cercanas al ambiente antes que el cilindro**, estabilizándose en un régimen más uniforme.

El cilindro, en contraste, presenta una superficie menor para el mismo volumen, por lo que **su enfriamiento es más lento y sostenido**. Esto implica que parte de la energía interna se disipa de manera más gradual, resultando en una respuesta térmica más extendida en el tiempo.

En conjunto, este gráfico evidencia cómo la geometría influye en la rapidez del enfriamiento superficial, destacando que la esfera se adapta más rápidamente al medio externo, mientras que el cilindro retiene su temperatura por más tiempo.

**Gráfico 3 – Tasa de pérdida de calor instantánea en el tiempo**

El Gráfico 3 muestra la evolución de la tasa de pérdida de calor instantánea hacia el ambiente para ambas geometrías. En las primeras horas, la esfera presenta una **tasa de pérdida de calor ligeramente mayor** que el cilindro, coherente con su mayor área superficial para un mismo volumen y, por lo tanto, con una transferencia convectiva más intensa en el inicio del enfriamiento.

A medida que avanza el tiempo, ambas curvas descienden y se aproximan  al equilibrio térmico. No obstante, la tasa de disipación de la esfera disminuye más rápidamente, mientras que el cilindro mantiene pérdidas más sostenidas en el tiempo. Este comportamiento indica que la esfera refleja una mayor eficiencia térmica. En contraste, el cilindro retiene temperaturas internas más bajas durante más tiempo, lo que implica en una disipación continua y más prolongada para el mismo volumen almacenado.

**Gráfico 4 – Sensibilidad de la pérdida de calor instantánea a distintas condiciones ambientales**

El Gráfico 4 muestra cómo varía la tasa de pérdida de calor instantánea cuando se modifican las condiciones ambientales para cada geometría:  
- **Caso base** con convección natural.  
- **Convección forzada leve** (h = 25 W/m²K), que incrementa el intercambio térmico.  
- **Ambiente frío** (T∞ = 0°C), que aumenta el gradiente térmico entre el estanque y el entorno.

En ambas geometrías, un mayor coeficiente convectivo o un ambiente más frío generan pérdidas de calor más elevadas, especialmente en las primeras horas. Sin embargo, el efecto relativo depende de la forma del cuerpo:

- **Cilindro:**  
  La convección forzada produce la mayor disipación inicial, mientras que el ambiente frío incrementa la pérdida de calor principalmente en etapas intermedias. La diferencia entre casos persiste por más tiempo debido a su mayor razón superficie/volumen efectiva en dirección radial.

- **Esfera:**  
  La convección forzada también genera la mayor pérdida instantánea inicial, pero la esfera reduce su tasa de disipación más rápidamente en todos los escenarios. El ambiente frío tiene un efecto más moderado que en el cilindro, lo que refleja nuevamente su **menor sensibilidad térmica geométrica**.

En conjunto, este gráfico evidencia que, bajo escenarios ambientales más exigentes, la esfera mantiene un desempeño más estable, mientras que el cilindro presenta variaciones más pronunciadas en su disipación térmica. Esto refuerza la importancia de considerar la geometría al diseñar estanques térmicos para aplicaciones reales sujetas a fluctuaciones climáticas.

## 5. Aplicación en Chile (IPCh)

El almacenamiento térmico es una tecnología particularmente relevante para Chile, cuya matriz energética se está expandiendo aceleradamente hacia fuentes renovables variables, especialmente la solar. En zonas como el Desierto de Atacama se registran algunos de los niveles de radiación más altos del mundo, lo que convierte al país en un candidato natural para sistemas termosolares de concentración (CSP) y almacenamiento de calor en estanques a alta temperatura.

Los resultados de este proyecto entregan criterios útiles para el diseño de estanques térmicos en este contexto. En particular, se observa que las pérdidas energéticas dependen fuertemente de la geometría del estanque, siendo la **esfera** más eficiente que el cilindro para un mismo volumen debido a su menor razón área/volumen. Esto implica que, para aplicaciones solares en campo chileno, las geometrías esféricas podrían contribuir a reducir pérdidas nocturnas, aumentar la estabilidad térmica y mejorar la eficiencia del sistema durante periodos de baja radiación.

Además, la sensibilidad del desempeño térmico frente a condiciones ambientales más exigentes, como temperaturas ambiente bajas o mayores coeficientes convectivos por viento, es relevante para zonas del norte chileno donde existen fuertes gradientes térmicos entre el día y la noche. En estos escenarios, la esfera demuestra una menor variabilidad térmica, lo que la hace más adecuada para plantas ubicadas en zonas desérticas con noches frías y vientos irregulares.

Finalmente, al mejorar el diseño de estanques de almacenamiento, esta tecnología puede apoyar los compromisos nacionales de descarbonización, permitiendo una operación más flexible de plantas solares térmicas y reduciendo la dependencia de respaldo fósil en el sistema eléctrico.

## 6. Conclusiones y Cierre

Este proyecto analizó cómo la geometría del estanque afecta la retención de energía térmica bajo un modelo transitorio de conducción radial con convección superficial, resuelto mediante el **Método de Líneas (MOL)**. Se compararon dos cuerpos de igual volumen —un cilindro y una esfera— evaluando perfiles de temperatura, curvas de enfriamiento y tasas de pérdida de calor bajo distintas condiciones ambientales.

Los resultados muestran de manera consistente que la **esfera presenta un comportamiento térmico superior**, manteniendo temperaturas internas más altas, disipando calor con mayor uniformidad y siendo menos sensible a variaciones del entorno. Su menor relación área/volumen la convierte en la opción más eficiente para reducir pérdidas energéticas en estanques de almacenamiento térmico.

El cilindro, aunque ampliamente utilizado en la industria por razones de manufactura y facilidad constructiva, presenta gradientes radiales más pronunciados y mayores pérdidas acumuladas. Esto sugiiere que, cuando el rendimiento térmico es prioritario, considerar alternativas geométricas puede mejorar sustancialmente la eficiencia global del sistema.

Desde el punto de vista metodológico, la implementación del **Método de Líneas** permitió transformar la ecuación de conducción transitoria en un sistema de EDOs numéricamente estable y flexible de resolver. Esta técnica facilitó analizar el comportamiento térmico con distintas geometrías y condiciones de borde, mostrando su capacidad para capturar efectos físicos relevantes como gradientes radiales, sensibilidad ambiental y regímenes transitorios, con una buena relación entre precisión y costo computacional. Su uso demuestra la importancia de las herramientas numéricas modernas en el diseño y evaluación de sistemas térmicos antes de su construcción física.

En conjunto, este estudio muestra que la geometría es un parámetro de diseño crítico para estanques térmicos y que la modelación numérica basada en MOL es una herramienta poderosa para evaluar su desempeño. Sus conclusiones pueden servir como base para el desarrollo de sistemas de almacenamiento energético más eficientes y sostenibles, especialmente relevantes para el contexto solar térmico chileno.

## Referencias

- Chile, G.d. Almacenamiento. 2024  28 de octubre de 2025]; Available from: https://generadoras.cl/almacenamiento/.
- (ETSAP), I.R.E.A.I.E.T.S.A.P., Thermal Energy Storage. 2013.
- Bird, R.B., W.E. Stewart, and E.N. Lightfoot, Transport phenomena. Revised second edition ed. 2007, New York: J. Wiley & Sons.
- Çengel, Y.A. & Ghajar, A.J. (2015). Heat and Mass Transfer: Fundamentals and Applications (5th ed.). McGraw-Hill Education.
- Keçebaş, A., Kayfeci, M., & Gedik, E. (2023). *Design and performance analysis of a spherical underground thermal energy storage tank for seasonal heat storage applications.* Journal of Energy Storage, 65, 107607.
- Dahash, A., Ochs, F., Tosatto, A., & Streicher, W. (2021). *Techno-economic and exergy analysis of tank and pit thermal energy storage for renewable district heating systems.* Renewable Energy, 180, 1358–1379. 