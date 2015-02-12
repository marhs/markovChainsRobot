#!/usr/local/bin/python3

# Trabajo AIA
# Localización de robots usando modelos ocultos de Markov. 
# Marco Herrero Serna - Grupo 2
# 

import copy, random

# Parte 1 - Implementacion de los algoritmos de muestreo, avance y Viterbi. 

class Markov:
  # Constructor
  def __init__(self, estados, observaciones, matrizTransicion, matrizInicial, matrizObservaciones):
    self.estados = estados
    self.observaciones = observaciones
    # Prob de ir desde a hasta b = self.probtransicion[a][b]
    self.probtransicion = matrizTransicion
    # Prob de empezar en a = self.probInicial[index(a)]
    self.probinicial = matrizInicial
    # Prob de una observacion o dado un estado e: matrizObservaciones[e][o]
    self.probobservaciones = matrizObservaciones
  
  # Algoritmo de muestreo
  # Devuelve una secuencia de estados y una secuencia de observaciones
  # [estados, observaciones]
  def muestreo(self, num):
    # Generacion del primer estado
    estados = []
    observaciones = []
    e_inicial = 0
    a, c = 0 , 0
    r = random.random()
    for n in self.probinicial:
      a = a + n
      if r < a:
        e_inicial = self.estados[c]
        break
      c = c + 1
    estados.append(e_inicial)

    # Generacion del resto de estados. 
    for n in range(1,num):
      r = random.random()
      estado_actual = self.estados.index(estados[n-1])
      b = 0
      l = len(estados)
      for m in range(len(self.estados)):
        b += self.probtransicion[estado_actual][m]
        if r < b:
          estados.append(self.estados[m])
          break
      if len(estados) == l:
        break


    # Generacion de observaciones
    for n in range(len(estados)):
      estado_actual = self.estados.index(estados[n])
      r = random.random()
      b = 0
      for m in range(len(self.observaciones)):
        b += self.probobservaciones[estado_actual][m]
        if r <= b:
          observaciones.append(self.observaciones[m])
          break

    return [estados,observaciones]
    

  # Algoritmo de avance (forward)
  def forward(self, observaciones):
    # Observaciones es una matriz de los INDICES de las observaciones obtenidas
    # Ej: [frio, lluvia] : [frio, frio, lluvia] = [0,0,1]
    a = []
    ob = self.obvToIndex(observaciones)
    # Alpha inicial
    for n in range(len(self.estados)):
      a.append(self.probobservaciones[n][ob[0]] * self.probinicial[n])

    # Iteraciones del algoritmo, calculo de los siguientes alphas. 
    for i in range(1, len(observaciones)):
      olda = copy.deepcopy(a)
      for k in range(len(self.estados)):
        x1 = self.probobservaciones[k][ob[i]]
        x2 = self.sumaColumnaAlpha(k,olda)
        a[k] = x1 * x2 
    return a
  
  def viterbi(self, observaciones):
    # De nuevo, observaciones es una matriz de los INDICES de las observaciones
    # [frio, lluvia] : [frio, frio, lluvia] = [0,0,1]

    res = [] # Aqui se almacenarán el resultado
    v = []   # [[v0(a),v0(b)],[v1(a),v1(b)]]
    pr = []  # Igual que el de arriba. 
    ob = self.obvToIndex(observaciones)

    # Primera iteración
    v.append([])
    pr.append([])
    for n in range(len(self.estados)):
      pr_sub = self.probobservaciones[n][ob[0]] * self.probinicial[n]
      v[0].append(pr_sub)
      pr[0].append(0)
    
    # Resto de iteraciones
    for i in range(1, len(ob)):
      v.append([])
      pr.append([])
      v_aux = []
      for j in range(len(self.estados)):
        aux = []
        for k in range(len(self.estados)):
          aux.append(v[i-1][k] * self.probtransicion[k][j])

        m = max(aux)
        v_aux.append(self.probobservaciones[j][ob[i]] * m)
        pr[i].append(aux.index(m))
      v[i] = v_aux
     
    # Paso final, salida de estados
    estado_final = v[len(v)-1].index(max(v[len(v)-1]))
    estados = [estado_final]
    pr.reverse()
    a = estado_final
    for n in range(len(pr)):
      a = pr[n][a]
      estados.append(a)
    pr.reverse()
    estados.reverse()
    return self.indexAEstados(estados[1:])

  # Calcula la suma de las columnas de la matriz de transicion*alpha (para algoritmo de forwarding)
  def sumaColumnaAlpha(self,k,alpha):
    res = 0
    for a in range(len(self.estados)):
      x1 = self.probtransicion[a][k] * alpha[a]
      res = res + x1
    return res

  # Devuelve una matriz de los indices de las observaciones, ya que trabajamos con indices. 
  def obvToIndex(self, mObs):
    res = []
    for n in mObs:
      res.append(self.observaciones.index(n))
    return res

  def indexAEstados(self,index):
    res = []
    for n in index:
      res.append(self.estados[n])
    return res


# Parte 2 - Modelo oculto de Markov para un problema de localización de robots. 

class MarkovRobot(Markov):
  
  # Dado un tablero de la forma 
  # [[0,1,0,1,1,0,...],
  #  [...
  #  
  # Siendo los 0 espacios libres y los 1 espacios ocupados. 

  def __init__(self, tablero, error):
    # Generación de los estados y observaciones
    estados = []
    for i in range(len(tablero)):
      for j in range(len(tablero[0])):
        if tablero[i][j] == 0:
          estados.append((i,j))
    observaciones = ['','N','S','E','O','NS','NE','NO','SE','SO','EO',
                          'NSE','NSO','NEO','SEO','NSEO']
    prob_inicial = 1/len(estados)

    # Generacion del vector de probabilidades iniciales
    probinicial = len(estados) * [prob_inicial]

    # Generacion de la matriz de transición
    probtransicion = []
    for n in estados:
      a = []
      for n in estados:
        a.append(0)
      probtransicion.append(a)

    for n in estados:
      vec = self.vecinos(n[0],n[1],tablero)
      if len(vec) != 0:
        prob = 1/len(vec)
      else:
        prob = 0
      for vecino in vec:
        x = estados.index(n)
        y = estados.index(vecino)
        probtransicion[x][y] = prob

    probobservaciones = []
    for estado in estados:
      a = []
      for observacion in observaciones:
        rep = self.prob_observacion(error, self.direcciones(estado[0],estado[1], tablero), observacion)
        a.append(rep)
      probobservaciones.append(a)

    Markov.__init__(self, estados, observaciones, probtransicion, probinicial, probobservaciones)

  # Mira si unas coordenadas están dentro del tablero. 
  def bordes(self, x, y, tablero):
    return not(x < 0 or y < 0 or x > len(tablero)-1 or y > len(tablero[0])-1)

  # Dada una casilla del tablero, devuelve los vecinos disponibles que hay. 
  def vecinos(self, x, y, tablero):
    res = []
    for m,n in [(-1,0),(0,-1),(0,1),(1,0)]:
      if self.bordes(x+m,y+n,tablero) and tablero[x+m][y+n] == 0:
        res.append((x+m,y+n))
    return res

  # Devuelve un string con las posibles direcciones
  def direcciones(self, x, y, tablero):
    res = []
    for m,n in [(-1,0),(0,-1),(0,1),(1,0)]:
      if self.bordes(x+m,y+n,tablero) and tablero[x+m][y+n] == 0:
        if (m,n) == (-1,0):
          res.append('N')
        elif (m,n) == (1,0):
          res.append('S')
        elif (m,n) == (0,-1):
          res.append('O')
        else:
          res.append('E')
    return res

  # Devuelve la probabilidad de una determinada observacion dada el error. 
  def prob_observacion(self, error, direcciones, observacion):
    res = 1
    for n in 'NSEO':
      if (n in observacion and n in direcciones) or ((not n in observacion) and (not n in direcciones)):
        res = res * (1-error)
      else:
        res = res * error
    return res
  
  # Hace viterbi en el robot
  def robot_viterbi(self, observaciones):
    return self.viterbi(observaciones)

  # Hace el algoritmo de avance en el robot (devolviendo la casilla, no una
  # matriz de alphas. 
  def robot_forward(self, observaciones):
    # Aqui usamos forward para calcular la casilla más probable
    res = self.forward(observaciones)
    return self.estados[res.index(max(res))]
  
  # Ejecuta muestreo de robots
  def robot_muestreo(self, number):
    return self.muestreo(number)

# Parte 3 - Experimentación con la implementación

# Distancia manhattan entre dos puntos
def dis_manhattan(a,b):
  return abs(a[0] - b[0]) + abs(a[1] - b[1])

# Devuelve el porcentaje de aciertos dada una lista respecto a otra
def porcentajeAciertos(secuenciaOrig, secuenciaNueva):
  res = 0.0
  for n in range(len(secuenciaOrig)):
    if secuenciaOrig[n] == secuenciaNueva[n]:
      res += 1.0
  return res/len(secuenciaOrig)

# Media aritmética
def media(lista):
  res = 0
  for n in lista:
    res += n
  return res/len(lista)

# Estadísticas para avance, dado un numero de experimentos
def stat_forward(tablero, error, inicio, final, experimentos):
  res = []
  rob_markov = MarkovRobot(tablero,error)
  for n in range(inicio, final+1):
    distancia = []
    for exp in range(experimentos):
      muest = rob_markov.muestreo(n)
      forw  = rob_markov.robot_forward(muest[1])
      distancia.append(dis_manhattan(muest[0][len(muest[0])-1], forw))
    res.append(media(distancia))
  return res

# Estadísticas para Viterbi, dado un número de experimentos
def stat_viterbi(tablero, error, inicio, final, experimentos):
  res = []
  rob_markov = MarkovRobot(tablero,error)
  for n in range(inicio, final+1):
    porcentaje = []
    for exp in range(experimentos):
      muest = rob_markov.muestreo(n)
      viterbi = rob_markov.robot_viterbi(muest[1])
      porcentaje.append(porcentajeAciertos(muest[0], viterbi))
    res.append(media(porcentaje))
  return res

# Estadísticas para avance, cambiando el error y con una longitud de observaciones
def stat_error_forward(tablero, observaciones, experimentos):
  res = []
  for n in range(0,11):
    rob_markov = MarkovRobot(tablero,n/10)
    distancia = []
    for exp in range(experimentos):
      muest = rob_markov.muestreo(observaciones)
      forw  = rob_markov.robot_forward(muest[1])
      distancia.append(dis_manhattan(muest[0][len(muest[0])-1], forw))
    res.append(media(distancia))
  return res

# Estadísticas para viterbi, cambiando el error y con una longitud de observaciones
def stat_error_viterbi(tablero, observaciones, experimentos):
  res = []
  for n in range(0,11):
    rob_markov = MarkovRobot(tablero,n/10)
    porcentaje = []
    for exp in range(experimentos):
      muest = rob_markov.muestreo(observaciones)
      viterbi = rob_markov.robot_viterbi(muest[1])
      porcentaje.append(porcentajeAciertos(muest[0], viterbi))
    res.append(media(porcentaje))
  return res



"""
Zona de ejecucion

Ejemplo 2 del tema 5 de teoría. 
"""
#matriz = [[0.7 , 0.3], 
          #[0.3 , 0.7]]
#observaciones = ['u','nou']
#matrizObservaciones = [[0.9,0.1],
                       #[0.2,0.8]]
#estados = ['l','nol']
#matrizInicial = [0.8, 0.2]

#tema5ejemplo2 = Markov(estados,observaciones, matriz, matrizInicial,matrizObservaciones)
#a = tema5ejemplo2.muestreo(20)[0]
## Cuadrícula del enunciado del trabajo. 

#robot =  [[0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0],
          #[1,1,0,0,1,0,1,1,0,1,0,1,0,1,1,1],
          #[1,0,0,0,1,0,1,1,0,0,0,0,0,1,1,0],
          #[0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0]]

#a = MarkovRobot(robot,1)
#print(a.muestreo(4))
 
