#from .base import BaseFilter
import numpy as np
import json
from itertools import combinations, product
import matplotlib.pyplot as plt
from scipy import signal

class chebyshevFilter():
    def __init__(self, order, fc, rp):
        self.order = order
        self.fc = fc
        self.rp = rp
        self.Sorder = self.get_order()
        self.ripple = self.get_repl()
        self.coeficientes = self.get_coeficientes()
        self.capacitor = [68e-12, 100e-12, 200e-12, 300e-12, 680e-12, 1e-9, 4.7e-9, 10e-9, 47e-9, 100e-9, 200e-9]
        self.resitor = [1, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.1, 5.6, 6.8, 8.2, 10, 12.0, 15.0, 18.0, 22.0, 27.0, 33.0, 39.0, 47.0, 51.0, 56.0, 68.0, 82.0, 100, 120.0, 150.0, 180.0, 220.0, 270.0, 330.0, 390.0, 470.0, 510, 560.0, 680.0, 820, 1000.0, 1200.0, 1500.0, 1800.0, 2200.0, 2700.0, 3300.0, 3900.0, 4700.0, 5100.0, 5600.0, 6800.0, 8200.0, 10000.0, 12000.0, 15000.0, 18000.0, 22000.0, 27000.0, 33000.0, 39000.0, 47000.0, 51000.0, 56000.0, 68000.0, 82000.0, 100000.0, 120000.0, 150000.0, 180000.0, 220000, 270000.0, 330000.0, 390000.0, 470000.0, 510000, 560000.0, 680000.0, 820000, 1000000.0, 1200000.0, 1500000.0, 1800000.0, 2200000.0, 2700000.0, 3300000.0, 3900000.0, 4700000.0, 5100000.0, 5600000.0, 6800000.0, 8200000]
        self.E12 = [1e-12, 1.2e-12, 1.5e-12, 1.8e-12, 2.2e-12, 2.7e-12, 3.3e-12, 3.9e-12, 4.7e-12, 5.6e-12, 6.8e-12, 8.2e-12, 1e-11, 1.2e-11, 1.5e-11, 1.8e-11, 2.2e-11, 2.7e-11, 3.3e-11, 3.9e-11, 4.7e-11, 5.6e-11, 6.8e-11, 8.2e-11, 1e-10, 1.2e-10, 1.5e-10, 1.8e-10, 2.2e-10, 2.7e-10, 3.3e-10, 3.9e-10, 4.7e-10, 5.6e-10, 6.8e-10, 8.2e-10, 1e-09, 1.2e-09, 1.5e-09, 1.8e-09, 2.2e-09, 2.7e-09, 3.3e-09, 3.9e-09, 4.7e-09, 5.6e-09, 6.8e-09, 8.2e-09, 1e-08, 1.2e-08, 1.5e-08, 1.8e-08, 2.2e-08, 2.7e-08, 3.3e-08, 3.9e-08, 4.8e-08, 5.6e-08, 6.8e-08, 8.2e-08, 1e-07, 1.2e-07, 1.5e-07, 1.8e-07, 2.2e-07, 2.7e-07, 3.3e-07, 3.9e-07, 4.7e-07, 5.6e-07, 6.8e-07, 8.2e-07, 1e-06, 1.2e-06, 1.5e-06, 1.8e-06, 2.2e-06, 2.7e-06, 3.3e-06, 3.9e-06, 4.7e-06, 5.6e-06, 6.8e-06, 8.2e-06]
        self.response = self.calculate_filters_values()

    def get_coeficientes(self):
        with open('filtros/table.json') as f:
            return json.load(f)

    def get_repl(self):
        if(self.rp == 0.5 or self.rp == 1 or self.rp == 2 or self.rp == 3):
            return f"Chebyshev_{self.rp}dB"

    def get_order(self):
        if(self.order < 10):
            return f"{self.order}"

    def calculate_filters_values(self):
        results = []
        for i in sorted(self.coeficientes[self.ripple][self.Sorder]):
            print(f"Etapa {i}:\n")
            if(self.coeficientes[self.ripple][self.Sorder][i]["bi"] <= 0.0):
                results.append(self.calculate_first_order(self.coeficientes[self.ripple][self.Sorder][i]["ai"]))    
            else:
                results.append(self.calculate_second_order(self.coeficientes[self.ripple][self.Sorder][i]["ai"],self.coeficientes[self.ripple][self.Sorder][i]["bi"]))
        return results
    def calculate_first_order(self, a): #retorna 
        response = []
        error = 0.004
        for i in self.capacitor:
            R1 = a / (2 * np.pi * self.fc * i)
            try:
                R = self.calculate_best_resistor(R1,error)
                if(i > 700e-12 and R[-1][2] > 500):
                    response.append((R1,i,R[-1][1],R[-1][2], R[-1][3], R[-1][4], R[-1][5],R[-1][6],1))
            except:
                print("Sin valor valido") 
        return response    

    def calculate_second_order(self, a, b): #retorna R1, R2, C1, C2, error en cap 2 ,C2real, C2_1, C2_2, C2_3, Arreglo
        response = []
        error = 0.1
        for i in self.capacitor:
            comp = i * 4 * (b/(a**2))
            C2_teoric = min([v for v in self.E12 if v > comp], default=None)
            try:
                C2 = self.calculate_best_capacitor(C2_teoric, error)
                r = (a * a * C2[-1][0] * C2[-1][0]) - (4 * b *  i * C2[-1][0])
                #print(f"{i}:{comp}:{r}:{C2[-1][1]}")
                if(i>700e-12 and r > 0):
                    R1 = (a * C2[-1][0] - np.sqrt(r) ) / (4 * np.pi * self.fc * i * C2[-1][0] )
                    R2 = (a * C2[-1][0] + np.sqrt(r) ) / (4 * np.pi * self.fc * i * C2[-1][0] )
                    response.append((R1, R2, i, C2[-1][0],C2[-1][1],C2[-1][2],C2[-1][3],C2[-1][4],C2[-1][5],C2[-1][6],2))
                    error = C2[-1][1]
            except:
                print("Valor no encontrado ")
        #print(response)
        return response

    def calculate_best_resistor(self, resis,tola = 0.004): #(resistecia deseada, capacitor, error, aproximacion, r1, r2, r3, modo) 0=serie, 1= paralelo
        response=[]
        i = resis
        for e in self.resitor:
            err = (abs(e - i))/i
            if(err < tola):
                tola = err
                response.append((i,err,e,e,0,0,0))

        for r1, r2 in combinations(self.resitor, 2):
            rt = r1 + r2   #serie
            err = (abs(rt - i))/i
            if(err < tola):
                tola = err
                response.append((i,err,rt,r1,r2,0,0))

            rt = (r1 * r2) / (r1 + r2)  #paralelo
            err = (abs(rt - i))/i
            if(err < tola):
                tola = err
                response.append((i,err,rt,r1,r2,0,1))

        for r1, r2, r3 in combinations(self.resitor, 3):
            rt = r1 + r2 + r3
            err = (abs(rt - i))/i
            if(err < tola):
                tola = err
                response.append((i,err,rt,r1,r2,r3,0))


            rt = 1 / (1/r1 + 1/r2 + 1/r3)
            err = (abs(rt - i))/i
            if(err < tola):
                tola = err
                response.append((i,err,rt,r1,r2,r3,1))

        return response


    def calculate_best_capacitor(self, cap, tola = 0.1):
        response=[]
        i = cap
        for e in self.capacitor:
            err = (abs(e - i))/i
            if(err < tola):
                tola = err
                response.append((i,err,e,e,0,0,0))

        for r1, r2 in combinations(self.capacitor, 2):
            rt = r1 + r2   #paralelo
            err = (abs(rt - i))/i
            if(err < tola):
                tola = err
                response.append((i,err,rt,r1,r2,0,1))

            rt = (r1 * r2) / (r1 + r2)  #serie
            err = (abs(rt - i))/i
            if(err < tola):
                tola = err
                response.append((i,err,rt,r1,r2,0,0))

        for r1, r2, r3 in combinations(self.capacitor, 3):
            rt = r1 + r2 + r3
            err = (abs(rt - i))/i
            if(err < tola):
                tola = err
                response.append((i,err,rt,r1,r2,r3,1))


            rt = 1 / (1/r1 + 1/r2 + 1/r3)
            err = (abs(rt - i))/i
            if(err < tola):
                tola = err
                response.append((i,err,rt,r1,r2,r3,0))

        return response


    def print_results():
        djhsk

