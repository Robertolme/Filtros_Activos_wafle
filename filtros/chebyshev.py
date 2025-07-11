#from .base import BaseFilter
import numpy as np
import json
from itertools import combinations, product

class chebyshevFilter():
    def __init__(self, order, fc, rp):
        self.order = order
        self.fc = fc
        self.rp = rp
        self.Sorder = self.get_order()
        self.ripple = self.get_repl()
        self.coeficientes = self.get_coeficientes()
        self.capacitor = [200e-9,100e-9, 47e-9, 10e-9, 4.7e-9, 1e-9, 680e-12, 300e-12, 200e-12, 100e-12, 68e-12]
        self.resitor = [1, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.1, 5.6, 6.8, 8.2, 10, 12.0, 15.0, 18.0, 22.0, 27.0, 33.0, 39.0, 47.0, 51.0, 56.0, 68.0, 82.0, 100, 120.0, 150.0, 180.0, 220.0, 270.0, 330.0, 390.0, 470.0, 510, 560.0, 680.0, 820, 1000.0, 1200.0, 1500.0, 1800.0, 2200.0, 2700.0, 3300.0, 3900.0, 4700.0, 5100.0, 5600.0, 6800.0, 8200.0, 10000.0, 12000.0, 15000.0, 18000.0, 22000.0, 27000.0, 33000.0, 39000.0, 47000.0, 51000.0, 56000.0, 68000.0, 82000.0, 100000.0, 120000.0, 150000.0, 180000.0, 220000, 270000.0, 330000.0, 390000.0, 470000.0, 510000, 560000.0, 680000.0, 820000, 1000000.0, 1200000.0, 1500000.0, 1800000.0, 2200000.0, 2700000.0, 3300000.0, 3900000.0, 4700000.0, 5100000.0, 5600000.0, 6800000.0, 8200000]
        self.counter = 1

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
        for i in sorted(self.coeficientes[self.ripple][self.Sorder]):
            print(f"Etapa {i}:\n")
            if(self.coeficientes[self.ripple][self.Sorder][i]["bi"] <= 0.0):
                response = self.calculate_best_resistor(self.calculate_first_order(self.coeficientes[self.ripple][self.Sorder][i]["ai"]))
                print(response)
            else:
                response = self.calculate_second_order(self.coeficientes[self.ripple][self.Sorder][i]["ai"],self.coeficientes[self.ripple][self.Sorder][i]["bi"])
                print(response)

    def calculate_first_order(self, a):
        response = []
        for i in self.capacitor:
            R1 = a / (2 * np.pi * self.fc * i)
            #print(f"C = {i}, R={R1}")
            response.append((R1,i))
        return response    

    def calculate_second_order(self, a, b):
        response_capacitor = []
        error = 0.01
        C1 = 300e-12
        comp = C1 * 4 * (b/(a**self.counter))
        self.counter += 1
        C2 = min([v for v in self.capacitor if v > comp], default=None)
        R1 = (a * C2 - np.sqrt((a**2 * C2 ** 2) - (4 * b *  C1 * C2)) ) / (4 * np.pi * self.fc * C1 * C2 )
        R2 = (a * C2 + np.sqrt((a**2 * C2 ** 2) - (4 * b *  C1 * C2)) ) / (4 * np.pi * self.fc * C1 * C2 )
        return (R1,R2,C1,C2)

    def calculate_best_resistor(self, resis): #(resistecia deseada, capacitor, error, aproximacion, r1, r2, r3, modo) 0=serie, 1= paralelo
        response=[]
        tola = 0.004
        for i in resis:
            #case 1
            for e in self.resitor:
                err = (abs(e - i[0]))/i[0]
                if(err < tola):
                    tola = err
                    response.append((i[0],i[1],err,rt,0,0,0,0))

            for r1, r2 in combinations(self.resitor, 2):
                rt = r1 + r2   #serie
                err = (abs(rt - i[0]))/i[0]
                if(err < tola):
                    tola = err
                    response.append((i[0],i[1],err,rt,r1,r2,0,0))

                rt = (r1 * r2) / (r1 + r2)  #paralelo
                err = (abs(rt - i[0]))/i[0]
                if(err < tola):
                    tola = err
                    response.append((i[0],i[1],err,rt,r1,r2,1))

            for r1, r2, r3 in combinations(self.resitor, 3):
                rt = r1 + r2 + r3
                err = (abs(rt - i[0]))/i[0]
                if(err < tola):
                    tola = err
                    response.append((i[0],i[1],err,rt,r1,r2,r3,0))


                rt = 1 / (1/r1 + 1/r2 + 1/r3)
                err = (abs(rt - i[0]))/i[0]
                if(err < tola):
                    tola = err
                    response.append((i[0],i[1],err,rt,r1,r2,r3,1))

            return response


    def calculate_best_capacitor(self, cap):
        response=[]
        tola = 0.01
        i = cap
        for e in self.capacitor:
            err = (abs(e - i))/i
            if(err < tola):
                tola = err
                response.append((i,err,e,0,0,0))

        for r1, r2 in combinations(self.capacitor, 2):
            rt = r1 + r2   #paralelo
            err = (abs(rt - i))/i
            if(err < tola):
                tola = err
                response.append((i,err,rt,r1,r2,1))

            rt = (r1 * r2) / (r1 + r2)  #serie
            err = (abs(rt - i))/i
            if(err < tola):
                tola = err
                response.append((i,err,rt,r1,r2,0))

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

