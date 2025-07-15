#from .base import BaseFilter
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.units import inch
from reportlab.lib import colors
from datetime import datetime

from filtros.base import BaseFilter

class chebyshevFilter(BaseFilter):
    def __init__(self, order, fc, rp):
        self.order = order
        self.fc = fc
        super().__init__(fc)
        self.rp = rp
        self.Sorder = self.get_order()
        self.ripple = self.get_repl()
        self.response = self.calculate_filters_values()

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

    def print_results(self):#R1, R2, C1, C2
        for i in sorted(self.coeficientes[self.ripple][self.Sorder]):
            s0 = 1
            if(self.response[int(i)-1][-1]["orden"] == 1):
                s1 = self.response[int(i)-1][-1]["R1"] * (self.response[int(i)-1][-1]["C1"])
                s2 = 0     
            else:
                s1 = self.response[int(i)-1][-1]["C1"] * ((self.response[int(i)-1][-1]["R1"])+ (self.response[int(i)-1][-1]["R2"]))
                s2 = self.response[int(i)-1][-1]["R1"] * (self.response[int(i)-1][-1]["R2"])*(self.response[int(i)-1][-1]["C1"])*(self.response[int(i)-1][-1]["C2"])

            num = [1]          # Numerador
            den = [s2, s1, s0]    # Denominador
            system = signal.TransferFunction(num, den)
            w, mag, phase = signal.bode(system, n=10000)

            plt.figure(figsize=(15, 12))
            plt.subplot(2, 1, 1)
            plt.semilogx(w/(2*np.pi), mag)
            plt.title('Bode Plot')
            plt.ylabel('Magnitude (dB)')
            plt.grid(True, which="both", ls="--")

            plt.subplot(2, 1, 2)
            plt.semilogx(w/(2*np.pi), phase)
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('Phase (degrees)')
            plt.grid(True, which="both", ls="--")

            plt.tight_layout()
            plt.savefig(f"plots/bode_plot{i}.png")   # Guarda la gráfica como PDF
            
        self.generate_pdf()

    def generate_pdf(self):
        doc = SimpleDocTemplate("salida.pdf", pagesize=A4)
        story = []
        styles = getSampleStyleSheet()

        # Título
        story.append(Paragraph("Reporte Automático de Filtros", styles['Title']))
        story.append(Paragraph(datetime.today().strftime("%d/%m/%Y"), styles['Normal']))
        story.append(Spacer(1, 12))

        # Sección: Resumen
        resumen = (
            f"En este reporte se presentan los resultados obtenidos en el diseño y simulación "
            f"de un filtro <b>{'Chebyshev'}</b> activo de <b>{self.Sorder}</b> orden, "
            f"con una frecuencia de corte de <b>{self.fc} Hz</b> y un ripple de <b>{self.ripple} dB</b>."
        )
        story.append(Paragraph("Resumen", styles['Heading2']))
        story.append(Paragraph(resumen, styles['Normal']))
        story.append(Spacer(1, 12))

        # Iterar sobre etapas
        for i in sorted(self.coeficientes[self.ripple][self.Sorder]):
            etapa = int(i)
            data = self.response[etapa - 1][-1]

            story.append(Paragraph(f"Etapa {etapa}", styles['Heading2']))

            if data["orden"] == 1:  # Template 2
                op = "+" if data["config"] == 0 else "||"
                info = f"""
                    R1: <b>{data["R1"]}</b><br/>
                    C1: <b>{data["C1"]}</b><br/><br/>
                    Arreglo resistivo recomendado:<br/>
                    R1 = {data["r1"]} {op} {data["r2"]} {op} {data["r3"]}
                """
            else:  # Template 3
                op = "+" if data["orden"] == 0 else "||"
                info = f"""
                    R1: <b>{data["R1"]}</b><br/>
                    R2: <b>{data["R2"]}</b><br/>
                    C1: <b>{data["C1"]}</b><br/>
                    C2: <b>{data["C2"]}</b><br/><br/>
                    Arreglo capacitivo recomendado:<br/>
                    C2 = {data["c1"]} {op} {data["c2"]} {op} {data["c3"]}
                """

            story.append(Paragraph(info, styles['Normal']))
            story.append(Spacer(1, 12))

            # Insertar imagen bode_plot{i}.png
            try:
                img = Image(f"plots/bode_plot{etapa}.png", width=7*inch, height=5*inch)
                story.append(img)
                story.append(Paragraph(f"Figura: Respuesta en frecuencia de etapa {etapa}", styles['Italic']))
                story.append(Spacer(1, 24))
            except Exception as e:
                story.append(Paragraph(f"[Imagen no disponible: bode_plot{etapa}.png]", styles['Normal']))

        # Generar PDF
        doc.build(story)







