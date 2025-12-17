import customtkinter as ctk
from tkinter import ttk, messagebox

from sympy import Function, dsolve, Eq, symbols, Integral, lambdify, sin, cos, exp, log, sqrt
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations, implicit_multiplication_application, convert_xor)

import matplotlib.pyplot as plt

ctk.set_appearance_mode("dark")

# - - - VARIABLES GLOBALES - - -

valores_x=[]
valores_y=[]
y_solucion=[]
error_abs=[]
error_rel=[]

# - - - PALETA COLORES Y FUENTE - - -

Fuente=("Leelawadee UI", 14)
FuenteTTK=("Leelawadee UI", 10)
FuenteTitulo=("Leelawadee UI", 22, "bold")
FuenteBoton=("Leelawadee UI", 14, "bold")
FuenteBotonTTK=("Leelawadee UI", 10, "bold")

Blanco="#F0F0F0"
Azul="#121729"
AzulClaro="#1B2347"
AzulFuerte="#151B36"

# - - - VENTANA PRINCIPAL Y ESTILOS - - -

ventana=ctk.CTk()
ventana.title("RK4")
ventana.geometry("1100x650")
ventana.config(bg=Azul)

#--filas--
ventana.grid_rowconfigure(0, weight=0)
ventana.grid_rowconfigure(1, weight=0)
ventana.grid_rowconfigure(2, weight=0)
ventana.grid_rowconfigure(3, weight=1)

#--columnas--
ventana.grid_columnconfigure(0, weight=1)
ventana.grid_columnconfigure(1, weight=1)

#--estilo treeview--
style=ttk.Style(ventana)
style.theme_use("clam")
style.configure("Treeview", background=AzulFuerte, foreground=Blanco, rowheight=24, fieldbackground=AzulFuerte, font=FuenteTTK)
style.configure("Treeview.Heading", background=AzulClaro, foreground=Blanco, font=FuenteBotonTTK, relief="flat")
style.map("Treeview.Heading", background=[("active", Azul)])
style.map("Treeview", background=[("selected", Azul)], foreground=[("selected", Blanco)])

# - - - FUNCIONES - - -

x, y=symbols("x y")
g=Function("g")

def btnGraficar():
    try:
        if not valores_x or not valores_y:
            messagebox.showwarning("Sin datos", "Debe calcular los valores antes de graficar.")
        else:
            plt.figure("Gráfica del método RK4")
            plt.plot(valores_x, valores_y, marker='o', color='cyan', label='Aproximación')
            plt.plot(valores_x, y_solucion, marker='o', color='red', label='Valor real')
            plt.title("Comparación Método RK4")
            plt.xlabel("x")
            plt.ylabel("y")
            plt.grid(True, linestyle='--', alpha=0.6)
            plt.legend()
            plt.show()
    except Exception as e:
        messagebox.showerror("Error al graficar", f"Ocurrió un error: {e}")

def btnCalcular():
    try:
        funcion=entry_funcion.get()
        x0_aux, y0_aux=entry_PVI.get().split(",")
        x0=float(x0_aux)
        y0=float(y0_aux)
        h=float(entry_h.get())
        xf=float(entry_xf.get())

        if h<0:
            messagebox.showerror("Error", "El tamaño de paso debe ser mayor que 0.")
        if xf<x0:
            messagebox.showerror("Error", "El valor de x final debe ser mayor al valor de x₀.")

        calcular(x0, y0, xf, h, funcion)
        
    except ValueError:
        messagebox.showerror("Error", "Ingrese los campos.")
    except SyntaxError:
        messagebox.showerror("Error", "Ingrese una función valida.")
    except ZeroDivisionError:
        messagebox.showerror("Error", "La ED no es continua en algun punto del intervalo ingresado.")
      
def calcular(x0, y0, xf, h, funcion):
    for widget in frame_tabla.winfo_children():
        widget.destroy()

    #--limpiar arreglos--
    valores_x.clear()
    valores_y.clear()
    y_solucion.clear()
    error_abs.clear()
    error_rel.clear()

    func=convertirFuncion(funcion)
    n=int((xf-x0)/h)

    #--resolver la edo--
    try:
        edo=Eq(g(x).diff(x), func.subs(y, g(x)))
        sol=dsolve(edo, g(x), ics={g(x0): y0})
        expr=sol.rhs

        tiene_solucion=not expr.has(Integral)
    except:
        tiene_solucion=False
  

    #--función simbólica--
    f=lambdify((x, y), func, "numpy") 

    #--función solución--
    if tiene_solucion:
        s=lambdify(x, expr, "numpy")
    else:
        s=None
        messagebox.showinfo("Aviso", "La ecuación diferencial no tiene solución analítica.\n" "Se mostrará solo la aproximación numérica.")

    valores_x.append(x0)
    valores_y.append(y0)
    y_solucion.append(y0)
    
    #--calcular aproximacion--
    for i in range(n):
        x_actual=x0 + h*(i+1)
        valores_x.append(x_actual)

        k1, k2, k3, k4 = 0, 0, 0, 0
        k1=f(valores_x[i], valores_y[i])
        k2=f(valores_x[i]+0.5*h, valores_y[i]+0.5*h*k1)
        k3=f(valores_x[i]+0.5*h, valores_y[i]+0.5*h*k2)
        k4=f(valores_x[i]+h, valores_y[i]+h*k3)
        y_actual=valores_y[i] + ((1/6)*h)*(k1 + 2*k2 + 2*k3 + k4)
        valores_y.append(y_actual)

        if s is not None:
            y_solucion.append(s(x_actual))
    
    if tiene_solucion:
        calcularErrores(n)
    
    #--crear tabla--
    tabla=ttk.Treeview(frame_tabla, columns=("iteracion","x", "y", "valor_real", "error_absoluto","error_relativo"), show="headings", height=8)
    tabla.heading("iteracion", text="n")
    tabla.heading("x", text="xₙ")
    tabla.heading("y", text="yₙ")
    tabla.heading("valor_real", text="Valor real")
    tabla.heading("error_absoluto", text="Error absoluto")
    tabla.heading("error_relativo", text="Error relativo %")

    tabla.column("iteracion", width=80, anchor="center")
    tabla.column("x", width=80, anchor="center")
    tabla.column("y", width=80, anchor="center")
    tabla.column("valor_real", width=100, anchor="center")
    tabla.column("error_absoluto", width=100, anchor="center")
    tabla.column("error_relativo", width=100, anchor="center")
    
    #--insertar resultados en la tabla--
    for i in range(n+1):
        tabla.insert("", "end", values=(i, round(valores_x[i], 4), round(valores_y[i], 4), (round(y_solucion[i], 4) if s is not None else "---"), (round(error_abs[i], 4) if s is not None else "---"), (round(error_rel[i], 4) if s is not None else "---")))
    
    tabla.pack(fill="both", expand=True)

def convertirFuncion(txt):
    trnsf=standard_transformations+(implicit_multiplication_application, convert_xor)
    expr=parse_expr(txt, transformations=trnsf)
    return expr

def calcularErrores(n):
    for i in range(n+1):
        error_abs.append(abs(y_solucion[i]-valores_y[i]))
        if y_solucion[i] == 0:
            error_rel.append(0)
        else:
            error_rel.append(abs((y_solucion[i]-valores_y[i])/y_solucion[i]))

# - - - COMPONENTES - - -

#--frame titulo--
frame_encabezado=ctk.CTkFrame(ventana, fg_color=AzulFuerte, bg_color=Azul, corner_radius=10, border_width=2, border_color=AzulFuerte)
frame_encabezado.grid(row=0, column=0, columnspan=2, sticky="nsew", pady=(20, 10), padx=20)

#--titulo--
titulo=ctk.CTkLabel(frame_encabezado, text="Método numerico Runge-Kutta de cuarto orden", font=FuenteTitulo, text_color=Blanco)
titulo.pack(pady=10)
descripcion1=ctk.CTkLabel(frame_encabezado, text="El método de Runge-Kutta de cuarto orden es probablemente uno de los procedimientos númericos más populares, así como preciso, usado para obeter soluciones aproximadas para un problema con valores iniciales y' = f(x, y), y(x₀) = y₀.", font=Fuente, text_color=Blanco, wraplength=800)
descripcion1.pack(pady=(0, 10))

#--informacion--
frame_izq=ctk.CTkFrame(ventana, fg_color=AzulFuerte, bg_color=Azul, corner_radius=10, border_width=2, border_color=AzulFuerte)
frame_izq.grid(row=1, column=0, pady=10, padx=(20, 10), sticky="nsew")
descripcion2=ctk.CTkLabel(frame_izq, text="El conjunto de valores usado con más frecuencia para realizar el cálculo es el siguiente: ", font=Fuente, text_color=Blanco, wraplength=300, justify="left")
descripcion2.pack(pady=(5, 2), padx=15, anchor="w")

formula_k1=ctk.CTkLabel(frame_izq,text="k₁ = f(xₙ, yₙ)", text_color=Blanco, font=Fuente)
formula_k1.pack(pady=2, padx=15, anchor="w")
formula_k2=ctk.CTkLabel(frame_izq,text="k₂ = f(xₙ + ½h, yₙ + ½hk₁)", text_color=Blanco, font=Fuente)
formula_k2.pack(pady=2, padx=15, anchor="w")
formula_k3=ctk.CTkLabel(frame_izq, text="k₃ = f(xₙ + ½h, yₙ + ½hk₂)", text_color=Blanco, font=Fuente)
formula_k3.pack(pady=2, padx=15, anchor="w")
formula_k4=ctk.CTkLabel(frame_izq,text="k₄ = f(xₙ + h, yₙ + hk₃)", text_color=Blanco, font=Fuente)
formula_k4.pack(pady=(2, 5), padx=15, anchor="w")

#--ingresar datos--
frame_der=ctk.CTkFrame(ventana, fg_color=AzulFuerte, bg_color=Azul, corner_radius=10, border_width=2, border_color=AzulFuerte)
frame_der.grid(row=1, column=1, pady=10, padx=(10, 20), sticky="nsew")

etiqueta_funcion=ctk.CTkLabel(frame_der, text="Ingrese la función dy/dx:", text_color=Blanco, font=Fuente)
etiqueta_funcion.grid(row=0, column=0, pady=5,padx=15,sticky="w")
etiqueta_PVI=ctk.CTkLabel(frame_der, text="Ingrese el PVI:", text_color=Blanco, font=Fuente) 
etiqueta_PVI.grid(row=1, column=0, pady=5, padx=15, sticky="w")
etiqueta_h=ctk.CTkLabel(frame_der, text="Ingrese el tamaño de paso:", text_color=Blanco, font=Fuente)
etiqueta_h.grid(row=2, column=0, pady=5, padx=15, sticky="w")
etiqueta_xf=ctk.CTkLabel(frame_der,text="Ingrese el punto de aproximación:", text_color=Blanco, font=Fuente)
etiqueta_xf.grid(row=3, column=0, pady=5, padx=15, sticky="w")

entry_funcion=ctk.CTkEntry(frame_der, fg_color=Azul, font=Fuente, width=400, border_width=0)
entry_funcion.grid(row=0, column=1, padx=15, pady=(15, 5), sticky="ew")
entry_PVI=ctk.CTkEntry(frame_der, fg_color=Azul, font=Fuente, width=400, border_width=0)
entry_PVI.grid(row=1, column=1, padx=15, pady=5, sticky="ew")
entry_h=ctk.CTkEntry(frame_der, fg_color=Azul, font=Fuente, width=400, border_width=0)
entry_h.grid(row=2, column=1, padx=15, pady=5, sticky="ew")
entry_xf=ctk.CTkEntry(frame_der, fg_color=Azul, font=Fuente, width=400, border_width=0)
entry_xf.grid(row=3, column=1, padx=15, pady=(5, 15), sticky="ew")

#--botones--
frame_botones=ctk.CTkFrame(ventana, fg_color=AzulFuerte, bg_color=Azul, corner_radius=10, border_width=2, border_color=AzulFuerte)
frame_botones.grid(row=2, column=0, columnspan=2, pady=10, padx=20, sticky="nsew")

boton_calcular=ctk.CTkButton(frame_botones, width=50, text="Calcular", text_color=Blanco, font=FuenteBoton, fg_color=Azul, hover_color=AzulClaro, command=btnCalcular)
boton_calcular.grid(row=0, column=1, padx=15, pady=10, ipady=2, sticky="ew")
boton_graficar=ctk.CTkButton(frame_botones, width=50, text="Graficar", text_color=Blanco, font=FuenteBoton, fg_color=Azul, hover_color=AzulClaro, command=btnGraficar)
boton_graficar.grid(row=0, column=2, padx=15, pady=10, ipady=2, sticky="ew")

frame_botones.grid_columnconfigure(0, weight=1)
frame_botones.grid_columnconfigure(1, weight=1)
frame_botones.grid_columnconfigure(2, weight=1)
frame_botones.grid_columnconfigure(3, weight=1)

#--frame datos--
frame_tabla=ctk.CTkFrame(ventana, fg_color=AzulFuerte, bg_color=Azul, corner_radius=10, border_width=2, border_color=AzulFuerte)
frame_tabla.grid(row=3, column=0, columnspan=2, pady=(10, 20), padx=20, ipadx=10, ipady=10, sticky="nsew")

#--tooltip funcion--
class ToolTip:
    def __init__(self, widget, text):
        self.widget=widget
        self.text=text
        self.tooltip=None

        widget.bind("<Enter>", self.mostrar)
        widget.bind("<Leave>", self.ocultar)

    def mostrar(self, event=None):
        if self.tooltip:
            return

        x=self.widget.winfo_rootx()+20
        y=self.widget.winfo_rooty()+self.widget.winfo_height()+5

        self.tooltip=ctk.CTkToplevel(self.widget)
        self.tooltip.overrideredirect(True)
        self.tooltip.geometry(f"+{x}+{y}")

        label=ctk.CTkLabel(
            self.tooltip,
            text=self.text,
            justify="left",
            fg_color="#1B2347",
            text_color="#F0F0F0",
            corner_radius=8,
            padx=10,
            pady=8,
            font=("Leelawadee UI", 11)
        )

        label.pack()

    def ocultar(self, event=None):
        if self.tooltip:
            self.tooltip.destroy()
            self.tooltip=None
            
mensaje_funcion=(
    "Use:\n"
    "  sin(x), cos(x), exp(x), log(x), sqrt(x)\n\n"
    "Se permite multiplicación implícita:\n"
    "  2xy, 3x(x+1)\n\n"
    "No use:\n"
    "  sen(x), ln(x), xsin(x)"
)

ToolTip(entry_funcion, mensaje_funcion)
ToolTip(entry_PVI, "Ingrese el valor de y(x₀) separado por comas. Ej: 0,1")
ToolTip(entry_h, "Tamaño de paso (h > 0)")
ToolTip(entry_xf, "Valor final de x")

ventana.mainloop()