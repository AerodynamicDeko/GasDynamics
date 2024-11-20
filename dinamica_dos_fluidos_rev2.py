import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
import gradio as gr

def mach_from_area_ratio(area_ratio, gama) -> float:
    def equation(M):
        return (1 / M) * ((2 / (gama + 1)) * (1 + (gama - 1) / 2 * M ** 2)) ** ((gama + 1) / (2 * (gama - 1))) - area_ratio
    M_initial_guess = 1.3
    M_solution, = fsolve(equation, M_initial_guess)
    return M_solution

def calculate_and_plot(gama, R, p, T, rt, raio_garganta, raio_saida, v):
    # Área da garganta e da saída
    A_star = math.pi * raio_garganta ** 2  
    A_saida = math.pi * raio_saida ** 2    

    max_area_ratio = A_saida / A_star  
    area_ratios = list(np.linspace(1.0, max_area_ratio, num=50)) 

    # Condições totais
    Tt = T * (1 + ((gama - 1) / 2))
    pt = p * (1 + (gama - 1) / 2) ** (gama / (gama - 1))

    # Calculo da velocidade do som e Mach na entrada
    a = math.sqrt(gama * R * T)
    M1 = v / a

    # Relações após onda de choque
    P1 = p
    P2 = P1 * ((gama + 1) / (2 * M1 ** 2)) / (1 + (gama - 1) / (2 * M1 ** 2))
    T2 = T * (1 + (gama - 1) / (2 * M1 ** 2)) / (1 + (gama - 1) / 2)
    ro2 = rt * (P2 / P1) ** (1 / gama)
    v2 = math.sqrt(2 * R * (Tt - T2))

    # Cálculo na seção divergente
    mach_numbers = []
    pressures = []
    temperatures = []
    densities = []
    velocities = []

    for area_ratio in area_ratios:
        M = mach_from_area_ratio(area_ratio, gama)
        T = Tt / (1 + (gama - 1) / 2 * M ** 2)
        p = pt / (1 + (gama - 1) / 2 * M ** 2) ** (gama / (gama - 1))
        rho = p / (R * T)
        v = M * math.sqrt(gama * R * T)

        mach_numbers.append(M)
        pressures.append(p)
        temperatures.append(T)
        densities.append(rho)
        velocities.append(v)

    # Cria DataFrame para armazenar resultados na seção divergente
    df_divergente = pd.DataFrame({
        "Área Relativa (A/A*)": area_ratios,
        "Mach": mach_numbers,
        "Pressão (Pa)": pressures,
        "Temperatura (K)": temperatures,
        "Densidade (kg/m³)": densities,
        "Velocidade (m/s)": velocities
    })

    # Plotando os resultados na seção divergente
    fig, axs = plt.subplots(3, 2, figsize=(10, 12))
    fig.suptitle("Parâmetros do Escoamento na Seção Divergente com Relação de Área")

    parameters = ["Mach", "Pressão (Pa)", "Temperatura (K)", "Densidade (kg/m³)", "Velocidade (m/s)"]
    for i, param in enumerate(parameters):
        ax = axs[i // 2, i % 2]
        ax.plot(df_divergente["Área Relativa (A/A*)"], df_divergente[param], label=param)
        ax.set_xlabel("Área Relativa (A/A*)")
        ax.set_ylabel(param)
        ax.set_ylim(df_divergente[param].min() * 0.95, df_divergente[param].max() * 1.05)
        ax.legend()

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    return fig, P2, T2, ro2, v2

# Definir a interface do Gradio com as novas sintaxes de inputs e outputs
interface = gr.Interface(
    fn=calculate_and_plot,
    inputs=[
        gr.Number(label="Gama"),
        gr.Number(label="R (Constante do gás)"),
        gr.Number(label="Pressão inicial (Pa)"),
        gr.Number(label="Temperatura inicial (K)"),
        gr.Number(label="Densidade inicial (kg/m³)"),
        gr.Number(label="Raio da garganta (m)"),
        gr.Number(label="Raio da saída (m)"),
        gr.Number(label="Velocidade inicial (m/s)")
    ],
    outputs=[
        gr.Plot(label="Gráficos dos Parâmetros na Seção Divergente"),
        gr.Textbox(label="Pressão após a onda de choque (Pa)"),
        gr.Textbox(label="Temperatura após a onda de choque (K)"),
        gr.Textbox(label="Densidade após a onda de choque (kg/m³)"),
        gr.Textbox(label="Velocidade na saída (m/s)")
    ],
    title="Análise de Escoamento em Tubeira",
    description="Digite os valores dos parâmetros iniciais para calcular os parâmetros na seção divergente de uma tubeira padrão."
)

interface.launch()(share=True)