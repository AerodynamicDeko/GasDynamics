import math
import numpy as np
import pandas as pd
import altair as alt
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def collect_data_from_user() -> dict: #essa seta serve para dizer qual o tipo de retorno da função
    """
    Function that collect data from users to perform calculation
    """
    print('Valores de saída de uma tubeira padrão')
    return {
        "gama": float(input("Insira o valor de Gama: ")),
        "R": float(input("Insira o valor de R: ")),
        "p": float(input("Insira o valor de p (Pa) antes da tubeira: ")),
        "T": float(input("Insira a temperatura (K) inicial do escoamento: ")),
        "rt": float(input("Insira a densidade antes da onda de choque (kg/m³): ")),
        "raio_garganta": float(input("Insira o raio da garganta (m): ")),
        "raio_saida": float(input("Insira o raio da saída (m): "))
    }

# Função para calcular o Mach em função da relação de área A/A*
def mach_from_area_ratio(area_ratio, gama) -> tuple:
    """
    Essa função, mach_from_area_ratio, calcula o número de Mach (M) correspondente 
    a uma relação de área dada (area_ratio = A/A*) e um valor específico de gama. 
    Ela usa a relação de área isentrópica para o escoamento compressível, e o método fsolve para 
    encontrar numericamente o valor de M que faz com que a equação seja satisfeita.
    Objetivo: Encontrar o número de Mach para uma dada relação de área e gama.
    Estratégia: Usar fsolve para resolver a equação transformada para f(M) = 0.
    Resultado: Retornar o número de Mach que corresponde à relação de área fornecida.
    """
    def equation(M):
        return (1 / M) * ((2 / (gama + 1)) * (1 + (gama - 1) / 2 * M ** 2)) ** ((gama + 1) / (2 * (gama - 1))) - area_ratio
    M_initial_guess = 1.3
    ''' O fsolve usa o método de Newton-Raphson para encontrar o valor de M a partir de um "chute inicial"'''
    return fsolve(equation, M_initial_guess)


def calculate_max_area_ratio() -> float:
    # Calcular áreas da garganta e saída
    # Área da garganta
    A_star = math.pi * data_from_users["raio_garganta"] ** 2
    # Área da saída
    A_saida = math.pi * data_from_users["raio_saida"] ** 2

    # Máxima relação A/A* baseada nos raios
    max_area_ratio = A_saida / A_star
    return max_area_ratio


def calculate_total_conditions(T_from_user, p_from_user, gama_) -> tuple:
    """

    """
    # Condições totais
    Tt = T_from_user * (1 + ((gama_ - 1) / 2))
    pt = p_from_user * (1 + (gama_ - 1) / 2) ** (gama_ / (gama_ - 1))
    return Tt, pt


def calculate_sound_speed_and_entry_mach(R_, T_, gama_):
    """

    """
    a = math.sqrt(gama_ * R_ * T_)
    v = float(input('Insira a velocidade inicial "entrada na tubeira": '))
    return v / a


def calculate_isentropic_relation_and_after_shock_wave(p_, gama_, M1_, T_, rt_, Tt_) -> None:
    P1 = p_
    P2 = P1 * ((gama_ + 1) / (2 * M1_ ** 2)) / (1 + (gama_ - 1) / (2 * M1_ ** 2))
    T2 = T_ * (1 + (gama_ - 1) / (2 * M1 ** 2)) / (1 + (gama_ - 1) / 2)
    ro2 = rt_ * (P2 / P1) ** (1 / gama_)
    v2 = math.sqrt(2 * R * (Tt_ - T2))

    print(f'Pressão após a onda de choque: {P2:.2f} Pa')
    print(f'Temperatura após a onda de choque: {T2:.2f} K')
    print(f'Densidade após a onda de choque: {ro2:.2f} kg/m³')
    print(f'Velocidade na saída: {v2:.2f} m/s')


def calculate_divergent_section(area_ratios_, gama_,  Tt_, pt_, R_) -> tuple:
    """

    """
    # Seção divergente com cálculo de Mach em função de A/A*
    mach_numbers = []
    pressures = []
    temperatures = []
    densities = []
    velocities = []

    for area_ratio in area_ratios_:
        M = mach_from_area_ratio(area_ratio, gama_)
        T = Tt_ / (1 + (gama_ - 1) / 2 * M ** 2)
        p = pt_ / (1 + (gama_ - 1) / 2 * M ** 2) ** (gama_ / (gama_ - 1))
        rho = p / (R_ * T)
        v = M * math.sqrt(gama_ * R_ * T)
        # Aqui vou usar a função append para adicionar os valores às listas que eu abri na linha 83
        mach_numbers.append(M)
        pressures.append(p)
        temperatures.append(T)
        densities.append(rho)
        velocities.append(v)

    return mach_numbers, pressures, temperatures, densities, velocities


def plot_data(area_ratios_, mach_numbers_, pressures_, temperatures_, densities_, velocities_) -> None:
    """

    """
    # Cria DataFrame para armazenar resultados na seção divergente
    df_divergente = pd.DataFrame({
        "Área Relativa (A/A*)": area_ratios_,
        "Mach": mach_numbers_,
        "Pressão (Pa)": pressures_,
        "Temperatura (K)": temperatures_,
        "Densidade (kg/m³)": densities_,
        "Velocidade (m/s)": velocities_
    })
    print("\nParâmetros na seção divergente:")
    print(df_divergente)

    # Plotando os resultados na seção divergente com ajuste automático de escala
    fig, axs = plt.subplots(3, 2, figsize=(10, 12))
    fig.suptitle("Parâmetros do Escoamento na Seção Divergente com Relação de Área")

    # Gráficos individuais para cada parâmetro com ajuste de limites
    parameters = ["Mach", "Pressão (Pa)", "Temperatura (K)", "Densidade (kg/m³)", "Velocidade (m/s)"]
    for i, param in enumerate(parameters):
        ax = axs[i // 2, i % 2]
        ax.plot(df_divergente["Área Relativa (A/A*)"], df_divergente[param], label=param)
        ax.set_xlabel("Área Relativa (A/A*)")
        ax.set_ylabel(param)
        ax.set_ylim(df_divergente[param].min() * 0.95, df_divergente[param].max() * 1.05)  # Ajuste para dar espaço nas extremidades
        ax.legend()

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    # plt.show()
    plt.savefig("plot.png")


if __name__ == '__main__':
    """
    [10:19, 12/11/2024] André Luiz: Infos para fluxo isentrópico
        Gama: 1.4
        R: 287
        Pressão na entrada: 101325
        Temperatura na entrada: 288.15
        Densidade na entrada: 1.225
        Raio da garganta: 0.05
        Raio da saída: 0.1
        Velocidade inicial na entrada: 300
    
    [10:19, 12/11/2024] André Luiz: para onda de choque:
        Gama: 1.66 helio
        R: 287
        (Pa) antes da tubeira: 101325
        (K) inicial do escoamento: 300
        densidade (kg/m³): 1.225
        Insira o raio da garganta (m): 0.03
        Insira o raio da saída (m): 0.95
        "entrada na tubeira": 300 m/s
    """

    #TODO - Padronizar o programa (ou tudo em ingles, ou em portugues - o default é inglês)
    #TODO - Padronizar os nomes de variáveis: variáveis declaradas com nomes de letras nao é uma boa prática
    #TODO - Documentar as funções e o tipo de dado de cada parâmetro na declaração da função

    data_from_users = collect_data_from_user()
    gama = data_from_users['gama']
    R = data_from_users['R']

    # Gerar uma série de valores de A/A* para a seção divergente, do A0 até o valor máximo
    # A função linspace cria 50 pontos entre as áreas de entrada e saída da tubeira
    area_ratios = list(np.linspace(1.0, calculate_max_area_ratio(), num=50))

    Tt, pt = calculate_total_conditions(data_from_users["T"], data_from_users["p"], gama)

    # Calculo da velocidade do som e Mach na entrada
    M1 = calculate_sound_speed_and_entry_mach(R_=R, T_=data_from_users['T'], gama_=gama)

    # Relações isentrópicas e propriedades após onda de choque
    calculate_isentropic_relation_and_after_shock_wave(p_=data_from_users["p"],
                                                       gama_=gama,
                                                       M1_=M1,
                                                       T_=data_from_users["T"],
                                                       rt_=data_from_users["rt"],
                                                       Tt_=Tt)
    
    # Seção divergente com cálculo de Mach em função de A/A*
    mach_numbers, pressures, temperatures, densities, velocities = calculate_divergent_section(area_ratios_=area_ratios,
                                                                                               gama_=gama,
                                                                                               Tt_=Tt,
                                                                                               pt_=pt,
                                                                                               R_=data_from_users["R"])
    
    plot_data(area_ratios_=area_ratios,
              mach_numbers_=mach_numbers,
              pressures_=pressures,
              temperatures_=temperatures,
              densities_=densities,
              velocities_=velocities)