import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


def calculate_order_parameter(num_cells, dat, dens, step, rng_seed):
    # para cada paso temporal calculamos el valor de los parametros de orden
    nematic_order = []
    polar_order = []
    nematic_2_order = []
    polar_2_order = []
    #last_tic = np.array(range(len(rng_seed)))
    for tic in range(0, dat, step): # Cambiar 500 a 0!
        # para cada paso. leemos el archivo
        # para cada paso de cada seed
        nematic_seed = np.array([])
        polar_seed = np.array([])
        nematic_2_seed = np.array([])
        polar_2_seed = np.array([])
        #i = 0
        for seed in rng_seed:
            # leemos el archivo correspondiente. Si no existe, entonces q lea el último existente
            # if os.path.exists(
            #     f"ejemplito_ovito/culture_nc={num_cells}_l={side}_rng_seed={seed}_step={tic:05}.dat"
            # ):
            #dat_actual = f"culture_nc={num_cells}_rho={dens}_rng_seed={seed}_step={tic:05}.dat"
            dat_actual = f"dat/{dens:.2f}/culture_initial_number_of_cells={num_cells}_density={dens}_force=Anisotropic_Grosmann_k=3.33_gamma=3_With_Noise_eta=0.033_With_Shrinking_rng_seed={seed}_step={tic:05}.dat"
            dat_5000 = f"dat/{dens:.2f}/culture_initial_number_of_cells={num_cells}_density={dens}_force=Anisotropic_Grosmann_k=3.33_gamma=3_With_Noise_eta=0.033_With_Shrinking_rng_seed={seed}_step=05000.dat"
            dat_2000 = f"dat/{dens:.2f}/culture_initial_number_of_cells={num_cells}_density={dens}_force=Anisotropic_Grosmann_k=3.33_gamma=3_With_Noise_eta=0.033_With_Shrinking_rng_seed={seed}_step=02000.dat"
            if os.path.exists(
                dat_actual
            ):
                df_tic = pd.read_csv(
                    dat_actual
                )
            elif os.path.exists(
                dat_5000
            ):
                df_tic = pd.read_csv(
                    dat_5000
                )
            elif os.path.exists(
                dat_2000
            ):
                df_tic = pd.read_csv(
                    dat_2000
                )
            else:
                print("Error")
                #last_tic[i] = tic
            # else:
            #     last_tic_i = last_tic[i]
            #     df_tic = pd.read_csv(
            #         f"ejemplito_ovito/culture_nc={num_cells}_l={side}_rng_seed={seed}_step={last_tic_i:05}.dat"
            #     )
            #i = i + 1
            # debemos calcular la cantidad de células elongadas
            # para eso creamos un nuevo dataframe df_elongated
            df_elongated = df_tic[df_tic["aspect_ratio"] != 1]
            num_elongated = len(df_elongated)
            # calculamos el sen(phi), cos(phi), sen(2phi), cos(2phi) para cada celula
            sin = np.sin(df_elongated["orientation"])
            cos = np.cos(df_elongated["orientation"])
            sin_2 = np.sin(2 * df_elongated["orientation"])
            cos_2 = np.cos(2 * df_elongated["orientation"])
            # y luego los sumamos
            sum_sin = sin.sum()
            sum_cos = cos.sum()
            sum_sin_2 = sin_2.sum()
            sum_cos_2 = cos_2.sum()
            # calculamos el parametro
            if num_elongated != 0:
                nematic = np.sqrt(sum_sin_2**2 + sum_cos_2**2) / num_elongated
                polar = np.sqrt(sum_sin**2 + sum_cos**2) / num_elongated
                nematic_2 = np.sqrt(sum_sin_2**2 + sum_cos_2**2) / num_cells
                polar_2 = np.sqrt(sum_sin**2 + sum_cos**2) / num_cells
            else:
                nematic = 0
                polar = 0
                nematic_2 = 0
                polar_2 = 0
            # y lo agregamos a la lista del seed
            nematic_seed = np.append(nematic_seed, nematic)
            polar_seed = np.append(polar_seed, polar)
            nematic_2_seed = np.append(nematic_2_seed, nematic_2)
            polar_2_seed = np.append(polar_2_seed, polar_2)
        # ahora calculamos el promedio de los parámetros para todas las seed
        nematic_mean = np.mean(nematic_seed)
        polar_mean = np.mean(polar_seed)
        nematic_2_mean = np.mean(nematic_2_seed)
        polar_2_mean = np.mean(polar_2_seed)

        nematic_order.append(nematic_mean)
        polar_order.append(polar_mean)
        nematic_2_order.append(nematic_2_mean)
        polar_2_order.append(polar_2_mean)
    return nematic_order, polar_order, nematic_2_order, polar_2_order


density = [0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90]
nc = 400
cell_area = np.pi
data = 30000 ###################
step = 100 #####################
delta_t = 0.05 ###################
# todos los seed
number_of_realizations=50

seed_1 = 0x87351080E25CB0FAD77A44A3BE03B491
rng_1 = np.random.default_rng(seed_1)

rng_seed_1 = rng_1.integers(
            low=2**20, high=2**50, size=number_of_realizations
        )

seed_2 = 1
rng_2 = np.random.default_rng(seed_2)

rng_seed_2 = rng_2.integers(
            low=2**20, high=2**50, size=number_of_realizations
        )

rng_seed = np.concatenate((rng_seed_1, rng_seed_2))

for dens in density:
    df = pd.DataFrame(
        {"tiempo[h]": np.array(list(range(0, data+step, step))) * delta_t} # cambiar data a data+1!
    )
    nematic_order, polar_order, nematic_2_order, polar_2_order = calculate_order_parameter(
        nc, data + 1, dens, step, rng_seed
    )
    df.insert(1, f"nematic_order_nc={nc}_rho={dens}", nematic_order)
    df.insert(2, f"polar_order_nc={nc}_rho={dens}", polar_order)
    df.insert(3, f"nematic_order_2_nc={nc}_rho={dens}", nematic_2_order)
    df.insert(4, f"polar_order_2_nc={nc}_rho={dens}", polar_2_order)


    # guardamos el dataframe
    #df.to_csv("order_parameters.csv", index=False)

    plt.figure()
    # hacemos el grafico con respecto a las horas
    plt.plot(
        df["tiempo[h]"],
        df[f"nematic_order_nc={nc}_rho={dens}"],
        label=f"Nematic order (Q)",
        color="blue",
        linewidth=1,
    )
    plt.plot(
        df["tiempo[h]"],
        df[f"polar_order_nc={nc}_rho={dens}"],
        color="green",
        label=f"Polar order (P)",
        linewidth=1,
    )
    plt.xlabel("System time")
    plt.ylabel("Order parameters")
    #plt.title(f"Evolution of the order parameters for a system with \u03C1={dens:.2f}") 
    plt.legend()
    plt.savefig(f"graphs/order parameters/order_parameters_rho={dens:.2f}.jpg", dpi=600)
    plt.close()

    plt.figure()
    # hacemos el grafico de los otros parametros de orden con respecto a las horas
    plt.plot(
        df["tiempo[h]"],
        df[f"nematic_order_2_nc={nc}_rho={dens}"],
        #label=f"Nematic order (Q)",
        label=r"Nematic order ($\hat{Q})$",
        color="blue",
        linewidth=1,
    )
    plt.plot(
        df["tiempo[h]"],
        df[f"polar_order_2_nc={nc}_rho={dens}"],
        #label=f"Polar order (P)",
        label=r"Polar order ($\hat{P})$",
        color="green",
        linewidth=1,
    )
    plt.xlabel("System time")
    plt.ylabel("Order parameters")
    #plt.title(f"Evolution of the order parameters for a system with \u03C1={dens:.2f}") 
    plt.legend()
    plt.savefig(f"graphs/order parameters/order_parameters_rho={dens:.2f}_2.jpg", dpi=600)
    plt.close()
#plt.show()
