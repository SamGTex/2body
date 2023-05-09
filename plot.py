import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def create_plot(h, df_r1, df_r2, path):
    # plot results of euler algorithm

    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(df_r1.iloc[:, 1].to_numpy(), df_r1.iloc[:, 2].to_numpy(), 'b.', label="mass 1")
    plt.plot(df_r2.iloc[:, 1].to_numpy(), df_r2.iloc[:, 2].to_numpy(), 'r.', label="mass 2")

    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.savefig(path)
    plt.clf()

def get_params(path):
    #get parameters from config.txt
    with open(path, "r") as f:
        for line in f:
            if "h=" in line:
                h = line.split("=")[1].strip()
            elif "algorithm=" in line:
                algorithm_int = int(line.split("=")[1].strip())
    
    # fix h: _ instad of . and 6 digits after comma
    h = h.replace(".", "_")
    h_floats = h.split("_")[1]
    if len(h_floats) < 6:
        h = h + (6-len(h_floats))*"0"

    # fix algorithm
    if algorithm_int==0:
        algorithm = "euler"
    elif algorithm_int==1: 
        algorithm = "verlet"
    else:
        print("algorithm not defined")
    
    return h, algorithm


if __name__ == "__main__":
    path_config = "config.txt"
    h, algorithm = get_params(path_config)

    path_write = f"build/A2_{algorithm}_h{h}.pdf"
    df_r1 = pd.read_csv(f"build/{algorithm}_r1_h{h}.csv")
    df_r2 = pd.read_csv(f"build/{algorithm}_r2_h{h}.csv")

    create_plot(h, df_r1, df_r2, path_write) #, df_verlet_r1, df_verlet_r2)