import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, FFMpegWriter
import pandas as pd
import time

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

def read_animate_save(h, forward_factor, algorithm, max_size=None, fps=60):
    # algorithm = euler, verlet
    # max_size of the file in mB
    # forward_factor = slow factor

    #read in data
    df1 = pd.read_csv(f"build/{algorithm}_r1_h{h}.csv")
    df2 = pd.read_csv(f"build/{algorithm}_r2_h{h}.csv")
    
    # reduce data
    if max_size is not None:
        factor = 8000/20
        while(True):
    
            if (df1.shape[0] > max_size*factor):
                df1 = df1.iloc[::2]
                
            elif df2.shape[0] > max_size*factor:
                df2 = df2.iloc[::2]
            else:
              break
        
    #coordinates
    time = df1.iloc[:,0].to_numpy()

    x1 = df1.iloc[:,1].to_numpy()
    y1 = df1.iloc[:,2].to_numpy()
    x2 = df2.iloc[:,1].to_numpy()
    y2 = df2.iloc[:,2].to_numpy()


    # get borders
    x_min = np.min(np.minimum.reduce([x1, x2]))
    x_max = np.max(np.maximum.reduce([x1, x2]))
    y_min = np.min(np.minimum.reduce([y1, y2]))
    y_max = np.max(np.maximum.reduce([y1, y2]))

    # create figure
    fig, ax = plt.subplots(dpi=300)
    plt.axis('off')

    line1, = ax.plot([], [], 'ro', markersize=20)
    line1_tr, = ax.plot([], [], 'r--', label='mass 1')

    line2, = ax.plot([], [], 'bo', markersize=40)
    line2_tr, = ax.plot([], [], 'b--', label='mass 2')

    # set x and y limits
    ax.set_xlim(x_min-0.5, x_max+0.5)
    ax.set_ylim(y_min-0.5, y_max+0.5)

    # update function for the animation
    def update(frame):
        line1_tr.set_data(x1[:frame], y1[:frame])
        line2_tr.set_data(x2[:frame], y2[:frame])
        line1.set_data(x1[frame], y1[frame])
        line2.set_data(x2[frame], y2[frame])

        return line1, line2, line1_tr, line2_tr

    # create the animation
    anim = FuncAnimation(fig, update, frames=len(time), interval=forward_factor, repeat=False)

    # FFMpegWriter to save the video
    writer = FFMpegWriter(fps=fps, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(f'videos/build/trajectory_{algorithm}_h{h}.mp4', writer=writer)

    # display animation
    plt.show()


if __name__ == "__main__":
    path_config = "config.txt"
    h, algorithm = get_params(path_config)

    time_start = time.time()
    read_animate_save(h=h, forward_factor=30, fps=60, algorithm=algorithm, max_size=3)
    time_end = time.time()
    print(f"Animation duration: {time_end-time_start:.1f}s")