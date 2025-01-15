'''
Problem 1
Consider a one-dimensional random walker that can move every second. With probability pl = 1/3 it moves
to the left, with probability pr = 1/3 it moves to right, and with probability ps = 1/3 it rests/stays and does
not move. Assuming at time t = 0, the random walker is at x = 0, plot the probability density function and
the cumulative probability function for t = 10, t = 100, and t = 1000 seconds. Make just two plots; each
showing all three time points.  Do you understand why these plots look different? The plots that you make 
should be designed well. For example, they should label curves, axes, etc.

Steps:
    1. create a random walk: variables will need to be adjustable
        a. pl = 1/3, pr = 1/3, ps = 1/3
    2. run simulations to get accurate statistics
    3. calculate probability density function
        a. @ t = 10
        b. @ t = 100
        c. @ t = 1000
        d. create plot showing all three time points
    4. calculate cumulative distribution function
        a. @ t = 10
        b. @ t = 100
        c. @ t = 1000
        d. create plot showing all three time points
    5. creating plots for the other pl
        a. pl = 0, pr = 1/2, ps = 1/2
'''

import numpy as np
import matplotlib.pyplot as plt

def simulate_random_walk(pl, pr, ps, max_time, num_simulations):
    #inputs: 
        #pl: probability of turning left
        #pr: prob of turning right
        #ps: prob of staying straight
        #max_time = highest timepoint desired
        #num_simulations = number of simulations desired
    # only three possible movements: -1 (left), 0 (stay), 1 (right); defining it this way simplifies plotting
    moves = [-1, 0, 1]
    #storing probabilities in a list
    probs = [pl, pr, ps]
    
    # store results of simulations at each of these timepoints
    timepoints = {10: [], 100: [], 1000: []}
    
    for _ in range(num_simulations): #for the number of simulations desired 
        position = 0  # start at 0, reset after each run
        for t in range(1, max_time + 1): # run up to the timepoint
            step = np.random.choice(moves, p=probs)
            position += step 
            if t in timepoints:
                timepoints[t].append(position)
    
    return timepoints


def plot_distributions(positions, time_points, title, filename_prefix=None):
    # specifics for plotting
    fig, ax = plt.subplots(1, 2, figsize=(14, 6))
    
    for t in time_points:
        # calculating the pdf
        counts, bins = np.histogram(positions[t], bins=range(min(positions[t]), max(positions[t])+1), density=True)
        centers = (bins[:-1] + bins[1:]) / 2
        ax[0].plot(centers, counts, label=f't = {t} seconds')
        
        # calculating the cdf
        cdf = np.cumsum(counts) * (bins[1] - bins[0])
        ax[1].plot(centers, cdf, label=f't = {t} seconds')
    
    # formatting the pdf plot with specifics
    ax[0].set_xlabel('Position')
    ax[0].set_ylabel('Probability Density')
    ax[0].set_title(f'{title} - Probability Density')
    ax[0].legend()
    
    # formatting the cdf plot with specifics
    ax[1].set_xlabel('Position')
    ax[1].set_ylabel('Cumulative Probability')
    ax[1].set_title(f'{title} - Cumulative Distribution')
    ax[1].legend()

    # saving the plots if name provided
    if filename_prefix:
        plt.savefig(f"{filename_prefix}.jpg", format="jpeg")
    
    plt.tight_layout() #auto-adjusts to prevent overlap
    plt.show()


# parameters
time_points = [10, 100, 1000]
num_simulations = 10000  # the higher the number, the tigther the curves; 1,000 worked well for me

# pl=1/3, pr=1/3, ps=1/3
positions1 = simulate_random_walk(pl=1/3, pr=1/3, ps=1/3, max_time=1000, num_simulations=num_simulations)
plot_distributions(positions1, time_points, title='Equal Probabilities (p_l=1/3, p_r=1/3, p_s=1/3)', filename_prefix="biased_probabilities_plot")

# pl=0, pr=1/2, ps=1/2
positions2 = simulate_random_walk(pl=0, pr=1/2, ps=1/2, max_time=1000, num_simulations=num_simulations)
plot_distributions(positions1, time_points, title='Equal Probabilities (p_l=0, p_r=1/2, p_s=1/2)', filename_prefix="biased_probabilities_plot")
