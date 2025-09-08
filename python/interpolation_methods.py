import pandas as pd
import numpy as np
import scipy
import copy
import sklearn
from sklearn import linear_model
import seaborn as sns
import matplotlib.pyplot as plt
import datetime


def knn_interpolation(samples, time_point_to_complete, K):
    """
    gets a list of all samples (as a dataframe. the columns names shoul be the timepoints), a time point to interpolate and K

    returns an interpolation of the data point using KNN kernel with Epanechnikov function
    """
    all_time_points = (list(samples.columns))
    all_time_points.pop(0)
    samples_dist = [float(time_point) - float(time_point_to_complete) for time_point in all_time_points]
    samples_dist = np.abs(samples_dist)
    samples_dist.sort()
    b = samples_dist[K]

    time_point_interpolated = None

    for i in range(samples.shape[1]):
        sample = float(samples.columns[i])
        if sample - b <= time_point_to_complete and time_point_to_complete <= sample + b:
            ker_norm = (time_point_to_complete - sample) / b
            kernel = 0.75 * (1 - (ker_norm ** 2))
            if type(time_point_interpolated) == np.ndarray:
                time_point_interpolated += kernel * np.array(samples.iloc[:, i])
            else:
                time_point_interpolated = kernel * np.array(samples.iloc[:, i])

    sum_time_point = sum(time_point_interpolated)
    #for i in range(time_point_interpolated.shape[0]):
        #time_point_interpolated[i] = time_point_interpolated[i] / sum_time_point
    return time_point_interpolated