import numpy as np
import pandas as pd

class rvcurve:
    def __init__(self, time, rv, rv_err = None):

        self.time = time
        self.rv_in_time = rv
        self.rv_err_in_time = rv_err

        self.rv_in_phase = None
        self.rv_err_in_phase = None


    def CALC_PHASES(self, porb, t_zero):

        phase = ((self.time - t_zero) / porb) % 1
        phase[phase > 1] -= 1
        phase[phase < 0] += 1

        self.phase = phase

        sorted_indices = np.argsort(self.phase)

        self.phase = self.phase[sorted_indices]
        self.rv_in_phase = self.rv_in_time[sorted_indices]
        if self.rv_err_in_time is not None:
            self.rv_err_in_phase = self.rv_err_in_time[sorted_indices]




class lightcurve:
    def __init__(self, time, flux, flux_err):

        self.time = []
        self.flux_in_time = []
        self.flux_err_in_time = []

        for i in range(len(time)):
            if (not np.isnan(flux[i])):
                self.time.append(time[i])
                self.flux_in_time.append(flux[i])
                self.flux_err_in_time.append(flux_err[i])
        self.time = np.array(self.time)
        self.flux_in_time = np.array(self.flux_in_time)
        self.flux_err_in_time = np.array(self.flux_err_in_time)

        self.flux_in_phase = None
        self.flux_err_in_phase = None

        self.binned_phase = None
        self.binned_flux = None
        self.binned_flux_err = None
        self.flux_var = None

    def CALC_PHASES(self, porb, t_zero):

        phase = ((self.time - t_zero) / porb) % 1
        phase[phase > 1] -= 1
        phase[phase < 0] += 1

        self.phase = phase

        sorted_indices = np.argsort(self.phase)

        self.phase = self.phase[sorted_indices]
        self.flux_in_phase = self.flux_in_time[sorted_indices]
        self.flux_err_in_phase = self.flux_err_in_time[sorted_indices]


    def BIN_LC(self, nbins):
        bin_edges = np.linspace(self.phase.min(), self.phase.max(), nbins+1)

        bin_medians = []
        bin_flux_err = []

        for i in range(nbins):
            mask = (self.phase >= bin_edges[i]) & (self.phase < bin_edges[i+1])
            num_points = np.sum(mask)

            if num_points > 0:
                median_value = np.median(self.flux_in_phase[mask])
                err_value = np.sqrt(np.sum(self.flux_err_in_phase[mask]**2) / num_points**2)
            else:
                median_value = np.nan
                err_value = np.nan

            bin_medians.append(median_value)
            bin_flux_err.append(err_value)



        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])


        self.binned_phase = bin_centers
        self.binned_flux = bin_medians
        self.binned_flux_err = bin_flux_err

        s = pd.Series(self.flux_in_phase)
        rolling_variance = s.rolling(10).var()
        self.flux_var  = np.mean(rolling_variance)


    def BIN_LC_NON_UNIFORM_old(self, nbins, phase1, phase2, oversampling_ratio, width):


        # Compute Gaussian distribution around each phase
        gaussian1 = np.exp(-(self.phase - phase1)**2 / (2 * width**2))
        gaussian2 = np.exp(-(self.phase - phase2)**2 / (2 * width**2))

        # Combine the two Gaussian distributions
        combined_gaussian = gaussian1 + gaussian2

        # Normalize combined_gaussian such that it has a peak value of oversampling_ratio
        combined_gaussian = combined_gaussian / np.max(combined_gaussian) * oversampling_ratio

        # Create a weighting function
        weights = 1 + combined_gaussian

        # Now, distribute the bins based on the weights
        cum_weights = np.cumsum(weights)
        bin_edges = np.interp(np.linspace(0, cum_weights[-1], nbins + 1), cum_weights, self.phase)

        bin_medians = []
        bin_flux_err = []

        for i in range(nbins):
            mask = (self.phase >= bin_edges[i]) & (self.phase < bin_edges[i+1])
            num_points = np.sum(mask)

            if num_points > 0:
                median_value = np.median(self.flux_in_phase[mask])
                err_value = np.sqrt(np.sum(self.flux_err_in_phase[mask]**2) / num_points**2)
            else:
                median_value = np.nan
                err_value = np.nan

            bin_medians.append(median_value)
            bin_flux_err.append(err_value)

        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        self.binned_phase = bin_centers
        self.binned_flux = bin_medians
        self.binned_flux_err = bin_flux_err

        s = pd.Series(self.flux_in_phase)
        rolling_variance = s.rolling(10).var()
        self.flux_var  = np.mean(rolling_variance)

    def BIN_LC_NON_UNIFORM(self, nbins, phase1, phase2, oversampling_ratio, width):
        # Helper function to compute the minimum cyclic distance
        def cyclic_distance(phase, center):
            direct_dist = np.abs(phase - center)
            cyclic_dist = 1 - direct_dist  # Considering wrap-around at 0 and 1
            return np.minimum(direct_dist, cyclic_dist)

        # Compute Gaussian distribution around each phase, considering cyclicity
        gaussian1 = np.exp(-cyclic_distance(self.phase, phase1)**2 / (2 * width**2))
        gaussian2 = np.exp(-cyclic_distance(self.phase, phase2)**2 / (2 * width**2))

        # Combine the two Gaussian distributions
        combined_gaussian = gaussian1 + gaussian2

        # Normalize combined_gaussian such that it has a peak value of oversampling_ratio
        combined_gaussian = combined_gaussian / np.max(combined_gaussian) * oversampling_ratio

        # Create a weighting function
        weights = 1 + combined_gaussian

        # Now, distribute the bins based on the weights
        cum_weights = np.cumsum(weights)
        bin_edges = np.interp(np.linspace(0, cum_weights[-1], nbins + 1), cum_weights, np.sort(self.phase))

        # Ensure bin edges wrap around correctly
        bin_edges = np.mod(bin_edges, 1)  # Wrap-around using modulo to ensure edges are within [0, 1]

        bin_medians = []
        bin_flux_err = []

        for i in range(nbins):
            # Use modulo to handle cyclic nature of phase space for mask
            mask = ((self.phase >= bin_edges[i]) & (self.phase < bin_edges[i+1])) | \
                   ((self.phase >= bin_edges[i] - 1) & (self.phase < bin_edges[i+1] - 1)) | \
                   ((self.phase >= bin_edges[i] + 1) & (self.phase < bin_edges[i+1] + 1))

            num_points = np.sum(mask)

            if num_points > 0:
                median_value = np.median(self.flux_in_phase[mask])
                err_value = np.sqrt(np.sum(self.flux_err_in_phase[mask]**2) / num_points**2)
            else:
                median_value = np.nan
                err_value = np.nan

            bin_medians.append(median_value)
            bin_flux_err.append(err_value)

        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        bin_centers = np.mod(bin_centers, 1)  # Ensure bin centers are within [0, 1]

        self.binned_phase = bin_centers
        self.binned_flux = bin_medians
        self.binned_flux_err = bin_flux_err

        s = pd.Series(self.flux_in_phase)
        rolling_variance = s.rolling(10).var()
        self.flux_var  = np.mean(rolling_variance)
