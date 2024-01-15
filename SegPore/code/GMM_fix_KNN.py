# -*- coding: UTF-8 -*-
"""
@Project : GMM 
@File    : myGMM.py
@Date    : 2023/2/21 14:12 
@Github  : https://github.com/guangzhaocs/SegPore
"""
import numbers
import numpy as np
from time import time
from scipy import linalg
from sklearn.mixture._gaussian_mixture import _compute_precision_cholesky, _estimate_log_gaussian_prob
from sklearn import cluster
from scipy.special._logsumexp import logsumexp


def check_random_state(seed):
    """Turn seed into a np.random.RandomState instance.

    Parameters
    ----------
    seed : None, int or instance of RandomState
        If seed is None, return the RandomState singleton used by np.random.
        If seed is an int, return a new RandomState instance seeded with seed.
        If seed is already a RandomState instance, return it.
        Otherwise raise ValueError.

    Returns
    -------
    :class:`numpy:numpy.random.RandomState`
        The random state object based on `seed` parameter.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, numbers.Integral):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError(
        "%r cannot be used to seed a numpy.random.RandomState instance" % seed
    )


def _estimate_gaussian_parameters(X, resp, reg_covar, covariance_type, fix_mu_flag, fix_mu, max_sigma):
    """Estimate the Gaussian distribution parameters.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        The input data array.

    resp : array-like of shape (n_samples, n_components)
        The responsibilities for each data sample in X.

    reg_covar : float
        The regularization added to the diagonal of the covariance matrices.

    covariance_type : {'full', 'tied', 'diag', 'spherical'}
        The type of precision matrices.

    Returns
    -------
    nk : array-like of shape (n_components,)
        The numbers of data samples in the current components.

    means : array-like of shape (n_components, n_features)
        The centers of the current components.

    covariances : array-like
        The covariance matrix of the current components.
        The shape depends of the covariance_type.
    """
    nk = resp.sum(axis=0) + 10 * np.finfo(resp.dtype).eps
    means = np.dot(resp.T, X) / nk[:, np.newaxis]
    if fix_mu_flag:
        means[0][0] = fix_mu
    covariances = {"full": _estimate_gaussian_covariances_full}[covariance_type](resp, X, nk, means, reg_covar, max_sigma)
    return nk, means, covariances


def _estimate_gaussian_covariances_full(resp, X, nk, means, reg_covar, max_sigma):
    """Estimate the full covariance matrices.

    Parameters
    ----------
    resp : array-like of shape (n_samples, n_components)

    X : array-like of shape (n_samples, n_features)

    nk : array-like of shape (n_components,)

    means : array-like of shape (n_components, n_features)

    reg_covar : float

    Returns
    -------
    covariances : array, shape (n_components, n_features, n_features)
        The covariance matrix of the current components.
    """
    MAX_COV = max_sigma ** 2
    n_components, n_features = means.shape
    covariances = np.empty((n_components, n_features, n_features))
    for k in range(n_components):
        diff = X - means[k]
        covariances[k] = np.dot(resp[:, k] * diff.T, diff) / nk[k]
        covariances[k].flat[::n_features + 1] += reg_covar
        covariances[k] = np.where(covariances[k] > MAX_COV, MAX_COV, covariances[k])
    return covariances


def trans_lable(_cluster_centers, _cluster_label, _fix_mu):
    _centers = _cluster_centers.flatten()
    _centers = np.abs(_centers - _fix_mu)
    _idx = np.argmin(_centers)
    if _idx == 0:
        return _cluster_label
    else:
        _new_label = np.zeros((len(_cluster_label),), dtype=int)
        for i in range(len(_cluster_label)):
            if _cluster_label[i] == _idx:
                _new_label[i] = 0
            elif _cluster_label[i] == 0:
                _new_label[i] = _idx
            else:
                _new_label[i] = _cluster_label[i]
    return _new_label


class GaussianMixtureModel:
    """

    tol :  EM iterations will stop when the lower bound average gain is below this threshold.
    reg_covar :  Non-negative regularization added to the diagonal of covariance. Allows to assure that the covariance
                 matrices are all positive.
    """
    def __init__(self, n_components=1, tol=1e-3, covariance_type='full',
                 reg_covar=1e-6, max_iter=100, n_init=10, init_params='kmeans',
                 random_state=None, verbose=0, verbose_interval=10,
                 fix_mu_flag=False, fix_mu=None, max_sigma=4.0):

        self.n_components = n_components
        self.tol = tol
        self.covariance_type = covariance_type
        self.reg_covar = reg_covar
        self.max_iter = max_iter
        self.n_init = n_init
        self.init_params = init_params
        self.random_state = random_state
        self.verbose = verbose
        self.verbose_interval = verbose_interval

        self.weights_init = None
        self.means_init = None
        self.precisions_init = None
        self.label = None

        self.weights_ = None
        self.means_ = None
        self.covariances_ = None
        self.lower_bound_ = None

        self.fix_mu_flag = fix_mu_flag
        self.fix_mu = fix_mu
        self.max_sigma = max_sigma

        if self.fix_mu_flag:
            assert self.fix_mu is not None

    def _n_parameters(self):
        """Return the number of free parameters in the model."""
        _, n_features = self.means_.shape
        if self.covariance_type == "full":
            cov_params = self.n_components * n_features * (n_features + 1) / 2.0
        elif self.covariance_type == "diag":
            cov_params = self.n_components * n_features
        elif self.covariance_type == "tied":
            cov_params = n_features * (n_features + 1) / 2.0
        elif self.covariance_type == "spherical":
            cov_params = self.n_components
        mean_params = n_features * self.n_components
        return int(cov_params + mean_params + self.n_components - 1)

    def _print_verbose_msg_init_beg(self, n_init):
        """Print verbose message on initialization."""
        if self.verbose == 1:
            print("Initialization %d" % n_init)
        elif self.verbose >= 2:
            print("Initialization %d" % n_init)
            self._init_prev_time = time()
            self._iter_prev_time = self._init_prev_time

    def _print_verbose_msg_init_end(self, ll):
        """Print verbose message on the end of iteration."""
        if self.verbose == 1:
            print("Initialization converged: %s" % self.converged_)
        elif self.verbose >= 2:
            print("Initialization converged: %s\t time lapse %.5fs\t ll %.5f" %
                  (self.converged_, time() - self._init_prev_time, ll))

    def _print_verbose_msg_iter_end(self, n_iter, diff_ll):
        """Print verbose message on initialization."""
        if n_iter % self.verbose_interval == 0:
            if self.verbose == 1:
                print("  Iteration %d" % n_iter)
            elif self.verbose >= 2:
                cur_time = time()
                print("  Iteration %d\t time lapse %.5fs\t ll change %.5f" % (
                    n_iter, cur_time - self._iter_prev_time, diff_ll))
                self._iter_prev_time = cur_time

    def _initialize_parameters(self, X, random_state):
        """Initialize the model parameters.

        Parameters
        ----------
        X : array-like of shape  (n_samples, n_features)

        random_state : RandomState
            A random number generator instance that controls the random seed
            used for the method chosen to initialize the parameters.
        """
        n_samples, _ = X.shape

        if self.init_params == 'kmeans':
            resp = np.zeros((n_samples, self.n_components))
            kmeans = cluster.KMeans(n_clusters=self.n_components, n_init=1,
                                    random_state=random_state).fit(X)
            if self.fix_mu_flag:
                label = trans_lable(kmeans.cluster_centers_, kmeans.labels_, self.fix_mu)
            else:
                label = kmeans.labels_
            self.label = label
            resp[np.arange(n_samples), label] = 1

        elif self.init_params == 'random':
            resp_init = random_state.rand(n_samples, self.n_components)
            resp = resp_init / resp_init.sum(axis=1)[:, np.newaxis]
        else:
            raise ValueError("Unimplemented initialization method '%s'"
                             % self.init_params)

        self._initialize(X, resp)

    def _initialize(self, X, resp):
        """Initialization of the Gaussian mixture parameters.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)

        resp : array-like of shape (n_samples, n_components)
        """
        n_samples, _ = X.shape

        weights, means, covariances = _estimate_gaussian_parameters(
            X, resp, self.reg_covar, self.covariance_type, self.fix_mu_flag, self.fix_mu, self.max_sigma)
        weights /= n_samples

        self.weights_ = (weights if self.weights_init is None
                         else self.weights_init)
        self.means_ = means if self.means_init is None else self.means_init
        self.label = np.argmax(resp, axis=1)

        if self.precisions_init is None:
            self.covariances_ = covariances
            self.precisions_cholesky_ = _compute_precision_cholesky(
                covariances, self.covariance_type)
        elif self.covariance_type == 'full':
            self.precisions_cholesky_ = np.array(
                [linalg.cholesky(prec_init, lower=True)
                 for prec_init in self.precisions_init])
        elif self.covariance_type == 'tied':
            self.precisions_cholesky_ = linalg.cholesky(self.precisions_init,
                                                        lower=True)
        else:
            self.precisions_cholesky_ = self.precisions_init

    def _e_step(self, X):
        """E step.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)

        Returns
        -------
        log_prob_norm : float
            Mean of the logarithms of the probabilities of each sample in X

        log_responsibility : array, shape (n_samples, n_components)
            Logarithm of the posterior probabilities (or responsibilities) of
            the point of each sample in X.
        """
        log_prob_norm, log_resp = self._estimate_log_prob_resp(X)
        return np.mean(log_prob_norm), log_resp

    def _m_step(self, X, log_resp):
        """M step.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)

        log_resp : array-like of shape (n_samples, n_components)
            Logarithm of the posterior probabilities (or responsibilities) of
            the point of each sample in X.
        """
        n_samples, _ = X.shape
        self.weights_, self.means_, self.covariances_ = (
            _estimate_gaussian_parameters(X, np.exp(log_resp), self.reg_covar,
                                          self.covariance_type, self.fix_mu_flag, self.fix_mu, self.max_sigma))
        self.label = np.argmax(log_resp, axis=1)
        self.weights_ /= n_samples
        self.precisions_cholesky_ = _compute_precision_cholesky(
            self.covariances_, self.covariance_type)

    def _estimate_log_prob(self, X):
        return _estimate_log_gaussian_prob(
            X, self.means_, self.precisions_cholesky_, self.covariance_type)

    def _estimate_log_weights(self):
        return np.log(self.weights_)

    def _compute_lower_bound(self, _, log_prob_norm):
        return log_prob_norm

    def _get_parameters(self):
        return (self.weights_, self.means_, self.covariances_,
                self.precisions_cholesky_)

    def _set_parameters(self, params):
        (self.weights_, self.means_, self.covariances_,
         self.precisions_cholesky_) = params

        # Attributes computation
        _, n_features = self.means_.shape

        if self.covariance_type == 'full':
            self.precisions_ = np.empty(self.precisions_cholesky_.shape)
            for k, prec_chol in enumerate(self.precisions_cholesky_):
                self.precisions_[k] = np.dot(prec_chol, prec_chol.T)

        elif self.covariance_type == 'tied':
            self.precisions_ = np.dot(self.precisions_cholesky_,
                                      self.precisions_cholesky_.T)
        else:
            self.precisions_ = self.precisions_cholesky_ ** 2

    def _estimate_weighted_log_prob(self, X):
        """Estimate the weighted log-probabilities, log P(X | Z) + log weights.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)

        Returns
        -------
        weighted_log_prob : array, shape (n_samples, n_component)
        """
        return self._estimate_log_prob(X) + self._estimate_log_weights()

    def _estimate_log_prob_resp(self, X):
        """Estimate log probabilities and responsibilities for each sample.

        Compute the log probabilities, weighted log probabilities per
        component and responsibilities for each sample in X with respect to
        the current state of the model.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)

        Returns
        -------
        log_prob_norm : array, shape (n_samples,)
            log p(X)

        log_responsibilities : array, shape (n_samples, n_components)
            logarithm of the responsibilities
        """
        weighted_log_prob = self._estimate_weighted_log_prob(X)
        log_prob_norm = logsumexp(weighted_log_prob, axis=1)
        with np.errstate(under='ignore'):
            # ignore underflow
            log_resp = weighted_log_prob - log_prob_norm[:, np.newaxis]
        return log_prob_norm, log_resp

    def fit(self, X):
        max_lower_bound = -np.infty
        self.converged_ = False
        n_samples, _ = X.shape

        for init in range(self.n_init):
            self._print_verbose_msg_init_beg(init)
            if init == 0:
                self._initialize_parameters(X, check_random_state(self.random_state))
            else:
                self._initialize_parameters(X, check_random_state(None))

            lower_bound = -np.infty
            for n_iter in range(1, self.max_iter + 1):
                prev_lower_bound = lower_bound

                log_prob_norm, log_resp = self._e_step(X)
                self._m_step(X, log_resp)

                lower_bound = self._compute_lower_bound(
                    log_resp, log_prob_norm)

                change = lower_bound - prev_lower_bound
                self._print_verbose_msg_iter_end(n_iter, change)

                if abs(change) < self.tol:
                    self.converged_ = True
                    break

            self._print_verbose_msg_init_end(lower_bound)

            if lower_bound > max_lower_bound:
                max_lower_bound = lower_bound
                best_params = self._get_parameters()
                best_n_iter = n_iter

        self._set_parameters(best_params)
        self.n_iter_ = best_n_iter
        self.lower_bound_ = max_lower_bound

        # Always do a final e-step to guarantee that the labels returned by
        # fit_predict(X) are always consistent with fit(X).predict(X)
        # for any value of max_iter and tol (and any random_state).
        _, log_resp = self._e_step(X)

        return log_resp.argmax(axis=1)

    def bic(self, X):
        """Bayesian information criterion for the current model on the input X.

        You can refer to this :ref:`mathematical section <aic_bic>` for more
        details regarding the formulation of the BIC used.

        Parameters
        ----------
        X : array of shape (n_samples, n_dimensions)
            The input samples.

        Returns
        -------
        bic : float
            The lower the better.
        """
        score_samples = logsumexp(self._estimate_weighted_log_prob(X), axis=1)
        score = score_samples.mean()
        return -2 * score * X.shape[0] + self._n_parameters() * np.log(
            X.shape[0]
        )