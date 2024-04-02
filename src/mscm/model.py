"""
TODO.
"""

import jax.numpy as jnp
import liesel.model as lsl
import tensorflow_probability.substrates.jax.distributions as tfd


class MSCMMultinomial(tfd.Multinomial):
    def _log_prob(self, counts):
        probs = jnp.where(counts > 0.1, self._probs, 1.0)
        return jnp.sum(jnp.log(probs) * counts, axis=-1)


def make_model(y, x, landscape=None):
    # N = y.shape[0]  # number of sites
    M = y.shape[1]  # number of species
    K = x.shape[1]  # number of covariates

    v_gamma = lsl.Param(
        value=jnp.zeros(M),
        distribution=lsl.Dist(tfd.Normal, loc=0.0, scale=10.0),
        name="gamma",
    )

    v_x = lsl.Obs(x, name="x")

    v_beta = lsl.Param(
        value=jnp.zeros(K),
        distribution=lsl.Dist(tfd.Normal, loc=0.0, scale=10.0),
        name="beta",
    )

    v_eta = lsl.Var(
        value=lsl.Calc(lambda x, beta: x @ beta, v_x, v_beta),
        name="eta",
    )

    def fn_psi(gamma, eta):
        return 1.0 / (1.0 + jnp.exp(-gamma - eta[:, jnp.newaxis]))

    v_psi = lsl.Var(
        value=lsl.Calc(fn_psi, v_gamma, v_eta),
        name="psi",
    )

    v_z = lsl.Param(
        value=(y > 0.1).astype("float32"),
        distribution=lsl.Dist(tfd.Bernoulli, probs=v_psi),
        name="z",
    )

    v_mu = lsl.Param(
        value=jnp.nan_to_num(y.mean(axis=0, where=y > 0.1), nan=1.0),
        distribution=lsl.Dist(tfd.HalfNormal, scale=10.0),
        name="mu",
    )

    v_lambda = lsl.Var(
        value=lsl.Calc(lambda z, mu: jnp.sum(z * mu, axis=1), v_z, v_mu),
        name="lambda",
    )

    v_n = lsl.Obs(
        value=y.sum(axis=1),
        distribution=lsl.Dist(tfd.Poisson, rate=v_lambda),
        name="n",
    )

    def fn_p(z, mu):
        rowsums = jnp.sum(z * mu, axis=1)
        rowsums = jnp.where(jnp.sum(z, axis=1) > 0.1, rowsums, 1.0)
        return (z * mu) / rowsums[:, jnp.newaxis]

    v_p = lsl.Var(
        value=lsl.Calc(fn_p, v_z, v_mu),
        name="p",
    )

    v_y = lsl.Obs(
        value=y,
        distribution=lsl.Dist(MSCMMultinomial, total_count=v_n, probs=v_p),
        name="y",
    )

    v_richness = lsl.Var(
        value=lsl.Calc(lambda z: jnp.sum(z > 0.1, axis=1), v_z),
        name="richness",
    )

    v_shannon = lsl.Var(
        # p * jnp.log(p) == nan sometimes
        value=lsl.Calc(lambda p: -jnp.sum(jnp.log(p**p), axis=1), v_p),
        name="shannon",
    )

    gb = lsl.GraphBuilder().add(v_n, v_y, v_richness, v_shannon)

    if landscape is not None:
        v_landscape = lsl.Var(landscape, name="landscape")

        def fn_richness_landscape(ls, z):
            return jnp.sum((ls @ z) > 0.1, axis=1)

        v_richness_landscape = lsl.Var(
            value=lsl.Calc(fn_richness_landscape, v_landscape, v_z),
            name="richness_landscape",
        )

        def fn_shannon_landscape(ls, p):
            p = ls @ p
            p = p / jnp.sum(p, axis=1)[:, jnp.newaxis]
            return -jnp.sum(jnp.log(p**p), axis=1)

        v_shannon_landscape = lsl.Var(
            value=lsl.Calc(fn_shannon_landscape, v_landscape, v_p),
            name="shannon_landscape",
        )

        gb.add(v_richness_landscape, v_shannon_landscape)

    return gb.build_model()
