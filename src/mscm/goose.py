"""
TODO.
"""

import jax.lax
import jax.numpy as jnp
import jax.random
import liesel.goose as gs


def make_z_gibbs_kernel(model):
    model = model._copy_computational_model()

    def transition(prng_key, model_state):
        model.state = model_state

        for node in model.nodes.values():
            node._outdated = False

        z = model.vars["z"].value
        y = model.vars["y"].value

        N = z.shape[0]
        M = z.shape[1]

        shift = model.log_prob
        keys = jax.random.split(prng_key, N * M)
        keys = jnp.reshape(keys, (N, M, 2))

        def transition_i(i, z):
            def transition_j(j, z):
                model.vars["z"].value = z.at[i, j].set(1.0)
                p1 = jnp.exp(model.log_prob - shift)

                model.vars["z"].value = z.at[i, j].set(0.0)
                p0 = jnp.exp(model.log_prob - shift)

                key = keys[i, j]
                p = p1 / (p0 + p1)
                z_ij = 1.0 * jax.random.bernoulli(key, p)
                z_ij = jax.lax.select(y[i, j] < 0.1, z_ij, 1.0)
                return z.at[i, j].set(z_ij)

            return jax.lax.fori_loop(0, M - 1, transition_j, z)

        z = jax.lax.fori_loop(0, N - 1, transition_i, z)
        return {model.vars["z"].value_node.name: z}

    return gs.GibbsKernel(["z"], transition)
