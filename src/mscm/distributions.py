"""
TODO.
"""

import jax
import jax.numpy as jnp
import jax.scipy.special
import tensorflow_probability.substrates.jax.bijectors as tfb
import tensorflow_probability.substrates.jax.distributions as tfd
from tensorflow_probability.python.internal.reparameterization import NOT_REPARAMETERIZED
from tensorflow_probability.substrates.jax.internal.parameter_properties import ParameterProperties


class NegBinLocScale(tfd.NegativeBinomial):
    def __init__(
        self,
        loc,
        scale,
        validate_args=False,
        allow_nan_stats=True,
        name="NegBinLocScale",
    ):
        parameters = dict(locals())
        total_count = 1.0 / scale

        super().__init__(
            total_count=total_count,
            probs=1.0 - total_count / (total_count + loc),
            validate_args=validate_args,
            allow_nan_stats=allow_nan_stats,
            require_integer_total_count=False,
            name=name,
        )

        self._parameters = parameters

    @classmethod
    def _parameter_properties(cls, dtype, num_classes=None):
        return {
            "loc": ParameterProperties(default_constraining_bijector_fn=lambda: tfb.Exp()),
            "scale": ParameterProperties(default_constraining_bijector_fn=lambda: tfb.Exp()),
        }


class Yule(tfd.Distribution):
    def __init__(
        self,
        loc,
        validate_args=False,
        allow_nan_stats=True,
        name="Yule",
    ):
        parameters = dict(locals())

        super().__init__(
            dtype=jnp.float32,
            reparameterization_type=NOT_REPARAMETERIZED,
            validate_args=validate_args,
            allow_nan_stats=allow_nan_stats,
            parameters=parameters,
            name=name,
        )

    def _log_prob(self, x):
        loc_inv = 1.0 / self._parameters["loc"]
        return jnp.log(loc_inv + 1.0) + jax.scipy.special.betaln(x + 1.0, loc_inv + 2.0)

    @classmethod
    def _parameter_properties(cls, dtype, num_classes=None):
        return {
            "loc": ParameterProperties(default_constraining_bijector_fn=lambda: tfb.Exp()),
        }
