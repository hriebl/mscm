import numpy as np
import pytest

import mscm

rng = np.random.default_rng(1337)


def test_model():
    y = np.zeros(shape=[40, 26], dtype=np.float32)
    x = rng.uniform(size=[40, 2]).astype(np.float32)
    model = mscm.make_model(y, x)

    model.vars["beta"].dist_node.kwinputs["scale"].value = 1.0
    model.vars["gamma"].dist_node.kwinputs["scale"].value = 1.0
    mscm.simulate_model(model, 42)

    assert model.log_prob == pytest.approx(-15301.772)
