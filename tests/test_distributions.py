import pytest

import mscm.distributions as dists


def approx(x):
    return pytest.approx(x, rel=1e-5 * abs(x))


def test_neg_bin_loc_scale():
    dist = dists.NegBinLocScale(1.0, 1.0)

    assert dist.log_prob(0.0) == approx(-0.6931472)
    assert dist.log_prob(1.0) == pytest.approx(-1.3862944)
    assert dist.log_prob(2.0) == pytest.approx(-2.0794415)

    dist = dists.NegBinLocScale(2.0, 1.0)

    assert dist.log_prob(0.0) == approx(-1.098612)
    assert dist.log_prob(1.0) == pytest.approx(-1.504077)
    assert dist.log_prob(2.0) == pytest.approx(-1.909543)

    dist = dists.NegBinLocScale(1.0, 2.0)

    assert dist.log_prob(0.0) == pytest.approx(-0.5493061)
    assert dist.log_prob(1.0) == pytest.approx(-1.6479184)
    assert dist.log_prob(2.0) == pytest.approx(-2.3410656)

    dist = dists.NegBinLocScale(2.0, 2.0)

    assert dist.log_prob(0.0) == pytest.approx(-0.804719)
    assert dist.log_prob(1.0) == pytest.approx(-1.721010)
    assert dist.log_prob(2.0) == pytest.approx(-2.231835)


def test_yule():
    dist = dists.Yule(1.0)

    assert dist.log_prob(0.0) == pytest.approx(-0.4054651)
    assert dist.log_prob(1.0) == pytest.approx(-1.7917595)
    assert dist.log_prob(2.0) == pytest.approx(-2.7080502)

    dist = dists.Yule(2.0)

    assert dist.log_prob(0.0) == pytest.approx(-0.5108256)
    assert dist.log_prob(1.0) == pytest.approx(-1.7635886)
    assert dist.log_prob(2.0) == pytest.approx(-2.5745188)

    dist = dists.Yule(3.0)

    assert dist.log_prob(0.0) == approx(-0.5596158)
    assert dist.log_prob(1.0) == pytest.approx(-1.7635886)
    assert dist.log_prob(2.0) == pytest.approx(-2.5367785)
