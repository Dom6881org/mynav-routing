from pathlib import Path

import numpy as np
import pytest

from mynav_routing import Polar

EX = Path(__file__).resolve().parent.parent / "examples" / "polaire-sugar.pol"

POL_SAMPLE = (
    "TWA/TWS, 4, 8, 12, 16\n"
    "0, 0.5, 1.0, 1.2, 1.3\n"
    "45, 2.5, 4.0, 5.0, 5.5\n"
    "90, 3.0, 5.0, 6.5, 7.0\n"
    "135, 2.8, 5.2, 6.8, 7.2\n"
    "180, 2.5, 4.5, 6.0, 6.5\n"
)


def test_parse_sample(tmp_path):
    p = tmp_path / "sample.pol"
    p.write_text(POL_SAMPLE)
    polar = Polar.from_file(str(p))
    assert polar.tws.tolist() == [4.0, 8.0, 12.0, 16.0]
    assert polar.twa.tolist() == [0.0, 45.0, 90.0, 135.0, 180.0]
    assert polar.speed(12, 90) == pytest.approx(6.5)


def test_parse_sugar_polar():
    polar = Polar.from_file(str(EX))
    assert polar.tws.tolist() == [4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 24.0]
    assert polar.twa[0] == 0.0
    assert polar.twa[-1] == 180.0
    # les lignes Up/Down (VMG) ne doivent pas créer de faux TWA
    assert len(polar.twa) == 16
    assert polar.speed(12, 90) == pytest.approx(6.4)
    # symétrie bâbord/tribord
    assert polar.speed(12, -90) == pytest.approx(6.4)


def test_interpolation():
    polar = Polar(
        [4, 8, 12],
        [0, 90, 180],
        np.array([[0.5, 1.0, 1.2], [3.0, 5.0, 6.5], [2.5, 4.5, 6.0]]),
    )
    assert polar.speed(10, 90) == pytest.approx(5.75)
    assert polar.speed(10, 45) == pytest.approx(3.425)
