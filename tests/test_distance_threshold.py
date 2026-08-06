import pytest

from compute_distance import is_within_distance_threshold


@pytest.mark.parametrize(
    ("distance", "expected"),
    [
        (0.00, True),
        (4.99, True),
        (5.00, True),
        (5.01, False),
        (31.74, False),
        (35.00, False),
    ],
)
def test_is_within_distance_threshold(
    distance: float,
    expected: bool,
) -> None:
    assert is_within_distance_threshold(distance) is expected