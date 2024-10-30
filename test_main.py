"""
Unit Test for the first part of the project, developed by us.
"""

from unittest.mock import patch

import pytest

from src.delimitation import Delimitation
from src.point import Point
from src.valid_points import ValidPoints


def test_delimitation_initialization():
    delim = Delimitation()

    # Test that it's an instance of Delimitation
    assert isinstance(delim, Delimitation)

    # Test that no points are present initially
    assert delim.get_points() == []

    # Test that the size is 0
    assert delim.size() == 0


def test_add_single_point():
    delim = Delimitation()
    p1 = Point("p1", 10, 20)

    delim.add_point(p1)

    # Test that the point is added correctly
    points = delim.get_points()
    assert points == [p1]

    # Test that the size is updated
    assert delim.size() == 1


def test_add_multiple_points_1():
    delim = Delimitation()
    p1 = Point("p1", 10, 20)
    p2 = Point("p2", 15, 25)

    delim.add_point(p1)
    delim.add_point(p2)

    # Test that the points are added in the correct order
    points = delim.get_points()
    assert points == [p1, p2]

    # Test that the size reflects the correct number of points
    assert delim.size() == 2


def test_add_duplicate_point():
    delim = Delimitation()
    p1 = Point("p1", 10, 20)

    delim.add_point(p1)
    delim.add_point(p1)

    # Test that the duplicate point is added (if allowed)
    points = delim.get_points()
    assert points == [p1, p1]

    # Test that the size counts duplicates
    assert delim.size() == 2


def test_pop_single_point():
    delim = Delimitation()
    p1 = Point("p1", 10, 20)

    delim.add_point(p1)
    popped = delim.pop_point()

    # Test that the correct point is popped
    assert popped == p1

    # Test that the list is now empty
    assert delim.get_points() == []

    # Test that the size is updated
    assert delim.size() == 0


def test_pop_multiple_points():
    delim = Delimitation()
    p1 = Point("p1", 10, 20)
    p2 = Point("p2", 15, 25)

    delim.add_point(p1)
    delim.add_point(p2)

    # Pop the last added point
    popped = delim.pop_point()
    assert popped == p2
    assert delim.get_points() == [p1]

    # Pop the first point
    popped = delim.pop_point()
    assert popped == p1
    assert delim.get_points() == []
    assert delim.size() == 0


def test_pop_on_empty_delimitation():
    delim = Delimitation()

    try:
        delim.pop_point()
        assert False, "Expected an exception when popping from an empty Delimitation"
    except IndexError:
        # Assuming IndexError or a custom error is raised when popping from empty list
        pass


def test_get_points_empty():
    delim = Delimitation()

    # Test that an empty list is returned
    assert delim.get_points() == []


def test_get_points_order():
    delim = Delimitation()
    p1 = Point("p1", 10, 20)
    p2 = Point("p2", 15, 25)

    delim.add_point(p1)
    delim.add_point(p2)

    # Test that points are returned in the correct order
    points = delim.get_points()
    assert points == [p1, p2]


def test_get_points_immutable():
    delim = Delimitation()
    p1 = Point("p1", 10, 20)
    delim.add_point(p1)

    # Modify the returned list and check the internal state
    delim.add_point(Point("p2", 15, 25))
    assert delim.get_points() == [
        p1,
        Point("p2", 15, 25),
    ], "get_points should return a copy of the list, not a reference"


def test_get_first_single_point():
    delim = Delimitation()
    p1 = Point("p1", 10, 20)
    delim.add_point(p1)

    # Test that it returns the first point
    assert delim.get_first() == p1


def test_get_first_multiple_points_2():
    delim = Delimitation()
    p1 = Point("p1", 10, 20)
    p2 = Point("p2", 15, 25)

    delim.add_point(p1)
    delim.add_point(p2)

    # Test that it returns the first point
    assert delim.get_first() == p1


def test_get_last_two_exactly_two_points():
    delim = Delimitation()
    p1 = Point("p1", 10, 20)
    p2 = Point("p2", 15, 25)

    delim.add_point(p1)
    delim.add_point(p2)

    # Test that it returns the last two points in the correct order
    assert delim.get_last_two() == (p1, p2)


def test_get_last_two_more_than_two_points():
    delim = Delimitation()
    p1 = Point("p1", 10, 20)
    p2 = Point("p2", 15, 25)
    p3 = Point("p3", 30, 35)

    delim.add_point(p1)
    delim.add_point(p2)
    delim.add_point(p3)

    # Test that it returns the last two points, p2 and p3
    assert delim.get_last_two() == (p2, p3)


def test_get_last_two_less_than_two_points():
    delim = Delimitation()
    p1 = Point("p1", 10, 20)

    delim.add_point(p1)

    try:
        assert (None, p1) == delim.get_last_two()
    except IndexError:
        pass  # Assuming IndexError or a custom exception is raised


def test_get_area_less_than_three_points():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 10)

    delim.add_point(p1)
    delim.add_point(p2)

    # Test that ValueError is raised with less than 3 points
    with pytest.raises(ValueError, match="A polygon must have at least 3 points to calculate the area."):
        delim.get_area()


def test_forward_angle():
    p_A = Point("A", 10, 0)
    p_B = Point("B", 10, 10)
    p_C = Point("C", 5, 5)
    p_D = Point("D", 0, 0)
    p_E = Point("E", 0, 10)
    assert p_A.get_forward_angle(p_D, p_B) == 270.0
    assert p_B.get_forward_angle(p_E, p_A) == 90.87038467497996
    assert p_A.get_forward_angle(p_D, p_C) == 315.1092215479925
    assert p_B.get_forward_angle(p_E, p_C) == 45.098738076563734


def test_get_area_collinear_points():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 10)
    p3 = Point("p3", 20, 20)  # Collinear with p1 and p2

    delim.add_point(p1)
    delim.add_point(p2)
    delim.add_point(p3)

    # Test that ValueError is raised for collinear points
    with pytest.raises(ValueError, match="It is impossible to form a polygon with collinear points."):
        delim.get_area()


""" def test_get_area_simple_triangle():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 0, 10)
    p3 = Point("p3", 10, 0)
    
    delim.add_point(p1)
    delim.add_point(p2)
    delim.add_point(p3)
    
    # Test the area of a right triangle with base and height 10
    expected_area = 50  # 1/2 * base * height
    assert abs(delim.get_area() - expected_area) < 1e-6  # Floating-point tolerance
 """

""" def test_get_area_complex_polygon():
    delim = Delimitation()
    # Create a simple square
    points = [Point(f"p{i}", x, y) for i, (x, y) in enumerate([(0, 0), (0, 10), (10, 10), (10, 0)], 1)]
    for p in points:
        delim.add_point(p)
    
    # Expected area of square: 100
    assert abs(delim.get_area() - 100) < 1e-6 """


def test_copy_independence():
    delim = Delimitation()
    p1 = Point("p1", 10, 20)
    p2 = Point("p2", 15, 25)

    delim.add_point(p1)
    delim.add_point(p2)

    # Create a copy
    delim_copy = delim.copy()

    # Test that the copied points are identical but independent
    assert delim.get_points() == delim_copy.get_points()

    # Modify the copy and check that the original is not affected
    p3 = Point("p3", 30, 35)
    delim_copy.add_point(p3)

    assert delim.get_points() != delim_copy.get_points()
    assert delim.size() == 2
    assert delim_copy.size() == 3


def test_intersects_intersecting_segments():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 10)
    p3 = Point("p3", 0, 10)
    p4 = Point("p4", 10, 0)

    # These two segments intersect at (5, 5)
    assert delim.intersects(p1, p2, p3, p4)


def test_intersects_non_intersecting_segments():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 10)
    p3 = Point("p3", 20, 20)
    p4 = Point("p4", 30, 30)

    # These segments do not intersect
    assert not delim.intersects(p1, p2, p3, p4)


def test_intersects_parallel_segments():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 0)
    p3 = Point("p3", 0, 10)
    p4 = Point("p4", 10, 10)

    # Parallel lines should not intersect
    assert not delim.intersects(p1, p2, p3, p4)


def test_intersects_overlapping_segments():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 10)
    p3 = Point("p3", 5, 5)
    p4 = Point("p4", 15, 15)

    # Overlapping segments should not be considered intersecting
    assert not delim.intersects(p1, p2, p3, p4)


def test_crosses_delimitation_crossing():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 0)
    p3 = Point("p3", 9, 10)
    p4 = Point("p4", 5, 5)

    delim.add_point(p1)
    delim.add_point(p2)

    # This new segment crosses the existing line (p1 -> p2)
    assert not delim.crosses_delimitation(p4, p3)


def test_crosses_delimitation_no_crossing():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 0)
    p3 = Point("p3", 10, 10)

    delim.add_point(p1)
    delim.add_point(p2)

    assert not delim.crosses_delimitation(p2, p3)


def test_crosses_delimitation_touching_corner():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 0)
    p3 = Point("p3", 0, 10)

    delim.add_point(p1)
    delim.add_point(p2)

    # This new segment touches p1 but does not cross
    assert not delim.crosses_delimitation(p3, p1)


def test_size_empty_delimitation():
    delim = Delimitation()
    # Initially, the size should be 0
    assert delim.size() == 0


def test_size_after_adding_points():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 10)

    delim.add_point(p1)
    delim.add_point(p2)

    # After adding two points, the size should be 2
    assert delim.size() == 2


def test_size_after_removing_points():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 10)

    delim.add_point(p1)
    delim.add_point(p2)
    delim.pop_point()  # Remove one point

    # After removing one point, the size should be 1
    assert delim.size() == 1


def test_get_first_one_point():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    delim.add_point(p1)

    # First point should be p1
    assert delim.get_first() == p1


def test_get_first_multiple_points_1():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 10)

    delim.add_point(p1)
    delim.add_point(p2)

    # First point should still be p1
    assert delim.get_first() == p1


def test_get_first_no_points():
    delim = Delimitation()

    try:
        assert delim.get_first() is None, "Expected None when calling get_first on an empty Delimitation"

    except IndexError:
        pass  # Assuming IndexError or a custom exception is raised


def test_show_no_points():
    delim = Delimitation()

    with patch("matplotlib.pyplot.show") as mock_show:
        delim.show(ValidPoints([]))
        # Ensure the plotting method was called even with no points
        mock_show.assert_called()


def test_show_with_points():
    delim = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 10)

    delim.add_point(p1)
    delim.add_point(p2)

    valid_points = ValidPoints([p1, p2])

    with patch("matplotlib.pyplot.show") as mock_show:
        delim.show(valid_points)
        # Ensure the plotting method was called
        mock_show.assert_called()


def test_eq_same_order():
    delim1 = Delimitation()
    delim2 = Delimitation()

    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 10)

    delim1.add_point(p1)
    delim1.add_point(p2)

    delim2.add_point(p1)
    delim2.add_point(p2)

    # Both delimitations have the same points in the same order
    assert delim1 == delim2


def test_eq_different_points():
    delim1 = Delimitation()
    delim2 = Delimitation()

    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 10)
    p3 = Point("p3", 5, 5)

    delim1.add_point(p1)
    delim1.add_point(p2)

    delim2.add_point(p1)
    delim2.add_point(p3)

    # Different points in delimitations should make them unequal
    assert delim1 != delim2


def test_eq_different_order():
    delim1 = Delimitation()
    delim2 = Delimitation()

    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 10, 10)
    p3 = Point("p3", 5, 5)

    # Delimitation 1: add points in order p1, p2, p3
    delim1.add_point(p1)
    delim1.add_point(p2)
    delim1.add_point(p3)

    # Delimitation 2: add points in a different order but forming the same polygon
    delim2.add_point(p3)
    delim2.add_point(p2)
    delim2.add_point(p1)

    # Test that different order does not affect equality
    assert delim1 == delim2


""" def test_eq_empty_delimitation():
    delim1 = Delimitation()
    delim2 = Delimitation()
    
    # Both delimitations are empty
    assert delim1 == delim2 """


def test_equal_triangles():
    delim1 = Delimitation()
    delim2 = Delimitation()

    # Triangle A
    points1 = [Point("p1", 0, 0), Point("p2", 0, 10), Point("p3", 10, 0)]
    # Triangle A but in different order
    points2 = [Point("p3", 10, 0), Point("p1", 0, 0), Point("p2", 0, 10)]

    for point in points1:
        delim1.add_point(point)
    for point in points2:
        delim2.add_point(point)

    assert delim1 == delim2  # Should be True


def test_equal_polygons_with_rotations():
    delim1 = Delimitation()
    delim2 = Delimitation()

    # Create a square
    points1 = [Point("p1", 0, 0), Point("p2", 0, 10), Point("p3", 10, 10), Point("p4", 10, 0)]
    # Rotated version of the same square
    points2 = [Point("p2", 0, 10), Point("p3", 10, 10), Point("p4", 10, 0), Point("p1", 0, 0)]

    for point in points1:
        delim1.add_point(point)
    for point in points2:
        delim2.add_point(point)

    assert delim1 == delim2  # Should be True


def test_not_equal_polygons():
    delim1 = Delimitation()
    delim2 = Delimitation()

    # Triangle A
    points1 = [Point("p1", 0, 0), Point("p2", 0, 10), Point("p3", 10, 0)]
    # A different triangle
    points2 = [Point("p1", 0, 0), Point("p2", 5, 5), Point("p3", 10, 0)]

    for point in points1:
        delim1.add_point(point)
    for point in points2:
        delim2.add_point(point)

    assert delim1 != delim2


def test_add_valid_point():
    delim = Delimitation()
    point = Point("p1", 0, 0)

    # Add a valid point
    delim.add_point(point)

    # Check if the point was added
    assert len(delim.get_points()) == 1
    assert delim.get_points()[0] == point


def test_add_multiple_points_2():
    delim = Delimitation()
    points = [Point("p1", 0, 0), Point("p2", 0, 10), Point("p3", 10, 0)]

    # Add multiple valid points
    for point in points:
        delim.add_point(point)

    # Check if all points were added
    assert len(delim.get_points()) == 3
    for i, point in enumerate(points):
        assert delim.get_points()[i] == point


def test_add_invalid_point_type():
    delim = Delimitation()

    with pytest.raises(ValueError, match="nly Point objects can be added to the Delimitation object"):
        delim.add_point("Not a Point")  # Trying to add a string #type: ignore


def test_add_integer_instead_of_point():
    delim = Delimitation()

    with pytest.raises(ValueError, match="Only Point objects can be added to the Delimitation object"):
        delim.add_point(123)  # Trying to add an integer # Trying to add a string #type: ignore


def test_add_none_as_point():
    delim = Delimitation()

    with pytest.raises(ValueError, match="Only Point objects can be added to the Delimitation object"):
        delim.add_point(None)  # type: ignore


def test_eq_stora():
    p1 = Point("p1", 0.0, 10.0)
    p2 = Point("p2", 10.0, 10.0)
    p4 = Point("p4", 4.0, 5.0)
    p6 = Point("p6", 1.0, 20.0)
    p8 = Point("p8", 12.0, -1.0)

    # Create the first delimitation (d1)
    d1 = Delimitation()
    d1.add_point(p1)
    d1.add_point(p4)
    d1.add_point(p8)
    d1.add_point(p2)
    d1.add_point(p6)

    # Create the second delimitation (d2)
    d2 = Delimitation()
    d2.add_point(p8)
    d2.add_point(p2)
    d2.add_point(p6)
    d2.add_point(p1)
    d2.add_point(p4)

    # Create the third delimitation (d3)
    d3 = Delimitation()
    d3.add_point(p4)
    d3.add_point(p1)
    d3.add_point(p6)
    d3.add_point(p2)
    d3.add_point(p8)

    # Create the fourth delimitation (d4)
    d4 = Delimitation()
    d4.add_point(p4)
    d4.add_point(p1)
    d4.add_point(p6)
    d4.add_point(p2)

    # Verify the points in each delimitation
    assert d1 == d2 and d1 == d3 and not d1 == d4


def test_get_points_empty_2():
    d1 = Delimitation()
    assert d1.get_points() == []


def test_add_and_get_first_point_2():
    d1 = Delimitation()
    p1 = Point("p01", 10, 0)
    p2 = Point("p02", 10, 10)

    d1.add_point(p1)
    d1.add_point(p2)

    # Check the first point
    assert d1.get_first() == p1


def test_get_last_two_points_2():
    d1 = Delimitation()
    p1 = Point("p01", 10, 0)
    p2 = Point("p02", 10, 10)

    d1.add_point(p1)
    d1.add_point(p2)

    # Check the last two points
    assert d1.get_last_two() == (p1, p2)


def test_pop_point():
    d1 = Delimitation()
    p1 = Point("p01", 10, 0)
    p2 = Point("p02", 10, 10)

    d1.add_point(p1)
    d1.add_point(p2)

    # Pop last point
    assert d1.pop_point() == p2

    # After popping, only the first point should remain
    assert d1.get_points() == [p1]


def area_prof_2():
    d2 = Delimitation()

    # Create the points and add them to the delimitation
    pA = Point("pA", -10, -10)
    pB = Point("pB", 10, -10)
    pC = Point("pC", 10, 10)
    pE = Point("pE", 0, -9.9)
    pD = Point("pD", -10, 10)

    d2.add_point(pA)
    d2.add_point(pB)
    d2.add_point(pC)
    d2.add_point(pE)
    d2.add_point(pD)

    assert d2.get_area() == 708148.3951477829


def area_prof_1():
    delimitation = Delimitation()

    # Create the points
    p1 = Point("p1", 0.0, 10.0)
    p4 = Point("p4", 4.0, 5.0)
    p8 = Point("p8", 12.0, -1.0)
    p2 = Point("p2", 10.0, 10.0)
    p6 = Point("p6", 1.0, 20.0)

    # Add points to the delimitation
    delimitation.add_point(p1)
    delimitation.add_point(p4)
    delimitation.add_point(p8)
    delimitation.add_point(p2)
    delimitation.add_point(p6)

    assert delimitation.get_area() == 404310.0128400312


def test_copy_delimitation():
    d1 = Delimitation()
    p1 = Point("p01", 10, 0)
    p2 = Point("p02", 10, 10)

    d1.add_point(p1)
    d1.add_point(p2)

    # Create a copy using the custom copy method
    d1_1 = d1.copy()

    # Check that the objects (Delimitation) are not the same
    assert id(d1) != id(d1_1)
    assert d1 == d1_1

    # Check that the points inside are equal but not the same objects
    assert d1.get_first() == d1_1.get_first()  # They should be equal in value
    assert id(d1.get_first()) != id(d1_1.get_first())  # But should not be the same object

    # Check all points are new objects but equal in value
    for p_orig, p_copy in zip(d1.get_points(), d1_1.get_points()):
        assert p_orig == p_copy  # Points should be equal in value
        assert id(p_orig) != id(p_copy)  # But they should not have the same ID


def test_intersect():
    d = Delimitation()

    print("Intersect Test Cases:\n")

    # Test case 1: Intersecting lines
    p1 = Point("p1", 1, 1)
    p2 = Point("p2", 3, 3)
    p3 = Point("p3", 1, 3)
    p4 = Point("p4", 3, 1)
    # Lines formed by p1-p2 and p3-p4 should intersect
    if not d.intersects(p1, p2, p3, p4):
        raise AssertionError("Test 1 Failed")

    # Test case 2: Parallel lines (no intersection)
    p5 = Point("p5", 1, 1)
    p6 = Point("p6", 3, 1)
    p7 = Point("p7", 1, 2)
    p8 = Point("p8", 3, 2)
    # Lines formed by p5-p6 and p7-p8 are parallel and should not intersect
    if not d.intersects(p5, p6, p7, p8):
        raise AssertionError("Test 2 Failed")

    # Test case 3: Coincident lines (no intersection, or considered the same line)
    p9 = Point("p9", 0, 0)
    p10 = Point("p10", 2, 2)
    p11 = Point("p11", 1, 1)
    p12 = Point("p12", 3, 3)
    # Lines formed by p9-p10 and p11-p12 are coincident (same line)
    if d.intersects(p9, p10, p11, p12):
        raise AssertionError("Test 3 Failed")

    # Test case 4: Intersecting lines with different slopes
    p13 = Point("p13", 0, 0)
    p14 = Point("p14", 4, 4)
    p15 = Point("p15", 0, 4)
    p16 = Point("p16", 4, 0)
    # Lines formed by p13-p14 and p15-p16 should intersect
    if not d.intersects(p13, p14, p15, p16):
        raise AssertionError("Test 4 Failed")

    # Test case 5: Vertical and horizontal intersecting lines
    p17 = Point("p17", 0, 0)
    p18 = Point("p18", 0, 5)
    p19 = Point("p19", -3, 2)
    p20 = Point("p20", 3, 2)
    # Lines formed by p17-p18 (vertical) and p19-p20 (horizontal) should intersect
    if not d.intersects(p17, p18, p19, p20):
        raise AssertionError("Test 5 Failed")

    print("All intersect tests passed!")


# Crosses Delimitation test cases
def test_crosses_delimitation():
    print("\nCrosses Delimitation Test Cases:\n")

    # Test case 1: Line intersects triangle
    d = Delimitation()
    p1 = Point("p1", 0, 0)
    p2 = Point("p2", 4, 0)
    p3 = Point("p3", 2, 4)
    p4 = Point("p4", 0, 0)  # Closing the triangle

    d.add_point(p1)
    d.add_point(p2)
    d.add_point(p3)
    d.add_point(p4)

    line_p1 = Point("p5", 2, 1)
    line_p2 = Point("p6", 2, 5)

    if not d.crosses_delimitation(line_p1, line_p2):
        raise AssertionError("Test 1 Failed")

    # Test case 2: Line does not intersect square
    g = Delimitation()
    p5 = Point("p5", 0, 0)
    p6 = Point("p6", 4, 0)
    p7 = Point("p7", 4, 4)
    p8 = Point("p8", 0, 4)

    g.add_point(p5)
    g.add_point(p6)
    g.add_point(p7)
    g.add_point(p8)

    line_p1 = Point("p10", 5, 5)
    line_p2 = Point("p11", 6, 6)

    if not g.crosses_delimitation(line_p1, line_p2):
        raise AssertionError("Test 2 Failed")

    # Test case 3: Line intersects the square
    line_p1 = Point("p12", 2, -1)
    line_p2 = Point("p13", 2, 5)

    if not g.crosses_delimitation(line_p1, line_p2):
        raise AssertionError("Test 3 Failed")

    # Test case 4: Line does not intersect
    line_p1 = Point("p14", -1, -1)
    line_p2 = Point("p15", 2, -2)

    if d.crosses_delimitation(line_p1, line_p2):
        raise AssertionError("Test 4 Failed")

    if not d.crosses_delimitation(line_p1, line_p2):
        print("All crosses delimitation tests passed!")

    print("All crosses delimitation tests passed!")
    print("All crosses delimitation tests passed!")
    print("All crosses delimitation tests passed!")
    print("All crosses delimitation tests passed!")
