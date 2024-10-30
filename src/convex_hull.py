"""
This module contains the ConvexHull class, which is responsible for calculating the convex hull of a set of points.
"""

from src.delimitation import Delimitation
from src.point import Point
from src.valid_points import ValidPoints


class ConvexHull:
    """
    A class to calculate the convex hull from a set of valid points.

    The convex hull is the smallest convex shape that encloses a set of points.
    This class uses the orientation method to find the hull by iterating over the points,
    ensuring that only those points that form part of the hull are retained.

    Attributes:
        __points (ValidPoints): The set of valid points from which the convex hull is calculated.

    Methods:
    find_delimitation():Returns the convex hull (delimitation) for the given set of points.
    """

    def __init__(self, valid_points: ValidPoints):
        """
        Initializes the ConvexHull class with a set of valid points.

        Args:
            valid_points (ValidPoints): An instance containing the points from which the convex hull is to be computed.

        Complexity: O(1)
        """
        self.__points = valid_points

    def find_delimitation(self) -> Delimitation:
        """
        Computes the convex hull (Delimitation) for the given set of points.
        This method identifies the starting point with the lowest longitude and
        iteratively adds points to the convex hull based on their angular orientation
        relative to the last two points added to the hull. The process continues until
        the starting point is reached again.

        Returns:
            Delimitation: An object containing the points that make up the convex hull.

        Complexity:
            O(p^2), where p is the number of points in the Delimitation object.
        """

        delimitation = Delimitation()
        points = self.__points.get_all_points()
        start = min(points, key=lambda p: p.get_longitude())
        delimitation.add_point(start)
        points.remove(start)

        while len(points) > 0:
            if delimitation.size() < 2:
                p2 = Point("ref", start.get_latitude() - 1, start.get_longitude() - 1)
                p3 = start
                next_point = min(points, key=lambda p: p3.get_forward_angle(p2, p))
                delimitation.add_point(next_point)
                points.remove(next_point)
                continue

            (p2, p3) = delimitation.get_last_two()

            if delimitation.size() < 3:
                next_point = min(points, key=lambda p: p3.get_forward_angle(p2, p))
            else:
                points.append(start)
                next_point = min(points, key=lambda p: p3.get_forward_angle(p2, p))

            if next_point == start:
                break

            delimitation.add_point(next_point)
            points.remove(next_point)

        return delimitation
