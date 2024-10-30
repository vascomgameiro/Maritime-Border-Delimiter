"""
This module contains the Delimitation class which models a polygon that can be constructed by sequentially adding points.
"""

import matplotlib.pyplot as plt

from src.point import Point
from src.valid_points import ValidPoints


class Delimitation:
    """
    The Delimitation ADT models a polygon that can be constructed by sequentially adding
    points.

    Attributes:
        __points (list): A List that stores all the points in the order they are added.

    Methods:
        add_point(Point): Adds a point to the Delimitation and updates relevant attributes.
        get_first(): Returns the first point added to the Delimitation, or None if no points exist.
        get_last_two(): Returns the last two points added as a tuple.
        get_points(): Returns the list of all points in the order they were added.
        pop_point(): Removes and returns the last point added to the Delimitation.
        get_area(): Calculates and returns the area of the polygon.
        copy(): Returns a deep copy of the current Delimitation object.
        intersects(Point, Point, Point, Point): Determines whether two line segments formed by the given points intersect.
        crosses_delimitation(Point, Point): Checks if a line segment formed by two points intersects with any existing segments in the Delimitation.
        show(ValidPoints): Displays the current Delimitation and points using a plotting tool.
    """

    def __init__(self):
        """
        This method initializes an empty Delimitation object.

        Complexity: O(1)
        """
        self.__points = []

    def get_points(self) -> list:
        """
        This method retrieves the Points that belong to the Delimitation object by the order
        they were added.

        Returns: list

        Complexity: O(1)
        """
        return self.__points

    def get_first(self) -> Point | None:
        """
        This method retrieves the first point that was added to the Delimitation object.

        Returns: Point

        Complexity: O(1)
        """
        return self.__points[0] if self.size() > 0 else None

    def get_last_two(self) -> tuple:
        """
        This method retrieves the last two Point objects that have been added to the
        Delimitation object,the rightmost being the last one added.

        Returns: tuple(Point, Point)

        Complexity: O(1)
        """
        if self.size() == 0:
            return (None, None)
        elif self.size() == 1:
            return (None, self.__points[-1])
        return (self.__points[-2], self.__points[-1])

    def get_area(self) -> float:
        """
        This method calculates the area of the polygon that is represented. If
        the Delimitation object has less than 3 points or collinear points,
        the method raises a ValueError.

        Returns: float

        Complexity: O(p), p being the number of points
        """
        vertices = self.get_points()
        n = len(vertices)

        if n < 3:
            raise ValueError("A polygon must have at least 3 points to calculate the area.")

        if self.__collinearity_test(vertices):
            raise ValueError("It is impossible to form a polygon with collinear points.")

        area = 0
        for i in range(n):
            x1, y1 = vertices[i].get_proj_x(), vertices[i].get_proj_y()
            x2, y2 = vertices[(i + 1) % n].get_proj_x(), vertices[(i + 1) % n].get_proj_y()
            area += x1 * y2 - x2 * y1

        return abs(area) / 2

    @staticmethod
    def __collinearity_test(points: list) -> bool:
        """
        This static method checks if a list of Point objects is collinear.

        Arguments: points(list)

        Returns: bool

        Complexity: O(p), p being the number of points
        """
        if len(points) < 3:
            return True

        x1, y1 = points[0].get_latitude(), points[0].get_longitude()
        x2, y2 = points[1].get_latitude(), points[1].get_longitude()

        for i in range(2, len(points)):
            x3, y3 = points[i].get_latitude(), points[i].get_longitude()
            if (x2 - x1) * (y3 - y1) != (y2 - y1) * (x3 - x1):
                return False

        return True

    def size(self) -> int:
        """
        This method calculates the number of points of the Delimitation object.

        Returns: int

        Complexity: O(1)
        """
        return len(self.get_points())

    def add_point(self, point: Point):
        """
        This method adds a Point object to the Delimitation object.

        Arguments: point(Point)

        Complexity: O(1)
        """
        if not isinstance(point, Point):
            raise ValueError("Only Point objects can be added to the Delimitation object")
        if point not in self.get_points():
            self.__points.append(point)

    def pop_point(self) -> Point:
        """
        This method removes the last Point object that had been added to the
        Delimitation object.

        Returns: Point

        Complexity: O(1)
        """
        if self.size() == 0:
            raise IndexError("Cannot remove point from an empty delimitation.")
        return self.__points.pop()

    def copy(self) -> "Delimitation":
        """
        This method performs a deep copy of the Delimitation object.

        Returns: Delimitation

        Complexity: O(p), p being the number of points
        """
        delimitation_copy = Delimitation()
        points = self.get_points()

        for i in range(self.size()):
            point_copy = Point(points[i].get_id(), points[i].get_latitude(), points[i].get_longitude())
            delimitation_copy.add_point(point_copy)

        return delimitation_copy

    @staticmethod
    def __min_max_lat_lon(p1: Point, p2: Point) -> tuple:
        """
        Helper function to return min and max latitudes and longitudes for a pair of points.

        This function takes two Point objects as arguments and returns a tuple of two tuples.
        The first tuple contains the minimum and maximum latitude, and the second tuple
        contains the minimum and maximum longitude.

        Args:
            p1 (Point): The first point.
            p2 (Point): The second point.

        Returns:
            tuple: A tuple of two tuples containing the min and max latitude, and min and max longitude.

        Complexity: O(1)
        """
        return (
            (min(p1.get_latitude(), p2.get_latitude()), max(p1.get_latitude(), p2.get_latitude())),
            (
                min(p1.get_longitude(), p2.get_longitude()),
                max(p1.get_longitude(), p2.get_longitude()),
            ),
        )

    def __is_within_segment(self, p1: Point, p2: Point, intersect_x: float, intersect_y: float) -> bool:
        """
        This static method determines whether a point is within a segment formed by two other points.

        Arguments: p1(Point), p2(Point), p3(Point)

        Returns: bool

        Complexity: O(1)
        """
        (x1_min, x1_max), (y1_min, y1_max) = self.__min_max_lat_lon(p1, p2)
        return (x1_min <= intersect_x <= x1_max) and (y1_min <= intersect_y <= y1_max)

    def intersects(self, p1: Point, p2: Point, p3: Point, p4: Point) -> bool:
        """
        This method determines whether the segment formed by the first two Point objects
        intersects the segment formed by the latter two Point objects.

        Arguments: p1(Point), p2(Point), p3(Point), p4(Point)

        Returns: bool

        Complexity: O(1)
        """
        if not all(isinstance(p, Point) for p in [p1, p2, p3, p4]):
            raise ValueError("The four points must be of type 'Point'.")

        intersection = self.__intersection_point(p1, p2, p3, p4)

        if intersection is None:
            return False

        intersect_x, intersect_y = intersection.get_latitude(), intersection.get_longitude()

        is_within_first_segment = self.__is_within_segment(p1, p2, intersect_x, intersect_y)
        is_within_second_segment = self.__is_within_segment(p3, p4, intersect_x, intersect_y)

        return is_within_first_segment and is_within_second_segment

    @staticmethod
    def __intersection_point(p1: Point, p2: Point, p3: Point, p4: Point) -> Point | None:
        """
        This method determines the intersection point, if it exists, between the segment formed by the
        first two Point objects and the segment formed by the latter two Point objects.

        Arguments: p1(Point), p2(Point), p3(Point), p4(Point)

        Returns: Point

        Complexity: O(1)
        """
        a1 = p2.get_longitude() - p1.get_longitude()
        b1 = p1.get_latitude() - p2.get_latitude()
        c1 = a1 * (p1.get_latitude()) + b1 * (p1.get_longitude())

        a2 = p4.get_longitude() - p3.get_longitude()
        b2 = p3.get_latitude() - p4.get_latitude()
        c2 = a2 * (p3.get_latitude()) + b2 * (p3.get_longitude())

        det = a1 * b2 - a2 * b1
        if abs(det) < 1e-9:
            return None

        x = (b2 * c1 - b1 * c2) / det
        y = (a1 * c2 - a2 * c1) / det
        return Point("intersection", x, y)

    def crosses_delimitation(self, p1: Point, p2: Point) -> bool:
        """
        This method determines whether the segment formed by the two points intersects the
        delimitation. If the new segment ends in the first point of the delimitation,
        while not intersecting any other segment belonging to the delimitation, this is not
        considered to be a crossing.

        Arguments: p1(Point), p2(Point)

        Returns: bool

        Complexity: O(p), p being the number of points
        """
        if not all(isinstance(p, Point) for p in [p1, p2]):
            raise ValueError("Both points must be of type 'Point'.")

        if self.size() == 0:
            raise ValueError("Delimitation object is empty.")

        points = self.get_points()
        for i in range(1, len(points)):
            p3 = points[i - 1]
            p4 = points[i]

            if self.intersects(p1=p3, p2=p4, p3=p1, p4=p2):
                intersection = self.__intersection_point(p1=p3, p2=p4, p3=p1, p4=p2)
                intersect_x, intersect_y = intersection.get_latitude(), intersection.get_longitude()  # type: ignore
                first_point = self.get_first()
                last_point = self.get_last_two()[1]

                if (intersect_x == first_point.get_latitude() and intersect_y == first_point.get_longitude()) or (  # type: ignore
                    intersect_x == last_point.get_latitude() and intersect_y == last_point.get_longitude()
                ):
                    continue
                return True

        return False

    def show(self, points: ValidPoints):
        """
        This method displays a window where all the stored points are displayed along
        with their identifier, and the complete delimitation formed by lines in between
        (possibly a subset of the) points is plotted.

        Arguments: points(ValidPoints)

        Complexity: O(p), p being the number of points in the Delimitation object
        """
        if not isinstance(points, ValidPoints):
            raise ValueError("The stored points must be of type 'ValidPoints'.")

        coords = []
        for point in points.get_all_points():
            coords.append((point.get_proj_x(), point.get_proj_y(), point.get_id()))

        _, ax = plt.subplots()

        if coords:
            latitudes, longitudes, point_ids = zip(*coords)
            ax.scatter(latitudes, longitudes, color="red")

            for i, point_id in enumerate(point_ids):
                ax.annotate(
                    point_id,
                    (latitudes[i], longitudes[i]),
                    textcoords="offset points",
                    xytext=(0, 5),
                    ha="center",
                )

        if self.size() > 1:
            delim_latitudes = [point.get_proj_x() for point in self.get_points()]
            delim_longitudes = [point.get_proj_y() for point in self.get_points()]
            ax.plot(delim_latitudes, delim_longitudes, color="blue", linewidth=1)

            if len(delim_latitudes) > 2:
                first_point = self.get_first()
                last_point = self.get_last_two()[1]
                ax.plot(
                    [first_point.get_proj_x(), last_point.get_proj_x()],  # type: ignore
                    [first_point.get_proj_y(), last_point.get_proj_y()],  # type: ignore
                    color="blue",
                    linewidth=1,
                )

        ax.set_title("Delimitation")
        ax.set_xlabel("Latitude")
        ax.set_ylabel("Longitude")
        plt.show()

    def __eq__(self, other: "Delimitation") -> bool:
        """
        This method checks if the Delimitation objects that are being compared
        represent the same complete polygon.

        Arguments: other(Delimitation)

        Returns: bool

        Complexity: O(p^2), p being the number of points in the Delimitation object
        """
        if not isinstance(other, Delimitation):
            return False

        points1 = self.get_points()
        points2 = other.get_points()
        set1 = set(points1)
        set2 = set(points2)
        if set1 == set2:
            return self.__check_rotations(points1, points2)

        return False

    @staticmethod
    def __check_rotations(list1: list, list2: list) -> bool:
        """
        This static method checks if list2 is a rotation of list1, considering both
        clockwise and counterclockwise rotations.

        Arguments: list1 (list), list2 (list)

        Returns: bool

        Complexity: O(p^2), p being the number of points in the lists
        """
        double_list1 = list1 + list1
        for i in range(len(list1)):
            if double_list1[i : i + len(list1)] == list2:
                return True

        reverse_list2 = list2[::-1]
        for i in range(len(list1)):
            if double_list1[i : i + len(list1)] == reverse_list2:
                return True

        return False

    def is_valid_delimitation(self, valid_points: ValidPoints, distance: int) -> bool:
        """
        Validates if the current delimitation formed by the points of the object is
        valid based on a set of conditions.

        Parameters: valid_points(ValidPoints), distance(int)

        Returns: bool

        Complexity: O(p^2), p being being the number of points in the Delimitation object
        """
        del_points = self.get_points()
        v_p_points = valid_points.get_all_points()
        if set(del_points).issubset(set(v_p_points)):
            if self.get_first().distance(self.get_last_two()[1]) > distance:  # type: ignore
                return False

            for i in range(self.size() - 1):
                if del_points[i].distance(del_points[i + 1]) <= distance:
                    continue
                else:
                    return False

            for i in range(len(del_points)):
                p1 = del_points[i]
                p2 = del_points[(i + 1) % len(del_points)]

                for j in range(len(del_points)):
                    p3 = del_points[j]
                    p4 = del_points[(j + 1) % len(del_points)]
                    if j == i or (j + 1) % len(del_points) == i or (j == (i + 1) % len(del_points)):
                        vertices = [(vertex.get_latitude(), vertex.get_longitude()) for vertex in (p1, p2, p3, p4)]
                        if (
                            self.intersects(p1, p2, p3, p4)
                            and (
                                self.__intersection_point(p1, p2, p3, p4).get_latitude(),  # type: ignore
                                self.__intersection_point(p1, p2, p3, p4).get_longitude(),  # type: ignore
                            )
                            not in vertices
                        ):
                            return False
                        continue
                    if self.intersects(p1, p2, p3, p4):
                        return False

            return True
        return False

    def is_point_inside(self, point: Point) -> bool:
        """
        Determines if a given point is inside the polygon formed by the delimitation points.

        Args:
            point (Point): The point to check.

        Returns:
            bool: True if the point is inside the polygon, False otherwise.

        Complexity: O(p), p being the number of points in the Delimitation object
        """
        n = len(self.get_points())
        if n < 3:
            return False

        x, y = point.get_longitude(), point.get_latitude()
        inside = False

        for i in range(n):
            p1 = self.get_points()[i]
            p2 = self.get_points()[(i + 1) % n]

            x1, y1 = p1.get_longitude(), p1.get_latitude()
            x2, y2 = p2.get_longitude(), p2.get_latitude()

            # Check if point is within y-range of edge and adjust ray crossing
            if (y1 > y) != (y2 > y):
                x_intersection = (x2 - x1) * (y - y1) / (y2 - y1) + x1
                if x < x_intersection:
                    inside = not inside

        return inside

    def __repr__(self) -> str:
        """
        This method returns a representation of the Delimitation object.

        Return: str

        Complexity: O(p), p being the number of points in the Delimitation object
        """
        points_repr = ", ".join([repr(point) for point in self.get_points()])
        return f"Delimitation([{points_repr}])"

    def __str__(self):
        """
        Returns a string representation of the object

        Return: str

        Complexity: O(p), p being the number of points in the Delimitation object
        """
        points_repr = ", ".join([repr(point) for point in self.get_points()])
        return f"[{points_repr}]"
        return f"[{points_repr}]"
